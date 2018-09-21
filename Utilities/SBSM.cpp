#include "SBSM.h"


SBSM::SBSM (size_t order, size_t nSamples, size_t nBreaks){
	baseMatrix = NULL;
	baseMatrix_t = NULL;
	btbMatrix = NULL;
	bsw = NULL;
	time = NULL;
	resize(order,nSamples,nBreaks);
}

void SBSM::resize(size_t order, size_t nSamples, size_t nBreaks){

	this->destroyBaseMatrix();
	this->order = order;
	this->nSamples = nSamples;

	time = gsl_vector_alloc(nSamples);

	for (size_t i =0; i<nSamples;i++)
	    gsl_vector_set(time, i, (double) i);
	if (nBreaks >= nSamples)
	      nBreaks = nSamples +order + 1;

	this->nBreaks = nBreaks;
	this->nSpline = nBreaks + order -2;

	createBaseMatrix();
}

void SBSM::createBaseMatrix(){

	cs *tmpbMatrix =cs_spalloc(nSpline, nSamples, 5*nSpline,1,1);
	cs *tmpbMatrix_t = cs_spalloc(nSamples, nSpline, 5*nSpline,1,1);

	int kstart = -1;
	int kend = -1;
	bsw = gsl_bspline_alloc(order, nBreaks);
	gsl_vector  *btdata = gsl_vector_alloc(nSpline);
	double c_sample = 0;
	double tol = 0;

	gsl_bspline_knots_uniform(0.0, (double)(nSamples-1), bsw);
	//Now create the matrix
	for (long i=0; i<(long)nSamples; i++){
	    gsl_bspline_eval(gsl_vector_get(time,i), btdata, bsw);
	    kstart = -1;
	    kend = -1;
	    for (int k=0; k<(int)nSpline; k++){
	      c_sample = gsl_vector_get(btdata, k);
	      if (kstart == -1 && c_sample > tol){
		kstart = k;
		continue;
	      }
	      if (kend ==-1 && kstart > -1 &&  c_sample <=tol )
		kend = k;
	    }
	    if (kend == -1 && kstart != -1)
	       kend = nSpline;
	    for (long k=kstart; k<kend; k++){
	      c_sample = gsl_vector_get(btdata, k);
	      cs_entry(tmpbMatrix, k, i , c_sample);
	      cs_entry(tmpbMatrix_t,i , k , c_sample);
	    }
	  }
	//Store it in sparse compress form

	baseMatrix = cs_compress(tmpbMatrix);
	baseMatrix_t = cs_compress (tmpbMatrix_t);
	btbMatrix = cs_multiply (baseMatrix,baseMatrix_t);
	cs_spfree(tmpbMatrix);
	cs_spfree(tmpbMatrix_t);

	gsl_vector_free(btdata);
}

VecDoub SBSM::fitData (double *dataVector, size_t nSamples){

	VecDoub tmpData(nSamples,0.0);
	double *tsi = new double [nSamples];
	double *coeff = new double [nSpline];

	for (size_t i=0; i<nSamples; i++){
	    	tmpData[i] = dataVector[i];
	    	tsi[i] = 0.0;
	}
	for (size_t i=0; i<nSpline; i++)
		coeff[i] = 0.0;
	cs_gaxpy (baseMatrix,tmpData.getData(),coeff);
	if (!cs_qrsol(3, btbMatrix, coeff)){
		cerr<<"SBSM::fitData(). Cannot create data B-Spline template."<<endl;
		exit(-1);
	}


	cs_gaxpy(baseMatrix_t, coeff, tsi);

	for (size_t i=0; i<nSamples; i++){
		tmpData[i]=tsi[i];
	}
	delete [] tsi;
	delete [] coeff;
	return tmpData;
}

VecDoub SBSM::fitData (VecDoub dataVector){
	return fitData(dataVector.getData(), dataVector.size());
}

void SBSM::destroyBaseMatrix(){
	if (baseMatrix)
		cs_spfree(baseMatrix);
	if (baseMatrix_t)
		cs_spfree(baseMatrix_t);
	if (btbMatrix)
		cs_spfree (btbMatrix);
	if (bsw)
		gsl_bspline_free(bsw);
	if (time)
		gsl_vector_free(time);
}

SBSM::~SBSM(){
	destroyBaseMatrix();
}
