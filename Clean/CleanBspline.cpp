
#include <iostream>
#include <cmath>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <omp.h>
#include <fftw3.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fit.h>
#include <assert.h>
#include <libgen.h>




using namespace std;

#include "CleanBspline.h"
#include "CleanPCA.h"
#include "CleanSelector.h"
#include "vector_utilities.h"
#include "sparseUtilities.h"
#include "AzElTemplateCalculator.h"

CleanBspline::CleanBspline(Array* dataArray, Telescope* tel) :
  Clean(dataArray,tel){
  baseMatrix = NULL;
  baseMatrix_t = NULL;
  bsw = NULL;
  time = NULL;
  baseMatrixSamples=0;
  baseMatrixDetectors=0;
  calibrated = 0;
  detPerHextant=NULL;
  nDetPerHextant.resize(0);
  fullMedianValues.resize(0);
  bright = 2.0/60.0;
  if (ap->getObservatory().compare("LMT")==0)
	  bright = 0.5/60.0;
  resample = ap->getResample();
  maxHex = 0;
  fixedData = NULL;
  totSamples = 0;
  scanStatus = NULL;
  currScan = 0;

#if defined(_OPENMP)
  tid = omp_get_thread_num() + 1;
#else
  tid = 0;
#endif
  if (resample <1){
	  cout<< "CleanBspline "<< tid<<" ). resample value is less than 1. Set to use all samples instead"<<endl;
	  resample = 1;
  }


}


void CleanBspline::removeBadBolos(){

	size_t nDetectors = dataArray->getNDetectors();
	size_t totSamples =  dataArray->detectors[0].hValues.size();

	int *di = dataArray->getDetectorIndices();

	VecDoub avg (totSamples,0.0);
	VecDoub detMean (nDetectors);
	VecInt count(totSamples, 0.0);
	for (size_t i=0; i<nDetectors; i++){
		detMean[i] = mean(&dataArray->detectors[di[i]].hValues[0], &dataArray->detectors[di[i]].hSampleFlags[0], totSamples);
		for (size_t j=0; j<totSamples; j++){
			if (dataArray->detectors[di[i]].hSampleFlags[j]){
				avg[j]+=(dataArray->detectors[di[i]].hValues[j]-detMean[i]);
				count[j]++;
			}
		}
	}
	for (size_t j=0; j<totSamples; j++)
		if (count[j]>0)
			avg[j]/=count[j];
		else
			avg[j]=0;


	VecDoub detPow(nDetectors,0.0);

	for (size_t i=0; i<nDetectors; i++){
		for (size_t j=0; j<totSamples; j++){
			if (dataArray->detectors[di[i]].hSampleFlags[j])
				detPow[i] += pow((dataArray->detectors[di[i]].hValues[j]-detMean[i])-avg[j], 2.0);
		}
	}

	double mPow = median(&detPow[0],nDetectors);
	size_t badBoloCount = 0;
	for (size_t i=0; i<nDetectors; i++)
		if (detPow[i]> 3.0*mPow){
			dataArray->detectors[di[i]].goodFlag = 0;
			badBoloCount++;
		}
	cout<<"CleanBspline()::removeBadBolo: Taking out "<<badBoloCount << "bolometers"<<endl;
	if (badBoloCount>0)
		dataArray->updateDetectorIndices();
}


bool CleanBspline::clean(){
	double cleanStripe = dataArray->getAp()->getCleanStripe();

	cleanStripe = abs(cleanStripe);

	calibrate();
//	removeBadBolos();
	dataArray->fakeAtmData(false);
	totSamples =  dataArray->detectors[0].hValues.size();
	atmTemplate = VecDoub (totSamples*dataArray->getNDetectors(),0.0);



	if (ap->getOrder() != 100)
		this->cleanScans();
	int *di = dataArray->getDetectorIndices();
	for (size_t i=0; i<(size_t)dataArray->getNDetectors(); i++)
		for (size_t j=0; j<totSamples; j++){
			if (ap->getOrder() != 100)
				dataArray->detectors[di[i]].hValues[j] -= atmTemplate[i*totSamples + j];
			else
				dataArray->detectors[di[i]].hValues[j] -= atmTemplate[i*totSamples + j];
			dataArray->detectors[di[i]].atmTemplate[j] -= atmTemplate[i*totSamples + j];
		}

//	this->removeLargeScaleResiduals(STRIPE_SCAN,0.1);
	if (cleanStripe != 0){
		this->removeLargeScaleResiduals(translateStripeMethod(),cleanStripe);
	}



	//delete gainFixer;
	//delete [] fixedData;
	atmTemplate.resize(0);
	return true;
}




bool CleanBspline::cleanScans (){

	Array *dataArray = this->dataArray;
	size_t nScans = telescope->scanIndex.ncols();
	size_t nDetectors=dataArray->getNDetectors();

	int *di = dataArray->getDetectorIndices();

	bool debug = false;

	if (!scanStatus){
		scanStatus = new bool[nScans];
		for (size_t k=0; k< nScans; k++)
			scanStatus[k] = true;
	}

//	VecInt boloGroup;
	double cleanPixelSize = ap->getCleanPixelSize()/3600.0;

	if (cleanPixelSize <= 0){
		cerr<<"CleanBspline(): Fatal Error cleanPixelSize is zero. Please check your ap xml file"<<endl;
		exit(-1);
	}

	cout<<"CleanBspline("<<tid<<")::clean(). Starting cleaning for "<< nDetectors<< "bolometers. This could take some time...."<<endl;

	for (size_t k =0; k<nScans; k++){
		size_t ei = telescope->scanIndex[1][k]+1;
		size_t si = telescope->scanIndex[0][k];
		size_t nSamples = ei-si;

//		size_t nDetectors = 0;
//		for (size_t ihex = 5; ihex <=6; ihex++){
//			char buff [10];
//			int *di;
//			sprintf(buff, "h%lu",ihex);
//			nDetectors = 0;
//			for (size_t ibb = 0; ibb<totDetectors; ibb++)
//				if (dataArray->detectors[diFull[ibb]].getName().find(buff) != string::npos)
//					nDetectors++;
//			di = new int [nDetectors];
//			size_t idi = 0;
//			for (size_t ibb = 0; ibb<totDetectors; ibb++)
//							if (dataArray->detectors[diFull[ibb]].getName().find(buff) != string::npos)
//								di[idi++]=diFull[ibb];

//		int *di;

		double increment = 1.0;
		size_t index;
		size_t oSamples=nSamples;

		if (resample > 1){
			long nnSamples = long(round(nSamples/resample));
			increment = double(nSamples-1)/double(nnSamples-1);
			nSamples =nnSamples;
		}

		size_t dataLen = nDetectors*nSamples;

		double *dataVector = new double [dataLen];
		gsl_vector *raVector = gsl_vector_alloc(dataLen);
		gsl_vector *decVector = gsl_vector_alloc(dataLen);
		bool *flags = new bool [dataLen];
		double *atmTemplate =NULL;
		double dist = 0.0;
		//cout<<"CleanBspline("<<tid<<")::clean(). Starting cleaning on scan: "<<" " <<k <<" of "<<nScans-1. <<" (Nbolo = "<<nDetectors<<")"<<endl;
		//Copy data to vectors

		for(size_t i=0; i<nDetectors; i++){
			for(size_t j=0; j<nSamples; j++){
				index = long(round(j*increment));
				dataVector[i*nSamples +j]=dataArray->detectors[di[i]].hValues[si+index];
				gsl_vector_set(raVector,i*nSamples+j, dataArray->detectors[di[i]].azElRaPhys[si+index]);
				gsl_vector_set(decVector,i*nSamples+j, dataArray->detectors[di[i]].azElDecPhys[si+index]);
				flags[i*nSamples+j] = (bool) dataArray->detectors[di[i]].hSampleFlags[si+index];
				dist = sqrt (pow(dataArray->detectors[di[i]].azElRaPhys[si+index],2.0)+pow(dataArray->detectors[di[i]].azElDecPhys[si+index],2.0));
				if (dist<this->bright)
					flags[i*nSamples+j] = false;
			}
		}

		if (debug){
			writeVecOut ("oData.txt", dataVector, nSamples*nDetectors);
		}

		if (fixFlags(dataVector,flags,(int)nDetectors,(int)nSamples)){
			atmTemplate = cottingham(dataVector, raVector, decVector, flags, nDetectors, nSamples, cleanPixelSize, 0);

			if (!atmTemplate){
				cout<<"CleanBspline("<<tid<<")::cleanScans(). Failed to produce and atmosphere template. Setting scan to 0.0"<<endl;
				for(size_t i=0; i<nDetectors; i++){
					for(size_t j=0; j<oSamples; j++){
						dataArray->detectors[di[i]].hValues[si+j]=0.0;
						dataArray->detectors[di[i]].hSampleFlags[si+j]=0.0;
					}
				}
				scanStatus[k] = false;
			}else{

				if (debug){
					writeVecOut ("residualData.txt", dataVector, dataLen);
					writeVecOut ("template1.txt", atmTemplate, dataLen);
					exit(-1);
				}
				for(size_t i=0; i<nDetectors; i++){
					//subtractTemplate(&dataArray->detectors[di[i]].hValues[si],oSamples,&atmTemplate[i*nSamples],nSamples,increment, NULL);
					subtractTemplate (&this->atmTemplate[i*totSamples + si], oSamples, &atmTemplate[i*nSamples], nSamples, increment, NULL, true);
				}

				delete [] atmTemplate;
			}
		}
		else{
			cout<<"CleanBspline("<<tid<<")::cleanScans(). Failed to produce and atmosphere template. Setting scan to 0.0"<<endl;
							for(size_t i=0; i<nDetectors; i++){
								for(size_t j=0; j<oSamples; j++){
									dataArray->detectors[di[i]].hValues[si+j]=0.0;
									dataArray->detectors[di[i]].hSampleFlags[si+j]=0.0;
								}
							}
		}

		//delete [] di;
		delete [] dataVector;
		delete [] flags;
		gsl_vector_free(raVector);
		gsl_vector_free(decVector);
		}



	return true;

}

CleanBsplineDestriping CleanBspline::translateStripeMethod(){

	if (ap->stripeMethod =="fft")
		return STRIPE_FFT;
	if (ap->stripeMethod =="pca")
		return STRIPE_PCA;
	if (ap->stripeMethod == "spline")
		return SUBARRAY_SPLINE;
	if (ap->stripeMethod == "scan")
		return STRIPE_SCAN;
	if (ap->stripeMethod == "azel" )
		return AZEL_TEMPLATE;

	return STRIPE_NONE;
}

double * CleanBspline::cottingham(double *dataVector, gsl_vector *raVector, gsl_vector *decVector, bool *flags, size_t nDetectors, size_t nSamples, double cleanPixelSize, int refBolo){

	size_t dataLen = nDetectors*nSamples;
	cs *ptMatrix = NULL;
	cs *ptMatrix_inverse = NULL;
	cs *pMatrix = NULL;
	cs *pMatrix_t = NULL;
	cs *tta = NULL;
	//Temporary matrixes
	cs *a1 =NULL;
	cs *a2 = NULL;
	cs *a3 =NULL;
	cs *bMatrix = NULL;
	//Temporary vector for data side
	double *v1 = NULL;
	double *v2 = NULL;
	//Projection side matrix
	cs *phi = NULL;
	//Data side vector
	double *tsi = NULL;

	size_t nP = 0;
	size_t nSp = 0;
	size_t rBoloIndex = 0;

	if (refBolo < 0)
		rBoloIndex = dataArray->getRefBoloIndex();
	else
		rBoloIndex = (size_t) refBolo;

	MatDoub relCoeff (nDetectors, 2);
	VecDoub energy (nDetectors,0.0);
	VecDoub meanTod (nSamples,0.0);
//	writeVecOut("oDataVector.txt", dataVector, nDetectors*nSamples);
	for (size_t id=0; id<nDetectors; id++){
		if (linfit_flags(&dataVector[rBoloIndex*nSamples], &flags[rBoloIndex*nSamples], &dataVector[id*nSamples], &flags[id*nSamples], nSamples,&relCoeff[id][0], &relCoeff[id][1],false)<0){
			cerr<<"Data on file "<<dataArray->getAp()->getDataFile()<< "has no valid data on this scan"<<endl;
			return NULL;
		}
		for (size_t is = 0; is <nSamples; is++){
			dataVector[id*nSamples+is] = (dataVector[id*nSamples+is]-relCoeff[id][0])/relCoeff[id][1];
			meanTod[is] += dataVector[id*nSamples+is]/nDetectors;
		}
	}

//		writeVecOut("newDataVector.txt", dataVector, nDetectors*nSamples);
//		writeVecOut("newFlagVector.txt", flags, nDetectors*nSamples);
		//exit(-1);
	for (size_t id=0; id<nDetectors; id++)
		for (size_t is = 0; is <nSamples; is++)
			energy[id]+= pow((dataVector[id*nSamples+is]-meanTod[is]),2.0);

//	int *di = dataArray->getDetectorIndices();
//	size_t *ide = new size_t [nDetectors];
//	gsl_sort_index(ide, &energy[0],1,nDetectors);
//
//	char buff [200];
//	sprintf(buff,"%s_badBolo.txt",ap->getMapFile().c_str());
//
//	ofstream dfile(buff,std::ios::app);
//	//dfile.open
//	cout <<"Write bad bolo information"<<endl;
//	for (size_t ibb=0; ibb<10; ibb++)
//		dfile<<dataArray->detectors[di[ide[nDetectors-ibb-1]]].getName()<<",";
//	dfile<<endl;
//	dfile.close();

//

//	exit(-1);



	this->createBaseMatrix(nSamples, nDetectors);
	nSp = this->baseMatrix->n;
	pMatrix = this->getPMatrix(raVector, decVector, cleanPixelSize);
	nP = pMatrix->n;

	if (!pMatrix){
		cerr<<"Problem creating sparse pointing matrix "<<endl;
		exit(-1);
	}
	pMatrix_t = cs_transpose(pMatrix,1);
	if (!pMatrix_t){
		cerr<<"Problem creating sparse pointing matrix transpose "<<endl;
		exit(-1);
	}
	ptMatrix = cs_multiply(pMatrix_t, pMatrix);
	if (!ptMatrix){
		cerr<<"Problem creating sparse projection matrix"<<endl;
		exit(-1);
	}
	ptMatrix_inverse = fast_sparse_invert(ptMatrix);
	if (!ptMatrix_inverse) {
		cerr<<"Problem creating sparse projection matrix inverse"<<endl;
		exit(-1);
	}
	//cs_droptol(ptMatrix_inverse, 1e-7);
	tta = cs_multiply(ptMatrix_inverse, pMatrix_t);
	if (!tta){
		cerr<<"Problem creating final sparse projection matrix "<<endl;
		exit(-1);
	}

	a1 = cs_multiply(baseMatrix_t, pMatrix);
	if (!a1){
		cerr<<"A1 matrix error"<<endl;
		exit(-1);
	}
	a2 = cs_multiply(tta, baseMatrix);
	if (!a2){
		cerr<<"a2 matrix error"<<endl;
	}
	bMatrix = cs_multiply(baseMatrix_t, baseMatrix);
	if (!bMatrix){
		cerr<<"bMatrix creation matrix error"<<endl;
		exit(-1);
	}
	a3 = cs_multiply(a1,a2);
	if (!a3){
		cerr<<"a3 creation matrix error"<<endl;
		exit(-1);
	}
	phi = cs_add(bMatrix, a3, 1.0, -1.0);
	if (!phi){
		cerr<<"phi creation matrix error"<<endl;
		exit(-1);
	}
	//cerr<<"Point Spline projection done. Starting data side"<<endl;

	cs_spfree(a2);
	cs_spfree(a3);
	cs_spfree(bMatrix);
	a2=NULL;
	a3=NULL;
	bMatrix = NULL;
	//Now compute data side matrix
	v1 = new double [nP];
	for (register size_t idata=0; idata < (size_t)nP; idata++)
		v1[idata] = 0.0;
	tsi = new double [nSp];
	for (register size_t idata=0; idata < (size_t)nSp; idata++)
		tsi[idata] = 0.0;
	v2 = new double [nSp];
	for (register size_t idata=0; idata < (size_t)nSp; idata++)
		v2[idata] = 0.0;
	if (!cs_gaxpy(tta, dataVector,v1)){
		cerr<<"Error creating v1 vector"<<endl;
		exit(-1);
	}
	if (!cs_gaxpy(baseMatrix_t, dataVector, tsi)){
		cerr<<"Error creating tsi vector"<<endl;
		exit(-1);
	}
	if (!cs_gaxpy(a1,v1,v2)){
		cerr<<"Error creating v2 vector"<<endl;
		exit(-1);
	}
	for (register size_t itsi = 0; itsi<nSp; itsi++){
		if (!finite(v2[itsi]) || !finite(tsi[itsi])){
			cerr<<"Nan detected on tsi vector"<<endl;
			exit(-1);
		}
		tsi[itsi]-=v2[itsi];
	}
	delete [] v1;
	delete [] v2;
	cs_spfree(a1);
	cs_spfree(tta);
	cs_spfree(pMatrix);
	cs_spfree(pMatrix_t);
	cs_spfree(ptMatrix);
	cs_spfree(ptMatrix_inverse);
	v1= NULL;
	v2= NULL;
	a1= NULL;
	tta = NULL;
	//Now solve linear system
	cs_dropnotfinite(phi);
	if (!cs_qrsol(3,phi,tsi)){
		cerr<<"CleanBspline(): Cannot solve Spline system for this observation"<<endl;
		exit(-1);
	}


	v1 = new double [dataLen];
	for (register size_t id = 0; id<dataLen;id++)
		v1[id] = 0.0;
	if (!cs_gaxpy(baseMatrix, tsi, v1)){
		cerr<<"CleanBspline():Cannot create atm template from Spline solution"<<endl;
		exit(-1);
	}


	delete [] tsi;
	cs_spfree(phi);
	tsi=NULL;
	phi = NULL;

	//Check for NaN and finite
	for (size_t i=0; i<nDetectors;i++)
		for(size_t j=0; j<nSamples; j++)
			if (!finite(v1[i*nSamples+j])){
				delete [] v1;
				return NULL;
			}


	double  *atmTemplate = new double[dataLen];
	bool  *cFlags = new bool [nSamples];
	double c0,c1;
	double pointDist;
	for (size_t i=0; i<nDetectors; i++){
		for (size_t j=0; j<nSamples; j++){
			pointDist = sqrt(pow(gsl_vector_get(raVector, i*nSamples+j),2)+pow(gsl_vector_get(decVector, i*nSamples+j),2));
			cFlags[j] = flags[i*nSamples+j]; // & (pointDist < bright);
		}
		linfit_flags(&dataVector[i*nSamples],cFlags, &v1[i*nSamples],cFlags, nSamples, &c0,&c1, true);
		//gsl_fit_linear(&dataVector[i*nSamples],1,&v1[i*nSamples],1,nSamples, &c0,&c1,&cov00,&cov01,&cov11,&chisq);

		for (size_t j=0; j<nSamples; j++){
			atmTemplate [i*nSamples+j]= (v1[i*nSamples+j]-c0)/c1;
			atmTemplate [i*nSamples+j] *= relCoeff[i][1];
			atmTemplate [i*nSamples+j] += relCoeff[i][0];
		}

	}
	delete [] cFlags;

//	writeVecOut("newAtmTemplate.txt", atmTemplate, nDetectors*nSamples);
//	exit(-1);

	delete []v1;
	return atmTemplate;

}


bool CleanBspline::removeScanPattern (double *dataVector, gsl_vector *azVector, gsl_vector *elVector,
		double *flags, size_t nDetectors, size_t nSamples){
	size_t ngood = 0;
	size_t nPar = 7;
	for (size_t i=0; i<nDetectors*nSamples;i++)
		if (flags[i])
			ngood++;

	gsl_matrix *S = gsl_matrix_alloc(ngood,nPar);
	gsl_matrix *SS = gsl_matrix_alloc(nPar,nPar);
	gsl_vector *d = gsl_vector_alloc(ngood);
	gsl_vector *Std = gsl_vector_alloc (nPar);
	gsl_vector *solution = gsl_vector_alloc(nPar);

	size_t igood = 0;
	size_t iPar;
	double az, el;
	double sr = dataArray->detectors[0].getSamplerate();
	for (size_t i =0; i< nDetectors; i++)
		for (size_t j=0; j<nSamples; j++){
			if (flags[i*nSamples+j]){
				iPar = 0;
				gsl_matrix_set(S,igood,iPar++, 1.0);
				gsl_matrix_set(S,igood, iPar++, double(j)/sr);
				az = gsl_vector_get(azVector,i*nSamples+j);
				el = gsl_vector_get (elVector,i*nSamples+j);
				gsl_matrix_set (S,igood, iPar++, az);
				gsl_matrix_set (S, igood, iPar++, el);
				gsl_matrix_set (S,igood, iPar++, az*az);
				gsl_matrix_set (S,igood, iPar++, az*el);
				gsl_matrix_set (S,igood, iPar++, el*el);
				gsl_vector_set (d,igood, dataVector[i*nSamples+j]);
				igood++;
			}
		}
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0, S,S,0.0, SS);
	gsl_blas_dgemv(CblasTrans, 1.0,S,d, 0.0, Std);
	gsl_permutation *p = gsl_permutation_alloc(nPar);
	int signum = 0;

	gsl_linalg_LU_decomp(SS, p, &signum);
	gsl_linalg_LU_solve(SS, p, Std, solution);
	gsl_matrix_free (SS);
	gsl_permutation_free(p);
	gsl_matrix_free (S);

	double res;
	for (size_t i =0; i< nDetectors; i++)
			for (size_t j=0; j<nSamples; j++){
				iPar = 0;
				gsl_vector_set(Std,iPar++, 1.0);
				gsl_vector_set(Std, iPar++, double(j)/sr);
				az = gsl_vector_get(azVector,i*nSamples+j);
				el = gsl_vector_get (elVector,i*nSamples+j);
				gsl_vector_set (Std, iPar++, az);
				gsl_vector_set (Std, iPar++, el);
				gsl_vector_set (Std, iPar++, az*az);
				gsl_vector_set (Std, iPar++, az*el);
				gsl_vector_set (Std, iPar++, el*el);
				gsl_blas_ddot(Std,solution,&res);
				dataVector[i*nSamples+j] -=res;
			}
	gsl_vector_free (Std);
	gsl_vector_free(d);
	gsl_vector_free (solution);

	return true;
}


void CleanBspline::removeLargeScaleResiduals (CleanBsplineDestriping method,double correlate){

	size_t nScans = telescope->scanIndex.ncols();
	size_t nDetectors = dataArray->getNDetectors();
	int* di = dataArray->getDetectorIndices();
	size_t si, ei, nSamples;

	double *dataVector = NULL;
	double *kVector = NULL;
	double *flagVector = NULL;
	double *atmVector = NULL;
	double *AzOffsets = new double [nDetectors];
	double *ElOffsets = new double [nDetectors];
	double *detSens = new double [nDetectors];
	gsl_vector *azVector= NULL;
	gsl_vector *elVector= NULL;

	for (size_t id = 0; id < nDetectors; id ++){
		AzOffsets[id]=dataArray->detectors[di[id]].azOffset;
		ElOffsets[id]=dataArray->detectors[di[id]].elOffset;
		detSens[id] = 1.0/dataArray->detectors[di[id]].getSensitivity();
	}

	double sampleRate = dataArray->detectors[di[0]].getSamplerate();

	//if (method != SUBARRAY_SPLINE){
		for(size_t k=0;k<nScans;k++){
			currScan = k+1;
			if (ap->getOrder() != 100)
				if (!scanStatus[k])
					continue;
			//cout<<"CleanBspline("<<tid<<")::clean().High pass clean on  scan: "<<k <<" of "<<nScans-1<<endl;
			si=telescope->scanIndex[0][k];
			ei=telescope->scanIndex[1][k]+1;
			nSamples = ei-si;

			dataVector = new double [nDetectors*nSamples];
			kVector = new double [nDetectors*nSamples];
			flagVector = new double [nDetectors*nSamples];
			atmVector = new double [nDetectors*nSamples];
			azVector = gsl_vector_alloc(nDetectors*nSamples);
			elVector = gsl_vector_alloc (nDetectors*nSamples);


			VecDoub fullFlags (nSamples,1.0);
			for(size_t i=0; i<nDetectors; i++){
				for(size_t j=0; j<nSamples; j++){
					dataVector[i*nSamples +j]=dataArray->detectors[di[i]].hValues[si+j];
					kVector[i*nSamples + j]=dataArray->detectors[di[i]].hKernel[si+j];
					atmVector[i*nSamples + j]=dataArray->detectors[di[i]].atmTemplate[si+j];
					flagVector[i*nSamples +j] = (double)dataArray->detectors[di[i]].hSampleFlags[si+j];
					if (flagVector[i*nSamples+j] ==0)
						fullFlags[j] =0;
					if (azVector!=NULL){
						gsl_vector_set(azVector, i*nSamples+j,dataArray->detectors[di[i]].azElRaPhys[si+j]);
						gsl_vector_set(elVector, i*nSamples+j,dataArray->detectors[di[i]].azElDecPhys[si+j]);
					}
					double dist = sqrt (pow(dataArray->detectors[di[i]].azElRaPhys[si+j],2.0)+pow(dataArray->detectors[di[i]].azElDecPhys[si+j],2.0));
					if (dist<this->bright)
						flagVector[i*nSamples+j] = false;
				}
			}

			switch (method){
			case STRIPE_PCA:
				removeCorrelations2(dataVector,nDetectors,nSamples, correlate, AzOffsets,ElOffsets, kVector, flagVector);
				break;
			case STRIPE_FFT:
				//cleanStripeFFT2(dataVector, kVector, nDetectors, nSamples,correlate, 0.00,sampleRate);
				cleanStripeFFT3(dataVector, kVector, nDetectors, nSamples,correlate,sampleRate);
				break;
			case STRIPE_SCAN:
				removeScanPattern(dataVector, azVector,elVector, flagVector, nDetectors, nSamples);
				break;
			case AZEL_TEMPLATE:
				//if (ap->getObservatory()=="LMT")
					azelResidualLMT(dataVector,flagVector, nDetectors, nSamples, azVector, elVector,detSens);
				//else
					//azelResidual(dataVector,flagVector, nDetectors, nSamples, AzOffsets, ElOffsets,detSens);
				break;
			case SUBARRAY_SPLINE:
				splineResidual(dataVector, flagVector, nDetectors, nSamples, azVector, elVector, correlate);
				break;
			default:
				cerr<<"Not a valid destriping method"<<endl;
				exit(-1);

			}

			for(size_t i=0; i<nDetectors; i++)
				for(size_t j=0; j<nSamples; j++){
					dataArray->detectors[di[i]].hValues[si+j]=dataVector[i*nSamples +j];
					dataArray->detectors[di[i]].hKernel[si+j]=kVector[i*nSamples +j];
				}

			delete [] dataVector;
			delete [] kVector;
			delete [] atmVector;
			delete [] flagVector;
			gsl_vector_free(azVector);
			gsl_vector_free (elVector);
			azVector=NULL;
			elVector= NULL;
		}
//	}else{
//		cleanScans(CORRELATE, correlate);
//	}

	delete [] AzOffsets;
	delete [] ElOffsets;
	delete [] detSens;
}



void CleanBspline::createBaseMatrix(int nSamples, int nDetectors){
	double tol = 0.0;
	if (baseMatrix != NULL){
		if (nSamples == baseMatrixSamples && nDetectors == baseMatrixDetectors)
			return;
		else
			this->destroyBaseMatrix();
	}

	int order = dataArray->getAp()->getOrder();
	double controlChunk = dataArray->getAp()->getControlChunk();
	double timeChunk = dataArray->getAp()->getTimeChunk();
	int nbreaks =0;
	int nSpline = 0;
	cs *tmpbMatrix=NULL;
	cs *tmpbMatrix_t=NULL;

	if (controlChunk <= 0.0 || timeChunk <= 0.0){
		nbreaks = nSamples +order +1;
	}else{
		nbreaks = round(timeChunk/controlChunk/this->resample);
		if (nbreaks >= nSamples){
			nbreaks = nSamples +order + 1;
		}
	}
	nSpline = nbreaks + order -2;
	//cerr<<"CleanBspline("<<tid<<"): Control Chunk is "<<controlChunk<< " Number of Splines: "<<nSpline<<endl;

	long dataLen = nSamples * nDetectors;
	gsl_vector* tmpB = NULL;
	bsw = gsl_bspline_alloc(order, nbreaks);
	if (bsw == NULL){
		cerr<<"CleanBspline::createBaseMatrix():Cannot allocate B-spline working space. Imploding"<< endl;
		exit(-1);
	}
	baseMatrixSamples = nSamples;
	baseMatrixDetectors = nDetectors;
	gsl_bspline_knots_uniform(0.0, (double)(baseMatrixSamples-1), bsw);
	time = gsl_vector_alloc(baseMatrixSamples);
	tmpB = gsl_vector_alloc(nSpline);
	for (int i=0; i<baseMatrixSamples; i++)
		gsl_vector_set(time, i, (double) i);

	tmpbMatrix = cs_spalloc(nSpline, dataLen, 5*nSpline,1,1);
	tmpbMatrix_t = cs_spalloc(dataLen, nSpline,5*nSpline,1,1);
	//cerr<<"CleanBspline::createBaseMatrix(): Base Matrix Allocated"<<endl;
	int kstart = 0;
	int kend = 0;
	double c_sample=0.0;
	for (long i=0; i<baseMatrixSamples; i++){
		gsl_bspline_eval(gsl_vector_get(time,i), tmpB, bsw);
		//gsl_matrix_set_col(tmpMatrix, i, tmpB);
		kstart = -1;
		kend = -1;
		for (int k=0; k<nSpline; k++){
			c_sample = gsl_vector_get(tmpB, k);
			if (!finite(c_sample)){
				cerr<<"Nan detected on BaseMatrix....Imploding"<<endl;
				exit(-1);
			}
			if (kstart == -1 && abs(c_sample) > tol){
				kstart = k;
				continue;
			}
			if (kend ==-1 && kstart > -1 &&  abs(c_sample) <=tol )
				kend = k;
		}
		if (kend == -1 && kstart != -1)
			kend = nSpline;
		for (long k=kstart; k<kend; k++){
			c_sample = gsl_vector_get(tmpB, k);
			if (c_sample == 0)
				continue;
			for (long j=0; j<nDetectors; j++){
				cs_entry(tmpbMatrix,  k, j*baseMatrixSamples + i , c_sample);
				cs_entry(tmpbMatrix_t,j*baseMatrixSamples + i , k , c_sample);
			}
		}
	}

	gsl_vector_free(tmpB);
	this->baseMatrix = cs_compress(tmpbMatrix_t);
	this->baseMatrix_t = cs_compress(tmpbMatrix);
	cs_spfree(tmpbMatrix_t);
	cs_spfree(tmpbMatrix);

	if (!this->baseMatrix || !this->baseMatrix_t){
		cerr<<"CleanBspline(): Cannot allocate Bspline base matrix. Imploding"<<endl;
		exit(-1);
	}

	return;
}

void CleanBspline::calibrate(){
	if (!this->calibrated){
		int nDetectors = dataArray->getNDetectors();
		//Calibrate signals before
		double avgTau=0;
		for(int i=0;i<nDetectors;i++){
			dataArray->detectors[i].estimateTau();
			avgTau += dataArray->detectors[i].estimatedTau;
		}
		avgTau /= nDetectors;
		for(int i=0;i<nDetectors;i++){
			dataArray->detectors[i].estimateExtinction(avgTau);
			dataArray->detectors[i].calibrate();
		}
		this->calibrated = true;
	}
}

void CleanBspline::downSample(gsl_vector** out,double *data, bool *flags, long nSamples, long downSample){

	int addOne =0;
	size_t newSample =0;
	if ((nSamples % downSample) != 0){
		cerr<<"CleanBspline::downSample() Warning number of samples is not exactly divisible by downsampling rate. Fixing" <<endl;
		addOne = 1;
	}
	newSample = nSamples/downSample + addOne;
	if (flags == NULL){
		flags = new bool[nSamples];
		for (long i =0; i<nSamples; i++)
			flags[i] = true;
	}
	if (*out == NULL)
		*out = gsl_vector_alloc(newSample);
	else{
		if ((*out)->size != newSample){
			cout<<"CleanBspline::downSample(). Warning. Output vector dimension mismatch. Deallocating"<<endl;
			gsl_vector_free(*out);
			*out=gsl_vector_alloc(newSample);
		}
	}

	for(size_t i =0, j=0; i<newSample-addOne; i++, j+=downSample)
		gsl_vector_set(*out, i, data[j]);

	if (addOne == 1)
		gsl_vector_set(*out, (*out)->size-1, data[nSamples-1]);

}


void CleanBspline::destroyBaseMatrix(){
	if (baseMatrix != NULL){
		cs_spfree(baseMatrix);
		cs_spfree(baseMatrix_t);
		gsl_bspline_free(bsw);
		gsl_vector_free(time);
		baseMatrixDetectors = 0;
		baseMatrixSamples = 0;
		baseMatrix = NULL;
		baseMatrix_t = NULL;
		bsw=NULL;
	}
}

bool CleanBspline::removeCorrelations (double *dataVector, double *flagVector, size_t nDetectors, size_t nSamples, double corrFactor, double *azoffset, double *eloffset, gsl_vector *azVector, gsl_vector *elVector){
	size_t totalSamples = nDetectors*nSamples;
	VecDoub outputVector (totalSamples);
	MatDoub corrMatrix (nDetectors,nDetectors,0.0);
	MatDoub signMatrix (nDetectors, nDetectors, 0.0);

	double dist;
	double minDist = 0.0*1.25*60.0/2.0;
	if (ap->getObservatory()=="LMT")
		minDist /=3.0;
	//Build correlation matrix
	for (size_t i=0; i<nDetectors;i++){
		corrMatrix[i][i] = 1.0;
		signMatrix[i][i] = 1.0;
		for (size_t j=i+1; j<nDetectors; j++){
			corrMatrix[i][j] = flagCorrelation(&(dataVector[i*nSamples]), &(flagVector[i*nSamples]), \
					&(dataVector[j*nSamples]), &(flagVector[j*nSamples]),\
					nSamples);
			dist = sqrt(pow(azoffset[i]-azoffset[j],2) + pow(eloffset[i]-eloffset[j],2));
			if (corrMatrix[i][j] < corrFactor || dist < minDist)
				corrMatrix[i][j]=0.0;
			corrMatrix[j][i]=corrMatrix[i][j];
			if (corrMatrix[i][j] >=0){
				signMatrix[i][j] = 1.0;
				signMatrix[j][i] = 1.0;
			}else{
				signMatrix[i][j] = signMatrix[j][i] = 1.0;
			}

		}
	}
	size_t nCorrBolo=0;
	bool corrRemoved = false;
	VecBool cleaned (nDetectors, false);
	for (size_t i=0; i<nDetectors; i++){
		nCorrBolo = 0;

		for (size_t j=0; j<nDetectors; j++)
			if (corrMatrix[i][j]!= 0  && cleaned[j] ==0)
				nCorrBolo++;
		//cout<<"Bolometer "<<i<< "has "<<nCorrBolo<<"correlations"<<endl;
		if (nCorrBolo >=2){
			corrRemoved = true;
			double *curData = new double [nCorrBolo*nSamples];
			bool *curFlags = new bool [nCorrBolo*nSamples];
			gsl_vector *curAz = gsl_vector_alloc (nCorrBolo*nSamples);
			gsl_vector *curEl = gsl_vector_alloc (nCorrBolo*nSamples);

			VecInt boloIndex (nCorrBolo,-1);
			double c0, c1;
			boloIndex[0]=i;
			size_t iIndex = 1;
//			double icov = 0;
			for (size_t ibolo = 0; ibolo<nDetectors; ibolo++)
				if (ibolo != i && corrMatrix[i][ibolo] !=0)
					boloIndex[iIndex++]=ibolo;

			for (size_t ibolo = 0; ibolo<nCorrBolo; ibolo++){
				linfit_flags(&dataVector[boloIndex[0]*nSamples],&flagVector[boloIndex[0]*nSamples] ,&dataVector[boloIndex[ibolo]*nSamples],&flagVector[boloIndex[ibolo]*nSamples] ,nSamples, &c0, &c1, false);
				for (size_t iSample = 0; iSample<nSamples; iSample++){
					curData[ibolo*nSamples+iSample]= (dataVector[boloIndex[ibolo]*nSamples + iSample]-c0)/c1;
					curFlags[ibolo*nSamples+iSample] = (bool) flagVector[boloIndex[ibolo]*nSamples + iSample];
					gsl_vector_set (curAz,ibolo*nSamples+iSample,gsl_vector_get(azVector,boloIndex[ibolo]*nSamples + iSample));
					gsl_vector_set (curEl,ibolo*nSamples+iSample, gsl_vector_get(elVector,boloIndex[ibolo]*nSamples + iSample));
				}
			}
//			writeVecOut ("idata.txt", curData, nCorrBolo*nSamples);
//			writeVecOut ("iaz.txt", curAz->data, nCorrBolo*nSamples);
//			writeVecOut ("iel.txt", curEl->data, nCorrBolo*nSamples);
//			cout <<"Number of Correlated Bolometers: "<< nCorrBolo;

			//double *atmTemplate = cottingham (curData, curAz, curEl, curFlags ,nCorrBolo, nSamples, cleanPixelSize);
			VecDoub atmTemp (nSamples,0.0);
			VecDoub count (nSamples, 0.0);
			double f;
			for (size_t i = 0; i< nCorrBolo; i++)
				for (size_t j = 0; j<nSamples; j++){
					f = curFlags[i*nSamples+j] ? 1.0 : 0.0;
					atmTemp[j] += curData[i*nSamples+j]*f;
					count[j] += f;
				}
			double badSample = 0.0;
			for (size_t j = 0; j<nSamples; j++){
				if (count[j] > 0.0)
					atmTemp[j]/=count[j];
				else{
					badSample ++;
					atmTemp[j]=0.0;
				}

			}
			if (badSample > 0.0)
				cout<<"Fraction of valid points"<<1.0-badSample/(double)nSamples<<endl;
			//if (atmTemplate){
				for (size_t k=0; k<nSamples; k++)
					outputVector[boloIndex[0]*nSamples +k] = atmTemp[k];
//				writeVecOut ("iAtmTemp.txt", atmTemplate, nSamples);
			//	delete [] atmTemplate;
//			}else{
//				for (size_t k=0; k<nSamples; k++)
//					outputVector[boloIndex[0]*nSamples +k] = 0.0;
//			}

			delete [] curData;
			delete [] curFlags;
			gsl_vector_free(curAz);
			gsl_vector_free(curEl);

		}else{
			cleaned[i] = -1;
			for (size_t k=0; k<nSamples; k++)
				outputVector[i*nSamples+k]=0.0;
		}


	}

	for (size_t i=0; i<totalSamples; i++)
		dataVector[i]-=outputVector[i];

	return corrRemoved;
}


bool CleanBspline::removeCorrelations2 (double *dataVector, size_t nDetectors, size_t nSamples, double corrFactor, double *azoffset, double *eloffset, double *kVector, double *flagsVector){

	MatDoub corrMatrix (nDetectors,nDetectors,0.0);
	MatDoub signMatrix (nDetectors, nDetectors, 0.0);
	VecDoub stddevs (nDetectors);
	double dist;
	double mn;
	double minDist = 1.0*60.0;

	//writeVecOut("inputData.txt",dataVector, nSamples*nDetectors);
	//Build correlation matrix
	for (size_t i=0; i<nDetectors;i++){
		corrMatrix[i][i] = 1.0;
		signMatrix[i][i] = 1.0;
		mn = mean(&dataVector[i*nSamples], nSamples);
		stddevs[i] = stddev(&dataVector[i*nSamples], nSamples, mn);
		for (size_t j=i+1; j<nDetectors; j++){
			corrMatrix[i][j] = flagCorrelation(&dataVector[i*nSamples], &flagsVector[i*nSamples], &dataVector[j*nSamples], &flagsVector[j*nSamples],nSamples);
			dist = sqrt(pow(azoffset[i]-azoffset[j],2) + pow(eloffset[i]-eloffset[j],2));
			if (round(corrMatrix[i][j]*100.0)/100.0 < corrFactor || dist < minDist)
				corrMatrix[i][j]=0.0;
			corrMatrix[j][i]=corrMatrix[i][j];// = abs(corrMatrix[i][j]);
			if (corrMatrix[i][j] >=0){
				signMatrix[i][j] = 1.0;
				signMatrix[j][i] = 1.0;
			}else{
				signMatrix[i][j] = signMatrix[j][i] = 1.0;
			}

		}
	}
//	writeMatOut("corrMatrix.txt", corrMatrix);
//	char buff [500];
	size_t nCorrBolo=0;
	bool corrRemoved = false;
	double *cleanedData = new double [nDetectors*nSamples];
	double *cleanedKernel = new double [nDetectors*nSamples];
	VecBool isCleaned (nDetectors, false);
	VecBool dummyFlags(nSamples,true);
	for (size_t i=0; i<nDetectors; i++){
		nCorrBolo = 0;
		for (size_t j=0; j<nDetectors; j++)
			if (corrMatrix[i][j]!= 0)
				nCorrBolo++;

		if (nCorrBolo >1){
			isCleaned[i] = true;
			corrRemoved = true;
			VecInt boloIndex (nCorrBolo,-1);
			VecDoub coeff(nCorrBolo,2);
			boloIndex[0]=i;
			size_t iIndex = 1;
			for (size_t ibolo = 0; ibolo<nDetectors; ibolo++)
				if (ibolo != i && corrMatrix[i][ibolo] !=0)
					boloIndex[iIndex++]=ibolo;

			gsl_matrix *det = gsl_matrix_alloc (nCorrBolo, nSamples);
			gsl_matrix *ket = gsl_matrix_alloc(nCorrBolo, nSamples);
			gsl_matrix *fl = gsl_matrix_alloc (nCorrBolo, nSamples);
			for (size_t ibolo = 0; ibolo< nCorrBolo; ibolo++)
				for (size_t isample=0; isample<nSamples; isample++){
					gsl_matrix_set(det,ibolo, isample, /*signMatrix[boloIndex[0]][boloIndex[ibolo]] * */dataVector[boloIndex[ibolo]*nSamples+isample]);
					gsl_matrix_set(ket,ibolo, isample, /*signMatrix[boloIndex[0]][boloIndex[ibolo]]* */dataVector[boloIndex[ibolo]*nSamples+isample]);
					gsl_matrix_set (fl, ibolo, isample , flagsVector[boloIndex[ibolo]*nSamples+isample]);
				}

			gsl_matrix *denom = gsl_matrix_alloc(nCorrBolo,nCorrBolo);
			gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.,fl,fl,0.,denom);
			gsl_matrix_add_constant(denom,-1.);

			gsl_matrix_mul_elements (det,fl);
			gsl_matrix_mul_elements (ket,fl);
			gsl_matrix *pcaCorr = gsl_matrix_alloc (nCorrBolo, nCorrBolo);
			gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.,det,det,0.,pcaCorr);
			gsl_matrix_div_elements(pcaCorr,denom);
			gsl_matrix_free(fl);
			gsl_matrix_free(denom);
			gsl_eigen_symmv_workspace* w=gsl_eigen_symmv_alloc(nCorrBolo);
			gsl_vector* eVals = gsl_vector_alloc(nCorrBolo);
			gsl_matrix* eVecs = gsl_matrix_alloc(nCorrBolo,nCorrBolo);
			gsl_eigen_symmv(pcaCorr,eVals,eVecs,w);
			gsl_eigen_symmv_sort(eVals,eVecs, GSL_EIGEN_SORT_ABS_DESC);
			gsl_eigen_symmv_free(w);
//			gsl_matrix_free(covMatrix);
			gsl_matrix_free (pcaCorr);



			gsl_matrix* eFunc = gsl_matrix_alloc(nCorrBolo,nSamples);
			gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,eVecs,det,0.,eFunc);
			gsl_matrix* kFunc = gsl_matrix_alloc (nCorrBolo, nSamples);
			gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,eVecs,ket,0.,kFunc);

			size_t n2cut = 2;
			if (n2cut >= nCorrBolo){
				cout<<"Wrong number of eigenvectors to cut..."<<n2cut<<","<<nCorrBolo<<endl;
				n2cut = 1;
			}
			for (size_t ibolo = 0.0; ibolo< n2cut; ibolo++)
				for (size_t isample=0; isample<nSamples; isample++){
					gsl_matrix_set (eFunc, ibolo, isample, 0.0);
					gsl_matrix_set (kFunc, ibolo, isample, 0.0);
				}
			gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,eVecs,eFunc,0.,det);
			gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,eVecs,kFunc,0.,ket);
			gsl_matrix_free(eFunc);
			gsl_matrix_free(kFunc);
			gsl_vector_free(eVals);
			gsl_matrix_free (eVecs);
			for (size_t k=0; k<nSamples; k++){
				cleanedData [i*nSamples+k] = gsl_matrix_get (det,0,k);
				cleanedKernel [i*nSamples+k] = gsl_matrix_get (ket,0,k);
			}
//			sprintf(buff,"outData_%03lu.txt", i);
//			writeVecOut(buff, &cleanedData [i*nSamples], nSamples);
			gsl_matrix_free (det);
			gsl_matrix_free (ket);
		}
//		else{
//			cout<<"Bolometer "<<i<< "has no correlations: "<<endl;
//		}
	}
//	static size_t quit =0;
//	if (quit==2)
//		exit(-1);
//	quit++;
	for (size_t i=0; i< nDetectors; i++)
		for (size_t j=0; j< nSamples; j++)
			if (isCleaned[i]){
				dataVector [i*nSamples+j] = cleanedData[i*nSamples+j];
				kVector [i*nSamples+j] = cleanedKernel[i*nSamples+j];
			}
	delete [] cleanedData;
	delete [] cleanedKernel;
	return corrRemoved;
}

bool CleanBspline::stripePCA(double *dataVector, double *kernelVector, double *flags, size_t nDetectors, size_t nSamples, size_t neigToCut, bool adaptive){

	gsl_matrix *dataMatrix = gsl_matrix_alloc (nDetectors, nSamples);
	gsl_matrix *kernelMatrix = gsl_matrix_alloc (nDetectors, nSamples);
	gsl_matrix *flaggedDataMatrix = gsl_matrix_alloc(nDetectors, nSamples);
	gsl_matrix *corrMatrix = gsl_matrix_alloc(nDetectors, nDetectors);
	double sum = 0;

	for (size_t i = 0; i<nDetectors; i++)
		for (size_t j=0; j<nSamples; j++){
			gsl_matrix_set(dataMatrix, i,j, dataVector[i*nSamples+j]);
			gsl_matrix_set (kernelMatrix, i,j, kernelVector[i*nSamples+j]);
			gsl_matrix_set (flaggedDataMatrix, i,j, dataVector[i*nSamples+j]* flags[i*nSamples+j]);
			sum += flags[i*nSamples+j];
		}
	double delta = double (sum)/double(nDetectors*nSamples);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans,1.0, flaggedDataMatrix, flaggedDataMatrix, 0.0, corrMatrix);
	gsl_matrix_free (flaggedDataMatrix);

	//This is a correction on covariance matrix estimator when missed data (flagged samples)
	//See arXiv 1201.2577v5
	double factor1 = (1.0/delta -1.0/pow(delta,2));
	double factor2 = 1.0/pow(delta,2);
	double entryij;
	for (size_t i = 0; i< nDetectors; i++)
		for (size_t j = 0; j< nDetectors; j++){
			entryij = gsl_matrix_get(corrMatrix, i,j);
			if (i==j)
				entryij = factor1*entryij + factor2*entryij;
			else
				entryij *= factor2;
			gsl_matrix_set (corrMatrix, i,j, entryij);
		}

	gsl_eigen_symmv_workspace* w=gsl_eigen_symmv_alloc(nDetectors);
	gsl_vector* eVals = gsl_vector_alloc(nDetectors);
	gsl_matrix* eVecs = gsl_matrix_alloc(nDetectors,nDetectors);
	gsl_eigen_symmv(corrMatrix,eVals,eVecs,w);
	gsl_eigen_symmv_sort(eVals,eVecs, GSL_EIGEN_SORT_VAL_DESC);
	gsl_eigen_symmv_free(w);
	gsl_matrix_free(corrMatrix);



	gsl_matrix* eFunc = gsl_matrix_alloc(nDetectors,nSamples);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,eVecs,dataMatrix,0.,eFunc);
	gsl_matrix* kFunc = gsl_matrix_alloc (nDetectors, nSamples);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,eVecs,kernelMatrix,0.,kFunc);

	//double damping = 0.0;

	if (adaptive){
		double medianVals, stdvals;
			robustMedian(eVals->data, nDetectors, neigToCut, &medianVals, &stdvals);
			double thresh = medianVals + neigToCut*stdvals;
		size_t nCut = 0;
		for (size_t ibolo = 0.0; ibolo< nDetectors; ibolo++)
			if (gsl_vector_get (eVals,ibolo) > thresh){
				nCut++;
				for (size_t isample=0; isample<nSamples; isample++){
					gsl_matrix_set (eFunc, ibolo, isample, 1.0/sqrt(abs(gsl_vector_get(eVals,ibolo)))*gsl_matrix_get (eFunc, ibolo,isample));
					gsl_matrix_set (kFunc, ibolo, isample, 1.0/sqrt(abs(gsl_vector_get(eVals,ibolo)))*gsl_matrix_get (kFunc, ibolo,isample));
				}
			}
		cout<<"CleanBspline("<<tid<<"): Removing "<< nCut << " eigenvectors"<<endl;
	}else{
		for (size_t ibolo = 0.0; ibolo< neigToCut; ibolo++)
			for (size_t isample=0; isample<nSamples; isample++){
				gsl_matrix_set (eFunc, ibolo, isample, 0.0);
				gsl_matrix_set (kFunc, ibolo, isample, 0.0);
			}
	}
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,eVecs,eFunc,0.,dataMatrix);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,eVecs,kFunc,0.,kernelMatrix);
	gsl_matrix_free(eFunc);
	gsl_matrix_free(kFunc);
	gsl_vector_free(eVals);
	gsl_matrix_free (eVecs);
	for (size_t i=0; i<nDetectors; i++)
		for (size_t k=0; k<nSamples; k++){
			dataVector [i*nSamples+k] = gsl_matrix_get (dataMatrix,i,k);
			kernelVector [i*nSamples+k] = gsl_matrix_get (kernelMatrix,i,k);
		}
	gsl_matrix_free (dataMatrix);
	gsl_matrix_free (kernelMatrix);

	return true;

}



MatDoub CleanBspline::linearGainCorrection(double *dataVector, size_t nDetectors, size_t nSamples, double *fakeAtm){

	//	double *average = new double [nSamples];
	//	VecDoub ithSample (nDetectors);

	//double relGain[nDetectors];
	//double relOffset[nDetectors];

	MatDoub relCoeff(2,nDetectors);    //[0] is relGain, [1] relOffset
	double cv0,cv1,cv2, chisq;


	//	for (size_t i=0; i<nSamples; i++){
	//		for (size_t j=0; j<nDetectors; j++)
	//			ithSample[j]=dataVector[j*nSamples + i];
	//		average[i] = gsl_stats_mean(ithSample.getData(),1,nDetectors);
	//	}

	for (size_t i=0; i<nDetectors; i++){
		gsl_fit_linear(&dataVector[0],1,&(dataVector[i*nSamples]),1,nSamples,&(relCoeff[1][i]),&(relCoeff[0][i]),&cv0,&cv1,&cv2,&chisq);
	}

	for (size_t i=0; i<nDetectors; i++)
		for (size_t j=0; j<nSamples; j++){
			dataVector[i*nSamples+j]= (dataVector[i*nSamples+j]-relCoeff[1][i])/relCoeff[0][i];
			if (fakeAtm)
				fakeAtm[i*nSamples+j]= (fakeAtm[i*nSamples+j]-relCoeff[1][i])/relCoeff[0][i];
		}
	//delete [] average;

	return relCoeff;
}







bool CleanBspline::fixFlags(double *dataVector, bool *flagsVector, int nDetectors, long nSamples){
	//  long prev=0;
	//  long next=0;
	//  double m=0;
	//  double data1=0.0;
	//  double data2=0.0;
	//  double value = 0.0;
	//  long fpos=0.0;
	double sum = 0.0;
	VecDoub sumPerDetector(nDetectors,0.0);
	//Check flag vector, if more than 50% of the data is flagged as invalid then dismiss the whole scan
	for (int i=0; i<nDetectors; i++)
		for (long j=0; j<nSamples; j++){
			sum += (double)flagsVector [i*nSamples + j];
			sumPerDetector[i]+=(double)flagsVector [i*nSamples + j];
		}
	if (sum <= 0.5*(double)nDetectors*(double)nSamples)
		return false;

//	if (sum == (double)nDetectors*(double)nSamples)
//		return true;
////	Now iterate in vector data to remove bad flags
//	VecDoub time (nSamples);
//	double iInter;
////	writeVecOut ("dVector.txt", dataVector, nSamples*nDetectors);
//	gsl_error_handler_t *hn =  gsl_set_error_handler_off ();
//	for (long j=0; j<nSamples;j++)
//		time[j]=j;
//	for (int i = 0; i<nDetectors; i++){
//		if (sumPerDetector[i] == nSamples)
//			continue;
//		VecDoub interpData (sumPerDetector[i]);
//		VecDoub interpTime (sumPerDetector[i]);
//		long ix = 0;
//		for (long j=0; j<nSamples;j++){
//			if (flagsVector[i*nSamples+j]){
//				interpData[ix]=dataVector[i*nSamples+j];
//				interpTime[ix++]= time[j];
//			}
//		}
//	    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
//	    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, sumPerDetector[i]);
//	    gsl_spline_init (spline, &interpTime[0], &interpData[0], sumPerDetector[i]);
//
//	    for (long j=0; j<nSamples;j++){
//	    			if (!flagsVector[i*nSamples+j]){
//	    				iInter  = gsl_spline_eval (spline, time[j], acc);
//	    				if (!gsl_finite(iInter))
//	    					return false;
//	    				dataVector[i*nSamples+j]=iInter;
//	    			}
//	    }
//	    gsl_spline_free (spline);
//	    gsl_interp_accel_free (acc);
//
//	}
//	gsl_set_error_handler(hn);
//	    writeVecOut ("fVector.txt", dataVector, nSamples*nDetectors);
//	    exit(-1);
	return true;
}

CleanBspline::~CleanBspline(){
	this->destroyBaseMatrix();
	if (scanStatus)
		delete [] scanStatus;
}

cs* CleanBspline::getPMatrix(gsl_vector *ra, gsl_vector *dec, double pixelSize){
	double minRa;
	double minDec;
	double maxRa;
	double maxDec;
	int sizeXX = 0;
	int sizeYY = 0;
	long dataLen = ra->size;
	long tPixels = 0;

	cs *matrix = NULL;
	cs *pMatrix = NULL;
	int *usedPixels = NULL;
	long nUsedPixels = 0;
	long *newPixelPosition = NULL;
	long *rowPos = NULL;
	long *colPos = NULL;

	gsl_vector_minmax(ra, &minRa, &maxRa);
	gsl_vector_minmax(dec, &minDec, &maxDec);

	sizeXX = ceil((maxRa-minRa)/pixelSize);
	sizeYY = ceil((maxDec-minDec)/pixelSize);



	tPixels = sizeXX*sizeYY;

	if (tPixels ==0){
		cerr<<"CleanBspline(): Fatal Error, no pixels in pointing matrix: "<<this->ap->getMapFile()<<endl;
		exit(-1);
	}
	//cerr <<"Creating Ponting Matrix of Nsamples: "<< dataLen<< " x Npixels"<<tPixels <<endl;
	//cerr <<"Ra (max,min) : ("<<maxRa<<","<<minRa<<") Dec (max,min) : ("<<maxDec<<","<<minDec<<")" <<endl;
	//cerr <<"Pixel Size: "<< pixelSize<<endl;
	usedPixels = new int[tPixels];
	rowPos = new long[dataLen];
	colPos = new long [dataLen];
	newPixelPosition = new long[tPixels];
	double xra = 0;
	double ydec = 0;
	long pos = 0;
	long i=0;
	long j=0;

	for (i=0; i< tPixels;i++)
		usedPixels[i] = 0;
	for (i=0; i< dataLen; i++){
		//Calculate pixel in xra and ydec of this samples
		xra =  round ((gsl_vector_get(ra, i) -minRa)/pixelSize);
		if (xra >= sizeXX)
			xra = sizeXX -1;
		ydec = round ((gsl_vector_get(dec, i) -minDec)/pixelSize);
		if (ydec >= sizeYY)
			ydec = sizeYY-1	;
		pos = xra * sizeYY + ydec;
		//pos = ydec * sizeXX  + xra;
		rowPos[i] = i;
		colPos[i] = pos;
		usedPixels[pos]++;
	}
	//trim pMatrix to the used pixels only
	//also redefine positions to used pixels only
	j=0;
	for (i=0; i<tPixels; i++){
		if (usedPixels[i] > 0){
			nUsedPixels++;
			newPixelPosition[i] =j++;
		}
		else{
			newPixelPosition[i] = -1;
		}
	}

	matrix = cs_spalloc(dataLen, nUsedPixels, dataLen,1,1);
	for (i=0; i<dataLen; i++){
		long npos = newPixelPosition[colPos[i]];
		if (npos >= 0){
			//cs_entry(matrix, i, npos, 1.0/usedPixels[colPos[i]]);
			cs_entry(matrix, i, npos, 1.0);
		}
		else{
			cerr<<"Wrong matrix position. Imploding"<<endl;
			exit(-1);
		}
	}

	pMatrix = cs_compress(matrix);
	cs_spfree(matrix);
	delete [] newPixelPosition;
	delete [] rowPos;
	delete [] colPos;
	delete [] usedPixels;
	return pMatrix;
}

void CleanBspline::subtractTemplate(double* detector, size_t nSamples,	double* aTemplate, size_t oSamples, double increment, char * outName, bool overwrite) {

	if (!aTemplate)
		return;

	if (nSamples==oSamples){
		for (size_t i=0; i<nSamples; i++){
			if (overwrite)
				detector[i]=aTemplate[i];
			else
				detector[i]-=aTemplate[i];
		}
		if (outName){
			writeVecOut(outName, aTemplate, nSamples);

		}
	}else{
		gsl_interp_accel *acc = gsl_interp_accel_alloc ();
		gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, oSamples);

		double *xTemp = new double [oSamples];
		double *iValues = new double [nSamples];

		for (size_t i=0; i<oSamples; i++)
			xTemp[i] = round(i*increment);
		gsl_spline_init (spline, xTemp, aTemplate, oSamples);

		for (size_t i=0; i<nSamples; i++){
			iValues[i] = gsl_spline_eval (spline, i, acc);
			if (overwrite)
				detector[i] = iValues[i];
			else
				detector[i]-= iValues[i];
		}

		if(outName)
			writeVecOut(outName,detector, nSamples);

		delete[] xTemp;
		delete [] iValues;
		gsl_spline_free (spline);
		gsl_interp_accel_free (acc);

	}

}


void CleanBspline::cleanStripeFFT2(double *data, double *kernel, int nDetectors, long nSamples, double cutOff, double cutLow, double sampleRate){
	double  *tmp = new double [nSamples];
	double  *ktmp = new double [nSamples];
	double  *freq = new double [nSamples];

	gsl_fft_real_wavetable * real = gsl_fft_real_wavetable_alloc (nSamples);
	gsl_fft_real_workspace * work = gsl_fft_real_workspace_alloc (nSamples);
	gsl_fft_halfcomplex_wavetable * hc = gsl_fft_halfcomplex_wavetable_alloc (nSamples);

	double maxf = 0.5*sampleRate;

	//Generate frequencies
	freq[0]= 0.0;
	for (long j=1; j<nSamples-1; (j+=2) ){
		freq[j] = freq[j+1] =  double(j)/(2.0*double(nSamples))*sampleRate;
	}

	if (nSamples % 2 == 0)
		freq[nSamples-1] = maxf;
	else
		freq[nSamples-1]= freq[nSamples-2] = maxf;

//	double factor [nSamples];
//	double spectrum [nSamples];


	for (int i = 0; i<nDetectors; i++){
		for (long j=0; j<nSamples; j++){
			tmp[j] = data [i*nSamples +j];
			ktmp[j] = kernel [i*nSamples +j];
		}
		gsl_fft_real_transform(tmp,1, nSamples,real,work);
		gsl_fft_real_transform(ktmp,1, nSamples,real,work);

		tmp[0]=0.0;
		ktmp[0]=0.0;
		for (long j=1; j<nSamples; j++){
			tmp[j] *= sqrt (1.0/ (1.0+ pow(cutOff/freq[j],10.0)));
			ktmp[j] *= sqrt (1.0/ (1.0+ pow(cutOff/freq[j],10.0)));
		}

		gsl_fft_halfcomplex_inverse (tmp, 1, nSamples, hc, work);
		gsl_fft_halfcomplex_inverse (ktmp, 1, nSamples, hc, work);
		for (long j=0; j<nSamples; j++){
			data[i*nSamples +j] =tmp[j];
			kernel[i*nSamples+j] =ktmp[j];
		}
	}

	gsl_fft_real_wavetable_free(real);
	gsl_fft_halfcomplex_wavetable_free (hc);
	gsl_fft_real_workspace_free(work);
	delete [] tmp;
	delete [] ktmp;
	delete [] freq;
}


void CleanBspline::cleanStripeFFT3(double *data, double *kernel, size_t nDetectors, size_t nSamples, double cutOff, double sampleRate){

  size_t newSamples = size_t(floor(pow(2,ceil(log(nSamples)/log(2)))));
  size_t pad = (newSamples-nSamples)/ 2;

  VecDoub freq (newSamples,0.0);
  VecDoub tmp (newSamples,0.0);
  VecDoub ktmp (newSamples,0.0);

  double maxf = 1.0/(2.0*newSamples)*sampleRate;
  //cout<<"CleanSplinne::cleanStripeFFT().Padding array to "<<newSamples<<endl;

  //Generate frequencies
  freq[0]= 0.0;
  for (size_t j=1; j<newSamples/2; j++){
	  freq[j] = freq[newSamples-j] =  j/(2.0*newSamples)*sampleRate;
  }
  freq [newSamples/2]= maxf;
  double m,b;
  double mk,bk;
  //Copy data into padded array
  for (size_t i = 0; i<nDetectors; i++){
	  for (size_t j=0; j<nSamples; j++){
		  tmp[j+pad] = data [i*nSamples +j];
		  ktmp[j+pad] = kernel [i*nSamples +j];
	  }
	  m = data [i*nSamples]/pad;
	  mk = kernel [i*nSamples]/pad;
	  for (size_t j=0; j<pad; j++){
		  tmp[j] = m*j;
		  ktmp[j] = mk*j;
	  }
      m = -data[(i+1)*nSamples-1]/(newSamples-nSamples-pad);
      mk = -kernel[(i+1)*nSamples-1]/(newSamples-nSamples-pad);
      b = data[(i+1)*nSamples-1];
      bk = kernel[(i+1)*nSamples-1];
	  for (size_t j=nSamples+pad; j<newSamples; j++){
		  tmp[j] = m*(j-nSamples-pad) + b ;
		  ktmp[j] = mk*(j-nSamples-pad) + bk;
	  }
	gsl_fft_real_radix2_transform(&tmp[0],1, newSamples);
	gsl_fft_real_radix2_transform(&ktmp[0],1, newSamples);


	tmp[0]=0.0;
	ktmp[0]=0.0;
	for (size_t j=1; j<newSamples; j++){
		tmp[j] *= sqrt (1.0/ (1.0+ pow(cutOff/freq[j],10.0)));
		ktmp[j] *= sqrt (1.0/ (1.0+ pow(cutOff/freq[j],10.0)));
	}
	gsl_fft_halfcomplex_radix2_inverse	 (&tmp[0], 1, newSamples);
	gsl_fft_halfcomplex_radix2_inverse (&ktmp[0], 1, newSamples);
    for (size_t j=0; j<nSamples; j++){
      data[i*nSamples +j] =tmp[j+pad];
      kernel[i*nSamples+j] =ktmp[j+pad];
    }
  }

}

void CleanBspline::splineResidual (double *data, double *flags, size_t nDetectors, size_t nSamples, gsl_vector  *azVector, gsl_vector  *elVector, double correlate){
	MatDoub corrMatrix (nDetectors,nDetectors,0.0);
	double dist;
	double mindist = 0.0/3600.0;

	if (dataArray->getAp()->getObservatory() == "LMT")
		mindist /=3.0;

	for (size_t i=0; i<nDetectors;i++){
		corrMatrix[i][i] = 1.0;

		//mn = mean(&data[i*nSamples], &flags[i*nSamples], nSamples);

		for (size_t j=i+1; j<nDetectors; j++){
			corrMatrix[i][j] = flagCorrelation(&data[i*nSamples], &flags[i*nSamples], &data[j*nSamples], &flags[j*nSamples],nSamples);
			dist = sqrt(pow(gsl_vector_get(azVector,i*nSamples)-gsl_vector_get(azVector,j*nSamples),2) + pow(gsl_vector_get(elVector,i*nSamples)-gsl_vector_get(elVector,j*nSamples),2));
			if (round(corrMatrix[i][j]*100.0)/100.0 < correlate || dist < mindist)
				corrMatrix[i][j]=0.0;
			corrMatrix[j][i]=corrMatrix[i][j];// = abs(corrMatrix[i][j]);
		}
	}
	MatDoub outData (nDetectors, nSamples, 0.0);
	VecBool isCleaned (nDetectors, false);

	size_t oSamples=nSamples;

	for (size_t ibolo = 0; ibolo < nDetectors; ibolo++){
		//Get the number of correlated bolometers
		size_t nBoloCorr = 0;
		for (size_t jbolo = 0; jbolo < nDetectors; jbolo++)
			if (corrMatrix [ibolo][jbolo]> 0.0)
				nBoloCorr++;
		//cerr<<"Number of correlated bolometers"<< ibolo << "->"<< nBoloCorr<<endl;
		if (nBoloCorr >= 2){

			double increment = 1.0;
			size_t index;


			if (resample > 1){
				long nnSamples = long(round(oSamples/resample));
				increment = double(oSamples-1)/double(nnSamples-1);
				nSamples =nnSamples;
			}

			size_t dataLen = nBoloCorr*nSamples;

			size_t boloIndex [nBoloCorr];
			boloIndex[0] = ibolo;
			size_t curBolo = 1;
			for (size_t jbolo = 0; jbolo < nDetectors; jbolo++)
				if (corrMatrix[ibolo][jbolo] > 0.0 && ibolo != jbolo){
					boloIndex[curBolo++]= jbolo;
					//cerr <<jbolo<<endl;
				}
			double *cBoloData = new double [dataLen];
			bool *cFlags = new bool [dataLen];
			gsl_vector *cAzVector = gsl_vector_alloc(dataLen);
			gsl_vector *cElVector = gsl_vector_alloc(dataLen);
			//copy data
			for (size_t jbolo = 0; jbolo < nBoloCorr; jbolo++){
				for (size_t jSample = 0; jSample < nSamples; jSample++){
					index = long(round(jSample*increment));
					cBoloData[jbolo*nSamples+jSample] = data [boloIndex[jbolo]*oSamples+index];
					cFlags [jbolo*nSamples+jSample] = flags [boloIndex[jbolo]*oSamples+index];
					gsl_vector_set (cAzVector,jbolo*nSamples+jSample,gsl_vector_get(azVector,boloIndex[jbolo]*oSamples+index));
					gsl_vector_set (cElVector, jbolo*nSamples+jSample, gsl_vector_get ( elVector ,boloIndex[jbolo]*oSamples+index));
				}
			}
			//writeVecOut ("osData", data, nDetectors*nSamples);
			//writeVecOut("sData.txt", cBoloData, nBoloCorr*nSamples);
			double *tmplt = cottingham(cBoloData, cAzVector, cElVector, cFlags, nBoloCorr, nSamples, this->ap->cleanPixelSize, 0);


			//writeVecOut("sTemplate.txt", tmplt, nBoloCorr*nSamples);
			//exit(-1);

			subtractTemplate (&outData[ibolo][0], oSamples, &tmplt[0], nSamples, increment, NULL, true);
			//for (size_t iSample = 0; iSample<nSamples; iSample++)
				//outData[ibolo][iSample] = tmplt [iSample];
			delete [] tmplt;
			delete [] cFlags;
			delete [] cBoloData;
			gsl_vector_free (cAzVector);
			gsl_vector_free (cElVector);
			isCleaned[ibolo]=true;

		}
	}

	nSamples = oSamples;
	for (size_t ibolo= 0; ibolo< nDetectors; ibolo++)
		if (isCleaned[ibolo])
			for (size_t jSample = 0; jSample<nSamples; jSample++)
				data[ibolo*nSamples+jSample] -= outData[ibolo][jSample];
}


void CleanBspline::azelResidualLMT (double *data, double *flags, size_t nDetectors, size_t nSamples, gsl_vector *azOffsets, gsl_vector *elOffsets, double *sens){
	gsl_matrix * cBoloData = gsl_matrix_alloc(nDetectors, nSamples);
	gsl_matrix * cBoloFlag = gsl_matrix_alloc (nDetectors, nSamples);
	gsl_matrix *az = gsl_matrix_alloc(nDetectors, nSamples);
	gsl_matrix *el = gsl_matrix_alloc(nDetectors, nSamples);
	gsl_vector *corrCoefs = gsl_vector_alloc (nDetectors);

	for (size_t ibolo = 0; ibolo < nDetectors; ibolo++){
			for (size_t jSample = 0; jSample < nSamples; jSample++){
				gsl_matrix_set (cBoloData,ibolo,jSample,data[ibolo*nSamples +jSample]);
				gsl_matrix_set (cBoloFlag, ibolo,jSample, flags [ibolo*nSamples+ jSample]);
				gsl_matrix_set (az, ibolo, jSample, gsl_vector_get(azOffsets,ibolo*nSamples+jSample));
				gsl_matrix_set (el, ibolo, jSample, gsl_vector_get(elOffsets,ibolo*nSamples+jSample));
			}

		gsl_vector_set (corrCoefs, ibolo, sens[ibolo]);

	}
	AzElTemplateCalculator aeTemp  (cBoloData, cBoloFlag,az,el);
//	aeTemp.overrideMode(LINEAR);
	aeTemp.calculateTemplate(corrCoefs);
	gsl_matrix *atmTemplate = aeTemp.getTemplate();
	for (size_t ibolo = 0; ibolo < nDetectors; ibolo++){
		for (size_t iSample = 0; iSample<nSamples; iSample++)
			data[ibolo*nSamples+iSample] -= gsl_matrix_get (atmTemplate, ibolo, iSample);
	}

	gsl_matrix_free(cBoloData);
	gsl_matrix_free(cBoloFlag);
	gsl_vector_free(corrCoefs);
	gsl_matrix_free(az);
	gsl_matrix_free(el);
}


void CleanBspline::azelResidual (double *data, double *flags, size_t nDetectors, size_t nSamples, double *azOffsets, double *elOffsets, double *sens){
	gsl_matrix * cBoloData = gsl_matrix_alloc(nDetectors, nSamples);
	gsl_matrix * cBoloFlag = gsl_matrix_alloc (nDetectors, nSamples);
	gsl_vector *az = gsl_vector_alloc(nDetectors);
	gsl_vector *el = gsl_vector_alloc(nDetectors);
	gsl_vector *corrCoefs = gsl_vector_alloc (nDetectors);

	for (size_t ibolo = 0; ibolo < nDetectors; ibolo++){
			for (size_t jSample = 0; jSample < nSamples; jSample++){
				gsl_matrix_set (cBoloData,ibolo,jSample,data[ibolo*nSamples +jSample]);
				gsl_matrix_set (cBoloFlag, ibolo,jSample, flags [ibolo*nSamples+ jSample]);
			}
		gsl_vector_set (az, ibolo, azOffsets[ibolo]);
		gsl_vector_set (el, ibolo, elOffsets[ibolo]);
		gsl_vector_set (corrCoefs, ibolo, sens[ibolo]);

	}
	AzElTemplateCalculator aeTemp  (cBoloData, cBoloFlag,az,el);
	//aeTemp.overrideMode(LINEAR);
	aeTemp.calculateTemplate(corrCoefs);
	gsl_matrix *atmTemplate = aeTemp.getTemplate();
	for (size_t ibolo = 0; ibolo < nDetectors; ibolo++){
		for (size_t iSample = 0; iSample<nSamples; iSample++)
			data[ibolo*nSamples+iSample] -= gsl_matrix_get (atmTemplate, ibolo, iSample);
	}
	//gsl_matrix_free(atmTemplate);
	gsl_matrix_free(cBoloData);
	gsl_matrix_free(cBoloFlag);
	gsl_vector_free(corrCoefs);
	gsl_vector_free(az);
	gsl_vector_free(el);
}


void CleanBspline::azelResidual (double *data, double *flags, size_t nDetectors, size_t nSamples, double *azOffsets, double *elOffsets, double correlate){
	MatDoub corrMatrix (nDetectors,nDetectors,0.0);
	double dist;
	double mindist = 0.0;

	for (size_t i=0; i<nDetectors;i++){
		corrMatrix[i][i] = 1.0;
		for (size_t j=i+1; j<nDetectors; j++){
			corrMatrix[i][j] = flagCorrelation(&data[i*nSamples], &flags[i*nSamples], &data[j*nSamples], &flags[j*nSamples],nSamples);
			dist = sqrt(pow(azOffsets[i]-azOffsets[j],2) + pow(elOffsets[i]-elOffsets[j],2));
			if (round(corrMatrix[i][j]*100.0)/100.0 < correlate || dist < mindist)
				corrMatrix[i][j]=0.0;
			corrMatrix[j][i]=corrMatrix[i][j];// = abs(corrMatrix[i][j]);
		}
	}
	MatDoub outData (nDetectors, nSamples, 0.0);
	VecBool isCleaned (nDetectors, false);

	for (size_t ibolo = 0; ibolo < nDetectors; ibolo++){
		//Get the number of correlated bolometers
		size_t nBoloCorr = 0;
		for (size_t jbolo = 0; jbolo < nDetectors; jbolo++)
			if (corrMatrix [ibolo][jbolo]> 0.0)
				nBoloCorr++;
		cout<<"Number of correlated bolometers"<< nBoloCorr<<endl;
		if (nBoloCorr > 2){
			size_t boloIndex [nDetectors];
			boloIndex[0] = ibolo;
			size_t curBolo = 1;
			for (size_t jbolo = 0; jbolo < nDetectors; jbolo++)
				if (corrMatrix[ibolo][jbolo] > 0.0)
					boloIndex[curBolo++]= jbolo;
			gsl_matrix * cBoloData = gsl_matrix_alloc(nBoloCorr, nSamples);
			gsl_matrix * cBoloFlag = gsl_matrix_alloc (nBoloCorr, nSamples);
			//copy data
			for (size_t jbolo = 0; jbolo < nBoloCorr; jbolo++){
				for (size_t jSample = 0; jSample < nSamples; jSample++){
					gsl_matrix_set (cBoloData,jbolo,jSample,data[boloIndex[jbolo]*nSamples +jSample]);
					gsl_matrix_set (cBoloFlag, jbolo,jSample, flags [boloIndex[jbolo]*nSamples+ jSample]);
				}
			}
			//Now copy coordinates
			gsl_vector *az = gsl_vector_alloc(nBoloCorr);
			gsl_vector *el = gsl_vector_alloc(nBoloCorr);
			gsl_vector *corrCoefs = gsl_vector_alloc (nBoloCorr);
			for (size_t jbolo = 0; jbolo < nBoloCorr; jbolo++){
				gsl_vector_set (az, jbolo, azOffsets[boloIndex[jbolo]]);
				gsl_vector_set (el, jbolo, elOffsets[boloIndex[jbolo]]);
				gsl_vector_set (corrCoefs, jbolo, corrMatrix[ibolo][boloIndex[jbolo]]);
			}

			AzElTemplateCalculator aeTemp  (cBoloData, cBoloFlag,az,el);
			aeTemp.calculateTemplate(corrCoefs);

//			//gsl_vector *btmp = gsl_vector_alloc(nSamples);
//			writeGslMatrix("dataTemp.txt",cBoloData);
//			//writeGslMatrix("azelTemplate.txt", tmplt);
			aeTemp.removeTemplate();
//			writeGslMatrix("azelResidual.txt", cBoloData);
//			gsl_matrix *tmplt = aeTemp.getTemplate();
			for (size_t iSample = 0; iSample<nSamples; iSample++)
				outData[ibolo][iSample] = gsl_matrix_get (cBoloData, 0, iSample);

			gsl_matrix_free(cBoloData);
			gsl_matrix_free(cBoloFlag);
			gsl_vector_free(az);
			gsl_vector_free(el);
			gsl_vector_free(corrCoefs);
			isCleaned[ibolo] = true;
//			exit(1);
		}
	}

	for (size_t ibolo= 0; ibolo< nDetectors; ibolo++)
		if (isCleaned[ibolo])
			for (size_t jSample = 0; jSample<nSamples; jSample++)
				data[ibolo*nSamples+jSample] = outData[ibolo][jSample];

}

//Getters

cs* CleanBspline::getBaseMatrix(){
	return baseMatrix;
}
