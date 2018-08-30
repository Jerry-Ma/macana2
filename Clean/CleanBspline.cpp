#include <omp.h>
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
#include <fftw3.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_spline.h>
#include <assert.h>
#include <libgen.h>

#include "CleanBspline.h"
#include "vector_utilities.h"
#include "sparseUtilities.h"
#include "AzElTemplateCalculator.h"


using namespace std;
#define USEMPFIT false
#define LINFIT true

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
  bright = ap->mask /3600.0;
  resample = ap->getResample();
  maxHex = 0;
  fixedData = NULL;
  totSamples = 0;
  scanStatus = NULL;
  currScan = 0;
  cleanKernel = false;

#if defined(_OPENMP)
  tid = omp_get_thread_num() + 1;
#else
  tid = 0;
#endif
  if (resample <1){
	  cout<< "CleanBspline ("<< tid<<"). resample value is less than 1. Set to use all samples instead"<<endl;
	  resample = 1;
  }

  if (bright>0.0){
	  cout<< "CleanBspline ("<< tid<<"). Masking out "<< ap->mask << " arcsec from the center of the map" <<endl;
  }


}


void CleanBspline::removeBadBolos(){

	size_t nDetectors = dataArray->getNDetectors();
	size_t totSamples =  dataArray->detectors[0].hValues.size();

	int *di = dataArray->getDetectorIndices();
	size_t refBolo = dataArray->getRefBoloIndex();

	VecDoub avg (totSamples,0.0);
	VecDoub detMean (nDetectors);
	VecInt count(totSamples, 0);
	MatDoub linCoeffs (nDetectors,2);
	for (size_t iBolo = 0; iBolo<nDetectors; iBolo++){

		linfit_flags(&dataArray->detectors[di[refBolo]].hValues[0], &dataArray->detectors[di[refBolo]].hSampleFlags[0], \
				&dataArray->detectors[di[iBolo]].hValues[0], &dataArray->detectors[di[iBolo]].hSampleFlags[0] , \
				totSamples, &linCoeffs[iBolo][0], &linCoeffs[iBolo][1], true);
		for (size_t iSample=0; iSample<totSamples; iSample++){
			if (dataArray->detectors[di[iBolo]].hSampleFlags[iSample]){
				avg[iSample]+= (dataArray->detectors[di[iBolo]].hValues[iSample]-linCoeffs[iBolo][0])/linCoeffs[iBolo][1];
				count[iSample]++;
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
				detPow[i] += pow(dataArray->detectors[di[i]].hValues[j]-(linCoeffs[i][1]*avg[j]+linCoeffs[i][0]), 2.0);
		}
	}

	double mPow = median(&detPow[0],nDetectors);
	size_t badBoloCount = 0;
	for (size_t i=0; i<nDetectors; i++)
		if (detPow[i]> 3.5*mPow){
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
		if (ap->getOrder() != 100)
			if (ap->getAzelMap() != 1)
					removeBadBolos();
		dataArray->fakeAtmData(this->cleanKernel);
		

	

	totSamples =  dataArray->detectors[0].hValues.size();
	atmTemplate = VecDoub (totSamples*dataArray->getNDetectors(),0.0);

	if (ap->getOrder() != 100)
		this->cleanScans();
	else{
		cerr<<"CleanBspline::Atmosphere subtraction is not done as requested by user (splineOrder=100)."<<endl;
		return true;
	}

	int *di = dataArray->getDetectorIndices();
	for (size_t i=0; i<(size_t)dataArray->getNDetectors(); i++)
		for (size_t j=0; j<totSamples; j++){
			dataArray->detectors[di[i]].atmTemplate[j]-= atmTemplate[i*totSamples + j];
			dataArray->detectors[di[i]].hValues[j] -= atmTemplate[i*totSamples + j];
			if (this->cleanKernel){
					dataArray->detectors[di[i]].hKernel[j] -= atmTemplate[i*totSamples + j];
			}
			

		}

	if (cleanStripe != 0){
		this->removeLargeScaleResiduals(translateStripeMethod(),cleanStripe);
	}
	atmTemplate.resize(0);
	this->destroyBaseMatrix();
	if (ap->getOrder() != 100)
		if (ap->getAzelMap() != 1)
				removeBadBolos();
	return true;
}



bool CleanBspline::cleanScans (){

	Array *dataArray = this->dataArray;
	size_t nScans = telescope->scanIndex.ncols();
	size_t nDetectors=dataArray->getNDetectors();

	int *di = dataArray->getDetectorIndices();
	int refBolo=-1;

	bool debug = false;

	if (!scanStatus){
		scanStatus = new bool[nScans];
		for (size_t k=0; k< nScans; k++)
			scanStatus[k] = true;
	}

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

		double increment = 1.0;
		size_t index;
		size_t oSamples=nSamples;

		if (resample > 1){
			long nnSamples = long(round(nSamples/resample));
			increment = double(nSamples-1)/double(nnSamples-1);
			nSamples =nnSamples;

			if (size_t (round(increment*(nSamples-1))) != oSamples-1){
				cerr<<"Cannot set the resample value to "<< resample << "please change it on your apFile"<<endl;
				exit(-1);
			}
		}

		size_t dataLen = nDetectors*nSamples;

		double *dataVector = new double [dataLen];
		gsl_vector *raVector = gsl_vector_alloc(dataLen);
		gsl_vector *decVector = gsl_vector_alloc(dataLen);
		bool *flags = new bool [dataLen];
		double *atmTemplate =NULL;
		double dist = 0.0;


		for(size_t i=0; i<nDetectors; i++){
			for(size_t j=0; j<nSamples; j++){
				index = long(round(j*increment));
				dataVector[i*nSamples +j]=dataArray->detectors[di[i]].hValues[si+index];
				gsl_vector_set(raVector,i*nSamples+j, dataArray->detectors[di[i]].hAzPhys[si+index]);
				gsl_vector_set(decVector,i*nSamples+j, dataArray->detectors[di[i]].hElPhys[si+index]);
				flags[i*nSamples+j] = true;// (bool) dataArray->detectors[di[i]].hSampleFlags[si+index];
				dist = sqrt (pow(dataArray->detectors[di[i]].hRa[si+index],2.0)+pow(dataArray->detectors[di[i]].hDec[si+index],2.0));
				dist*=180.0/M_PI;
				if (dist<this->bright)
					flags[i*nSamples+j] = false;
			}
		}

		if (debug){
			char buffer [1000];
			strcpy(buffer, this->ap->getDataFile());
			string filename = basename(buffer);
			filename = filename.substr(0,filename.find(".nc"));
			sprintf(buffer, "observed_%s_ps%02d_cc%03d.txt", filename.c_str(), int(this->ap->getCleanPixelSize()), int (1.0/this->ap->getControlChunk()));
			writeVecOut (buffer, dataVector, nSamples*nDetectors);
		}

		if (fixFlags(dataVector,flags,(int)nDetectors,(int)nSamples)){
			atmTemplate = cottingham(dataVector, raVector, decVector, flags, nDetectors, nSamples, cleanPixelSize, refBolo);

			if (!atmTemplate){
				cout<<"CleanBspline("<<tid<<")::cleanScans(). Failed to produce and atmosphere template. on scan "<<k <<"for file"<<ap->getDataFile()<<endl;
				for(size_t i=0; i<nDetectors; i++){
					for(size_t j=0; j<oSamples; j++){
						dataArray->detectors[di[i]].hValues[si+j]=0.0;
						dataArray->detectors[di[i]].hSampleFlags[si+j]=0.0;
					}
				}
				scanStatus[k] = false;
			}else{

				if (debug){
					char buffer[1000];
					strcpy(buffer, this->ap->getDataFile());
					string filename = basename(buffer);
					filename = filename.substr(0,filename.find(".nc"));
					sprintf(buffer, "residual_%s_ps%02d_cc%03d.txt", filename.c_str(), int(this->ap->getCleanPixelSize()), int (1.0/this->ap->getControlChunk()));
					writeVecOut (buffer, dataVector, nSamples*nDetectors);
					//writeVecOut ("residualData.txt", dataVector, dataLen);
					sprintf(buffer, "atmtemplate_%s_ps%02d_cc%03d.txt", filename.c_str(), int(this->ap->getCleanPixelSize()), int (1.0/this->ap->getControlChunk()));
					writeVecOut (buffer, atmTemplate, dataLen);
					exit(-1);
				}
				for(size_t i=0; i<nDetectors; i++){
					//subtractTemplate(&dataArray->detectors[di[i]].hValues[si],oSamples,&atmTemplate[i*nSamples],nSamples,increment, NULL);
					subtractTemplate (&this->atmTemplate[i*totSamples + si], oSamples, &atmTemplate[i*nSamples], nSamples, increment, NULL, true);
				}

				if (debug){
					writeVecOut ("templateInterp.txt", &this->atmTemplate[si], oSamples);
					writeVecOut ("fullOdata.txt", &dataArray->detectors[di[0]].hValues[si], oSamples);
					exit(-1);
				}

				delete [] atmTemplate;
			}
		}
		else{
			cout<<"CleanBspline("<<tid<<")::cleanScans(). Not enough valid samples  to produce and atmosphere template. on scan "<<k <<"for file"<<ap->getDataFile()<<endl;
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
	if (ap->stripeMethod == "azel" ){
		cout <<"CleanBspline ("<<tid<<")Using Az/El polynomial template of order:"<<ap->cleanStripe<<endl;
		return AZEL_TEMPLATE;
	}
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
	VecDoub AzOff(nDetectors);
	VecDoub ElOff(nDetectors);
	int * di = dataArray->getDetectorIndices();
	//VecDoub energy (nDetectors,0.0);
	//VecDoub meanTod (nSamples,0.0);

	for (size_t id=0; id<nDetectors; id++){
		AzOff[id]= dataArray->detectors[di[id]].azOffset;
		ElOff[id]= dataArray->detectors[di[id]].elOffset;

		VecBool testFlags (nSamples, true);

#ifdef LINFIT
			if (linfit_flags(&dataVector[rBoloIndex*nSamples], &flags[rBoloIndex*nSamples], &dataVector[id*nSamples], &flags[id*nSamples], nSamples,&relCoeff[id][0], &relCoeff[id][1],USEMPFIT)<0){
			//if (linfit_flags(&dataVector[rBoloIndex*nSamples], &testFlags, &dataVector[id*nSamples],  &testFlags, nSamples,&relCoeff[id][0], &relCoeff[id][1],USEMPFIT)<0){
				cerr<<"Data on file "<<dataArray->getAp()->getDataFile()<< "has no valid data on this scan"<<endl;
				return NULL;
			}
			for (size_t is = 0; is <nSamples; is++)
				dataVector[id*nSamples+is] = (dataVector[id*nSamples+is]-relCoeff[id][0])/relCoeff[id][1];
#else
		relCoeff[id][0]= mean ( &dataVector[id*nSamples], &flags[id*nSamples], nSamples);
		relCoeff[id][1]=1.0;

		for (size_t is = 0; is <nSamples; is++)
			dataVector[id*nSamples+is]-=relCoeff[id][0];
#endif
	}

	this->createBaseMatrix(nSamples, nDetectors);//, 2,AzOff,ElOff);
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
	bool success = true;
	for (register size_t itsi = 0; itsi<nSp; itsi++){
		if (!finite(v2[itsi]) || !finite(tsi[itsi])){
			cerr<<"Nan detected on tsi vector"<<endl;
			success = false;
			break;
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
	if (!success){
		delete [] tsi;
		cs_spfree(phi);
		return NULL;
	}
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
	for (size_t i=0; i<nDetectors; i++){
		for (size_t j=0; j<nSamples; j++){
			cFlags[j] = flags[i*nSamples+j];
			//cFlags[j] = true;
		}
		linfit_flags(&dataVector[i*nSamples],cFlags, &v1[i*nSamples],cFlags, nSamples, &c0,&c1, USEMPFIT);

		for (size_t j=0; j<nSamples; j++){
			atmTemplate [i*nSamples+j]= (v1[i*nSamples+j]-c0)/c1;
			atmTemplate [i*nSamples+j] *= relCoeff[i][1];
			atmTemplate [i*nSamples+j] += relCoeff[i][0];
		}
	}

	delete [] cFlags;

	delete []v1;
	return atmTemplate;

}


void CleanBspline::removeLargeScaleResiduals (CleanBsplineDestriping method,double correlate){

	size_t nScans = telescope->scanIndex.ncols();
	size_t nDetectors = dataArray->getNDetectors();
	int* di = dataArray->getDetectorIndices();
	size_t si, ei, nSamples;

	double *dataVector = NULL;
	double *kVector = NULL;
	double *flagVector = NULL;
	double *flagVectorDist = NULL;
	double *atmVector = NULL;
	double *AzOffsets = new double [nDetectors];
	double *ElOffsets = new double [nDetectors];
	double *detSens = new double [nDetectors];
	gsl_vector *azVector= NULL;
	gsl_vector *elVector= NULL;

	for (size_t id = 0; id < nDetectors; id ++){
		AzOffsets[id]=dataArray->detectors[di[id]].azOffset;
		ElOffsets[id]=dataArray->detectors[di[id]].elOffset;
		detSens[id] = dataArray->detectors[di[id]].getSensitivity()*1e-3; //Set value to Jy
	}

	double sampleRate = dataArray->detectors[di[0]].getSamplerate();

		for(size_t k=0;k<nScans;k++){
			currScan = k+1;
			if (ap->getOrder() != 100)
				if (!scanStatus[k])
					continue;
			si=telescope->scanIndex[0][k];
			ei=telescope->scanIndex[1][k]+1;
			nSamples = ei-si;

			dataVector = new double [nDetectors*nSamples];
			kVector = new double [nDetectors*nSamples];
			flagVector = new double [nDetectors*nSamples];
			flagVectorDist = new double [nDetectors*nSamples];
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
					flagVectorDist[i*nSamples+j] = 1.;
					if (flagVector[i*nSamples+j] ==0)
						fullFlags[j] =0;
					if (azVector!=NULL){
						gsl_vector_set(azVector, i*nSamples+j,dataArray->detectors[di[i]].hAzPhys[si+j]);
						gsl_vector_set(elVector, i*nSamples+j,dataArray->detectors[di[i]].hElPhys[si+j]);
					}
					double dist = sqrt (pow(dataArray->detectors[di[i]].hRa[si+j],2.0)+pow(dataArray->detectors[di[i]].hDec[si+j],2.0));
					dist *= 180.0/M_PI;
					if (dist<this->bright)
						flagVectorDist[i*nSamples+j] = 0;
				}

				double ngood=0;
				detSens[i] = stddev(dataVector, flagVectorDist, i*nSamples,(i+1)*nSamples-1, &ngood);
				if (detSens[i]!=detSens[i]){
					cerr<<"Error, no  good data on Az/El template calculation on file "<<ap->getDataFile()<< "valid points are: " <<ngood<<endl;
					exit(-1);
				}
			}

			switch (method){
			case STRIPE_PCA:
				pca (dataVector, flagVector, nSamples, nDetectors, ap->getCleanStripe());
				break;
			case STRIPE_FFT:
				cleanStripeFFT(dataVector, kVector, nDetectors, nSamples,correlate,sampleRate);
				break;
			case AZEL_TEMPLATE:
					azelResidualLMT(dataVector,flagVector,flagVectorDist, kVector, nDetectors, nSamples, azVector, elVector,detSens);
				break;
			case STRIPE_NONE:
			default:
				cerr<<"CleanSpline("<<tid<<") Fatal Error: not a valid destriping method. Check you Analysis Parameter file. Program Terminated"<<endl;
				exit(-1);

			}

			for(size_t i=0; i<nDetectors; i++)
				for(size_t j=0; j<nSamples; j++){
					dataArray->detectors[di[i]].hValues[si+j]=dataVector[i*nSamples +j];
					dataArray->detectors[di[i]].hKernel[si+j]=kVector[i*nSamples +j];
					dataArray->detectors[di[i]].hSampleFlags[si+j]=flagVector[i*nSamples+j];
				}

			delete [] dataVector;
			delete [] kVector;
			delete [] atmVector;
			delete [] flagVector;
			delete [] flagVectorDist;
			gsl_vector_free(azVector);
			gsl_vector_free (elVector);
			azVector=NULL;
			elVector= NULL;
		}

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
	cerr<<"CleanBspline("<<tid<<"): Control Chunk is "<<controlChunk<< " Number of Splines: "<<nSpline<<endl;

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


void CleanBspline::createBaseMatrix(int nSamples, int nDetectors, int nCells, VecDoub azOff, VecDoub elOff){
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
	//Deal with book-keeping
	double maxx, minx, maxy, miny, resx,resy;
	int xi, yi;
	maxmin(azOff,&maxx,&minx);
	maxmin(elOff,&maxy,&miny);
	resx = (maxx-minx)/nCells;
	resy = (maxy-miny)/nCells;
	size_t celli = pow(nCells,2);
	VecInt bpos (nDetectors,0);

	for (size_t ibolo = 0; ibolo< (size_t)nDetectors; ibolo++){
		xi = floor((azOff[ibolo]-minx)/resx);
		yi = floor((elOff[ibolo]-miny)/resy);
		if (xi == nCells)
			xi--;
		if (yi == nCells)
			yi--;
		bpos[ibolo] = xi*nCells +yi;
	}


	tmpbMatrix = cs_spalloc(celli*nSpline, dataLen, celli*nSpline,1,1);
	tmpbMatrix_t = cs_spalloc(dataLen, celli*nSpline,celli*nSpline,1,1);

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
				cs_entry(tmpbMatrix, bpos[j]*nSpline + k, j*baseMatrixSamples + i , c_sample);
				cs_entry(tmpbMatrix_t,j*baseMatrixSamples + i , bpos[j]*nSpline + k , c_sample);
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
	cerr<<"Downsample data: "<< nSamples <<","<<downSample<<","<< addOne<<endl;
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


bool CleanBspline::pca (double *dataVector, double *flags, size_t nSamples, size_t nDetectors, size_t neig2cut){

	gsl_matrix *det = gsl_matrix_alloc (nDetectors, nSamples);
	//gsl_matrix *ket = gsl_matrix_alloc(nDetectors, nSamples);
	gsl_matrix *fl = gsl_matrix_alloc (nDetectors, nSamples);
	gsl_matrix *denom = gsl_matrix_alloc(nDetectors,nDetectors);
	gsl_matrix *pcaCorr = gsl_matrix_alloc(nDetectors,nDetectors);
	gsl_eigen_symmv_workspace* w=gsl_eigen_symmv_alloc(nDetectors);
	gsl_vector* eVals = gsl_vector_alloc(nDetectors);
	gsl_matrix* eVecs = gsl_matrix_alloc(nDetectors,nDetectors);
	gsl_matrix* eFunc = gsl_matrix_alloc(nDetectors,nSamples);
	double c1,c0;
	VecDoub flagt(nSamples);
	VecDoub datat(nSamples);

	//Copy data

	for (size_t ibolo = 0; ibolo< nDetectors; ibolo++)
		for (size_t isample=0; isample<nSamples; isample++){
			gsl_matrix_set (det, ibolo, isample, dataVector[ibolo*nSamples+isample]);
			gsl_matrix_set (fl, ibolo, isample, flags[ibolo*nSamples+isample]);
		}


	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.,fl,fl,0.,denom);
	gsl_matrix_add_constant(denom,-1.);
	gsl_matrix_mul_elements (det,fl);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.,det,det,0.,pcaCorr);
	gsl_matrix_div_elements(pcaCorr,denom);

	gsl_matrix_free (fl);

	gsl_eigen_symmv(pcaCorr,eVals,eVecs,w);
	gsl_eigen_symmv_sort(eVals,eVecs, GSL_EIGEN_SORT_ABS_DESC);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,eVecs,det,0.,eFunc);

	for (size_t ibolo = neig2cut; ibolo< nDetectors; ibolo++)
		for (size_t isample=0; isample<nSamples; isample++){
			gsl_matrix_set (eFunc, ibolo, isample, 0.0);
			//gsl_matrix_set (kFunc, ibolo, isample, 0.0);
	}

	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,eVecs,eFunc,0.,det);

	for (size_t ibolo = 0; ibolo< nDetectors; ibolo++){
		for (size_t isample=0; isample<nSamples; isample++)
			datat[isample]=gsl_matrix_get(det,ibolo,isample);
		linfit_flags(&dataVector[ibolo*nSamples], &flags[ibolo*nSamples], &datat[0], &flags[ibolo*nSamples], nSamples, &c0, &c1, USEMPFIT);
		for (size_t isample=0; isample<nSamples; isample++)
			dataVector[ibolo*nSamples+isample]-= datat[isample]/c1;
	}
	gsl_matrix_free (det);
	gsl_eigen_symmv_free(w);
	gsl_matrix_free(pcaCorr);
	gsl_matrix_free (denom);
	gsl_vector_free(eVals);
	gsl_matrix_free (eVecs);
	gsl_matrix_free(eFunc);

	return true;
}


MatDoub CleanBspline::linearGainCorrection(double *dataVector, size_t nDetectors, size_t nSamples, double *fakeAtm){


	MatDoub relCoeff(2,nDetectors);    //[0] is relGain, [1] relOffset
	double cv0,cv1,cv2, chisq;


	for (size_t i=0; i<nDetectors; i++){
		gsl_fit_linear(&dataVector[0],1,&(dataVector[i*nSamples]),1,nSamples,&(relCoeff[1][i]),&(relCoeff[0][i]),&cv0,&cv1,&cv2,&chisq);
	}

	for (size_t i=0; i<nDetectors; i++)
		for (size_t j=0; j<nSamples; j++){
			dataVector[i*nSamples+j]= (dataVector[i*nSamples+j]-relCoeff[1][i])/relCoeff[0][i];
			if (fakeAtm)
				fakeAtm[i*nSamples+j]= (fakeAtm[i*nSamples+j]-relCoeff[1][i])/relCoeff[0][i];
		}

	return relCoeff;
}







bool CleanBspline::fixFlags(double *dataVector, bool *flagsVector, int nDetectors, long nSamples){

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

		for (size_t i=0; i<oSamples; i++){
			if (round(i*increment)>=nSamples){
				cerr<<"bad interpolation value"<<endl;
				exit(-1);
			}
			xTemp[i] = round(i*increment);
		}

		if (xTemp[oSamples-1] != nSamples-1){
			cerr<<"bad interpolation value"<<endl;
			exit(-1);
		}
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


void CleanBspline::cleanStripeFFT(double *data, double *kernel, size_t nDetectors, size_t nSamples, double cutOff, double sampleRate){

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
//
//void CleanBspline::splineResidual (double *data, double *flags, size_t nDetectors, size_t nSamples, gsl_vector  *azVector, gsl_vector  *elVector, double correlate){
//	VecDoub correlation (nDetectors);
//	VecDoub signs (nDetectors);
//	VecBool bflags (nDetectors*nSamples, true);
//	VecDoub stddevs (nDetectors);
//	double minCor = ap->getCleanStripe();
//	double iCor, maxCor;
//	size_t maxIt =100, curIt = 0;
//	size_t indexMaxStd=0;
//	double maxStd = 0., std;
//	double increment=1.;
//	size_t rSamples=nSamples;
//	size_t index;
//
//
//	if (resample > 1){
//		long nnSamples = long(round(nSamples/resample));
//		increment = double(nSamples-1)/double(nnSamples-1);
//		rSamples =nnSamples;
//	}
//
//	VecDoub rdata (nDetectors*rSamples);
//	VecBool rflags (nDetectors*rSamples);
//	gsl_vector *rAz = gsl_vector_alloc (nDetectors*rSamples);
//	gsl_vector *rEl = gsl_vector_alloc (nDetectors*rSamples);
//
//	do{
//		//Calculate maximum std bolomter
//
//		maxStd = 0.;
//		maxCor = 0.0;
//		for (size_t ibolo = 0; ibolo<nDetectors; ibolo++){
//			stddevs[ibolo] = stddev(data, flags,ibolo*nSamples, (ibolo+1)*nSamples);
//		}
//
//		maxStd = median(stddevs);
//
//		for (size_t ibolo = 0; ibolo<nDetectors; ibolo++){
//			if (maxStd == stddevs[ibolo]){
//				indexMaxStd = ibolo;
//				break;
//			}
//		}
//
//		for (size_t ibolo = 0; ibolo<nDetectors; ibolo++){
//			 iCor = flagCorrelation(&data[ibolo*nSamples], &flags[ibolo*nSamples], &data[indexMaxStd*nSamples], &flags[indexMaxStd*nSamples],nSamples);
//			 signs[ibolo] = iCor >= 0?1:-1;
//			 if (ibolo != indexMaxStd && iCor >maxCor){
//				 maxCor = iCor;
//			 }
//		}
//		cout<<"CleanBspline::splineResidual("<<tid<<"). Maximum correlation on iteration "<<curIt<< " is: "<<maxCor<<endl;
//		if (maxCor < minCor)
//			break;
//
//		for (size_t jbolo = 0; jbolo < nDetectors; jbolo++){
//			for (size_t jSample = 0; jSample < rSamples; jSample++){
//				index = long(round(jSample*increment));
//				rdata[jbolo*rSamples+jSample] = signs[jbolo]* data [jbolo*nSamples+index];
//				rflags [jbolo*rSamples+jSample] = (bool) flags [jbolo*nSamples+index];
//				gsl_vector_set (rAz,jbolo*rSamples+jSample,gsl_vector_get(azVector,jbolo*nSamples+index));
//				gsl_vector_set (rEl, jbolo*rSamples+jSample, gsl_vector_get ( elVector ,jbolo*nSamples+index));
//			}
//		}
//
//		//writeVecOut("rData.txt", &rdata[0], nDetectors*rSamples );
//
//		double *tmplt = cottingham(&rdata[0], rAz, rEl, &rflags[0], nDetectors, rSamples, this->ap->cleanPixelSize, (int)indexMaxStd);
//
//		//writeVecOut("lData.txt", &rdata[0], nDetectors*rSamples );
//		//writeVecOut("sTemplate.txt", tmplt, rSamples*nDetectors);
//
//		//exit(-1);
//
//
//		for (size_t ibolo = 0; ibolo<nDetectors; ibolo++){
//			subtractTemplate (&data[ibolo*nSamples], nSamples, &tmplt[ibolo*nSamples], rSamples, increment, NULL, true);
//			for (size_t jSample = 0; jSample < rSamples; jSample++){
//				 data [ibolo*nSamples+jSample] *= signs[ibolo];
//			}
//		}
//
//	}while (++curIt <maxIt);
//
//	gsl_vector_free(rAz);
//	gsl_vector_free(rEl);
//
////	for (size_t i=0; i<nDetectors;i++){
////		corrMatrix[i][i] = 1.0;
////		for (size_t j=i+1; j<nDetectors; j++){
////			corrMatrix[i][j] = flagCorrelation(&data[i*nSamples], &flags[i*nSamples], &data[j*nSamples], &flags[j*nSamples],nSamples);
////			dist = sqrt(pow(gsl_vector_get(azVector,i*nSamples)-gsl_vector_get(azVector,j*nSamples),2) + pow(gsl_vector_get(elVector,i*nSamples)-gsl_vector_get(elVector,j*nSamples),2));
////			if (round(corrMatrix[i][j]*100.0)/100.0 < correlate || dist < mindist)
////				corrMatrix[i][j]=0.0;
////			corrMatrix[j][i]=corrMatrix[i][j];// = abs(corrMatrix[i][j]);
////		}
////	}
////	MatDoub outData (nDetectors, nSamples, 0.0);
////	VecBool isCleaned (nDetectors, false);
////
////	size_t oSamples=nSamples;
////
////	for (size_t ibolo = 0; ibolo < nDetectors; ibolo++){
////		//Get the number of correlated bolometers
////		size_t nBoloCorr = 0;
////		for (size_t jbolo = 0; jbolo < nDetectors; jbolo++)
////			if (corrMatrix [ibolo][jbolo]> 0.0)
////				nBoloCorr++;
////		//cerr<<"Number of correlated bolometers"<< ibolo << "->"<< nBoloCorr<<endl;
////		if (nBoloCorr >= 2){
////
////			double increment = 1.0;
////			size_t index;
////
////
////			if (resample > 1){
////				long nnSamples = long(round(oSamples/resample));
////				increment = double(oSamples-1)/double(nnSamples-1);
////				nSamples =nnSamples;
////			}
////
////			size_t dataLen = nBoloCorr*nSamples;
////
////			size_t boloIndex [nBoloCorr];
////			boloIndex[0] = ibolo;
////			size_t curBolo = 1;
////			for (size_t jbolo = 0; jbolo < nDetectors; jbolo++)
////				if (corrMatrix[ibolo][jbolo] > 0.0 && ibolo != jbolo){
////					boloIndex[curBolo++]= jbolo;
////					//cerr <<jbolo<<endl;
////				}
////			double *cBoloData = new double [dataLen];
////			bool *cFlags = new bool [dataLen];
////			gsl_vector *cAzVector = gsl_vector_alloc(dataLen);
////			gsl_vector *cElVector = gsl_vector_alloc(dataLen);
////			//copy data
////			for (size_t jbolo = 0; jbolo < nBoloCorr; jbolo++){
////				for (size_t jSample = 0; jSample < nSamples; jSample++){
////					index = long(round(jSample*increment));
////					cBoloData[jbolo*nSamples+jSample] = data [boloIndex[jbolo]*oSamples+index];
////					cFlags [jbolo*nSamples+jSample] = flags [boloIndex[jbolo]*oSamples+index];
////					gsl_vector_set (cAzVector,jbolo*nSamples+jSample,gsl_vector_get(azVector,boloIndex[jbolo]*oSamples+index));
////					gsl_vector_set (cElVector, jbolo*nSamples+jSample, gsl_vector_get ( elVector ,boloIndex[jbolo]*oSamples+index));
////				}
////			}
////			//writeVecOut ("osData", data, nDetectors*nSamples);
////			//writeVecOut("sData.txt", cBoloData, nBoloCorr*nSamples);
////			double *tmplt = cottingham(cBoloData, cAzVector, cElVector, cFlags, nBoloCorr, nSamples, this->ap->cleanPixelSize, 0);
////
////
////			//writeVecOut("sTemplate.txt", tmplt, nBoloCorr*nSamples);
////			//exit(-1);
////
////			subtractTemplate (&outData[ibolo][0], oSamples, &tmplt[0], nSamples, increment, NULL, true);
////			//for (size_t iSample = 0; iSample<nSamples; iSample++)
////				//outData[ibolo][iSample] = tmplt [iSample];
////			delete [] tmplt;
////			delete [] cFlags;
////			delete [] cBoloData;
////			gsl_vector_free (cAzVector);
////			gsl_vector_free (cElVector);
////			isCleaned[ibolo]=true;
////
////		}
////	}
////
////	nSamples = oSamples;
////	for (size_t ibolo= 0; ibolo< nDetectors; ibolo++)
////		if (isCleaned[ibolo])
////			for (size_t jSample = 0; jSample<nSamples; jSample++)
////				data[ibolo*nSamples+jSample] -= outData[ibolo][jSample];
//}
//
//
//
//void CleanBspline::splineIndividual (double *data, double *flags, size_t nDetectors, size_t nSamples, size_t nbreaks)
//{
//
//    size_t order =4;
//    size_t nSpline = nbreaks + order -2;
//    //cerr<<"CleanBspline("<<tid<<"): Control Chunk is "<<controlChunk<< " Number of Splines: "<<nSpline<<endl;
//
//    gsl_bspline_workspace *bsw = gsl_bspline_alloc(order, nbreaks);
//    if (bsw == NULL){
//        cerr<<"CleanBspline::splineIndividual():Cannot allocate B-spline working space. Imploding"<< endl;
//        exit(-1);
//    }
//    gsl_matrix *baseMatrix = gsl_matrix_alloc(nSamples, nSpline);
//    gsl_bspline_knots_uniform(0.0, (double)(nSamples-1), bsw);
//    gsl_vector *tmpD = gsl_vector_alloc(nSamples);
//    gsl_vector *tmpT= gsl_vector_alloc(nSamples);
//    gsl_vector *tmpB = gsl_vector_alloc(nSpline);
//    gsl_vector *coeff = gsl_vector_alloc(nSpline);
//
//    for (size_t i = 0; i < nSamples; i++)
//    {
//      /* compute B_j(xi) for all j */
//      gsl_bspline_eval(i, tmpB, bsw);
//
//      /* fill in row i of X */
//      for (size_t j = 0; j < nSpline; j++)
//        {
//          gsl_matrix_set(baseMatrix, i, j, gsl_vector_get(tmpB, j));
//        }
//    }
//    double chisq;
//    gsl_matrix *cov = gsl_matrix_alloc(nSpline,nSpline);
//    gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(nSamples, nSpline);
//    for (size_t i=0; i<nDetectors; i++){
//        for (size_t j=0; j<nSamples; j++)
//            gsl_vector_set(tmpD,j,data[i*nSamples+j]);
//        gsl_multifit_linear(baseMatrix, tmpD, coeff, cov, &chisq, mw);
//        gsl_blas_dgemv(CblasNoTrans,1.0,baseMatrix,coeff,0.0,tmpT);
//        gsl_vector_sub(tmpD,tmpT);
//        for (size_t j=0; j<nSamples; j++)
//            data[i*nSamples+j] = gsl_vector_get(tmpD,j);
//    }
//
//    gsl_vector_free(tmpD);
//    gsl_vector_free(tmpT);
//    gsl_vector_free(tmpB);
//    gsl_vector_free(coeff);
//    gsl_matrix_free(cov);
//    gsl_matrix_free(baseMatrix);
//    gsl_multifit_linear_free(mw);
//    gsl_bspline_free(bsw);
//}




void CleanBspline::azelResidualLMT (double *data, double *flags, double *flagdis, double *kvector, size_t nDetectors, size_t nSamples, gsl_vector *azOffsets, gsl_vector *elOffsets, double *sens){
	gsl_matrix * cBoloData = gsl_matrix_alloc(nDetectors, nSamples);
	gsl_matrix * cBoloFlag = gsl_matrix_alloc (nDetectors, nSamples);
	gsl_matrix *az = gsl_matrix_alloc(nDetectors, nSamples);
	gsl_matrix *el = gsl_matrix_alloc(nDetectors, nSamples);
	gsl_vector *corrCoefs = gsl_vector_alloc (nDetectors);

	for (size_t ibolo = 0; ibolo < nDetectors; ibolo++){
			for (size_t jSample = 0; jSample < nSamples; jSample++){
				gsl_matrix_set (cBoloData,ibolo,jSample,data[ibolo*nSamples +jSample]);
				gsl_matrix_set (cBoloFlag, ibolo,jSample, (int)flags [ibolo*nSamples+ jSample] & (int)flagdis[ibolo*nSamples+ jSample]);
				gsl_matrix_set (az, ibolo, jSample, gsl_vector_get(azOffsets,ibolo*nSamples+jSample));
				gsl_matrix_set (el, ibolo, jSample, gsl_vector_get(elOffsets,ibolo*nSamples+jSample));
			}

		gsl_vector_set (corrCoefs, ibolo, 1.0/pow(sens[ibolo],2.0));

	}
	AzElTemplateCalculator aeTemp  (cBoloData, cBoloFlag,az,el);
	aeTemp.overrideMode(int(ap->getCleanStripe()));
	aeTemp.calculateTemplate2(corrCoefs);
//	aeTemp.decorrelateCoeffs(1);
//	aeTemp.updateTemplate();
	double atTmp;
	gsl_matrix *atmTemplate = aeTemp.getTemplate();
	for (size_t ibolo = 0; ibolo < nDetectors; ibolo++){
		for (size_t iSample = 0; iSample<nSamples; iSample++){
			atTmp = gsl_matrix_get (atmTemplate, ibolo, iSample);
			if (atTmp!=atTmp){
				data[ibolo*nSamples+iSample] = 0.0;
				flags[ibolo*nSamples+iSample] = 0;
				cerr<<"CleanBSpline("<<tid<<")Az/El Template Object did not produce any valid data on file: " << ap->getDataFile()<<endl;
				//exit(-1);
			}
			else{
				data[ibolo*nSamples+iSample] -= atTmp;
				if (this->cleanKernel)
					data[ibolo*nSamples+iSample] -= atTmp;
				gsl_matrix_set(cBoloData, ibolo, iSample,data[ibolo*nSamples+iSample]);
			}
		}
	}

	gsl_matrix_free(cBoloData);
	gsl_matrix_free(cBoloFlag);
	gsl_vector_free(corrCoefs);
	gsl_matrix_free(az);
	gsl_matrix_free(el);
	
}

//Getters

cs* CleanBspline::getBaseMatrix(){
	return baseMatrix;
}
