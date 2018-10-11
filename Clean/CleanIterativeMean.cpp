#include "nr3.h"
#include "vector_utilities.h"
#include "CleanIterativeMean.h"

CleanIterativeMean::CleanIterativeMean(Array* dataArray, Telescope *telescope):
	Clean(dataArray,telescope){
	this->maxDelta=1e-4;
	this->maxIteration = 100;
}

void CleanIterativeMean::signMultiply(VecDoub signs) {

	size_t nBolo = dataArray->getNDetectors();
	int *di = dataArray->getDetectorIndices();
	size_t nSamples = dataArray->getNSamples();


	for (size_t ibolo=0; ibolo<nBolo; ibolo++)
		if (signs[ibolo]<0){
			for (size_t isample=0; isample<nSamples; isample++)
				dataArray->detectors[di[ibolo]].hValues[isample]*= -1;
		}
}



bool CleanIterativeMean::clean(){

	size_t nBolo = dataArray->getNDetectors();
	int *di = dataArray->getDetectorIndices();
	size_t nSamples = dataArray->getNSamples();
	double m, prevStd, newStd,delta,startStd,maxCorr;
	size_t curIt = 0;
	VecDoub signs;

	for (size_t ibolo = 0; ibolo<nBolo; ibolo++){
		//m = mean(dataArray->detectors[di[ibolo]].hValues,dataArray->detectors[di[ibolo]].hSampleFlags);
		m = median(dataArray->detectors[di[ibolo]].hValues);
		for (size_t jsample=0; jsample<nSamples; jsample++)
			dataArray->detectors[di[ibolo]].hValues[jsample]-=m;
	}
	size_t ndata = dataArray->detectors[di[0]].hValues.size();
	double medi=0;
	VecDoub isample (nBolo);
	for (size_t iss = 0; iss<ndata; iss++){
		for (size_t ibolo = 0; ibolo<nBolo; ibolo++){
			isample[ibolo]= dataArray->detectors[di[ibolo]].hValues[iss];
		}
		medi = median(isample);
		for (size_t ibolo = 0; ibolo<nBolo; ibolo++){
			dataArray->detectors[di[ibolo]].hValues[iss]-=medi;
		}
	}

//	startStd =prevStd = linearMeanFit();
//
//	cout <<"Starting mean Deviation: "<<prevStd<<endl;
//	do{
//		signs = getSigns(&maxCorr);
//		if (maxCorr <0.6)
//			break;
//		signMultiply(signs);
//		newStd = linearMeanFit();
//		signMultiply(signs);
//		delta = abs(newStd)/startStd;
//		prevStd = newStd;
//
//	}while(++curIt<maxIteration);

	return true;
}


VecDoub CleanIterativeMean::getSigns(double *maxCorr){
	size_t nBolo = dataArray->getNDetectors();
	int *di = dataArray->getDetectorIndices();
	size_t nSamples = dataArray->getNSamples();
	size_t mIndex = 0;
	double maxStd = 0.0;

	VecDoub result (nBolo, 0.0);

	*maxCorr = 0;

	for (size_t ibolo = 0; ibolo<nBolo; ibolo++)
		if (stddev(dataArray->detectors[di[ibolo]].hValues, dataArray->detectors[di[ibolo]].hSampleFlags) > maxStd)
			mIndex = ibolo;

	for (size_t ibolo = 0; ibolo<nBolo; ibolo++){
		result[ibolo]= flagCorrelation(&dataArray->detectors[di[ibolo]].hValues[0], &dataArray->detectors[di[ibolo]].hSampleFlags[0],\
									   &dataArray->detectors[di[mIndex]].hValues[0], &dataArray->detectors[di[mIndex]].hSampleFlags[0],\
									   nSamples);
		if (abs(result[ibolo])>*maxCorr && ibolo!=mIndex)
			*maxCorr = abs(result[ibolo]);
	}
	return result;
}

double CleanIterativeMean::linearMeanFit(){
	double meanStd=0.0;
	size_t nBolo = dataArray->getNDetectors();
	int *di = dataArray->getDetectorIndices();
	size_t nSamples = dataArray->getNSamples();
	VecDoub boloMean (nSamples,0.);
	VecDoub boloSamples (nSamples, 0.);
	VecBool boloMeanFlags (nSamples, false);

	for (size_t ibolo = 0; ibolo<nBolo; ibolo++)
		for (size_t jsample = 0; jsample<nSamples; jsample++){
			boloMean[jsample]+=dataArray->detectors[di[ibolo]].hValues[jsample]*
								dataArray->detectors[di[ibolo]].hSampleFlags[jsample];
			boloSamples[jsample] +=	dataArray->detectors[di[ibolo]].hSampleFlags[jsample]?1.:0.;
		}
	for (size_t jsample=0; jsample<nSamples; jsample++)
		if (boloSamples[jsample]>0){
			boloMeanFlags[jsample] = true;
			boloMean[jsample]/=boloSamples[jsample];
		}


	meanStd = stddev(boloMean,boloMeanFlags);

	double c0,c1;
	for (size_t ibolo = 0; ibolo<nBolo; ibolo++){
		linfit_flags(&boloMean[0], &boloMeanFlags[0], \
				&dataArray->detectors[di[ibolo]].hValues[0], &dataArray->detectors[di[ibolo]].hSampleFlags[0] , \
				nSamples, &c0, &c1, false);
		for (size_t jsample=0; jsample<nSamples; jsample++){
			dataArray->detectors[di[ibolo]].hValues[jsample]-= c1*boloMean[jsample]+ c0;
		}
	}
	return meanStd;
}

CleanIterativeMean::~CleanIterativeMean() {
}
