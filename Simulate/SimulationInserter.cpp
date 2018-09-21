#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "Detector.h"
#include "SimulatorInserter.h"
#include "SimParams.h"
#include "SBSM.h"
#include "vector_utilities.h"

SimulatorInserter::SimulatorInserter(MapNcFile *map, SimParams *spar){

	this->map = map;
	this->r= NULL;
	if (spar != NULL){

		this->sigOnly = !spar->getAddSignal();
		this->atmFreq = spar->getAtmFreq();
		this->fluxFactor = -1.0 * spar->getFluxFactor();
		this->noiseChunk = spar->getNoiseChunk();

		//if (this->fluxFactor == 0)
			//this->fluxFactor = -1.0;

		//Noise

		if (this->noiseChunk != 0)
			createNoiseGenerator();
		isTemp = false;
	} else{
		isTemp = true;
	}
}

SimulatorInserter::SimulatorInserter(MapNcFile *map){
	this->map = map;
	this->r=NULL;
	this->sigOnly = false;
	this->atmFreq = 0;
	this->fluxFactor = -1.0;
	this->noiseChunk = 0;
	isTemp = false;
}

bool SimulatorInserter::insertIntoArray(Array *arr){
	int *di = NULL;
	size_t nbolo = 0;
	size_t np = 0;
	VecDoub interp(0);
	VecDoub noise (0);
	arr->updateDetectorIndices();
	di = arr->getDetectorIndices();
	nbolo = arr->getNDetectors();

	SBSM *bspline = NULL;
	bool sameAtm = false;
	VecDoub atmTemplate (0);

	if (!this->isTemp){
	  if (this->atmFreq !=0){
	    
	    double atmFreq =this->atmFreq;
	    if (this->atmFreq < 0){
	      sameAtm = true;
	      atmFreq*=-1;
	    }
	    double sr = arr->detectors[di[0]].getSamplerate();
	    double minFreq = sr/(2.0*arr->detectors[di[0]].getNSamples());
	    if (atmFreq > sr/2.0){
	      cerr<<"SimulatorInserter()::insertIntoArray(). Requested Atm Freq exceed sampling frequency. Clipping to " <<sr/2.0<<"Hz"<<endl;
	      atmFreq = sr/2.0;
	    }
	    if (atmFreq < minFreq){
	      atmFreq = minFreq;
	      cerr<<"SimulatorInserter()::insertIntoArray(). Requested Atm Freq is less than minimum frequency."<<"Clipping to " <<minFreq<<"Hz"<<endl;
	    }
	    size_t nBreaks = size_t(round(arr->detectors[di[0]].getNSamples()/sr*atmFreq));
	    cout<<"SimulatorInserter()::insertIntoArray(). Creating atmosphere signals with frequency:" << atmFreq << "Hz"<<endl;
	    bspline = new SBSM(4,arr->detectors[di[0]].getNSamples(),nBreaks);
	    if (!bspline){
	      cerr<<"SimulatorInserter()::insertIntoArray(). Not enough memory to create atmosphere template"<<endl;
	      exit (-1);
	    }
	  }
	  double mBolo;
	  double calFactor;
	  for (size_t i=0; i<nbolo; i++){
	    np = arr->detectors[di[i]].hValues.size();
	    mBolo = median(&arr->detectors[di[i]].hValues[0], np);
	    calFactor = arr->detectors[di[i]].getCalibrationFactor();
	    interp.resize(np);
	    interp = this->map->fastMapSignal(arr->detectors[di[i]].hRa, arr->detectors[di[i]].hDec);
	    for (size_t is = 0; is < np; is++)
	      if (!finite(interp[is])){
		char buff [200];
		sprintf(buff, "Interp_%ul.txt", (unsigned int)i);
		writeVecOut(buff, atmTemplate.getData(), atmTemplate.size());
		exit(-1);
	      }
	    
	    if (bspline){
	      
	      if (sameAtm && i==0){
		atmTemplate = bspline->fitData(arr->detectors[di[0]].hValues);
		//writeVecOut("atmTemp.txt", atmTemplate.getData(), atmTemplate.size());
	      }else if (!sameAtm){
		atmTemplate = bspline->fitData(arr->detectors[di[i]].hValues);
		for (size_t is = 0; is < np; is++)
		  if (!finite(atmTemplate[is])){
		    char buff [200];
		    sprintf(buff, "atmTemp_%ul.txt", (unsigned int)i);
		    writeVecOut(buff, atmTemplate.getData(), atmTemplate.size());
		    sprintf(buff, "oData_%ul.txt", (unsigned int)i);
		    writeVecOut(buff, &arr->detectors[di[i]].hValues[0], np);
		    exit(-1);
		  }
		
	      }
	    }
	    if (!sigOnly){
	      for (size_t j=0; j<np; j++)
		arr->detectors[di[i]].hValues[j] += this->fluxFactor*interp[j]/calFactor;
	      
	    }else{
	      
	      for (size_t j=0; j<np; j++)
		arr->detectors[di[i]].hValues[j] = this->fluxFactor*interp[j]/calFactor;
	    }
	    //arr->detectors[di[i]].setAtmTemplate(interp);
	    if (this->noiseChunk != 0){
	      noise = this->createNoiseSignal(arr->detectors[di[i]]);
	      for (size_t j=0; j<np; j++)
		arr->detectors[di[i]].hValues[j] += noise[j];
	    }
	    if (bspline){
	      for (size_t j=0; j<np; j++)
		arr->detectors[di[i]].hValues[j] += atmTemplate[j];
	    } else{
	      for (size_t j=0; j<np; j++)
		arr->detectors[di[i]].hValues[j] += 0.0;
	    }
	    
	    
	//		if (i==17){
	//			writeVecOut ("interp_fast.txt",arr->detectors[di[i]].hValues.getData(), np);
	//			exit(-1);
	//		}

	  }
	  //exit(-1);
	  noise.resize(0);
	}else{
	  cerr<<"Creating template signals"<<endl;
	  for (size_t i=0; i<nbolo; i++){
	    arr->detectors[di[i]].atmTemplate = this->map->fastMapSignal(arr->detectors[di[i]].hRa, arr->detectors[di[i]].hDec);
	  }
	  
	}
	if (bspline)
	  delete bspline;
	return true;
}

VecDoub SimulatorInserter::createNoiseSignal(Detector det){

	//Number of samples to get the stddev from
	double chunkSample = det.getSamplerate()*this->noiseChunk;
	size_t  nSamples = det.getNSamples();
	VecDoub noiseSignal (nSamples);
	double scanSdev=0;
	size_t si=floor (chunkSample);
	size_t se=ceil (2*chunkSample);
	scanSdev = stddev(det.hValues,si,se);
	for (register size_t j=0; j<nSamples; j++)
			noiseSignal[j]= gsl_ran_gaussian(this->r,scanSdev);
	return noiseSignal;
}


void SimulatorInserter::createNoiseGenerator(){
	const gsl_rng_type * T=gsl_rng_default;
	this->r= gsl_rng_alloc (T);
}

void SimulatorInserter::deleteNoiseGenerator(){
	if (this->r != NULL)
		gsl_rng_free (this->r);
}

SimulatorInserter::~SimulatorInserter(){
	deleteNoiseGenerator();
}

