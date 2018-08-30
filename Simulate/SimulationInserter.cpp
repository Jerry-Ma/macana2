#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_trig.h>
#include <cmath>

#include "Detector.h"
#include "SimulatorInserter.h"
#include "SimParams.h"
#include "SBSM.h"
#include "AtmosphereRandom.h"
#include "vector_utilities.h"

SimulatorInserter::SimulatorInserter(MapNcFile *map, SimParams *spar){

	this->map = map;
	this->r= NULL;
	if (spar != NULL){

		this->sigOnly = !spar->getAddSignal();
		this->atmFreq = spar->getAtmFreq();
		this->atmSeed = spar->getAtmSeed();
		this->fluxFactor = -1.0 * spar->getFluxFactor();
		this->noiseChunk = spar->getNoiseChunk();
		this->resample = spar->getResample();

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
	this->atmSeed = 0;
	this->fluxFactor = -1.0;
	this->noiseChunk = 0;
	this->resample=0;
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
	AtmosphereRandom *atmrng=NULL;
	bool sameAtm = false;
	VecDoub atmTemplate (0);

	double sr = arr->detectors[di[0]].getSamplerate();
	double minFreq = sr/(2.0*arr->detectors[di[0]].getNSamples());

	if (!this->isTemp){
		if (this->atmFreq !=0 || this->atmSeed !=0){
//					if (this->atmFreq !=0. && this->atmSeed !=0){
//						cerr<<"Simulator Inserter Error. You must select either atmFreq != 0 or atmFreq !=0. Check your AP File"<<endl;
//						exit(-1);
//					}
					double atmFreq =this->atmFreq;
					double seed = this->atmSeed;
					if (this->atmFreq < 0){
						sameAtm = true;
						atmFreq*=-1;
					}

					if (atmFreq > sr/2.0){
						cerr<<"SimulatorInserter()::insertIntoArray(). Requested Atm Freq exceed sampling frequency. Clipping to " <<sr/2.0<<"Hz"<<endl;
						atmFreq = sr/2.0;
					}
					if (atmFreq < minFreq && atmFreq > 0){
						atmFreq = minFreq;
						cerr<<"SimulatorInserter()::insertIntoArray(). Requested Atm Freq is less than minimum frequency."<<"Clipping to " <<minFreq<<"Hz"<<endl;
					}
					size_t nBreaks = size_t(round(arr->detectors[di[0]].getNSamples()/sr*atmFreq));
					if (seed < 0)
						seed = -1;
					if (atmFreq !=0){
						bspline = new SBSM(4,arr->detectors[di[0]].getNSamples(),nBreaks);
						cout<<"SimulatorInserter()::insertIntoArray(). Creating atmosphere signals with frequency:" << atmFreq << "Hz"<<endl;
						if (!bspline){
							cerr<<"SimulatorInserter()::insertIntoArray(). Not enough memory to create atmosphere template"<<endl;
							exit (-1);
						}
					}
					else if (atmSeed > 0){
						atmrng = new AtmosphereRandom(arr->detectors[di[0]].getNSamples(), seed);
						cout<<"SimulatorInserter()::insertIntoArray(). Creating random 1/f atmosphere signals with seed:" << atmrng->getSeed()<<endl;
					}
		}
		double calFactor;
		for (size_t i=0; i<nbolo; i++){
			np = arr->detectors[di[i]].hValues.size();
			calFactor = arr->detectors[di[i]].getCalibrationFactor();
			interp.resize(np);
			interp = this->map->fastMapSignal(arr->detectors[di[i]].hRa, arr->detectors[di[i]].hDec);
			if (bspline){
						if (sameAtm && i==0){
							atmTemplate = bspline->fitData(arr->detectors[di[0]].hValues);
						}else if (!sameAtm){
							atmTemplate = bspline->fitData(arr->detectors[di[i]].hValues);
						}
			}
			if (atmrng){

						if (sameAtm && i==0){
							atmTemplate = atmrng->getSimulatedData(arr->detectors[di[0]].hValues);
							
							//writeVecOut("atmTemp.txt", atmTemplate.getData(), atmTemplate.size());
						}else if (!sameAtm){
							atmTemplate = atmrng->getSimulatedData(arr->detectors[di[i]].hValues);
						}
						VecDoub gradient = atmrng->getSimulatedAirmassGradient (arr->detectors[di[i]].hEl, arr->getAvgTau());
						//cout<<"SimulationInserter. Generated gradient with average opacity: " << arr->getAvgTau()<<endl;
						for (size_t ig = 0; ig<np; ig++)
							atmTemplate [ig] -= gradient[ig];
			}
			if (!sigOnly){
				for (size_t j=0; j<np; j++)
					arr->detectors[di[i]].hValues[j] += this->fluxFactor*interp[j]/calFactor;

			}else{

				for (size_t j=0; j<np; j++)
					arr->detectors[di[i]].hValues[j] = this->fluxFactor*interp[j]/calFactor;
			}
			if (this->noiseChunk != 0){
				noise = this->createNoiseSignal(arr->detectors[di[i]]);
				for (size_t j=0; j<np; j++)
					arr->detectors[di[i]].hValues[j] += noise[j];
			}
			if (atmTemplate.size()>0)
				for (size_t j=0; j<np; j++)
									arr->detectors[di[i]].hValues[j] += atmTemplate[j];
		}
		noise.resize(0);
	}else{
		cerr<<"Creating template signals"<<endl;
		for (size_t i=0; i<nbolo; i++){
				arr->detectors[di[i]].atmTemplate = this->map->fastMapSignal(arr->detectors[di[i]].hRa, arr->detectors[di[i]].hDec);
		}

	}
	if (bspline)
		delete bspline;
	if (atmrng)
		delete atmrng;

	if (resample > 0){
		if (resample > sr/2.0){
				cerr<<"SimulatorInserter()::insertIntoArray(). Requested Resample exceed sampling frequency. Clipping to " <<sr/2.0<<"Hz"<<endl;
				resample = sr/2.0;
			}
			if (resample < minFreq && resample > 0){
				resample = minFreq;
				cerr<<"SimulatorInserter()::insertIntoArray(). Requested Resample is less than minimum frequency."<<"Clipping to " <<minFreq<<"Hz"<<endl;
			}
			size_t nBreaks = size_t(round(arr->detectors[di[0]].getNSamples()/sr*resample));
			bspline = new SBSM(4,arr->detectors[di[0]].getNSamples(),nBreaks);
			cout<<"SimulatorInserter()::insertIntoArray(). Resampling bolometer signals with frequency:" << resample << "Hz"<<endl;
			if (!bspline){
				cerr<<"SimulatorInserter()::insertIntoArray(). Not enough memory to create resampler"<<endl;
				exit (-1);
			}

			for (size_t i =0; i<nbolo; i++){
				arr->detectors[di[i]].hValues = bspline->fitData(arr->detectors[di[i]].hValues);
			}

			delete bspline;
	}

	return true;
}

VecDoub SimulatorInserter::createNoiseSignal(Detector det){

	//Number of samples to get the stddev from
	double chunkSample = det.getSamplerate()*this->noiseChunk;
	size_t  nSamples = det.getNSamples();
	VecDoub noiseSignal (nSamples);
	double nScans = ceil(nSamples/chunkSample);
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

