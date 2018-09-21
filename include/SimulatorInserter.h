#ifndef SIMULATOR_INSERTER_H
#define SIMULATOR_INSERTER_H


#include <gsl/gsl_rng.h>

#include "MapNcFile.h"

#include "SimParams.h"
#include "Detector.h"
#include "Array.h"


class SimulatorInserter{
	private:
		void init();
		void createNoiseGenerator();
		void deleteNoiseGenerator();
		VecDoub createNoiseSignal(Detector det);
	protected:
		double atmFreq;
		double noiseChunk;
		double fluxFactor;
		bool sigOnly;
		Array *arr;
		MapNcFile *map;
		//Noise
		gsl_rng *r;
		bool isTemp;
		bool interpolateMapSignals(Detector *di);

	public:
		SimulatorInserter(MapNcFile *map,SimParams *spar);
		SimulatorInserter(MapNcFile *map);
		bool insertIntoArray (Array *arr);
		~SimulatorInserter();
};

#endif
