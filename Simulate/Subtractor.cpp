#include "Subtractor.h"
#include "vector_utilities.h"


 Subtractor::Subtractor(MapNcFile* map):SimulatorInserter(map) {
	 this->map = map;
	 atmFreq = 0.0;
	 seed=0.0;
	 noiseChunk = 0;
	 fluxFactor = 1.0;
	 sigOnly = false;
	 isTemp = false;
	 resample=0;
}

void Subtractor::subtract(Array* dataArray) {
	insertIntoArray(dataArray);
}

Subtractor::~Subtractor(){

}
