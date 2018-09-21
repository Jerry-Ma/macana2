#include "Subtractor.h"
#include "vector_utilities.h"


 Subtractor::Subtractor(MapNcFile* map):SimulatorInserter(map) {
	 this->map = map;
	 atmFreq = 0.0;
	 noiseChunk = 0;
	 fluxFactor = 1.0;
	 sigOnly = false;
	 isTemp = false;
}

void Subtractor::subtract(Array* dataArray) {
	insertIntoArray(dataArray);
	//int *di = dataArray->getDetectorIndices();
}

Subtractor::~Subtractor(){

}
