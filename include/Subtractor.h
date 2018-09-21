#ifndef SUBTRACTOR_H
#define SUBTRACTOR_H

#include "SimulatorInserter.h"
#include "Array.h"
#include "MapNcFile.h"


class Subtractor : public SimulatorInserter{
	public:
		Subtractor (MapNcFile *map);
		void subtract (Array *dataArray);
		~Subtractor ();
};
#endif
