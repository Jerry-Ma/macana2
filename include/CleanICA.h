#ifndef _CLEAN_ICA_H
#define _CLEAN_ICA_H
#include <omp.h>
#include <itpp/itbase.h>
#include "Clean.h"


class CleanICA: public Clean{
private:
	size_t neig2Cut;				////Number of ICA components to estimate
	size_t tid;						////Thread id for output formatting
	bool averageType;
	void componentFit (double *rawData, size_t nSamples, itpp::mat components, size_t ncomp);		////Use a linear combination of the decomposed components to estimate an atmosphere template for each bolometer
	void removeArrayAverage (size_t si, size_t nSamples, int refBolo=-1);			////Removed average signal from all bolometers on a scan basis
	void removeFullAverage (int refBolo=-1);			////Removed average signal from all bolometers without taking into account the scan limits
	VecBool createFlags (Detector det, size_t si, size_t nSamples);					////Update flags for average signal estimation
	

public:
	CleanICA(Array* dataArray, Telescope *telescope);
    ~CleanICA();
	bool clean();																	////Start cleaning process
	void setNeig2Cut(int neig2Cut);													////Set the number of component
};


#endif
