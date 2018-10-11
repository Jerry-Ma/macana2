#ifndef _CLEAN_APCA_H
#define _CLEAN_APCA_H


#include "Clean.h"

#include <gsl/gsl_matrix.h>

class CleanAPCA: public Clean{
private:

	size_t cutStd;
	double fcut;
	double sampleRate;
	size_t tid;

	VecDoub generateHighPassData(double *rawData, size_t nDetectors, size_t nSamples);
	VecDoub generateLowPassData (double *rawData, size_t nDetectors, size_t nSamples);
	size_t applyPCA (gsl_matrix *det, gsl_matrix *ket);
	VecDoub filterData (double *rawData, size_t nDetectors, size_t nSamples, bool type = false);

public:
	CleanAPCA(Array* dataArray, Telescope *telescope);
    ~CleanAPCA();
	bool clean();
};


#endif
