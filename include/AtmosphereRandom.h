#ifndef _ATMOSPHERE_RANDOM_H
#define _ATMOSPHERE_RANDOM_H

#include "nr3.h"
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

class AtmosphereRandom{
private:
	size_t nSamples;
	long seed;
	double fluxGradient;
	double *randomSequence;
	double *spectrum;
	double *copyData;
	gsl_fft_real_wavetable *real;
	gsl_fft_real_workspace *work;
	gsl_fft_halfcomplex_wavetable *hc;

	void updateSpectrum();
public:
	AtmosphereRandom(size_t nSamples, long seed = -1, double F0=20e-3);
	VecDoub getSimulatedData(VecDoub data);
	VecDoub getSimulatedAirmassGradient (VecDoub hEl, double avgTau, double tau0=0.056);
	~AtmosphereRandom();
	long getSeed();

};


#endif
