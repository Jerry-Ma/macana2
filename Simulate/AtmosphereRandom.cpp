#include <cstring>
#include <cmath>
#include <gsl/gsl_sf_trig.h>

#include "AtmosphereRandom.h"
#include "GslRandom.h"
#include "vector_utilities.h"

AtmosphereRandom::AtmosphereRandom(size_t nSamples, long seed, double F0){
	this->nSamples = nSamples;
	this->randomSequence = new double[nSamples];
	this->spectrum = new double[nSamples];
	this->copyData = new double[nSamples];
	GslRandom r;

	if (seed >0){
		this->seed = seed;
		r.reseed(this->seed);
	}else
		this->seed = r.getSeed();
	for (size_t i = 0; i<nSamples; i++)
		randomSequence[i] = r.gaussDeviate();
	real = gsl_fft_real_wavetable_alloc (nSamples);
	work = gsl_fft_real_workspace_alloc (nSamples);
	hc = gsl_fft_halfcomplex_wavetable_alloc (nSamples);
	gsl_fft_real_transform(randomSequence,1, nSamples,real,work);
	//Gradient of atmosphere
	this->fluxGradient = F0* (1.0 + 0.1*r.gaussDeviate());
	cout<<"Initialized AtmosphereRandom Object with zenith gradient flux: "<<this->fluxGradient<<endl;
}

VecDoub AtmosphereRandom::getSimulatedData(VecDoub data){
	VecDoub result (nSamples);
	double dmean = mean(data);
	double dstd = stddev(data);
	double rmean, rstd;
	std::memcpy (copyData, &data[0], sizeof(double)*nSamples);
	updateSpectrum();

	for (size_t i=0; i<nSamples; i++)
		result[i]= randomSequence[i]*spectrum[i];
	gsl_fft_halfcomplex_inverse (&result[0], 1, nSamples, hc, work);
	rmean= mean(result);
	rstd = stddev(result);
	for (size_t i=0; i<nSamples; i++)
		result[i]= (result[i]-rmean)*dstd/rstd +dmean;
	return result;
}


void AtmosphereRandom::updateSpectrum(){
	gsl_fft_real_transform(copyData,1, nSamples,real,work);
	size_t limit = nSamples -1;
	spectrum[0]= copyData[0];
	if (nSamples % 2==0){
		limit-=1;
		spectrum[nSamples-1]= copyData[nSamples-1];
	}
	for (size_t j=1; j<limit; (j+=2) ){
		spectrum[j] = spectrum[j+1] =  sqrt(pow(copyData[j],2.0) + pow(copyData[j+1],2.0));
	}

}

long AtmosphereRandom::getSeed(){
	return seed;
}

VecDoub AtmosphereRandom::getSimulatedAirmassGradient (VecDoub hEl, double avgTau, double tau0){
	VecDoub grad (nSamples,0.0);
	for (size_t i=0; i<nSamples; i++)
		grad[i] = fluxGradient*exp(avgTau/tau0)/gsl_sf_cos((90.0-hEl[i])*!M_PI/180.);
	return grad;
}

AtmosphereRandom::~AtmosphereRandom(){
	delete [] randomSequence;
	delete [] copyData;
	delete [] spectrum;
	gsl_fft_real_wavetable_free(real);
	gsl_fft_halfcomplex_wavetable_free (hc);
	gsl_fft_real_workspace_free(work);
}


