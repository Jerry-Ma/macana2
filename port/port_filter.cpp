//Filter Include
#include "port_filter.h"

namespace filter {

//Filter for Macana
bool digiFilter(const double fLow, const double fHigh, const double aGibbs,
                   const int nTerms, VecDoub &coefOut)
{
  // computes Kaiser weights W(N,K) for digital filters.
  // W = COEF = returned array of Kaiser weights
  // N = value of N in W(N,K), i.e. number of terms
  // A = Size of gibbs phenomenon wiggles in -DB.

  double alpha;
  if(aGibbs <= 21.){
    alpha = 0.;
  }
  if(aGibbs >= 50){
    alpha = 0.1102*(aGibbs-8.7);
  }
  if(aGibbs > 21. && aGibbs < 50){
    alpha = 0.5842*pow(aGibbs-21.,0.4) + 0.07886*(aGibbs-21.);
  }

  //Band stop?
  double fStop = (fHigh < fLow) ? 1. : 0.;

  //Arg
  VecDoub arg(nTerms);
  for(double i=0;i<nTerms;i++) arg[i]=(i+1.)/nTerms;

  //Coef
  VecDoub coef(nTerms);
  for(int i=0;i<nTerms;i++){
    coef[i] = gsl_sf_bessel_I0(alpha*sqrt(1.-pow(arg[i],2))) /
      gsl_sf_bessel_I0(alpha);
  }

  //t
  VecDoub t(nTerms);
  for(double i=0;i<nTerms;i++) t[i]=(i+1.)*M_PI;

  //here we go
  for(int i=0;i<nTerms;i++){
    coef[i] *= (sin(t[i]*fHigh)-sin(t[i]*fLow))/t[i];
  }

  //build the uncentered version
  for(int i=0;i<nTerms;i++) coefOut[i] = coef[nTerms-i-1];
  coefOut[nTerms] = fHigh-fLow-fStop;
  for(int i=0;i<nTerms;i++) coefOut[i+nTerms+1] = coef[i];

  //renormalize
  double area=0;
  //for(int i=0;i<2*nTerms+1;i++) area += coefOut[i];
  //for(int i=0;i<2*nTerms+1;i++) coefOut[i] /= area;

  return 1;

}


/*----------------------------------------------------------------------------------------------*/

//Function for Macana Filter Convolution

void MacanaFilter(VecDoub &hValues, const VecDoub &digFiltTerms, VecDoub &storage,
                  const int nSamples, const int nTerms, const int nCoef){

    double tmp;
    int th=nSamples-nTerms;
    int tt;

    //Convolve it!
    for(int t=nTerms;t<th;t++){
        tmp=0.;
        tt = t-nTerms;
        for(int i = 0;i<nCoef;i++){
            tmp += hValues[tt+i]*digFiltTerms[i];
        }
        storage[t] = tmp;
    }

    //Replace old values
    for(int i=nTerms;i<th;i++) {
      hValues[i] = storage[i];
    }

}

} //namespaace

