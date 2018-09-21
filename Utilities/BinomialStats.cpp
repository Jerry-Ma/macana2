#include <iostream>
#include <cmath>
#include <gsl/gsl_randist.h>
using namespace std;
#include "nr3.h"
#include "BinomialStats.h"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>

///BinomialStats constructor
/**Just initializes some members**/
BinomialStats::BinomialStats(int nTrials, double nSuccess, double t)
{
  n = nTrials;
  k = nSuccess;
  thresh = t;
}


//----------------------------- o ---------------------------------------

//calculates binomial confidence intervals following IDL
//aztec_calc_binomial_error.pro form and nomenclature
bool BinomialStats::calcIntervals()
{
  //zoom in on where the binomial matters
  double npt = 1000.0;
  if(n/10.0 > 1000.0){
    npt = double(ceil(n/10.));
  }

  vector<double> pt(npt);
  for(int i=0;i<npt;i++){
    pt[i] = i/npt;
  }

  vector<double> ptb(npt);
  for(int i=0;i<npt;i++){
    ptb[i] = gsl_cdf_binomial_Q(k, pt[i], n);
  }

  //a decrease would be a sign of a problem, so
  //anything below zero difference is bad but make
  //sure it's not just bit error
  vector<int> wtest;
  for(int i=0;i<npt-1;i++){
    if((ptb[i+1] - ptb[i]) < -1.0e-10){
      wtest.push_back(i);
    }
  }
  if(wtest.size() > 0){
    vector<int> wfull;
    for(int i=0;i<npt;i++){
      if(ptb[i] >= 1.0){
        wfull.push_back(i);
      }
    }
    if(wfull.size() == 0){
      cerr << "binoConfidence::calcIntervals(): severe underflow";
      cerr << "problem, attempting to compensate" << endl;
      for(int i=wtest[0]+1;i<npt;i++){
	ptb[i] = 1.0;
      }
    }
    else{
      for(int i=wfull[0];i<npt;i++) ptb[i] = 1.0;
    }
    wfull.clear();
  }
  wtest.clear();

  //get rid of useless low end
  vector<int> wz;
  double total = 0.0;
  for(int i=0;i<npt;i++){
    if(ptb[i] < (1.0/pow(n, 2))){
      wz.push_back(i);
    }
    total += ptb[i];
  }
  double lowpt = (wz.size() > 0) ? pt[wz.back()] : 0.0;
  //protect against total underflow
  if(total == 0) lowpt = 0.0;

  //get rid of useless high end
  //test against float, not double to avoid bit error
  vector<int> wo;
  for(int i=0;i<npt;i++){
    if(float(ptb[i]) == 1.0){
      wo.push_back(i);
    }
  }
  //if 1.0 at the first point, force it to the next
  //one or higher
  double highpt=-9999.;
  if(wo.size() == 0){
    highpt = 1.0;
  }
  else{
    if(wo[0] == 0){
      if(n >= 10) highpt = 10.0*npt/n*pt[1];
      if(n < 5) highpt = 1.0;
    }
    else{
      highpt = pt[wo[0]];
    }
  }
  if(highpt == -9999.){
    cerr << "BinomialStats: something terribly wrong with highpt." << endl;
    exit(1);
  }

  //array of gridded underlying probabilities of TRUE
  //must purposely ignore p=1, b/c that gives no options
  long np = (n/10.0 > long(5000)) ? long(np/10.0) : 5000;
  vector<double> p(np);
  for(int i=0;i<np;i++){
    p[i] = (double(i)/np)*(highpt - lowpt) + lowpt;
  }

  //get exact probability of k TRUEs for n tries at each p
  double binsize = p[1] - p[0];
  gsl_vector* f = gsl_vector_alloc(np);
  for(int i=0;i<np;i++){
    gsl_vector_set(f, i, binsize*n*(gsl_cdf_binomial_Q(k, p[i], n) -
				    gsl_cdf_binomial_Q(k+1, p[i],n)));
  }

  int wm = 0;
  double max = gsl_vector_get(f, 0);
  for(int i=1;i<np;i++){
    if(gsl_vector_get(f, i) > max){
      max = gsl_vector_get(f, i);
      wm = i;
    }
  }
  gsl_permutation* a = gsl_permutation_alloc(np);
  gsl_sort_vector_index(a, f);
  gsl_permutation_reverse(a);

  //should be gaussian-like, or at least continuously
  //downward from peak so we can sort them and do
  //easier search to find "fastest to thresh
  //confidence", where "thresh" might be 68% or something
  vector<float> sums(np);
  for(int i=0;i<np;i++){
    double total = 0.0;
    for(int j=0;j<i;j++){
      total += gsl_vector_get(f, gsl_permutation_get(a, j));
    }
    sums[i] = total;
  }

  //store the probabilities
  probU.resize(np);
  probK.resize(np);
  for(int i=0;i<np;i++){
    probU[i] = p[i];
    probK[i] = gsl_vector_get(f, i);
  }

  //find how far out must go to get "thresh" (eg 68%)
  //of total probability
  double min = abs(sums[0] - thresh);
  for(int i=1;i<np;i++){
    if(abs(sums[i] - thresh) < min){
      min = abs(sums[i] - thresh);
    }
  }
  vector<int> wv;
  for(int i=0;i<np;i++){
    if(abs(sums[i] - thresh) == min){
      wv.push_back(i);
    }
  }
  if((wv.size() == 0) || (wv.size() == unsigned(np))){
    cerr << "binoConfidence::calcIntervals(): Serious problem";
    cerr << "calculating binomial error range..." << endl;
    cerr << "binoConfidence::calcIntervals(): inputs were" << endl;
    cerr << "n = " << n << endl;
    cerr << "k = " << k << endl;
    cerr << "thresh = " << thresh << endl;
  }
  unsigned int minInd = gsl_permutation_get(a, 0);
  unsigned int maxInd = gsl_permutation_get(a, 0);
  for(int i=1;i<wv[0]+1;i++){
    if(gsl_permutation_get(a, i) < minInd){
      minInd = gsl_permutation_get(a, i);
    }
    if(gsl_permutation_get(a, i) > maxInd){
      maxInd = gsl_permutation_get(a, i);
    }
  }
  errorLow = float(wm - minInd)*binsize*n;
  errorHigh = float(maxInd - wm)*binsize*n;

 //done with f and a
  gsl_vector_free(f);
  gsl_permutation_free(a);

  return 1;
}


//----------------------------- o ---------------------------------------

float BinomialStats::getErrLow()
{
  return errorLow;
}


//----------------------------- o ---------------------------------------

float BinomialStats::getErrHigh()
{
  return errorHigh;
}


//----------------------------- o ---------------------------------------

double BinomialStats::getProbU(int i)
{
  return probU[i];
}


//----------------------------- o ---------------------------------------

double BinomialStats::getProbK(int i)
{
  return probK[i];
}


//----------------------------- o ---------------------------------------

BinomialStats::~BinomialStats()
{

}
