#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <omp.h>
#include <random>
using namespace std;

#include "GslRandom.h"

//constructor
GslRandom::GslRandom()
{
  //create the generator
  r = gsl_rng_alloc(gsl_rng_mt19937);  //using the MT19937 generator

  //gsl_rng * g = gsl_rng_alloc (gsl_rng_taus);

  //set the random number seed based on the current system time
  //which is seconds since a long time ago
//  time_t t;
//  time(&t);
//  seed = long(t);
  long tid;
#if defined(_OPENMP)
	  tid = omp_get_thread_num() + 1;
#else
	  tid = 1;
#endif
//  seed += tid;
  //set the generator
  seed = std::random_device()();
  gsl_rng_set(r,seed);
}

void GslRandom::reseed(long seed){
	this->seed = seed;
	gsl_rng_set (r,seed);
}


long GslRandom::getSeed(){
	return seed;
}
//return a gaussian deviate with variance=1
double GslRandom::gaussDeviate()
{
  return gsl_ran_ugaussian(r);
}

//returns a uniform deviate on the range [a,b]
double GslRandom::uniformDeviate(double a, double b)
{
  return gsl_ran_flat(r, a, b);
}

GslRandom::~GslRandom()
{
  gsl_rng_free(r);
}

  
