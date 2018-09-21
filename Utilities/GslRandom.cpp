#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>
using namespace std;

#include "GslRandom.h"

//constructor
GslRandom::GslRandom()
{
  //create the generator
  r = gsl_rng_alloc(gsl_rng_mt19937);  //using the MT19937 generator

  //set the random number seed based on the current system time
  //which is seconds since a long time ago
  //GW - 17Nov15 - this turns out to have been a bad idea as 
  //when running in parallel mode this constructor could be 
  //called several times a second.  The following increases
  //the resolution to ns.
  timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  seed = ts.tv_sec*1e9 + ts.tv_nsec;
  gsl_rng_set(r,seed);
}

//constructor with a given seed
GslRandom::GslRandom(long seed)
{
  //create the generator
  r = gsl_rng_alloc(gsl_rng_mt19937);  //using the MT19937 generator
  gsl_rng_set(r,seed);
}


//return a gaussian deviate with variance=1
double GslRandom::gaussDeviate()
{
  double ret;
#pragma omp critical
  ret = gsl_ran_ugaussian(r);
  return ret;
}

//returns a uniform deviate on the range [a,b]
double GslRandom::uniformDeviate(double a, double b)
{
  double ret;
#pragma omp critical
  ret=gsl_ran_flat(r, a, b);
  return ret;
}

GslRandom::~GslRandom()
{
  gsl_rng_free(r);
}

  
