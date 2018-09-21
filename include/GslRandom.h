#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <stdio.h>
#include <gsl/gsl_rng.h>

//uses the gsl random number generator suite
class GslRandom
{

 private:
  long seed;                //random number generator seed
  gsl_rng * r;              //the generator

 public:
  GslRandom();
  GslRandom(long seed);
  double gaussDeviate();    //returns a gaussian deviate with variance 1
  double uniformDeviate(double a, double b);  //uniform deviate on [a,b]
  ~GslRandom();
};


#endif
