#ifndef _BinomialStats_h_
#define _BinomialStats_h_

#include "nr3.h"

class BinomialStats
{
 protected:
  int n;
  double k;
  double thresh;

  //confidence intervals
  double errorLow;
  double errorHigh;

  //probabilties
  VecDoub probU;
  VecDoub probK;

 public:
  //methods
  BinomialStats(int nTrials, double nSuccess, double t);
  bool calcIntervals();
  float getErrLow();
  float getErrHigh();
  double getProbU(int i);
  double getProbK(int i);
  ~BinomialStats();
};

#endif
