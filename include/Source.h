#ifndef _SOURCE_H_
#define _SOURCE_H_

#include "AnalParams.h"
#include "TimePlace.h"

///Source - what was observed
/**The Source class contains all the information about the source that was
   observed.  Of particular importance are the source coordinates.
   Note that since the LMT and ASTE dealt with sources somewhat 
   differently, different methods had to be adopted to deal with each
   of them.
**/
class Source
{
protected:
  AnalParams* ap;

  const char* dataFile;             ///<netcdf data file for raw data
  int nSamples;                     ///<number of data samples
  TimePlace* timePlace;             ///<access to the time

  //protected methods
  bool getLMTSourceData();
  bool getASTESourceData();
  bool parseRaDecString(string ra, double* t1, double* t2, double* t3);
  bool getSourceValues(const char* sigName, double *tData);
  bool signalCondition0();
  bool signalCondition1();
  bool fixRollover();
  bool alignWithDetectors();
  bool testNovas();

public:
  string name;                     ///<source name

  //source pointing signals
  VecDoub hSourceAz;               ///<source azimuth vs time (radians)
  VecDoub hSourceEl;               ///<source elevation vs time (radians)
  VecDoub hSourceRa;               ///<source ra vs time (radians)
  VecDoub hSourceDec;              ///<source dec vs time (radians)

  //constructor
  Source(AnalParams* ap, TimePlace* tP);
  ~Source();
};

#endif
