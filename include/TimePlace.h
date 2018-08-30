#ifndef _TIMEPLACE_H_
#define _TIMEPLACE_H_

#include "AnalParams.h"

///TimePlace - keeps track of where and when.
/** The TimePlace object simply keeps the time and telescope location.
    \todo At some point we should add the UT correction predictions so
    that the ASTE ra/dec to az/el (and visa versa) transformations are
    more accurate.
 **/
class TimePlace
{
protected:
  //analysis parameters
  AnalParams* ap;              ///<pointer to our analysis parameters

  const char* dataFile;        ///<netcdf data file for raw data
  int nSamples;                ///<number of data samples

  //protected methods
  bool getLMTTiming();
  bool getASTETiming();
  bool getTimePlaceValues(const char* sigName, double *tData);
  bool signalCondition0();
  bool signalCondition1();
  bool fixRollover();
  bool alignWithDetectors();
  bool findUniqueTelUtc();


public:
  VecDoub telUtc;               ///<UTC as reported by telescope (IN HOURS!)
  VecDoub telUtcUnique;         ///<unique values of telUtc
  VecInt locUnique;             ///<locations of unique values of telUtc
  int nUnique;                  ///<number of unique values in telUtc
  VecDoub detUtc;               ///<UTC as reported by detector backend
  VecDoub lst;                  ///<LST
  double  longitude;            ///<telescope longitude [rad]
  double  latitude;             ///<telescope latitude [rad]
  double  elevation;            ///<telescope elevation [m]
  double  timeOffset;           ///<offset between tel and backend time [s]
  double  utDate;               ///<obs date in decimal years (this is epoch)
  double  julDate;              ///<julian date calculated from utDate
  int     year;                 ///<year of data file
  int     month;                ///<month of data file
  int     day;                  ///<day of data file

  //constructor
  TimePlace(AnalParams* ap);
  ~TimePlace();
};

#endif
