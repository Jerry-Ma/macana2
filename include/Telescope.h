#ifndef _TELESCOPE_H_
#define _TELESCOPE_H_

#include "AnalParams.h"
#include "Source.h"
#include "TimePlace.h"


///Telescope - a class describing how the telescope boresight moves
/** The Telescope class describes the movement of the
    telescope.
**/
class Telescope
{
protected:
  //analysis parameters
  AnalParams* ap;              ///<pointer to our analysis parameters

  const char* dataFile;        ///<netcdf data file for raw data
  int nSamples;                ///<number of data samples
  int nScans;                  ///<number of scans in observation

  //need access to the time
  TimePlace* timePlace;       ///<pointer to provide access to the time
  double utDate;              ///<the utDate storage only

  //some header info
  string projectID;           ///<The project ID assigned by the observatory
  double M2XReq;              ///<The commanded x offset of M2
  double M2YReq;              ///<The commanded y offset of M2
  double M2ZReq;              ///<The commanded z offset of M2
  double M1ZernikeC0;         ///<The 0th zernike coefficient of M1
  double M1ZernikeC1;         ///<The 1st zernike coefficient of M1
  double M1ZernikeC2;         ///<The 2nd zernike coefficient of M1
  double azUserOff;           ///<Azimuth user offset
  double elUserOff;           ///<Elevation user offset

  // obsNum
  int obsNum;
  
  //protected methods
  bool getTelescopeValues(const char* sigName, double *tData);
  bool signalCondition0();
  bool signalCondition1();
  bool fixRollover();
  bool alignWithDetectors();

public:
  //telescope pointing signals on host CPU
  VecDoub hTelRa;              ///<telescope boresight ra (radians)
  VecDoub hTelDec;             ///<telescope boresight dec (radians)
  VecDoub hTelRaPhys;          ///<tel ra in tangential projection (rad)
  VecDoub hTelDecPhys;         ///<tel dec in tangential projection (rad)
  VecDoub hTelAzPhys;          ///<tel az in tangential projection (rad)
  VecDoub hTelElPhys;          ///<tel el in tangential projection (rad)
  VecDoub hTelAzAct;           ///<encoder-based azimuth (rad)
  VecDoub hTelElAct;           ///<encoder-based elevation (rad)
  VecDoub hTelAzDes;           ///<target azimuth (rad)
  VecDoub hTelElDes;           ///<target elivation (rad)
  VecDoub hTelAzCor;           ///<various azimuth corrections (rad)
  VecDoub hTelElCor;           ///<various elevation corrections (rad)

  //Information about scanning
  MatInt scanIndex;            ///<arrays of scan indices [start, finish]

  //parallactic angle of each pointing sample
  VecDoub paraAngle;           ///<angle between curent azimuth and J2000 Ra

  //public methods
  Telescope(AnalParams* ap, TimePlace* timePlace, Source* source);
  bool getLMTPointing();
  bool getASTEPointing();
  bool defineScans(AnalParams* ap, double timeChunk);
  void checkScanLengths(VecBool obsflags, int samplerate);
  bool makeTelRaDec(TimePlace* tp, Source* source);
  bool absToPhysEqPointing();
  bool absToPhysHorPointing(Source* source);
  bool calcParallacticAngle();
  bool getM2Offsets();
  bool getM1Zernike();
  bool getObsNumFromFile();
  bool getProjectIDFromFile();
  double getM2XReq();
  double getM2YReq();
  double getM2ZReq();
  double getM1ZernikeC0();
  double getM1ZernikeC1();
  double getM1ZernikeC2();
  int getObsNum();
  double getUTDate();
  double getAzUserOff();
  double getElUserOff();
  string getProjectID();
  ~Telescope();
};

#endif
