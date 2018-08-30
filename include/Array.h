#ifndef _ARRAY_H_
#define _ARRAY_H_

#include "nr3.h"
#include "Detector.h"
#include "Telescope.h"
#include "AnalParams.h"

///Array - a collection of detectors
/** The Array class is responsible for everything that is common
    across all the detectors and for all methods that work on
    the detector array as a whole (such as PCA cleaning, for 
    example).
**/
class Array
{
protected:
  ///analysis parameters
  AnalParams* ap;

  ///where to find the data
  const char* dataFile;        ///<netcdf data file for raw data
  const char* bolostatsFile;   ///<xml file with ancillary data

  ///detector stuff
  string* detectorIdsStr;      ///<list of the detectors (string)
  VecInt detectorIds;          ///<list of the detectors (int)
  VecInt detectorInd;          ///<indexes of the detectors with goodFlag=1
  int nDetectors;              ///<number of detectors in the array (goodflag=1)
  int nSamples;                ///<number of samples in the time stream

  ///pointing related
  NcToken refPix;              ///<the reference pixel name
  int refPixId;                ///<the ref pixel's detector id
  double boresightOffsets[2];  ///<Az/El delta between tel bs and ref pix
  double maxX;                 ///<maximum value of Ra or Az coordinate
  double minX;                 ///<minimum value of Ra or Az coordinate
  double maxY;                 ///<maximum value of Dec or El coordinate
  double minY;                 ///<minimum value of Dec or El coordinate

  size_t refBoloIndex;			///<Index of the Reference Bolometer (h2b2) in the detectorInd array

  void updateRefBoloIndex();

public:
  ///detectors and values
  Detector* detectors;         ///<the set of detectors
  int nFiltTerms;              ///<the number of digital filter terms
  VecDoub digFiltTerms;        ///<the dig. filt. coef. to lowpass the det data
  MatDoub eVectors;            ///<the saved PCA eigenvectors

  ///calibration
  double tau;                  ///<the array-averaged value of tau

  ///methods
  Array(AnalParams* ap);
  bool getLMTAddons();
  bool getASTEAddons();
  bool populate();
  bool updateDetectorIndices();
  int* getDetectorIndices();
  int getNDetectors();
  int getNSamples();
  double getMaxX();
  double getMinX();
  double getMaxY();
  double getMinY();
  double getAvgTau();
  bool findMinMaxXY();
  bool fakeFlaggedData(Telescope* telescope);
  bool fakeAtmData(bool addToKernel = true);
  size_t getRefBoloIndex();
  AnalParams *getAp();
  ~Array();
};

#endif
