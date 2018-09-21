#ifndef _DETECTOR_H_
#define _DETECTOR_H_

#include <netcdfcpp.h>
#include "Telescope.h"
#include "AnalParams.h"
#include "Source.h"

///Detector - the base element of an array.
/**This class contains everything that a detector
   can do on its own.  It includes some elements of
   despiking, calibration (except for calculating the
   mean extinction, nonlinearity correction, etc.
**/
class Detector
{
protected:
  //items common to all detectors in analysis
  const char* dataFile;              ///<netcdf data file for raw data
  const char* bolostatsFile;         ///<text file with beammap data
  NcToken ref_pix;                   ///<the reference pixel
  int rawSamplerate;                 ///<sample rate for raw data
  AnalParams* ap;                    ///<pointer to our AnalParams driver

  //identification and status
  string name;                      ///<detector name in data file
  int  id;                           ///<a more convenient detector ID
  int ncdfLocation;                 ///<variable number in dataFile

  //psf attributes
  double beamSigAz;                  ///<psf std deviation in azimuth
  double beamSigEl;                  ///<psf std deviation in el
  double fcf;                        ///<flux conversion factor

  //time stream characteristics
  int nSamples;                     ///<number of timestream samples
  double samplerate;                 ///<current sample rate

  //calibration
  double dc;                         ///<bolo dc offset and slope
  double dc2tau[3];                  ///<conversion [offset,slope]
  double dc2tauErr[3];               ///<error in conversion params
  double dc2responsivity[2];         ///<conversion [offset,slope]
  double dc2responsivityErr[2];      ///<error in conversion params
  double voltsToJanskies;            ///<volts to janskies convertion
  double extinction;                 ///<inferred atmosphere extinction
  double fecGain;                    ///<electronics gain
  double sensitivity;                ///<sensitivity from beammaps mJy-rts

  //despiking
  int despikeWindow;                ///<size of region around spike to discard

  //pointing
  double maxX;                       ///<maximum value of Ra or Az coordinate
  double minX;                       ///<minimum value of Ra or Az coordinate
  double maxY;                       ///<maximum value of Dec or El coordinate
  double minY;                       ///<minimum value of Dec or El coordinate

  //analysis status
  bool isDespiked;  
  bool isCleaned;  
  bool isLowpassed;  
  bool isCalibrated;  
  bool isDownsampled;  
  bool isPointingGenerated;

  //private methods
  bool estimateResponsivity();
  double cmdToGain(int cmd);


public:
  bool goodFlag;                     ///<is this bolometer worthwhile to include
  double responsivity;               ///<responsivity (volts per watt)

  //time stream characteristics
  VecDoub hValues;                   ///<pointer to det values on cpu host
  VecBool hSampleFlags;              ///<pointer to sampleflags on host
  VecDoub hRa;                       ///<pointer to ra values on host
  VecDoub hDec;                      ///<pointer to dec values on host
  VecDoub hKernel;                   ///<pointer to the kernel signal values
  VecDoub dValues;                   ///<pointer to det values on gpu device
  VecBool dSampleFlags;              ///<pointer to sampleflags on device
  VecDoub dRa;                       ///<pointer to ra values on device
  VecDoub dDec;                      ///<pointer to dec values on device
  VecDoub azElRa;
  VecDoub azElDec;
  VecDoub azElRaPhys;
  VecDoub azElDecPhys;
  VecDoub atmTemplate;

  //calibration
  double estimatedTau;               ///<tau based on bolodc conversion
  VecDoub scanWeight;                ///<variance^-1 of the samples in a scan

  double azOffset;                   ///<az offset from reference pixel
  double elOffset;                   ///<el offset from reference pixel
  //public methods
  Detector();
  void initialize(AnalParams* analParams, const char* dId);
  int getNSamples();
  double getSamplerate();
  bool getBoloValues(double *bData);
  bool despike(double nSigmaSpikes);
  bool despike(int nSigmaSpikes, int gpuId);
  bool lowpass(double* digFilterTerms, int nTerms);
  bool lowpass(int gpuId);
  bool downsample(double desiredSamplerate);
  bool estimateTau();
  bool estimateExtinction(double tau);
  bool calibrate();
  bool calculateScanWeight(Telescope* tel);
  void setFcf(double fcf);
  double calculateSensitivity(Telescope* tel);
  void calibrateSensitivity();
  bool getPointing(Telescope* tel, TimePlace* tp, Source* source);
  bool getAzElPointing (Telescope *tel);
  bool makeKernelTimestream(Telescope* tel);
  void addGaussian(MatDoub &params, int detNumber);
  bool setAtmTemplate(VecDoub temp);
  void throwXmlError(string p);
  double getMaxX();
  double getMinX();
  double getMaxY();
  double getMinY();
  double getFcf();
  double getSensitivity();
  double getCalibrationFactor();
  double getExtinction();
  double getFecGain();
  int getId();
  string getName();
  ~Detector();
};

#endif /* _DETECTOR_H_ */
