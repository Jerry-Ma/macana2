#ifndef _ANALPARAMS_H_
#define _ANALPARAMS_H_

#include <string>
#include "nr3.h"

using namespace std;
///AnalParams - All of the analysis parameters and directions
/** The AnalParams class runs the show.  It contains all of the
    information needed to do the entire analysis, including
    input/output file names and paths, analysis parameters for
    mapmaking, etc.  All switches in the analysis should be
    included here.
    \todo Make a menu of analysis options to choose from and
    propagate those choices into macana.cpp.
**/

class AnalParams
{
 protected:
  ///Analysis Steps
  bool mapIndividualObservations;
  bool coaddObservations;
  bool fitCoadditionToGaussian;
  bool produceNoiseMaps;
  bool applyWienerFilter;
  bool postReductionAnalysis;

  ///Observation files and paths
  string apXml;                     ///<the xml file for the analysis
  int nFiles;                       ///<the number of files to be reduced
  string rawDataPath;               ///<the path for the raw data
  string bsPath;                    ///<the path for the corresponding bstats
  string mapPath;                   ///<the path for the output map nc files
  string* fileList;                 ///<the list of files to be reduced
  string* bstatList;                ///<the list of bolostat files
  string* mapFileList;              ///<the output map file names
  MatDoub bsOffsetList;             ///<the corresponding boresight offsets

  const char* dataFile;             ///<the current netcdf file to be reduced
  const char* bolostatsFile;        ///<the current bolostats file
  const char* mapFile;              ///<the current output map file name
  string sourceName;                ///<the source name from the data file    

  ///coaddition
  string coaddOutPath;              ///<path of output coadded maps nc file
  string coaddOutFile;              ///<filename of output coadded map nc file

  ///noise realizations
  int nNoiseMapsPerObs;             ///<number of noise maps for each obs
  int nRealizations;                 ///<number of noise realizto generate
  string noisePath;                  ///<path for the noise realizations
  string avgNoisePsdFile;            ///<filename for avg noise psd nc file
  string avgNoiseHistFile;           ///<filename for avg noise hist nc file

  ///coverage threshold for norms
  double coverageThreshold;         ///<coverage cut for noise normalizations

  ///despiking
  double despikeSigma;              ///<despiking threshold

  ///weights
  bool approximateWeights;

  ///lowpassing and subsampling
  double lowpassFilterKnee;
  double subsampleFrequency;

  ///scan definition
  double timeOffset;
  double timeChunk;
  bool cutTurnArounds;
  int cutSamplesAtEndOfScans;

  //azel map?
  int azelMap;
  
  ///pointing
  double bsOffset[2];
  double* masterGridJ2000;
  double* masterGridJ2000_init;
  double pixelSize;

  ///Observatory specific items
  string observatory;
  string timeVarName;

  ///Wiener Filter options
  bool wienerFilter;                  ///wiener filter the maps?
  bool gaussianTemplate;              ///use a gaussian instead of kernel?
  double gaussianTemplateFWHM;        ///the gaussian template's fwhm [radians]
  bool lowpassOnly;                   ///only lowpass the maps?
  bool highpassOnly;                  ///only highpass the maps?
  bool normalizeErrors;               ///normalize the weight matrix?
 
  ///Threaded operation parameters
  int nThreads;

  bool saveTimeStreams;

  ///Source finding parameters and switches
  bool findSources;                    ///switch to turn on source finding
  double beamSize;                     ///beam size in radians
  double covCut;                       ///fraction of max coverage considered
  double snglSourceWin;                ///RADIUS masked out around source (rad)
  double sourceSigma;                  ///"number of sigmas" considered a source
  bool negativeToo;                    ///also identify negative parts of map
  bool mapNegative;                    ///search negative of signal map
  bool sfCentroidSources;              ///centroid sources after finding?
  bool sfFitGaussians;                 ///fit Gaussians to sources after finding
  int psSideLength;                    ///# pixels on source postage stamp side

  ///Completeness calculation parameters and switches
  bool calcCompleteness;               ///switch to turn on completeness calc
  int nFluxBins;                       ///# of flux bins
  int nSynthSources;                   ///# of synthetic sources in flux bins
  double minFlux;                      ///min flux bin used in mJy
  double maxFlux;                      ///max flux bin used in mJy
  double recovS2N;                     ///S/N thresh for recovering synth source
  double recovRadius;                  ///dist. from synth source for recovery



 public:

  ///cleaning PCA parameters
  double cutStd;
  int neigToCut;
  int icaMode;

  ///cleanning Cuttingham method
  double cleanPixelSize;		///Pix size for pointing mat in arcsec
  int order;					///Spline order
  double cleanStripe;			///Apply large scale resid sustraction
  string stripeMethod;
  double controlChunk;			///Time sep of Bspline control points
  double resample;
  double mask;
  ///cleaning High Order Atm Template
  int tOrder;

  AnalParams(string apXmlFile);
  AnalParams(AnalParams *ap);
  bool setDataFile(int index);
  bool determineObservatory();
  void throwXmlError(string p);
  string getObservatory();
  string getTimeVarName();
  const char* getDataFile();
  const char* getBolostatsFile();
  bool getMapIndividualObservations();
  bool getCoaddObservations();
  bool getFitCoadditionToGaussian();
  bool getProduceNoiseMaps();
  bool getApplyWienerFilter();
  bool getWienerFilter();
  double getDespikeSigma();
  bool getApproximateWeights();
  int getAzelMap();
  double getLowpassFilterKnee();
  double getSubsampleFrequency();
  double getTimeOffset();
  double getTimeChunk();
  bool getCutTurnArounds();
  int getCutSamplesAtEndOfScans();
  double getCutStd();
  int getNeigToCut();
  double getCleanPixelSize();
  void setCleanPixelSize(double pixSize);
  int getOrder();
  double getCleanStripe();
  double getControlChunk();
  void setControlChunk(double control);
  int getNThreads();
  bool getSaveTimestreams();
  double* getBsOffset();
  double* getMasterGridJ2000();
  bool   setMasterGridJ2000(double ra, double dec);
  double getPixelSize();
  bool getGaussianTemplate();
  double getGaussianTemplateFWHM();
  bool getLowpassOnly();
  bool getHighpassOnly();
  bool getNormalizeErrors();
  double getCoverageThreshold();

  bool getPostReductionAnalysis();
  double getBeamSize();
  double getCovCut();
  double getSnglSourceWin();
  double getSourceSigma();
  bool getNegativeToo();
  bool getMapNegative();
  bool getSFCentroidSources();
  bool getSFFitGaussians();
  int getPsSideLength();
  int getNFluxBins();
  int getNSynthSources();
  double getMinFlux();
  double getMaxFlux();
  double getRecovS2N();
  double getRecovRadius();
  int getResample();
  double getFcut();
  double getFalpha();


  void setOrder (int Order);
  void setNeig2cut (int n2cut);

  int getTOrder ();

  int getNFiles();
  string getMapFile();
  string getMapFileList(int i);
  string getSourceName();
  bool setSourceName(string name);
  string getCoaddOutFile();
  string getNoisePath();
  int getNRealizations();
  string getAvgNoiseHistFile();
  string getAvgNoisePsdFile();
  int getNNoiseMapsPerObs();
  bool getFindSources();
  bool getCalcCompleteness();
  ~AnalParams();
};

#endif
