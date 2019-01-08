#ifndef _OBSERVATION_H_
#define _OBSERVATION_H_

#include "nr3.h"
#include "AnalParams.h"
#include "Array.h"
#include "Map.h"
#include "Telescope.h"

///Observation - a single observation (mapping) of one area of sky.
/** The Observation class manages the mapmaking for a single observation
    of the sky (that is, a single netcdf file containing raw data).
    The entities of this class are analagous to what comes out of
    aztec_quickmap.pro in the idl utilities.
**/
class Observation
{
 protected:
  AnalParams* ap = nullptr;              ///<pointer to our analysis parameters

  //carried over from maps
 //int nrows;                   ///<number of rows in the image
 //int ncols;                   ///<number of columns in the image
 //int nPixels;                 ///<total number of pixels in the image
 // VecDoub rowCoordsPhys;       ///<row coordinates in tangential projection
 // VecDoub colCoordsPhys;       ///<column coordinates in tangential projection

  //file output
  string ncdfFile;             ///<the output ncdf files containing everything
  bool saveTimestreams = false;

 public:
  //dimensions
  int nrows = 0;                   ///<number of rows in the image
  int ncols = 0;                   ///<number of columns in the image
  int npixels = 0;                 ///<total number of pixels in the image

  //the maps
  VecDoub rowCoordsPhys;       ///<row coordinates in tangential projection
  VecDoub colCoordsPhys;       ///<column coordinates in tangential projection

  Map* atmTemplate = nullptr;
  Map* signal = nullptr;                 ///<the signal map values
  Map* weight = nullptr;                 ///<the weight map (1/error^2)
  Map* kernel = nullptr;                 ///<kernel map
  Map* inttime = nullptr;                ///<inttime map (in seconds)
  MatDoub noiseMaps;           ///<array of noise map for each observaiton

  MatDoub beammapSignal;       ///<array of detector beammap signal values
  MatDoub beammapWeight;       ///<array of detector beammap weight values
  MatDoub beammapSigma;      ///<array of detector beammap inttime values
  MatDoub beammapIntTime;      ///<array of detector beammap inttime values

  //absolute coordinates
  double pixelSize=0;           ///<map pixel size in radians
  VecDoub masterGrid;         ///<map center coordinates (sky tangent point)
  MatDoub xCoordsAbs;         ///<matrix of absolute on-sky x-coordinates
  MatDoub yCoordsAbs;         ///<matrix of absolute on-sky y-coordinates

  Array * array = nullptr;
  Telescope *tel = nullptr;

  //methods
  Observation(AnalParams* ap);
  MatDoub calculateWeights(Array* a, Telescope* tel);
  bool generateBeammaps(Array* a, Telescope* tel);
  bool generateMaps(Array* a, Telescope* tel);
  bool histogramSignal();
  bool histogramSignal(double cc);
  bool histogramSignal(int nbins, double cc);
  void mapGenerationPrep(Array* a);
  bool signalMapPsd();
  bool signalMapPsd(double cc);
  bool writeBeammapsToFits(string filename);
  bool writeBeammapsToNcdf(string filename);
  bool writeBeammapsToBmp(string bmpFilename, int nDetectors);
  bool writeFitParamsToNcdf(string ncdfFilename, MatDoub &fitParams, Array* a);
  bool writeSensToNcdf(string ncdfFilename, Array* a);
  bool writeObservationToNcdf(string ncdfFilename);
  bool writeObservationToFits(string fitsFilename);
  bool writeObservationToBmp(string bmpFilename);
  ~Observation();  
};


#endif

