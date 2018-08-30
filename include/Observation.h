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
  AnalParams* ap;              ///<pointer to our analysis parameters

  //carried over from maps
  int nrows;                   ///<number of rows in the image
  int ncols;                   ///<number of columns in the image
  int nPixels;                 ///<total number of pixels in the image
 // VecDoub rowCoordsPhys;       ///<row coordinates in tangential projection
 // VecDoub colCoordsPhys;       ///<column coordinates in tangential projection

  //file output
  string ncdfFile;             ///<the output ncdf files containing everything
  bool saveTimestreams;

 public:
  //the maps
  VecDoub rowCoordsPhys;       ///<row coordinates in tangential projection
  VecDoub colCoordsPhys;       ///<column coordinates in tangential projection

  Map* signal;                 ///<the signal map values
  Map* weight;                 ///<the weight map (1/error^2)
  Map* kernel;                 ///<kernel map
  Map* inttime;                ///<inttime map (in seconds)
  MatDoub noiseMaps;           ///<array of noise map for each observaiton
  Map* atmTemplate;

  //absolute coordinates
  double pixelSize;           ///<map pixel size in radians
  VecDoub masterGrid;         ///<map center coordinates (sky tangent point)
  MatDoub xCoordsAbs;         ///<matrix of absolute on-sky x-coordinates
  MatDoub yCoordsAbs;         ///<matrix of absolute on-sky y-coordinates

  Array * array;
  Telescope *tel;

  //methods
  Observation(AnalParams* ap);
  bool generateMaps(Array* a, Telescope* tel);
  bool histogramSignal();
  bool histogramSignal(double cc);
  bool histogramSignal(int nbins, double cc);
  bool signalMapPsd();
  bool signalMapPsd(double cc);
  bool writeObservationToNcdf(string ncdfFilename);
  bool writeObservationToFits(string fitsFilename);
  bool writeObservationToBmp(string bmpFilename);
  ~Observation();  
};


#endif

