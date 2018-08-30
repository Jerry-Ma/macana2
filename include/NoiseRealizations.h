#ifndef _NOISEREALIZATIONS_H_
#define _NOISEREALIZATIONS_H_

#include "nr3.h"
#include "AnalParams.h"
#include "Array.h"
#include "Map.h"
#include "Telescope.h"

///NoiseRealizations - noise realizations generated from Observations
/** A NoiseRealizations object is a set of noise realizations generated
    from the noise maps produced from a set of observations.  Because
    there can be many noise realizations (typically 100), most of this
    class is file-based and so io speed is particularly important.
**/
class NoiseRealizations
{
 protected:
  AnalParams* ap;            ///<pointer to our analysis parameters

  //map carryovers
  int nrows;                 ///<number of rows in each noise realization
  int ncols;                 ///<number of columns in each noise realization
  int nPixels;
  double pixelSize;          ///<pixel size in radians
  VecDoub rowCoordsPhys;     ///<row coordinates in tangential projection
  VecDoub colCoordsPhys;     ///<column coordinates in tangential projection
  MatDoub xCoordsAbs;        ///<matrix of sphere coordinates in ra/az
  MatDoub yCoordsAbs;        ///<matrix of sphere coordinates in dec/el
  VecDoub masterGrid;        ///<tangential point on sphere

 public:
  //the noise files
  int nNoiseFiles;           ///<number of noise realizations to make
  string noisePath;          ///<the full path to the directory of output files
  string* noiseFiles;        ///<the noise realizations output filenames

  //a noise map (really storage)
  Map* noise;                ///<storage for each noise realization

  //average noise histogram
  VecDoub histBins;          ///<lower bin edges for average noise histogram
  VecDoub histVals;          ///<corresponding bin histogram values

  //average noise psd
  VecDoub psd;               ///<the average psd of the noise realizations
  VecDoub psdFreq;           ///<the corresponding frequencies
  //2d psd
  MatDoub psd2d;
  MatDoub psd2dFreq;

  //average filtered rms    
  double averageFilteredRms; ///<the average rms of the filtered noise realizations

  NoiseRealizations(AnalParams* ap);
  bool generateNoiseRealizations(Coaddition* cmap);
  bool writeNoiseMapToNcdf(string thisNoiseFile, Map* noise);
  bool makeAverageHistogram(bool filtered);
  bool makeAveragePsd();
  bool calculateAverageFilteredRms(double cov);
  bool normalizeErrors(double cov);
  bool calcFilteredNoiseMapsHistogram(double cov);
  ~NoiseRealizations();
};

#endif
