#ifndef _MAP_H_
#define _MAP_H_
#include "nr3.h"
#include "AnalParams.h"
#include "Array.h"
#include "Telescope.h"

///Map - the base class of imaging products
/** A Map is everything associated with an image of a field.  Map
    objects are used in the Observation, Coaddition, and
    NoiseRealization classes.  To save memory I have not included the
    absolute sky matrices (ra and dec) in this class since they are
    shared amongst many maps.  
    \todo The weight map is common amongst many maps and so
    as-implemented with a direct copy this is very wastful in terms of
    memory.  Figure out how to pass this weight map by reference
    isntead.
**/
class Map
{
 protected:
  //counting
  int nPixels;                ///<total number of pixels in map
  int nrows;                  ///<number of rows (changing ra or az)
  int ncols;                  ///<number of columns (changing dec or el)

  //coordinates
  double pixelSize;           ///<pixel size in radians
  VecDoub rowCoordsPhys;      ///<vector of row delta-source coordinates
  VecDoub colCoordsPhys;      ///<vector of column delta-source coordinates
  bool isAzEl;                ///<is the coordinate system Az/El?
  
  //for coverage cut
  double coverageCut;         ///<the coverage cut as defined in idl utils
  double weightCut;           ///<coverage threshold in weight
  VecInt cutXRange;           ///<vector [low column index, high column index]
  VecInt cutYRange;           ///<vector [low row index, high row index]

 public:
  //psd
  VecDoub psd;                ///<the psd of the map
  VecDoub psdFreq;            ///<the corresponding frequencies
  //2dpsd
  MatDoub psd2d;
  MatDoub psd2dFreq;

  //histogram
  VecDoub histBins;           ///<lower bin edges for signal histogram
  VecDoub histVals;           ///<corresponding bin histogram values

  //storage
  string mapName;             ///<name associated with this map in ncdf file
  string mapFile;             ///<nc file containing the map

  //the map
  MatDoub image;              ///<the actual map
  MatDoub weight;             ///<the corresponding weight map
  MatBool coverageBool;       ///<boolean map indicating good coverage region

  //methods
  Map(string mapName, int nr, int nc, double pixsz, 
      MatDoub &wt, VecDoub &rcp, VecDoub &ccp);
  bool calcMapPsd(double covCut);
  bool calcMapHistogram(int nbins, double covCut);
  bool calcMapHistogram();
  double fitToGaussianMasked(VecDoub &params, VecInt &fixme, VecDoub &fixVals, double *iguess=nullptr, int deg=40);
  double fitToGaussian(VecDoub &params, VecInt &fixme, VecDoub &fixVals, double *iguess=nullptr, int deg=40);
  double fitToGaussian(VecDoub &params);
  double fitToGaussian(VecDoub &params, int deg);
  double fitToGaussian();
  bool writeMapToNcdf();
  bool writeMapsToNcdf(string mapFilename);
  int getNrows();
  int getNcols();
  int getNPixels();
  double getPixelSize();
  double getRowCoordsPhys(int i);
  double getColCoordsPhys(int i);
  double getCoverageCut();
  bool setCoverageCut(double cc);
  bool setCoverageCutRanges();
  bool findWeightThresh();
  bool makeCovBoolMap();
  bool raDecPhysToIndex(double ra, double dec, int* irow, int* icol);
  bool raDecAbsToIndex(double ra, double dec, int* irow, int* icol);
  bool indexToRaDecPhys(int index, double* ra, double* dec);
  bool indexToRaDecAbs(int index, double* ra, double* dec);
  double getXCoordsAbs(int i, int j);
  double getYCoordsAbs(int i, int j);
  int getCutXRangeLow();
  int getCutXRangeHigh();
  int getCutYRangeLow();
  int getCutYRangeHigh();
	double select(vector<double> input, int index);
  ~Map();
};

#endif
