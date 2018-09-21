#ifndef _COADDITION_H_
#define _COADDITION_H_

#include "nr3.h"
#include "AnalParams.h"
#include "Array.h"
#include "Map.h"
#include "Telescope.h"
#include "PointSource.h"

///Coaddition - a coaddition of many maps to form a single map set.
/** The Coaddition class takes the output of many Observations 
    to make a set of coadded maps.  Currently the coaddition is done
    very naively and so perhaps it can be improved.
    \todo Implement wiener filtering of the coadded maps
    \todo Review coaddition procedure to make sure it matches idl utils.
    \todo Carefully compare to idl utilities for several cases.
    \todo Check this for large maps.
**/
class Coaddition
{
 protected:
  AnalParams* ap;              ///<pointer to our analysis parameters

  //map carryovers
  int nrows;                   ///<number of rows in coadded maps
  int ncols;                   ///<number of columns in coadded maps
  int nPixels;                 ///<total number of map pixels
  double pixelSize;            ///<pixel size in radians
  VecDoub rowCoordsPhys;       ///<row coordinates in tangential projection
  VecDoub colCoordsPhys;       ///<column coordinates in tangential projection
  MatDoub xCoordsAbs;          ///<matrix of sphere coordinates in ra/az
  MatDoub yCoordsAbs;          ///<matrix of sphere coordinates in dec/el
  VecDoub masterGrid;          ///<tangential point on sphere

  //source finding
  double snglSourceWin;        ///<RADIUS masked out around source in radians
  double sourceSigma;          ///<"number of sigmas" considered a source
  int nSources;                ///<number of sources found
  bool synth;                  ///<synthetic source finding flag

  //tau information
  VecDoub individualMapsTau;
  double coaddAvgTau;

 public:
  Map* signal;                 ///<coadded unfiltered signal map
  Map* kernel;                 ///<coadded unfiltered kernel map
  Map* weight;                 ///<coadded unfiltered weight map
  Map* inttime;                ///<coadded unfiltered inttime map
  Map* filteredSignal;         ///<coadded filtered signal map
  Map* filteredKernel;         ///<coadded filtered kernel map
  Map* filteredWeight;         ///<coadded filtered weight map
  Map* tWeight;				   ///<coadded filtered weight map calculation from all maps clean timestreams
  Map* tSignal;				   ///<coadded filtered signal map calculation from all maps clean timestreams
  Map* tKernel;				   ///<coadded filtered kernel map calculation from all maps clean timestreams
  PointSource* sources;        ///<pointer to array of PointSources found

  //storage
  string ncdfFile;             ///<output ncdf filename

  //needs to be removed once Wiener filter is done
  double maxPreS2N;            ///<normalization constant for making fake S/N

  //methods
  Coaddition(AnalParams* ap);
  bool coaddMaps();
  bool coaddTimeStreams();
  bool histogramSignal();
  bool histogramSignal(double cc);
  bool histogramSignal(int nbins, double cc);
  bool histogramFilteredSignal();
  bool histogramFilteredSignal(double cc);
  bool histogramFilteredSignal(int nbins, double cc);
  int getNrows();
  int getNcols();
  double getPixelSize();
  double getRowCoordsPhys(int i);
  double getColCoordsPhys(int i);
  double getXCoordsAbs(int i, int j);
  double getYCoordsAbs(int i, int j);
  bool writeCoadditionToNcdf();
  bool writeCoadditionToFits(string fitsFilename);
  bool writeFilteredMapsToNcdf();
  bool findSources();
  bool findSources(double* signalMap, double* s2nMap, Coaddition* realCoadd,
                   double maxS2N); ///<last argument needs to go when Wiener
                                   ///<filter is done
  int getNSources();
  bool normalizeErrors(double noiseRms, double cov);
  ~Coaddition();
};

#endif


