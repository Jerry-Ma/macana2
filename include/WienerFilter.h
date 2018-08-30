#ifndef _WIENERFILTER_H_
#define _WIENERFILTER_H_

#include "nr3.h"
#include "AnalParams.h"
#include "Coaddition.h"
#include "NoiseRealizations.h"

///WienerFilter - a collection of detectors
/** The Array class is responsible for everything that is common
    across all the detectors and for all methods that work on
    the detector array as a whole (such as PCA cleaning, for 
    example).
**/
class WienerFilter
{
protected:
  ///analysis parameters
  AnalParams* ap;

  ///storage for current and future wf
  MatDoub Denom;               ///<denominator to filter (calc once)
  MatDoub Nume;                ///<numerator to the filter (changes)
  MatDoub rr;
  MatDoub vvq;
  MatDoub tplate;

  bool uniformWeight;          ///<set to force uniform weighting

  int nrows;                   ///<number of rows in coadded maps
  int ncols;                   ///<number of columns in coadded maps
  int nx;                      ///<xdim of wiener filtered maps
  int ny;                      ///<ydim of wiener filtered maps
  double pixelSize;            ///<pixel size in radians
  double diffx;
  double diffy;

public:
  ///methods
  WienerFilter(AnalParams* ap, Coaddition* cmap);
  bool prepareTemplate(Coaddition *cmap);
  bool prepareGaussianTemplate(Coaddition *cmap);
  bool filterCoaddition(Coaddition* cmap);
  bool filterNoiseMaps(NoiseRealizations *nr);
  bool calcRr(Coaddition* cmap);
  bool calcVvq();
  bool calcNumerator(MatDoub &mflt);
  bool calcDenominator();
  void simpleWienerFilter2d(Coaddition *cmap);
//  void gaussFilter(Coaddition *cmap);
  ~WienerFilter();
};

#endif
