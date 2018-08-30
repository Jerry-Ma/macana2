#ifndef _COMPLETENESSSIM_H_
#define _COMPLETENESSSIM_H_

#include "nr3.h"
#include "Coaddition.h"

///CompletenessSim - class for calculating map completeness
/**The CompletenessSims class contains all information regarding
   the completeness of a coadded map along with the methods used
   to calculate and save that information.
**/
class CompletenessSim
{
 protected:
  AnalParams* ap;                  ///<analysis parameters pointer
  Coaddition* cmap;                ///<coadded map pointer

  int nBins;                       ///<number of flux bins to calculate at
  int nSources;                    ///<number of sources per bin to inject
  VecDoub inputFluxes;             ///<fluxes completeness calculated at (mJy)

  VecInt synthRowInput;            ///<synthetic sources row input indices
  VecInt synthColInput;            ///<synthetic sources column input indices
  MatDoub synthRaInput;            ///<synthetic sources RA input coords
  MatDoub synthDecInput;           ///<synthetic sources Dec input coords

  MatInt statusKey;                ///<synthetic source recovery statuses
  MatDoub synthFluxes;             ///<synthetic sources recovered fluxes
  MatDoub synthS2N;                ///<synthetic sources recovered S/N
  MatDoub synthNoises;             ///<synthetic sources recovered noises
  MatDoub synthRas;                ///<synthetic sources recovered RAs
  MatDoub synthDecs;               ///<synthetic sources recovered Decs

  VecDoub completeness;            ///<completeness at each flux bin
  VecDoub compErr;                 ///<normal dist approx to binomial error
  VecDoub compErrLow;              ///<low underlying binomial dist error
  VecDoub compErrHigh;             ///<high underlying binomial dist error
  VecInt nConclusive;              ///<number conclusive results of simulation
                                   /// i.e. statusKey = 0 or 1

 public:

  //methods
  CompletenessSim(AnalParams* ap, Coaddition* cmap);
  bool calcCompleteness();
  double gCirc(double ra1, double dec1,
	       double ra2, double dec2);
  int getNBins();
  int getNSources();
  double getInputFlux(int i);
  int getStatusKey(int i, int j);
  double getCompleteness(int i);
  double getCompErr(int i);
  double getCompErrL(int i);
  double getCompErrH(int i);
  int getNConclusive(int i);
  bool addCompletenessToNcdf();
  ~CompletenessSim();
};

#endif
