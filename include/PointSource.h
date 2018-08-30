#ifndef _POINTSOURCE_H_
#define _POINTSOURCE_H_

#include "nr3.h"
#include "AnalParams.h"
#include "Map.h"

///PointSource - what SourceLocate::findSources() found
/**The PointSource class contains all information about the
   source that was found using SourceLocate::findSources()
   on a coadded map.
**/
class PointSource
{
public:
  //source identification
  int sID;                         ///<source ID
  int nSourcesParentMap;           ///<number of sources in parent map

  //parent map members
  AnalParams* ap;                  ///<analysis parameters pointer
  double* mapSignal;               ///<parent signal map pointer
  double* mapWeight;               ///<parent weight map pointer
  double* mapRowCoordsPhys;        ///<parent map phys row coords pointer
  double* mapColCoordsPhys;        ///<parent map phys column coords pointer
  double* mapXCoordsAbs;           ///<parent map absolute x coords pointer
  double* mapYCoordsAbs;           ///<parent map absolute y coords pointer
  int mapNRows;                    ///<number rows in parent map
  int mapNCols;                    ///<number columns in parent map

  //locations
  double centerRaAbs;               ///<absolute Ra of peak pixel
  double centerDecAbs;              ///<absoulte Dec of peak pixel
  double centerRaPhys;              ///<delta-source Ra of peak pixel
  double centerDecPhys;             ///<delta-source Dec of peak pixel
  int centerXPos;                   ///<x coordinate index of peak pixel
  int centerYPos;                   ///<y coordinate index of peak pixel
  double raCentroid;                ///<centroided absolute Ra
  double decCentroid;               ///<centroided absolute Dec
  double raPhysCentroid;            ///<centroided physical Ra
  double decPhysCentroid;           ///<centroided physiclal Dec
  double xCentroid;                 ///<centroided x coordinate index
  double yCentroid;                 ///<centroided y coordinate index

  //source properties
  double centerFlux;                ///<flux of source peak pixel (Jy)
  double centerNoise;               ///<niose of source peak pixel (Jy)
  double centerS2N;                 ///<S/N of source peak pixel

  //postage stamp
  Map* postageStamp;                ///<map object containing the postage stamp
  MatDoub xCoordsAbs;               ///<row absolute coordinates of stamp
  MatDoub yCoordsAbs;               ///<col absolute coordinates of stamp
  VecDoub gaussParams;              ///<Gaussian fit parameters of source

  //storage
  string parentMapFile;             ///<ncdf filename source found in
  double centroidWin;               ///<source centroid window (radians)
  double pixelSize;

  //fit output
  double redchisq;

  //methods
  bool initialize(AnalParams* ap,
		  double* mapSignal,
                  double *mapWeight,
		  double *mapRowCoordsPhys,
                  double* mapColCoordsPhys,
		  double* mapXCoordsAbs,
                  double* mapYCoordsAbs,
		  int mapNRows, int mapNCols,
		  double pixelSize, std::string parentMapFile);
  bool makePostageStamp();
  bool centroidSource();
  bool fitGaussianToSource(double fwhmx=0, double fwhmy=0);
  bool addSourceToNCDF();

  //constructor
  PointSource();
  ~PointSource();
};

#endif
