#ifndef _SIMPARAMS_H_
#define _SIMPARAMS_H_

#include <string.h>
#include <nr3.h>

#include "tinyxml2.h"



class SimParams{


 protected:
  ///Observation files and paths
  string apXml;                     ///<the xml file for the analysis
  int nFiles;                       ///<the number of files to be reduced
  string rawDataPath;               ///<the path for the raw data
  string bsPath;                    ///<the path for the corresponding bstats
  string mapPath;                   ///<the path for the input simulated files
  string mapFile;					///<input simulated file


  //parameters
  bool addSignal;       			//True add signal to current streams, False simulated signal only-->
  double fluxFactor;
  double atmFreq;
  double noiseChunk;

  void throwXmlError(string p);

 public:
  SimParams(tinyxml2::XMLElement *apXml);
  SimParams(SimParams *sp);
  string getMapFile();
  ~SimParams();
  bool getAddSignal();
  double getFluxFactor();
  double getAtmFreq();
  double getNoiseChunk();

};

#endif
