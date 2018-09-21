#include <netcdfcpp.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <cstdio>
#include <sys/stat.h>
#include <iostream>
using namespace std;

#include "nr3.h"
#include "AnalParams.h"
#include "tinyxml2.h"
#include "SimParams.h"

///AnalParams constructor
/** The AnalParams constructor pulls analysis parameters from an xml
    file (using tinyxml2) and stores them in the object for other
    classes to pull from.  
**/
AnalParams::AnalParams(string apXmlFile)
{
  beammapping = 0;
  //start out by generating the "global" random number generator
  macanaRandom = new GslRandom();

  //get all analysis parameters and file list from xml file
  apXml = apXmlFile;
  tinyxml2::XMLDocument xml;

  struct stat buf;

  if(stat(apXml.c_str(),&buf) == -1){
	  cerr<<"AnalParams()::Fatal Error. Analysis parameter file not found: "<<
			  apXml.c_str()<< endl;
	  exit (-1);
  }


  xml.LoadFile(apXml.c_str());

  this->simParams = NULL;
  //pack up the analysis parameters and steps first
  tinyxml2::XMLElement* xAnalysis;
  tinyxml2::XMLElement* xSteps;
  tinyxml2::XMLElement* xParameters;
  tinyxml2::XMLElement* xtmp;
  xAnalysis = xml.FirstChildElement("analysis");

  if(!xAnalysis){
    cerr << "XML error in " << apXml << ": ";
    cerr << "<analysis> not found." << endl;
    exit(1);
  }

  //the analysis steps
  xSteps = xAnalysis->FirstChildElement("analysisSteps");
  if(!xSteps){
    cerr << "XML error in " << apXml << ": ";
    cerr << "<analysisSteps> not found." << endl;
    exit(1);
  }
  xtmp = xSteps->FirstChildElement("mapIndividualObservations");
  if(!xtmp) throwXmlError("analysisSteps:mapIndividualObservations not found.");
  mapIndividualObservations = atoi(xtmp->GetText());

  xtmp = xSteps->FirstChildElement("coaddObservations");
  if(!xtmp) throwXmlError("analysisSteps:coaddObservations not found.");
  coaddObservations = atoi(xtmp->GetText());

  xtmp = xSteps->FirstChildElement("fitCoadditionToGaussian");
  if(!xtmp) throwXmlError("analysisSteps:fitCoadditionToGaussian not found.");
  fitCoadditionToGaussian = atoi(xtmp->GetText());

  xtmp = xSteps->FirstChildElement("produceNoiseMaps");
  if(!xtmp) throwXmlError("analysisSteps:produceNoiseMaps not found.");
  produceNoiseMaps = atoi(xtmp->GetText());

  xtmp = xSteps->FirstChildElement("applyWienerFilter");
  if(!xtmp) throwXmlError("analysisSteps:applyWienerFilter not found.");
  applyWienerFilter = atoi(xtmp->GetText());

  //sanity check analysis steps
  if(!coaddObservations){
    if(fitCoadditionToGaussian){
      cerr << endl;
      cerr << "XML error in " << apXml << ": ";
      cerr << "In order to fit the coadded map to a Gaussian \n";
      cerr << "you must coadd the individual observations ";
      cerr << "in the same run. \n (i.e., <coaddObservations> 1 </coaddObservations>)";
      cerr << endl << endl;
      exit(1);
    }
    if(produceNoiseMaps){
      cerr << endl;
      cerr << "XML error in " << apXml << ": ";
      cerr << "In order to produce noise maps \n";
      cerr << "you must coadd the individual observations ";
      cerr << "in the same run. \n (i.e., <coaddObservations> 1 </coaddObservations>)";
      cerr << endl << endl;
      exit(1);
    }
    if(applyWienerFilter){
      cerr << endl;
      cerr << "XML error in " << apXml << ": ";
      cerr << "In order to apply the wiener filter \n";
      cerr << "you must coadd the individual observations ";
      cerr << "in the same run. \n (i.e., <coaddObservations> 1 </coaddObservations>)";
      cerr << endl;
      cerr << "(Don't forget that you ALSO need to make the noise maps " << endl;
      cerr << "(i.e., <produceNoiseMaps> 1 </produceNoiseMaps>)";
      cerr << endl << endl;
      exit(1);
    }
  }
  if(applyWienerFilter && !produceNoiseMaps){
    cerr << endl;
    cerr << "XML error in " << apXml << ": ";
    cerr << "In order to apply the wiener filter you must produce noise \n maps ";
    cerr << "in the same run. (i.e., <produceNoiseMaps> 1 </produceNoiseMaps>)";
    cerr << endl << endl;
    exit(1);
  }


  //the parameters section
  xParameters = xAnalysis->FirstChildElement("parameters");
  if(!xParameters){
    cerr << "<parameters> not found in input xml file.";
    cerr << "These MUST be set." << endl;
    exit(1);
  }
  xtmp = xParameters->FirstChildElement("despikeSigma");
  if(!xtmp) throwXmlError("despikeSigma not found.");
  despikeSigma = atof(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("lowpassFilterKnee");
  if(!xtmp) throwXmlError("lowpassfilterKnee not found.");
  lowpassFilterKnee = atof(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("timeOffset");
  if(!xtmp) throwXmlError("timeOffset not found.");
  timeOffset = atof(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("timeChunk");
  if(!xtmp) throwXmlError("timeChunk not found.");
  timeChunk = atof(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("cutTurnArounds");
  if(!xtmp){
    cutTurnArounds=1;
  } else {
    cutTurnArounds = bool(atoi(xtmp->GetText()));
  }

  xtmp = xParameters->FirstChildElement("cutSamplesAtEndOfScans");
  if(!xtmp){
    cutSamplesAtEndOfScans=20;
  } else {
    cutSamplesAtEndOfScans = atoi(xtmp->GetText());
  }

  xtmp = xParameters->FirstChildElement("coverageThreshold");
  if(!xtmp){
    coverageThreshold=0.75;
  } else {
    coverageThreshold = atof(xtmp->GetText());
  }

  xtmp = xParameters->FirstChildElement("neigToCut");
  if(!xtmp) throwXmlError("neigToCut not found.");
  neigToCut = atoi(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("cutStd");
  if(!xtmp) throwXmlError("cutStd not found.");
  cutStd = atof(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("noiseMapsPerObs");
  if(!xtmp){
	  nNoiseMapsPerObs = 5;
  }else nNoiseMapsPerObs = atoi(xtmp->GetText());
  
  //++Cottingham method stuff
  xtmp = xParameters->FirstChildElement("splineOrder");
  if(!xtmp) throwXmlError("splineOrder not found.");
  order = atoi(xtmp->GetText());
  
  xtmp = xParameters->FirstChildElement("cleanPixelSize");
  if(!xtmp) throwXmlError("cleanPixelSize not found.");
  cleanPixelSize = atof(xtmp->GetText());
  
  xtmp = xParameters->FirstChildElement("cleanStripe");
  if(!xtmp) throwXmlError("cleanStripe not found.");
  cleanStripe = atof(xtmp->GetText());
  const char *strMethod = xtmp->Attribute("method");
  stripeMethod.assign("");
  if (strMethod)
	  stripeMethod.assign(strMethod);
  
  xtmp = xParameters->FirstChildElement("controlChunk");
  if(!xtmp) throwXmlError("controlChunk not found.");
  controlChunk = atof(xtmp->GetText());
  
  xtmp = xParameters->FirstChildElement("resample");
  if(!xtmp){
	  resample = 1;
  }else
	  resample = atoi(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("tOrder");
  if(!xtmp)
	  tOrder = 0;
  else
	  tOrder = atoi(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("threadNumber");
  if(!xtmp){
    nThreads = 1;
  } else nThreads = atoi(xtmp->GetText());
  if (nThreads < 1) nThreads = 1;

  saveTimeStreams=false;
  xtmp = xParameters->FirstChildElement("saveTimeStreams");
  if(xtmp){
	  if (atoi(xtmp->GetText())!=0)
	  	  saveTimeStreams=true;
  }


  xtmp = xParameters->FirstChildElement("pixelSize");
  if(!xtmp) throwXmlError("pixelSize not found.");
  pixelSize = atof(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("approximateWeights");
  if(!xtmp) throwXmlError("approximateWeights not found.");
  approximateWeights = atoi(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("azelMap");
  if(!xtmp){
    azelMap = 0;
  } else {
    azelMap = atoi(xtmp->GetText());
  }

  //the master grid - set to zero to have Source determine it
  masterGridJ2000 = new double [2];
  masterGridJ2000_init = new double [2];
  xtmp = xParameters->FirstChildElement("masterGridJ2000_0");
  if(!xtmp) throwXmlError("masterGridJ2000_0 not found.");
  masterGridJ2000_init[0] = atof(xtmp->GetText());
  xtmp = xParameters->FirstChildElement("masterGridJ2000_1");
  if(!xtmp) throwXmlError("masterGridJ2000_1 not found.");
  masterGridJ2000_init[1] = atof(xtmp->GetText());
  masterGridJ2000[0] = masterGridJ2000_init[0];
  masterGridJ2000[1] = masterGridJ2000_init[1];


  //print these out to see that we got it right
  cerr << endl;
  cerr << "------------------ o -----------------------" << endl;
  cerr << "mapIndividualObservations: " << mapIndividualObservations << endl;
  cerr << "coaddObservations: " << coaddObservations << endl;
  cerr << "fitCoadditionToGaussian: " << fitCoadditionToGaussian << endl;
  cerr << "produceNoiseMaps: " << produceNoiseMaps << endl;
  cerr << "applyWienerFilter: " << applyWienerFilter << endl;
  cerr << "------------------ o -----------------------" << endl;
  cerr << "despikeSigma: " << despikeSigma << endl;
  cerr << "lowpassfilterKnee: " << lowpassFilterKnee << endl;
  cerr << "timeOffset: " << timeOffset << endl;
  cerr << "timeChunk: " << timeChunk << endl;
  cerr << "neigToCut: " << neigToCut << endl;
  cerr << "cutStd: " << cutStd << endl;
  cerr << "pixelSize: " << pixelSize << endl;
  cerr << "azelMap: " << azelMap << endl;
  cerr << "approximateWeights: " << approximateWeights << endl;
  cerr << "splineOrder: " << order << endl;
  cerr << "cleanPixelSize: " << cleanPixelSize << endl;
  cerr << "controlChunk: " << controlChunk<<endl;
  cerr << "cleanStripe: " <<cleanStripe<<endl;
  cerr << "resample: "<< resample <<endl;
  cerr << "ThreadNumber: " <<nThreads<<endl;
  cerr << "initial Mastergrid: [" << masterGridJ2000[0];
  cerr << "," << masterGridJ2000[1] << "]" << endl;

  
  //get and check the filenames and paths given for existance


  //the coaddition path and filenames
  if(coaddObservations){
    tinyxml2::XMLElement* xCoadd;
    xCoadd = xAnalysis->FirstChildElement("coaddition");
    if(!xCoadd){
      cerr << "AnalParams: ap: coaddition settings not found." << endl;
      exit(1);
    } else {
      xtmp = xCoadd->FirstChildElement("mapPath");
      if(!xtmp) throwXmlError("mapPath not found.");
      coaddOutPath = xtmp->GetText();
      if(stat(coaddOutPath.c_str(),&buf) == -1){
	cerr << "AnalParams: ap: mapPath not found: " << coaddOutPath << endl;
	exit(1);
      }
      xtmp = xCoadd->FirstChildElement("mapFile");
      if(!xtmp) throwXmlError("mapFile not found.");
      coaddOutFile = coaddOutPath;
      coaddOutFile.append(xtmp->GetText());
    }
  }


  //the noise realization values
  if(produceNoiseMaps){
    tinyxml2::XMLElement* xNoise;
    xNoise = xAnalysis->FirstChildElement("noiseRealization");
    if(!xNoise){
      cerr << "AnalParams: ap: noiseRealization settings not found." << endl;
      exit(1);
    } else {
      xtmp = xNoise->FirstChildElement("nRealizations");
      if(!xtmp) throwXmlError("noiseRealization->nRealizations not found.");
      nRealizations = atoi(xtmp->GetText());
      if(nRealizations == 0) {
	cerr << "noiseRealizations->nRealizations = 0. Set to non-zero ";
	cerr << "value or set produceNoiseMaps=0." << endl;
      }
      
      xtmp = xNoise->FirstChildElement("noisePath");
      if(!xtmp) throwXmlError("noisePath not found.");
      noisePath = xtmp->GetText(); 
      if(stat(noisePath.c_str(),&buf) == -1){
	cerr << "AnalParams: ap: noisePath not found: " << noisePath << endl;
	exit(1);
      }
      
      avgNoisePsdFile=noisePath;
      xtmp = xNoise->FirstChildElement("avgNoisePsdFile");
      if(!xtmp) throwXmlError("avgNoisePsdFile name not set.");
      avgNoisePsdFile.append(xtmp->GetText());
      
      avgNoiseHistFile = noisePath;
      xtmp = xNoise->FirstChildElement("avgNoiseHistFile");
      if(!xtmp) throwXmlError("avgNoiseHistFile name not set.");
      avgNoiseHistFile.append(xtmp->GetText());
    }
  }


  //are we to wiener filter the maps?
  if(applyWienerFilter){
    tinyxml2::XMLElement* xWienerFilter;
    xWienerFilter = xAnalysis->FirstChildElement("wienerFilter");
    if(xWienerFilter){
      xtmp = xWienerFilter->FirstChildElement("gaussianTemplate");
      if(!xtmp) throwXmlError("gaussianTemplate not set");
      gaussianTemplate = (atof(xtmp->GetText()) != 0);
      if(gaussianTemplate){
	xtmp = xWienerFilter->FirstChildElement("gaussianTemplateFWHM");
	if(!xtmp) throwXmlError("Wiener Filter gaussian template requested but \ngaussianTemplateFWHM not set.");
	gaussianTemplateFWHM = (atof(xtmp->GetText()));
      } else gaussianTemplateFWHM = 0.;
      xtmp = xWienerFilter->FirstChildElement("lowpassOnly");
      if(!xtmp) throwXmlError("lowpassOnly not set");
      lowpassOnly = (atof(xtmp->GetText()) != 0);
      
      xtmp = xWienerFilter->FirstChildElement("highpassOnly");
      if(!xtmp) throwXmlError("highpassOnly not set");
      highpassOnly = (atof(xtmp->GetText()) != 0);
      
      xtmp = xWienerFilter->FirstChildElement("normalizeErrors");
      if(!xtmp) throwXmlError("normalizeErrors not set");
      normalizeErrors = (atof(xtmp->GetText()) != 0);
    } else {
      cerr << "AnalParams: ap: wienerFilter settings not found." << endl;
      exit(1);
    }
  }


  //check to see if the altKernel parameters are set
  //These are parameters for choosing an alternative kernel shape
  //to the normal Gaussian beam shape
  tinyxml2::XMLElement* xAltKernel;
    xAltKernel = xAnalysis->FirstChildElement("altKernel");
    if(xAltKernel){
      altKernel = 1;
      xtmp = xAltKernel->FirstChildElement("kernelName");
      if(!xtmp) throwXmlError("kernelName not set");
      kernelName = xtmp->GetText();
      if(kernelName == "core"){
	xtmp = xAltKernel->FirstChildElement("coreR0");
	if(!xtmp) throwXmlError("altKernel=core but coreR0 not set");
	coreR0 = atof(xtmp->GetText());
	xtmp = xAltKernel->FirstChildElement("coreP");
	if(!xtmp) throwXmlError("altKernel=core but coreP not set");
	coreP = atof(xtmp->GetText());
	xtmp = xAltKernel->FirstChildElement("coreAxisRatio");
	if(!xtmp) throwXmlError("altKernel=core but coreAxisRatio not set");
	coreAxisRatio = atof(xtmp->GetText());
      } else {
	cerr << "Unknown altKernel name supplied." 
	     << "  Currently support 'core' only." << endl;
	exit(1);
      }
    }


  //the observations
  tinyxml2::XMLElement* xObs;
  xObs = xAnalysis->FirstChildElement("observations");
  xtmp = xObs->FirstChildElement("rawDataPath");
  if(!xtmp) throwXmlError("rawDataPath not found.");
  rawDataPath = xtmp->GetText();
  if(stat(rawDataPath.c_str(),&buf) == -1){
    cerr << "AnalParams: ap: rawDataPath not found: " << rawDataPath << endl;
    exit(1);
  }

  xtmp = xObs->FirstChildElement("bsPath");
  if(!xtmp) throwXmlError("bsPath not found.");
  bsPath = xtmp->GetText();
  if(stat(bsPath.c_str(),&buf) == -1){
    cerr << "AnalParams: ap: bsPath not found: " << bsPath << endl;
    exit(1);
  }

  xtmp = xObs->FirstChildElement("mapPath");
  if(!xtmp) throwXmlError("mapPath not found.");
  mapPath = xtmp->GetText();
  if(stat(mapPath.c_str(),&buf) == -1){
    cerr << "AnalParams: ap: mapPath not found: " << mapPath << endl;
    exit(1);
  }

  xtmp = xObs->FirstChildElement("nFiles");
  if(!xtmp) throwXmlError("nFiles not found.");
  nFiles = atoi(xtmp->GetText());
  fileList = new string [nFiles];
  bstatList = new string [nFiles];
  mapFileList = new string [nFiles];
  bsOffsetList.resize(2,nFiles);
  char fId[20];
  //int n;
  string ferror;
  for(int i=0;i<nFiles;i++){
    //n = sprintf(fId, "f%d",i);
	sprintf(fId, "f%d",i);
    ferror = fId;
    xtmp = xObs->FirstChildElement(fId)->FirstChildElement("fileName");
    if(!xtmp) throwXmlError(ferror.append(": fileName not found."));
    fileList[i].assign(rawDataPath);
    fileList[i].append(xtmp->GetText());
    if(stat(fileList[i].c_str(),&buf) == -1){
    cerr << "AnalParams: ap: data file not found: " << fileList[i] << endl;
    exit(1);
    }

    xtmp = xObs->FirstChildElement(fId)->FirstChildElement("bsName");
    if(!xtmp) throwXmlError(ferror.append(": bsName not found."));
    bstatList[i].assign(bsPath);
    bstatList[i].append(xtmp->GetText());
    if(stat(bstatList[i].c_str(),&buf) == -1){
    cerr << "AnalParams: ap: bstat not found: " << bstatList[i] << endl;
    exit(1);
    }

    xtmp = xObs->FirstChildElement(fId)->FirstChildElement("mapName");
    if(!xtmp) throwXmlError(ferror.append(": mapName not found."));
    mapFileList[i].assign(mapPath);
    mapFileList[i].append(xtmp->GetText());

    xtmp = xObs->FirstChildElement(fId)->FirstChildElement("bsOffset_0");
    if(!xtmp) throwXmlError(ferror.append(": bsOffset_0 not found."));
    bsOffsetList[0][i] = atof(xtmp->GetText());

    xtmp = xObs->FirstChildElement(fId)->FirstChildElement("bsOffset_1");
    if(!xtmp) throwXmlError(ferror.append(": bsOffset_1 not found."));
    bsOffsetList[1][i] = atof(xtmp->GetText());
  }

  cerr << "Planning to analyze " << nFiles << " files." << endl;

  xtmp = xAnalysis->FirstChildElement("simulate");
  if (xtmp != NULL)
	  this->simParams = new SimParams(xtmp);

  doSubtract = false;
  subtractFirst =false;
  xtmp = xAnalysis->FirstChildElement("subtract");
  if (xtmp != NULL){
    struct stat st;
    tinyxml2::XMLElement *stmp = xtmp->FirstChildElement("subPath");
    if (stmp != NULL){
      subtractPath.assign(stmp->GetText());
      if(stat(subtractPath.c_str(),&st) == -1){
        cerr << "AnalParams: subtractPath not such file/directory: " << subtractPath << endl;
        exit(1);
      }
      tinyxml2::XMLElement *sfiletmp = xtmp->FirstChildElement("subFile");
      if (sfiletmp !=NULL){
        subtractFile.assign(sfiletmp->GetText());
        if(stat((subtractPath+subtractFile).c_str(),&st) == -1){
          cerr << "AnalParams: subtractFile not  such file/directory: " << subtractFile << endl;
        }
        doSubtract = true;
      }
      tinyxml2::XMLElement *sFirst = xtmp->FirstChildElement("subFirst");
      if (sFirst)
      subtractFirst = true;
    }
  }



  //do we have post reduction information?
  tinyxml2::XMLElement* xPostReduction;
  xPostReduction = xAnalysis->FirstChildElement("postReductionAnalysis");
  if(xPostReduction){
    cerr << "Planning to run post-reduction analysis." << endl;
    postReductionAnalysis = true;
    tinyxml2::XMLElement* xSourceFinding;
    //source finding values
    xSourceFinding = xPostReduction->FirstChildElement("sourceFinding");
    if(xSourceFinding){
      cerr << "   - Source Finding " << endl;
      findSources = true;
      xtmp = xSourceFinding->FirstChildElement("beamSize");
      if(!xtmp) throwXmlError("source finding beamSize not found.");
      beamSize = atof(xtmp->GetText());

      xtmp = xSourceFinding->FirstChildElement("covCut");
      if(!xtmp) throwXmlError("source finding covCut not found.");
      covCut = atof(xtmp->GetText());

      xtmp = xSourceFinding->FirstChildElement("snglSourceWin");
      if(!xtmp) throwXmlError("source finding sngleSourceWin not found.");
      snglSourceWin = atof(xtmp->GetText());
      
      xtmp = xSourceFinding->FirstChildElement("sourceSigma");
      if(!xtmp) throwXmlError("source finding sourceSigma not found.");
      sourceSigma = atof(xtmp->GetText());

      xtmp = xSourceFinding->FirstChildElement("negativeToo");
      if(!xtmp){
	negativeToo = 0;
      } else {
	negativeToo = atoi(xtmp->GetText());
      }

      xtmp = xSourceFinding->FirstChildElement("mapNegative");
      if(!xtmp){
	mapNegative = 0;
      } else {
	mapNegative = atoi(xtmp->GetText());
      }

      xtmp = xSourceFinding->FirstChildElement("centroidSources");
      if(!xtmp) throwXmlError("Source finding centroidSources not found.");
      sfCentroidSources = atoi(xtmp->GetText());

      xtmp = xSourceFinding->FirstChildElement("fitGaussians");
      if(!xtmp) throwXmlError("Source finding fitGaussians not found.");
      sfFitGaussians = atoi(xtmp->GetText());
      
      xtmp = xSourceFinding->FirstChildElement("psSideLength");
      if(!xtmp) throwXmlError("source finding psSideLength not found.");
      psSideLength = atoi(xtmp->GetText());
    }

    //completeness calculation values
    tinyxml2::XMLElement* xCompletenessCalc;
    //source finding values
    xCompletenessCalc = xPostReduction->FirstChildElement("calcCompleteness");
    if(xCompletenessCalc){
      cerr << "   - Completeness Calculation " << endl;
      xtmp = xCompletenessCalc->FirstChildElement("nFluxBins");
      if(!xtmp) throwXmlError("calc completeness nFluxBins not found.");
      nFluxBins = atoi(xtmp->GetText());
      
      xtmp = xCompletenessCalc->FirstChildElement("nSynthSources");
      if(!xtmp) throwXmlError("calc completeness nSynthSources not found.");
      nSynthSources = atoi(xtmp->GetText());
      
      xtmp = xCompletenessCalc->FirstChildElement("minFlux");
      if(!xtmp) throwXmlError("calc completeness minFlux not found.");
      minFlux = atof(xtmp->GetText());
      
      xtmp = xCompletenessCalc->FirstChildElement("maxFlux");
      if(!xtmp) throwXmlError("calc completeness maxFlux not found.");
      maxFlux = atof(xtmp->GetText());

      xtmp = xCompletenessCalc->FirstChildElement("recovS2N");
      if(!xtmp) throwXmlError("calc completeness recovS2N not found.");
      recovS2N = atof(xtmp->GetText());
      
      xtmp = xCompletenessCalc->FirstChildElement("recovRadius");
      if(!xtmp) throwXmlError("calc completeness recovRadius not found.");
      recovRadius = atof(xtmp->GetText());
    }
  }

  cerr << "------------------ o -----------------------" << endl;
}

//----------------------------- o ---------------------------------------

///AnalParams constructor for beammapping
/* This AnalParams constructor requires fewer parameters. 
Requires additional BeammapSourceFlux and output paths for beammapping
*/

AnalParams::AnalParams(string apXmlFile, int beammap)
{
  beammapping = 1;
  macanaRandom = new GslRandom();

  //get all analysis parameters and file list from xml file
  apXml = apXmlFile;
  tinyxml2::XMLDocument xml;
  struct stat buf;

  if(stat(apXml.c_str(),&buf) == -1){
    cerr<< "AnalParams()::Fatal Error. Analysis parameter file not found: " << apXml.c_str() << endl;
    exit(-1);
  }
  xml.LoadFile(apXml.c_str());


  this->simParams = NULL;
  //pack up the analysis parameters and steps first
  tinyxml2::XMLElement* xAnalysis;
  tinyxml2::XMLElement* xParameters;
  tinyxml2::XMLElement* xtmp;
  xAnalysis = xml.FirstChildElement("analysis");

  if(!xAnalysis){
    cerr << "XML error in " << apXml << ": ";
    cerr << "<analysis> not found." << endl;
    exit(1);
  }

  //the parameters section
  xParameters = xAnalysis->FirstChildElement("parameters");
  if(!xParameters){
    cerr << "<parameters> not found in input xml file.";
    cerr << " These MUST be set." << endl;
    exit(1);
  }
  xtmp = xParameters->FirstChildElement("despikeSigma");
  if(!xtmp) throwXmlError("despikeSigma not found.");
  despikeSigma = atof(xtmp->GetText());
  xtmp = xParameters->FirstChildElement("lowpassFilterKnee");
  if(!xtmp) throwXmlError("lowpassFilterKnee not found.");
  lowpassFilterKnee = atof(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("timeOffset");
  if(!xtmp) throwXmlError("timeOffset not found.");
  timeOffset = atof(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("timeChunk");
  if(!xtmp) throwXmlError("timeChunk not found.");
  timeChunk = atof(xtmp->GetText());

	xtmp = xParameters->FirstChildElement("cutTurnArounds");
  if(!xtmp){
    cutTurnArounds=1;
  } else {
    cutTurnArounds = bool(atoi(xtmp->GetText()));
  }

  xtmp = xParameters->FirstChildElement("cutSamplesAtEndOfScans");
  if(!xtmp){
    cutSamplesAtEndOfScans=20;
  } else {
    cutSamplesAtEndOfScans = atoi(xtmp->GetText());
  }

  xtmp = xParameters->FirstChildElement("writeBeammapToNcdf");
  if(!xtmp){
    writeBeammapToNcdf=0;
  } else {
    writeBeammapToNcdf = bool(atoi(xtmp->GetText()));
  }

  xtmp = xParameters->FirstChildElement("cleanIterationCap");
  if(!xtmp){
    cerr << "default clean iteration cap set to 10" << endl;
    cleanIterationCap = 10;
  } else {
    cleanIterationCap = atoi(xtmp->GetText());
  }

  xtmp = xParameters->FirstChildElement("cleanIterationCutoff");
  if(!xtmp){
    cerr << "default clean iteration cutoff set to 0.05" << endl;
  } else {
    cleanIterationCutoff = atof(xtmp->GetText());
  }

  xtmp = xParameters->FirstChildElement("coverageThreshold");
  if(!xtmp){
    coverageThreshold=0.75;
  } else {
    coverageThreshold = atof(xtmp->GetText());
  }

  xtmp = xParameters->FirstChildElement("neigToCut");
  if(!xtmp) throwXmlError("neigToCut not found.");
  neigToCut = atoi(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("cutStd");
  if(!xtmp) throwXmlError("cutStd not found.");
  cutStd = atof(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("noiseMapsPerObs");
  if(!xtmp){
    nNoiseMapsPerObs = 5;
  }else nNoiseMapsPerObs = atoi(xtmp->GetText());

  //++Cottingham method stuff
  xtmp = xParameters->FirstChildElement("splineOrder");
  if(!xtmp) throwXmlError("splineOrder not found.");
  order = atoi(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("cleanPixelSize");
  if(!xtmp) throwXmlError("cleanPixelSize not found.");
  cleanPixelSize = atof(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("cleanStripe");
  if(!xtmp) throwXmlError("cleanStripe not found.");
  cleanStripe = atof(xtmp->GetText());
  const char *strMethod = xtmp->Attribute("method");
  stripeMethod.assign("");
  if (strMethod)
    stripeMethod.assign(strMethod);

	xtmp = xParameters->FirstChildElement("controlChunk");
  if(!xtmp) throwXmlError("controlChunk not found.");
  controlChunk = atof(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("resample");
  if(!xtmp){
    resample = 1;
  }else
    resample = atoi(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("tOrder");
  if(!xtmp)
    tOrder = 0;
  else
    tOrder = atoi(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("threadNumber");
  if(!xtmp){
    nThreads = 1;
  } else nThreads = atoi(xtmp->GetText());
  if (nThreads < 1) nThreads = 1;

  saveTimeStreams=false;
  xtmp = xParameters->FirstChildElement("saveTimeStreams");
  if(xtmp){
    if (atoi(xtmp->GetText())!=0)
        saveTimeStreams=true;
  }


  xtmp = xParameters->FirstChildElement("pixelSize");
  if(!xtmp) throwXmlError("pixelSize not found.");
  pixelSize = atof(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("approximateWeights");
  if(!xtmp) throwXmlError("approximateWeights not found.");
  approximateWeights = atoi(xtmp->GetText());

  xtmp = xParameters->FirstChildElement("azelMap");
  if(!xtmp){
    azelMap = 0;
  } else {
    azelMap = atoi(xtmp->GetText());
  }

	//the master grid - set to zero to have Source determine it
  masterGridJ2000 = new double [2];
  masterGridJ2000_init = new double [2];
  xtmp = xParameters->FirstChildElement("masterGridJ2000_0");
  if(!xtmp) throwXmlError("masterGridJ2000_0 not found.");
  masterGridJ2000_init[0] = atof(xtmp->GetText());
  xtmp = xParameters->FirstChildElement("masterGridJ2000_1");
  if(!xtmp) throwXmlError("masterGridJ2000_1 not found.");
  masterGridJ2000_init[1] = atof(xtmp->GetText());
  masterGridJ2000[0] = masterGridJ2000_init[0];
  masterGridJ2000[1] = masterGridJ2000_init[1];

  //print these out to see that we got it right
  cerr << endl;
  cerr << "------------------ o -----------------------" << endl;
  cerr << "mapIndividualObservations: " << mapIndividualObservations << endl;
  cerr << "coaddObservations: " << coaddObservations << endl;
  cerr << "fitCoadditionToGaussian: " << fitCoadditionToGaussian << endl;
  cerr << "produceNoiseMaps: " << produceNoiseMaps << endl;
  cerr << "applyWienerFilter: " << applyWienerFilter << endl;
  cerr << "------------------ o -----------------------" << endl;
  cerr << "despikeSigma: " << despikeSigma << endl;
  cerr << "lowpassfilterKnee: " << lowpassFilterKnee << endl;
  cerr << "timeOffset: " << timeOffset << endl;
  cerr << "timeChunk: " << timeChunk << endl;
  cerr << "neigToCut: " << neigToCut << endl;
  cerr << "cutStd: " << cutStd << endl;
  cerr << "pixelSize: " << pixelSize << endl;
  cerr << "azelMap: " << azelMap << endl;
  cerr << "approximateWeights: " << approximateWeights << endl;
  cerr << "splineOrder: " << order << endl;
  cerr << "cleanPixelSize: " << cleanPixelSize << endl;
  cerr << "controlChunk: " << controlChunk<<endl;
  cerr << "cleanStripe: " <<cleanStripe<<endl;
  cerr << "resample: "<< resample <<endl;
  cerr << "ThreadNumber: " <<nThreads<<endl;
  cerr << "initial Mastergrid: [" << masterGridJ2000[0];
  cerr << "," << masterGridJ2000[1] << "]" << endl;

  //the observations
  tinyxml2::XMLElement* xObs;
  xObs = xAnalysis->FirstChildElement("observations");
  xtmp = xObs->FirstChildElement("rawDataPath");
  if(!xtmp) throwXmlError("rawDataPath not found.");
  rawDataPath = xtmp->GetText();
  if(stat(rawDataPath.c_str(),&buf) == -1){
    cerr << "AnalParams: ap: rawDataPath not found: " << rawDataPath << endl;
    exit(1);
  }

	xtmp = xObs->FirstChildElement("bsPath");
  if(!xtmp) throwXmlError("bsPath not found.");
  bsPath = xtmp->GetText();
  if(stat(bsPath.c_str(),&buf) == -1){
    cerr << "AnalParams: ap: bsPath not found: " << bsPath << endl;
    exit(1);
  }

  xtmp = xObs->FirstChildElement("outBeammapInfoPath");
  if(!xtmp) throwXmlError("outBeammapInfoPath not found.");
  outBeammapInfoPath = xtmp->GetText();
  if(stat(outBeammapInfoPath.c_str(),&buf) == -1){
    cerr << "AnalParams: ap: outBeammapInfoPath not found: " << outBeammapInfoPath << endl;
    exit(1);
  }

  xtmp = xObs->FirstChildElement("outBeammapNcdfPath");
  if(!xtmp) throwXmlError("outBeammapNcdfPath not found.");
  outBeammapNcdfPath = xtmp->GetText();
  if(stat(outBeammapNcdfPath.c_str(),&buf) == -1){
    cerr << "AnalParams: ap: outBeammapNcdfPath nout found: " << outBeammapNcdfPath << endl;
    exit(1);
  }

  xtmp = xObs->FirstChildElement("nFiles");
  if(!xtmp) throwXmlError("nFiles not found.");
  nFiles = atoi(xtmp->GetText());
  outBeammapInfoList = new string [nFiles];
  outBeammapNcdfList = new string [nFiles];
  fileList = new string [nFiles];
  bstatList = new string [nFiles];
  mapFileList = new string [nFiles];
  bsOffsetList.resize(2,nFiles);
  beammapSourceFluxList.resize(nFiles);
  char fId[20];
  //int n;
  string ferror;
  for(int i=0;i<nFiles;i++){
    //n = sprintf(fId, "f%d",i);
    sprintf(fId, "f%d",i);
    ferror = fId;
    xtmp = xObs->FirstChildElement(fId)->FirstChildElement("fileName");
    if(!xtmp) throwXmlError(ferror.append(": fileName not found."));
    fileList[i].assign(rawDataPath);
    fileList[i].append(xtmp->GetText());
    if(stat(fileList[i].c_str(),&buf) == -1){
    cerr << "AnalParams: ap: data file not found: " << fileList[i] << endl;
    exit(1);
    }

    xtmp = xObs->FirstChildElement(fId)->FirstChildElement("bsName");
    if(!xtmp) throwXmlError(ferror.append(": bsName not found."));
    bstatList[i].assign(bsPath);
    bstatList[i].append(xtmp->GetText());
    if(stat(bstatList[i].c_str(),&buf) == -1){
    cerr << "AnalParams: ap: bstat not found: " << bstatList[i] << endl;
    exit(1);
    }

    xtmp = xObs->FirstChildElement(fId)->FirstChildElement("outBeammapInfoName");
    if(!xtmp) throwXmlError(ferror.append(": outBeammapInfoName not found."));
    outBeammapInfoList[i].assign(outBeammapInfoPath);
    outBeammapInfoList[i].append(xtmp->GetText());
    
    if(writeBeammapToNcdf){
      xtmp = xObs->FirstChildElement(fId)->FirstChildElement("outBeammapNcdfName");
      if(!xtmp) throwXmlError(ferror.append(": outBeammapNcdfName not found."));
      outBeammapNcdfList[i].assign(outBeammapNcdfPath);
      outBeammapNcdfList[i].append(xtmp->GetText());
    }

    xtmp = xObs->FirstChildElement(fId)->FirstChildElement("beammapSourceFlux");
    if(!xtmp) throwXmlError(ferror.append(": beammapSourceFlux not found."));
    beammapSourceFluxList[i] = atof(xtmp->GetText());
  }

  cerr << "------------------ o -----------------------" << endl;

}


///AnalParams constructor
/** The AnalParams copy  constructor pulls analysis parameters from 
 *  another AnalParms object. Intended for multhithread operation
**/
//----------------------------- o ---------------------------------------
AnalParams::AnalParams(AnalParams *ap){
  this->apXml = ap->apXml;
  this->nFiles = ap->nFiles;
  this->rawDataPath=ap->rawDataPath;               ///<the path for the raw data
  this->bsPath=ap->bsPath;                    ///<the path for the corresponding bstats
  this->mapPath = ap->mapPath;                   ///<the path for the output map nc files
  this->bsOffsetList=ap->bsOffsetList;
  this->dataFile = ap->dataFile;
  this->bolostatsFile = ap->bolostatsFile;
  this->mapFile = ap->mapFile;
  this->mapIndividualObservations = ap->mapIndividualObservations;
  this->coaddObservations = ap->coaddObservations;
  this->fitCoadditionToGaussian = ap->fitCoadditionToGaussian;
  this->produceNoiseMaps = ap->produceNoiseMaps;
  this->applyWienerFilter = ap->applyWienerFilter;
  this->sourceName = ap->sourceName;
  this->coaddOutPath = ap->coaddOutPath;
  this->coaddOutFile = ap->coaddOutFile;
  this->nRealizations = ap->nRealizations;
  this->noisePath = ap->noisePath;
  this->avgNoiseHistFile = ap->avgNoiseHistFile;
  this->avgNoisePsdFile = ap->avgNoisePsdFile;
  this->nNoiseMapsPerObs = ap->nNoiseMapsPerObs;
  this->despikeSigma = ap->despikeSigma;
  this->approximateWeights = ap->approximateWeights;
  this->azelMap = ap->azelMap;
  this->lowpassFilterKnee = ap->lowpassFilterKnee;
  this->subsampleFrequency = ap->subsampleFrequency;
  this->timeOffset = ap->timeOffset;
  this->timeChunk = ap->timeChunk;
  this->cutTurnArounds = ap->cutTurnArounds;
  this->cutSamplesAtEndOfScans = ap->cutSamplesAtEndOfScans;
  this->cutStd = ap->cutStd;
  this->neigToCut=ap->neigToCut;
  this->cleanPixelSize = ap->cleanPixelSize;
  this->order = ap->order;
  this->cleanStripe = ap->cleanStripe;
  this->stripeMethod.assign(ap->stripeMethod);
  this->controlChunk = ap->controlChunk;
  this->resample = ap->resample;
  this->bsOffset[0] = ap->bsOffset[0];
  this->bsOffset[1] = ap->bsOffset[1];
  this->pixelSize = ap->pixelSize;
  this->observatory = ap->observatory;
  this->timeVarName = ap->timeVarName;
  this->nThreads = ap-> nThreads;
  this->saveTimeStreams = ap->saveTimeStreams;
  this->tOrder = ap->tOrder;
  if (ap->simParams != NULL)
	  this->simParams = new SimParams(ap->simParams);
  else
	  this->simParams = NULL;
  this->doSubtract = ap->doSubtract;
  this->subtractFile = ap->subtractFile;
  this->subtractPath = ap->subtractPath;
  this->templateFile = ap->templateFile;
  this->subtractFirst = ap->subtractFirst;

  this->kernelName = ap->kernelName;
  this->coreR0 = ap->coreR0;
  this->coreP = ap->coreP;
  this->coreAxisRatio = ap->coreAxisRatio;


  ///post reduction values (added 2-7-13 NB)
  this->beamSize = ap->beamSize;
  this->covCut = ap->covCut;
  this->snglSourceWin = ap->snglSourceWin;
  this->sourceSigma = ap->sourceSigma;
  this->negativeToo = ap->negativeToo;
  this->mapNegative = ap->mapNegative;
  this->sfCentroidSources = ap->sfCentroidSources;
  this->sfFitGaussians = ap->sfFitGaussians;
  this->psSideLength = ap->psSideLength;
  this->nFluxBins = ap->nFluxBins;
  this->nSynthSources = ap->nSynthSources;
  this->minFlux = ap->minFlux;
  this->maxFlux = ap->maxFlux;
  this->recovS2N = ap->recovS2N;
  this->recovRadius = ap->recovRadius;
  this->coverageThreshold = ap->coverageThreshold;

  //Post-analysis stuff
  this->postReductionAnalysis = ap->postReductionAnalysis;
  this->calcCompleteness = ap->calcCompleteness;
  this->gaussianTemplate = ap->gaussianTemplate;
  this->highpassOnly = ap->highpassOnly;
  this->findSources = ap->findSources;
  this->gaussianTemplateFWHM = ap->gaussianTemplateFWHM;
  this->lowpassOnly = ap->lowpassOnly;
  this->normalizeErrors = ap->normalizeErrors;
  this->wienerFilter = ap->wienerFilter;
  
  if (this->nFiles <=0){
    cerr<<"AnalParams(): Copy constructor cannot continue without analyzing files."<<endl<< "Program ending"<<endl;
    exit(-1);
  }
  
  this->fileList = new string [nFiles];
  this->bstatList = new string [nFiles];
  this->mapFileList = new string [nFiles];
  for (int fix=0; fix<this->nFiles; fix++){
    this->fileList[fix]=ap->fileList[fix];
    this->bstatList[fix]= ap->bstatList[fix];
    this->mapFileList[fix]=ap->mapFileList[fix];
  }
  
  this->masterGridJ2000 = new double [2];
  this->masterGridJ2000_init = new double [2];
  
  this->masterGridJ2000[0] = ap->masterGridJ2000[0];
  this->masterGridJ2000[1] = ap->masterGridJ2000[1];
  this->masterGridJ2000_init[0] = ap->masterGridJ2000_init[0];
  this->masterGridJ2000_init[1] = ap->masterGridJ2000_init[1];

  //each thread should have its own generator seeded in its own way.
  //I'll use the original thread's generator to see the new ones.
  long seed = ap->macanaRandom->uniformDeviate(0,1e8);
  this->macanaRandom = new GslRandom(seed);
  
}
//----------------------------- o ---------------------------------------

///Resets the relevant AnalParams values for new data file
/** Here we reset all relevant AnalParams values for each raw data
    file that we are reducing.  It's vital that this be reviewed
    periodically so that we don't forget to reset any values.
**/
bool AnalParams::setDataFile(int index)
{
  dataFile = fileList[index].c_str();
  bolostatsFile = bstatList[index].c_str();
  mapFile = mapFileList[index].c_str();
  if(beammapping == 1){
    outBeammapInfo = outBeammapInfoList[index].c_str();
    if(writeBeammapToNcdf){
      outBeammapNcdf = outBeammapNcdfList[index].c_str();
    }
    beammapSourceFlux = beammapSourceFluxList[index]; 
  }
  bsOffset[0] = bsOffsetList[0][index];
  bsOffset[1] = bsOffsetList[1][index];

  //reset the mastergrid
  setMasterGridJ2000(masterGridJ2000_init[0], masterGridJ2000_init[1]);  

  cerr << endl;
  cerr << "Now starting: " << dataFile << "." << endl;

  //determine the observatory
  determineObservatory();
  cerr << "Observatory is " << observatory << endl;
  
  return 1;
}


//----------------------------- o ---------------------------------------


bool AnalParams::determineObservatory()
{
  //is it the LMT?
#pragma omp critical (dataio)
{
  NcError ncerror(NcError::silent_nonfatal);
  NcFile ncfid(dataFile, NcFile::ReadOnly);
  NcVar* timeVar = ncfid.get_var("Data.AztecBackend.time");
  if(timeVar){
    cerr << "AnalParams::determinObservatory(): ";
    cerr << "This is an LMT data file." << endl;
    observatory.assign("LMT");
    timeVarName.assign("Data.AztecBackend.time");
  }
  else{
  //must check to see if it's ASTE or JCMT
	  NcVar* randomAsteSig = ncfid.get_var("aste_windtime");
	  if(randomAsteSig){
		observatory.assign("ASTE");
		timeVarName.assign("time");
	  }else{
		observatory.assign("JCMT");
		timeVarName.assign("time");
	  }
  }
  ncfid.close();
}
  return 1;
}


//----------------------------- o ---------------------------------------


void AnalParams::throwXmlError(string p)
{
  cerr << "XML error in " << apXml << ": " << p << endl;
  exit(-1);
}


//----------------------------- o ---------------------------------------


string AnalParams::getObservatory()
{
  return observatory;
}


//----------------------------- o ---------------------------------------


bool AnalParams::getMapIndividualObservations()
{
  return mapIndividualObservations;
}


//----------------------------- o ---------------------------------------


bool AnalParams::getCoaddObservations()
{
  return coaddObservations;
}


//----------------------------- o ---------------------------------------


bool AnalParams::getFitCoadditionToGaussian()
{
  return fitCoadditionToGaussian;
}


//----------------------------- o ---------------------------------------


bool AnalParams::getProduceNoiseMaps()
{
  return produceNoiseMaps;
}


//----------------------------- o ---------------------------------------


bool AnalParams:: getApplyWienerFilter()
{
  return applyWienerFilter;
}


//----------------------------- o ---------------------------------------


string AnalParams::getTimeVarName()
{
  return timeVarName;
}


//----------------------------- o ---------------------------------------


const char* AnalParams::getDataFile()
{
  return dataFile;
}


//----------------------------- o ---------------------------------------


const char* AnalParams::getBolostatsFile()
{
  return bolostatsFile;
}


//----------------------------- o ---------------------------------------


const char* AnalParams::getOutBeammapInfo()
{
  return outBeammapInfo;
}


//----------------------------- o ---------------------------------------


const char* AnalParams::getOutBeammapNcdf()
{
  return outBeammapNcdf;
}


//----------------------------- o ---------------------------------------


double AnalParams::getDespikeSigma()
{
  return despikeSigma;
}


//----------------------------- o ---------------------------------------


bool AnalParams::getApproximateWeights()
{
  return approximateWeights;
}


//----------------------------- o ---------------------------------------


int AnalParams::getAzelMap()
{
  return azelMap;
}


//----------------------------- o ---------------------------------------


int AnalParams::getWriteBeammapToNcdf()
{
  return writeBeammapToNcdf;
}


//----------------------------- o ---------------------------------------

int AnalParams::getCleanIterationCap()
{
  return cleanIterationCap;
}

//----------------------------- o ---------------------------------------

double AnalParams::getCleanIterationCutoff()
{
  return cleanIterationCutoff;
}

//----------------------------- o ---------------------------------------


double AnalParams::getLowpassFilterKnee()
{
  return lowpassFilterKnee;
}


//----------------------------- o ---------------------------------------


double AnalParams::getSubsampleFrequency()
{
  return subsampleFrequency;
}


//----------------------------- o ---------------------------------------


double AnalParams::getTimeOffset()
{
  return timeOffset;
}

//----------------------------- o ---------------------------------------


double AnalParams::getTimeChunk()
{
  return timeChunk;
}

//----------------------------- o ---------------------------------------

bool AnalParams::getCutTurnArounds()
{
  return cutTurnArounds;
}

//----------------------------- o ---------------------------------------

int AnalParams::getCutSamplesAtEndOfScans()
{
  return cutSamplesAtEndOfScans;
}
//----------------------------- o ---------------------------------------


double AnalParams::getCutStd()
{
  return cutStd;
}

//----------------------------- o ---------------------------------------


int AnalParams::getNeigToCut()
{
  return neigToCut;
}

//----------------------------- o ---------------------------------------

double AnalParams::getCleanPixelSize()
{
  return cleanPixelSize;
}

//----------------------------- o ---------------------------------------

int AnalParams::getOrder()
{
  return order;
}

//----------------------------- o ---------------------------------------

double AnalParams::getCleanStripe()
{
  return cleanStripe;
}

//----------------------------- o ---------------------------------------

double AnalParams::getControlChunk()
{
  return controlChunk;
}

//----------------------------- o ---------------------------------------
int AnalParams::getNThreads()
{
  return nThreads;
}

//----------------------------- o ---------------------------------------

bool AnalParams::getSaveTimestreams()
{
 return saveTimeStreams;
}


//----------------------------- o ---------------------------------------

double* AnalParams::getBsOffset()
{
  return bsOffset;
}

//----------------------------- o ---------------------------------------

double* AnalParams::getMasterGridJ2000()
{
  return masterGridJ2000;
}


//----------------------------- o ---------------------------------------

bool AnalParams::setMasterGridJ2000(double ra, double dec)
{
  masterGridJ2000[0] = ra;
  masterGridJ2000[1] = dec;
  return 1;
}


//----------------------------- o ---------------------------------------

double AnalParams::getPixelSize()
{
  return pixelSize;
}


//----------------------------- o ---------------------------------------

double AnalParams::getCoverageThreshold()
{
  return coverageThreshold;
}


//----------------------------- o ---------------------------------------

int AnalParams::getNFiles()
{
  return nFiles;
}


//----------------------------- o ---------------------------------------

string AnalParams::getMapFile()
{
  string tmp(mapFile);
  return tmp;
}


//----------------------------- o ---------------------------------------

string AnalParams::getMapFileList(int i)
{
  return mapFileList[i];
}


//----------------------------- o ---------------------------------------

string AnalParams::getSourceName()
{
  return sourceName;
}


//----------------------------- o ---------------------------------------

bool AnalParams::setSourceName(string name)
{
  sourceName.assign(name);
  return 1;
}


//----------------------------- o ---------------------------------------


string AnalParams::getCoaddOutFile()
{
  return coaddOutFile;
}


//----------------------------- o ---------------------------------------


string AnalParams::getNoisePath()
{
  return noisePath;
}


//----------------------------- o ---------------------------------------



int AnalParams::getNRealizations()
{
  return nRealizations;
}

//----------------------------- o ---------------------------------------



string AnalParams::getAvgNoiseHistFile()
{
  return avgNoiseHistFile;
}


//----------------------------- o ---------------------------------------



string AnalParams::getAvgNoisePsdFile()
{
  return avgNoisePsdFile;
}


//----------------------------- o ---------------------------------------


double AnalParams::getBeamSize()
{
  return beamSize;
}


//----------------------------- o ---------------------------------------


double AnalParams::getCovCut()
{
  return covCut;
}


//----------------------------- o ---------------------------------------


double AnalParams::getSnglSourceWin()
{
  return snglSourceWin;
}


//----------------------------- o ---------------------------------------


double AnalParams::getSourceSigma()
{
  return sourceSigma;
}


//----------------------------- o ---------------------------------------


bool AnalParams::getNegativeToo()
{
  return negativeToo;
}


//----------------------------- o ---------------------------------------


bool AnalParams::getMapNegative()
{
  return mapNegative;
}


//----------------------------- o ---------------------------------------


bool AnalParams::getSFCentroidSources()
{
  return sfCentroidSources;
}


//----------------------------- o ---------------------------------------


bool AnalParams::getSFFitGaussians()
{
  return sfFitGaussians;
}


//----------------------------- o ---------------------------------------


int AnalParams::getPsSideLength()
{
  return psSideLength;
}

//----------------------------- o ---------------------------------------


int AnalParams::getNFluxBins()
{
  return nFluxBins;
}


//----------------------------- o ---------------------------------------


int AnalParams::getNSynthSources()
{
  return nSynthSources;
}


//----------------------------- o ---------------------------------------


double AnalParams::getMinFlux()
{
  return minFlux;
}


//----------------------------- o ---------------------------------------


double AnalParams::getMaxFlux()
{
  return maxFlux;
}


//----------------------------- o ---------------------------------------


double AnalParams::getRecovS2N()
{
  return recovS2N;
}


//----------------------------- o ---------------------------------------


double AnalParams::getRecovRadius()
{
  return recovRadius;
}


//----------------------------- o ---------------------------------------


double AnalParams::getBeammapSourceFlux()
{
  return beammapSourceFlux;
}


//----------------------------- o ---------------------------------------


bool AnalParams::getPostReductionAnalysis()
{
  return postReductionAnalysis;
}


//----------------------------- o ---------------------------------------


bool AnalParams::getGaussianTemplate()
{
  return gaussianTemplate;
}


//----------------------------- o ---------------------------------------


double AnalParams::getGaussianTemplateFWHM()
{
  return gaussianTemplateFWHM;
}


//----------------------------- o ---------------------------------------


bool AnalParams::getLowpassOnly()
{
  return lowpassOnly;
}


//----------------------------- o ---------------------------------------


bool AnalParams::getHighpassOnly()
{
  return highpassOnly;
}


//----------------------------- o ---------------------------------------


bool AnalParams::getNormalizeErrors()
{
  return normalizeErrors;
}

void AnalParams::setOrder (int Order){
	this->order=Order;
}
void AnalParams::setNeig2cut (int n2cut){
	this->neigToCut=n2cut;
}

int AnalParams::getResample(){
	return resample;
}

int AnalParams::getTOrder(){
	return tOrder;
}

SimParams* AnalParams::getSimParams(){
	return this->simParams;
}

string AnalParams::getTemplateFile(){
	return this->templateFile;
}

string AnalParams::getSubtractFile() {
	return subtractPath+"/"+subtractFile;
}

bool AnalParams::getDoSubtract() {
	return doSubtract;
}

int AnalParams::getNNoiseMapsPerObs(){
  return nNoiseMapsPerObs;
}

bool AnalParams::getFindSources(){
  return findSources;
}

bool AnalParams::getCalcCompleteness(){
  return calcCompleteness;
}

bool AnalParams::getAltKernel(){
  return altKernel;
}

string AnalParams::getKernelName(){
  return kernelName;
}

double AnalParams::getCoreR0(){
  return coreR0;
}

double AnalParams::getCoreP(){
  return coreP;
}

double AnalParams::getCoreAxisRatio(){
  return coreAxisRatio;
}

//----------------------------- o ---------------------------------------

AnalParams::~AnalParams()
{
  delete [] masterGridJ2000;
  delete [] masterGridJ2000_init;
  delete [] fileList;
  delete [] mapFileList;
  delete macanaRandom;
  if (simParams !=NULL)
	  delete simParams;
}
