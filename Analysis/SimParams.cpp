#include <sys/stat.h>

#include "SimParams.h"
#include "tinyxml2.h"

SimParams::SimParams(tinyxml2::XMLElement *apXml){

	  tinyxml2::XMLElement *xSim = apXml;
	  tinyxml2::XMLElement *xMap;
	  tinyxml2::XMLElement *xtmp;
	  struct stat buf;

	  xtmp = xSim->FirstChildElement("addSignal");
	  if(!xtmp) throwXmlError("SimParams()::addSignal not found.");
	  this->addSignal = atoi(xtmp->GetText());

	  xtmp = xSim->FirstChildElement("atmFreq");
	  if(!xtmp)
		  this->atmFreq = 0.;
	  else
		  this->atmFreq = atof(xtmp->GetText());

	  xtmp = xSim->FirstChildElement("resample");
	  if(!xtmp)
		  this->resample = 0.;
	  else{
		  this->resample= atof(xtmp->GetText());
		  cout<<"SimParams(): Requested to resample model signal with frequency: "<<this->resample<<"Hz"<<endl;
	  }

	  xtmp = xSim->FirstChildElement("atmSeed");
	  if(!xtmp)
		  this->atmSeed = 0;
	  else
		  this->atmSeed= atol(xtmp->GetText());

	  if (atmFreq !=0. && atmSeed != 0)
		  throwXmlError("SimParams Error. Options atmFreq and atmSeed are mutually exclusive. Please check your Analysis parameter file.");

	  xtmp = xSim->FirstChildElement("noiseChunk");
	  if(!xtmp)
		  this->noiseChunk = 0.;
	  else
		  this->noiseChunk = atof(xtmp->GetText());

	  xtmp = xSim->FirstChildElement("fluxFactor");
	  if(!xtmp) throwXmlError("SimParams()::fluxFactor addSignal not found.");
	  this->fluxFactor = atof(xtmp->GetText());


	  xtmp = xSim->FirstChildElement("simPath");
	  if(!xtmp) throwXmlError("SimParams()::simPath not found.");
	  this->mapPath = xtmp->GetText();
	  if(stat(this->mapPath.c_str(),&buf) == -1){
		  cerr << "SimParams: rawDataPath not file/directory: " << rawDataPath << endl;
		  exit(1);
	  }

	  xtmp = xSim->FirstChildElement("simFile");
	  if(!xtmp) throwXmlError("SimParams()::simFile not found.");
	  this->mapFile = xtmp->GetText();
	  if(stat(this->getMapFile().c_str(),&buf) == -1){
		  cerr << "SimParams: mapFile not file/directory: " << rawDataPath << endl;
		  exit(1);
	  }

}

SimParams::SimParams(SimParams *sp){
	apXml = sp->apXml;
	nFiles = sp->nFiles;                    ///<the number of files to be reduced
	rawDataPath = sp->rawDataPath;               ///<the path for the raw data
	bsPath=sp->bsPath;                    ///<the path for the corresponding bstats
	mapPath=sp->mapPath;                   ///<the path for the input simulated files
	mapFile=sp->mapFile;
	addSignal =sp->addSignal;       			//True add signal to current streams, False simulated signal only-->
	fluxFactor = sp->fluxFactor;
	atmFreq = sp->atmFreq;
	noiseChunk = sp->noiseChunk;
	atmSeed = sp->atmSeed;
	resample= sp->resample;
}


bool SimParams::getAddSignal()
{
	return addSignal;
}
double SimParams::getFluxFactor()
{
	return fluxFactor;
}
double SimParams::getAtmFreq()
{
	return atmFreq;
}
double SimParams::getNoiseChunk(){
	return noiseChunk;
}

long SimParams::getAtmSeed(){
	return atmSeed;
}

double SimParams::getResample(){
	return resample;
}

void SimParams::updateSeed(long fileNumber){
	if (atmSeed >0){
		atmSeed+=fileNumber;
	}
}

string SimParams::getMapFile(){
	string fullpath = this->mapPath + this->mapFile;
	return fullpath;
}

void SimParams::throwXmlError(string p)
{
  cerr << "XML error in " << apXml << ": " << p << endl;
  cerr.flush();
  exit(-1);
}


SimParams::~SimParams(){
}
