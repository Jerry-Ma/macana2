#include <netcdfcpp.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;

#include "nr3.h"
#include "tinyxml2.h"
#include "Source.h"
#include "TimePlace.h"
#include "astron_utilities.h"
#include "vector_utilities.h"

///Source constructor
/**Sets up the observatory (LMT or ASTE) and other initialization.
 **/
Source::Source(AnalParams* anal, TimePlace* tp)
{
  ap = anal;
  dataFile = ap->getDataFile();
  string timeVarName = ap->getTimeVarName();
  string observatory = ap->getObservatory();
  timePlace = tp;
#pragma omp critical (dataio)
  {
	  //set the number of available samples
	  NcFile ncfid(dataFile, NcFile::ReadOnly);
	  NcVar* timeVar = ncfid.get_var(timeVarName.c_str());
	  long int* edges = timeVar->edges();
	  nSamples = edges[0]*edges[1];
	  //cleanup
	  delete [] edges;



	  ncfid.close();
  }
  //get the source data
  if(observatory.compare("LMT") == 0){
    getLMTSourceData();
  } else if(observatory.compare("ASTE") == 0){
    getASTESourceData();
  } else {
    cerr << "JCMT files not yet supported." << endl;
    exit(1);
  }
  //set the source name in ap for future use
  ap->setSourceName(name);

}


//----------------------------- o ---------------------------------------

///Fetches LMT source data from file and does some signal conditioning.
/**Gets the LMT source data from the raw data file and performs
   signal conditioning to remove dropouts, NaNs and spikes. 
   Alignment with the detectors also ocurrs here.
**/
bool Source::getLMTSourceData()
{
  //open the datafile
#pragma omp critical (dataio)
	{
		NcFile ncfid(dataFile, NcFile::ReadOnly);

		//get the source Name
		NcVar* sourceNameVar = ncfid.get_var("Header.Source.SourceName");
		char* tmpStr;
		tmpStr = sourceNameVar->as_string(0);
		name.assign(tmpStr);
		delete [] tmpStr;
		name.erase(name.find(" "));
		cerr << "Source:: source name is " << name << endl;

		double* mg = ap->getMasterGridJ2000();
		if(mg[0] == 0. && mg[1] == 0.){
			//mastergrid to be set by source data from datafile
			//J2000 data comes from netcdf file header
			double sra, sdec;
			NcVar* sraVar=ncfid.get_var("Header.Source.Ra");
			sraVar->get(&sra,1,0,0,0,0);
			NcVar* sdecVar=ncfid.get_var("Header.Source.Dec");
			sdecVar->get(&sdec,1,0,0,0,0);
			ap->setMasterGridJ2000(sra,sdec);
		}
		ncfid.close();
	}
		//get the rest of the source signals
	hSourceAz.resize(nSamples);
	getSourceValues("Data.AztecBackend.SourceAz", &hSourceAz[0]);
	hSourceEl.resize(nSamples);
	getSourceValues("Data.AztecBackend.SourceEl", &hSourceEl[0]);
	hSourceRa.resize(nSamples);
	getSourceValues("Data.AztecBackend.SourceRa", &hSourceRa[0]);
	hSourceDec.resize(nSamples);
	getSourceValues("Data.AztecBackend.SourceDec", &hSourceDec[0]);



  //apply levels 0 signal conditioning, this removes NaNs and dropouts
  signalCondition0();
  //TODO::
  //now take out spikes that are not NaNs
  signalCondition1();

  //apply the rollover correction
  fixRollover();

  //and finally, do the alignment with detector signals
  alignWithDetectors();

  //set the AnalParams.mastergrid if it is zeros

  
  //close the dataFile
  return 1;
}


//----------------------------- o ---------------------------------------

///Fetches ASTE source data from file and does some signal conditioning.
/**Gets the ASTE source data from the raw data file and performs
   signal conditioning to remove dropouts, NaNs and spikes. 
   Alignment with the detectors also ocurrs here.
   Unlike the LMT version, the source Az and El must be calcualted
   here.
   \todo - update the time estimates for Novas.
**/
bool Source::getASTESourceData()
{
	//open the datafile
#pragma omp critical (dataio)
	{
		NcFile ncfid(dataFile, NcFile::ReadOnly);

		//get the source Name
		//source_name = "\"1310+323\""
		NcAtt* nameAtt = ncfid.get_att("source_name");
		char *tmpstr = NULL;
		tmpstr = nameAtt->as_string(0);
		name.assign(tmpstr);
		delete [] tmpstr;
		//name = name.substr(2,name.length()-3);
		cerr << "Source:: source name is " << name << endl;

		//use ra/dec from datafile attribute to generate
		//source ra/dec
		NcAtt* raAtt = ncfid.get_att("aste_header_SrcPosX");
		tmpstr = raAtt->as_string(0);
		string rastr = tmpstr;
		delete [] tmpstr;
		rastr = rastr.substr(0,rastr.length());
		double t1, t2, t3;
		parseRaDecString(rastr, &t1, &t2, &t3);
		double sra = (t1+t2/60.+t3/3600.)/24.*TWO_PI;
		NcAtt* decAtt = ncfid.get_att("aste_header_SrcPosY");
		tmpstr = decAtt->as_string(0);
		string decstr = tmpstr;
		delete [] tmpstr;
		decstr = decstr.substr(0,decstr.length());
		parseRaDecString(decstr, &t1, &t2, &t3);
		double st1 = (t1 < 0) ? -1. : 1;
		double sdec = st1*(abs(t1)+t2/60.+t3/3600.) / DEG_RAD;

		hSourceRa.resize(nSamples);
		hSourceDec.resize(nSamples);
		for(int i=0;i<nSamples;i++){
			hSourceRa[i] = sra;
			hSourceDec[i] = sdec;
		}

		//now use Novas to calculate Az/El for source
		hSourceAz.resize(nSamples);
		hSourceEl.resize(nSamples);
		raDecToAzEl(timePlace, &hSourceRa[0], &hSourceDec[0],
				&hSourceAz[0], &hSourceEl[0], nSamples);

		//set the AnalParams.mastergrid if it is zeros
		double* mg = ap->getMasterGridJ2000();
		if(mg[0] == 0. && mg[1] == 0.){
			//mastergrid to be set by source data from datafile
			//J2000 data comes from netcdf file header
			ap->setMasterGridJ2000(sra,sdec);
		}

		delete nameAtt;
		delete raAtt;
		delete decAtt;

		//close the dataFile
		ncfid.close();
	}
	return 1;
}


//----------------------------- o ---------------------------------------


///parses ra/dec string into component parts
bool Source::parseRaDecString(string ra, double* t1, double* t2, double* t3)
{
  string tmp;
  int found;
  found=ra.find_first_of(":",0);
  tmp = ra.substr(0,found);
  ra = ra.substr(found+1,ra.length());
  *t1 = atof(tmp.c_str());

  found=ra.find_first_of(":",0);
  tmp = ra.substr(0,found);
  ra = ra.substr(found+1,ra.length());
  *t2 = atof(tmp.c_str());

  tmp = ra.substr(0,ra.length());
  *t3 = atof(tmp.c_str());

  return 1;
}


//----------------------------- o ---------------------------------------


///returns vector of Source values
bool Source::getSourceValues(const char* sigName, double *tData)
{
	//open the netcdf file
#pragma omp critical (dataio)
	{
		NcFile ncfid(dataFile, NcFile::ReadOnly);
		NcVar* tele=ncfid.get_var(sigName);

		//need size of detector array
		long int* edges = tele->edges();

		//here's a vector for the data
		MatDoub hm(*edges,*(edges+1));
		Doub *phm = &hm[0][0];

		//fetch the data from the file and put it in something pointing at hm
		NcBool test;
		test = tele->get(phm,edges[0],edges[1],0,0,0);

		//repack the data into vector form
		int i;
		int j;
		for(i=0;i<edges[0];i++){
			for(j=0;j<edges[1];j++){
				tData[edges[1]*i+j]=hm[i][j];
			}
		}

		//cleanup
		delete [] edges;
	}
	return 1;
}


//----------------------------- o ---------------------------------------


///take out NaNs and dropouts from dataset
bool Source::signalCondition0()
{
  //need to remove bad values a la aztec_makeLMTPointing.pro
  removeDropouts(&hSourceAz[0],nSamples);
  removeDropouts(&hSourceEl[0],nSamples);
  removeDropouts(&hSourceRa[0],nSamples);
  removeDropouts(&hSourceDec[0],nSamples);
  deNanSignal(&hSourceAz[0],nSamples);
  deNanSignal(&hSourceEl[0],nSamples);
  deNanSignal(&hSourceRa[0],nSamples);
  deNanSignal(&hSourceDec[0],nSamples);
  return 1;
}


//----------------------------- o ---------------------------------------

///resets signals that are out of bounds to the mean of adjacent samples
bool Source::signalCondition1()
{
  double low=1.e-9;
  double high=5.*PI;
  correctOutOfBounds(&hSourceAz[0], low, high, nSamples);
  correctOutOfBounds(&hSourceEl[0], low, high, nSamples);
  correctOutOfBounds(&hSourceRa[0], low, high, nSamples);
  correctOutOfBounds(&hSourceDec[0], low, high, nSamples);
  return 1;
}


//----------------------------- o ---------------------------------------

///fixes rollover of signals with periodic boundary conditions
bool Source::fixRollover()
{
  correctRollover(&hSourceAz[0], PI, 1.99*PI, 2.0*PI, nSamples);
  correctRollover(&hSourceEl[0], PI, 1.99*PI, 2.0*PI, nSamples);
  correctRollover(&hSourceRa[0], PI, 1.99*PI, 2.0*PI, nSamples);
  correctRollover(&hSourceDec[0], PI, 1.99*PI, 2.0*PI, nSamples);
  return 1;
}


//----------------------------- o ---------------------------------------


///interpolates the pointing signals to align with detector signals
bool Source::alignWithDetectors()
{
  //pointing signals start out referenced to telUtc
  //here we interpolate them to detUtc to align them
  //with the detector signals

  //here's the model for this
  VecDoub xx(timePlace->nUnique),yy(timePlace->nUnique);   //original signals

  cerr << "Source::alignWithDetectors() with " << nSamples;
  cerr << " samples." << endl;

  //deal with the xx vector first
  for(int i=0;i<timePlace->nUnique;i++)
    xx[i] = timePlace->telUtcUnique[i]+timePlace->timeOffset;

  //interpolate hSourceAz
  for(int i=0;i<timePlace->nUnique;i++) 
    yy[i]=hSourceAz[timePlace->locUnique[i]];
  interpolateLinear(xx, yy, timePlace->nUnique,
		    &timePlace->detUtc[0], &hSourceAz[0], nSamples);

  //interpolate hSourceEl
  for(int i=0;i<timePlace->nUnique;i++) 
    yy[i]=hSourceEl[timePlace->locUnique[i]];
  interpolateLinear(xx, yy, timePlace->nUnique,
		    &timePlace->detUtc[0], &hSourceEl[0], nSamples);

  //interpolate hSourceRa
  for(int i=0;i<timePlace->nUnique;i++) 
    yy[i]=hSourceRa[timePlace->locUnique[i]];
  interpolateLinear(xx, yy, timePlace->nUnique,
		    &timePlace->detUtc[0], &hSourceRa[0], nSamples);

  //interpolate hSourceDec
  for(int i=0;i<timePlace->nUnique;i++) 
    yy[i]=hSourceDec[timePlace->locUnique[i]];
  interpolateLinear(xx, yy, timePlace->nUnique,
		    &timePlace->detUtc[0], &hSourceDec[0], nSamples);

  cerr << "Source::alignWithDetectors(): done." << endl;
  return 1;

}


//----------------------------- o ---------------------------------------


Source::~Source()
{

}
