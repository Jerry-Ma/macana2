#include <netcdfcpp.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;

#include "nr3.h"
#include "tinyxml2.h"
#include "AnalParams.h"
#include "TimePlace.h"
#include "Telescope.h"
#include "astron_utilities.h"
#include "vector_utilities.h"


///TimePlace constructor
/**
   Determines observatory and then collects timing signals.
 **/
TimePlace::TimePlace(AnalParams* anal)
{
  //preliminaries
  ap = anal;
  dataFile = ap->getDataFile();
  string timeVarName = ap->getTimeVarName();
  string observatory = ap->getObservatory();

  //open the file and get dimensions 
#pragma omp critical (dataio)
{
  NcFile ncfid(dataFile, NcFile::ReadOnly);
  NcVar* timeVar = ncfid.get_var(timeVarName.c_str());
  long int* edges = timeVar->edges();
  nSamples = edges[0]*edges[1];
  ncfid.close();
  delete [] edges;
}
  //get the timing signals
  if(observatory.compare("LMT") == 0){
    getLMTTiming();
  } else if(observatory.compare("ASTE") == 0){
    getASTETiming();
  } else {
    cerr << "JCMT files not yet supported." << endl;
    exit(1);
  }

}


//----------------------------- o ---------------------------------------

///collects LMT timing data and signal conditions it
/**Collects the LMT timing data and then applies the standard
   signal conditioning.  Also applies the time offset which is a fixed
   latency (at least it's assumed to be fixed!) between the detector
   signals and the telescope signals.  Much of this code mimics what
   is found in aztec_makeLMTPointing.pro.
**/
bool TimePlace::getLMTTiming()
{
  //get the single-valued signals
#pragma omp critical (dataio)
{
  NcFile ncfid(dataFile, NcFile::ReadOnly);
  NcVar* longVar = ncfid.get_var("Header.TimePlace.ObsLongitude");
  longVar->get(&longitude,1,0,0,0,0);
  NcVar* latVar = ncfid.get_var("Header.TimePlace.ObsLatitude");
  latVar->get(&latitude,1,0,0,0,0);
  NcVar* elevVar = ncfid.get_var("Header.TimePlace.ObsElevation");
  elevVar->get(&elevation,1,0,0,0,0);
  NcVar* utdVar = ncfid.get_var("Header.TimePlace.UTDate");
  utdVar->get(&utDate,1,0,0,0,0);
  ncfid.close();
}
  //now run through variables and fetch them
  telUtc.resize(nSamples);
  getTimePlaceValues("Data.AztecBackend.TelUtc", &telUtc[0]);
  detUtc.resize(nSamples);
  getTimePlaceValues("Data.AztecBackend.AztecUtc", &detUtc[0]);
  lst.resize(nSamples);
  getTimePlaceValues("Data.AztecBackend.TelLst", &lst[0]);

  //apply levels 0 signal conditioning
  signalCondition0();

  //LMT's UTC is in radians so convert to hours
  for(int i=0;i<nSamples;i++){
    telUtc[i] = telUtc[i]*24./2./PI;
  }

  //now take out spikes that are not NaNs
  signalCondition1();

  //HACK Alert!!!!
  //this should not be necessary, but it is
  //assume aztec and lmt utc signals are offset from
  //each other by a fixed amount and apply a correction
  double meanUtc=0.;
  for(int i=0;i<nSamples;i++) meanUtc += (detUtc[i]-telUtc[i]);
  meanUtc /= nSamples;
  for(int i=0;i<nSamples;i++) detUtc[i] = detUtc[i]-meanUtc;

  //apply the rollover correction to the utc and lst signals
  fixRollover();

  //find the unique values of telUtc
  findUniqueTelUtc();

  //if tOff is nonzero then use it for timeOffset, otherwise use the
  //default value of 0.125s
  double tOff = ap->getTimeOffset();
  if(tOff != 0.) timeOffset = tOff/3600.; else timeOffset = 0.125/3600.;
  cerr << "TimePlace(): Using timeOffset=" << timeOffset*3600. << "s" << endl;

  //do the alignment of the lst signal
  alignWithDetectors();

  //Julian date calculation
  double modJulDate;
  juldayFromUTDate(utDate, &julDate, &modJulDate);  

  //go backwards and get the calendar date from the juldate
  caldateFromJulDate(julDate, &year, &month, &day);

  return 1;
}


//----------------------------- o ---------------------------------------

///collects ASTE timing data and signal conditions it
/**Collects the ASTE timing data and then applies the standard
   signal conditioning.  Also applies the time offset which is a fixed
   latency (at least it's assumed to be fixed!) between the detector
   signals and the telescope signals.  Much of this code mimics what
   is found in aztec_makeASTEPointing.pro.
**/
bool TimePlace::getASTETiming()
{
  //get the single-valued signals
  longitude = -67.703304 / DEG_RAD;
  latitude = -22.971657 / DEG_RAD;
  elevation = 4861.9;
  string date;
#pragma omp critical (dataio)
{
  NcFile ncfid(dataFile, NcFile::ReadOnly);
  NcAtt* dateAtt = ncfid.get_att("data_file");
  char *tmpDate = dateAtt->as_string(0);
  date.assign(tmpDate);
  delete [] tmpDate;
  delete dateAtt;
}
  int hour, minute, second;
  string smonth, sday, syear, shour, sminute, ssecond;
  int found=date.find_first_of("_",1);
  smonth = date.substr(1,found-1);
  date = date.substr(found+1,date.length());
  month = atoi(smonth.c_str());

  found=date.find_first_of("_",1);
  sday = date.substr(0,found);
  date = date.substr(found+1,date.length());
  day = atoi(sday.c_str());

  found=date.find_first_of("_",1);
  syear = date.substr(0,found);
  date = date.substr(found+1,date.length());
  year = atoi(syear.c_str());

  found=date.find_first_of("_",1);
  shour = date.substr(0,found);
  date = date.substr(found+1,date.length());
  hour = atoi(shour.c_str());

  found=date.find_first_of("_",1);
  sminute = date.substr(0,found);
  date = date.substr(found+1,date.length());
  minute = atoi(sminute.c_str());

  found=date.find_first_of("_",1);
  ssecond = date.substr(0,date.length()-4);
  date = date.substr(found+1,date.length());
  second = atoi(ssecond.c_str());

  //set the UT Date
  utDate = epoch(month, day, year);
  //now run through variables and fetch them
  telUtc.resize(nSamples);
  getTimePlaceValues("aste_utc", &telUtc[0]);
  detUtc.resize(nSamples);
  getTimePlaceValues("aztec_utc", &detUtc[0]);
  lst.resize(nSamples);
  getTimePlaceValues("aste_lst", &lst[0]);
  for(int i=0;i<nSamples;i++) lst[i] = lst[i]/24.*TWO_PI;

  //apply levels 0 signal conditioning
  signalCondition0();

  //now take out spikes that are not NaNs
  signalCondition1();

  //apply the rollover correction to the utc and lst signals
  fixRollover();
  //find the unique values of telUtc
  findUniqueTelUtc();

  //if tOff is nonzero then use it for timeOffset, otherwise use the
  //default value of 0.125s
  double tOff = ap->getTimeOffset();
  if(tOff != 0.) timeOffset = tOff; else timeOffset = 0.167/3600.;
  cerr << "TimePlace(): Using timeOffset=" << timeOffset*3600. << "s" << endl;

  //do the alignment of the lst signal
  alignWithDetectors();

  //Julian date calculation
  double modJulDate;
  juldayFromUTDate(utDate, &julDate, &modJulDate);  

  return 1;
}


//----------------------------- o ---------------------------------------


///take out NaNs and dropouts from dataset
bool TimePlace::signalCondition0()
{
  //need to remove bad values a la aztec_makeLMTPointing.pro
  removeDropouts(&telUtc[0],nSamples);
  deNanSignal(&telUtc[0],nSamples);
  removeDropouts(&detUtc[0],nSamples);
  deNanSignal(&detUtc[0],nSamples);
  removeDropouts(&lst[0],nSamples);
  deNanSignal(&lst[0],nSamples);
  return 1;
}


//----------------------------- o ---------------------------------------


///correct signals that appear out of bounds (ie, totally unreasonable)
bool TimePlace::signalCondition1()
{
  double low=1.e-9;
  double high=30;
  correctOutOfBounds(&telUtc[0], low, high, nSamples);
  return 1;
}


//----------------------------- o ---------------------------------------

///apply rollover correction to signals with periodic boundary conditions
bool TimePlace::fixRollover()
{
  correctRollover(&telUtc[0], 10., 23., 24., nSamples);
  correctRollover(&detUtc[0], 10., 23., 24., nSamples);
  correctRollover(&lst[0], PI, 1.99*PI, TWO_PI, nSamples);
  return 1;
}


//----------------------------- o ---------------------------------------


///returns vector of TimePlace values
bool TimePlace::getTimePlaceValues(const char* sigName, double *tData)
{
  //open the netcdf file

#pragma omp critical (dataio)
{
  NcFile ncfid(dataFile, NcFile::ReadOnly);
  NcVar* timer=ncfid.get_var(sigName);
  
  //need size of detector array
  long int* edges = timer->edges();
  
  //here's a vector for the data
  MatDoub hm(*edges,*(edges+1));
  Doub *phm = &hm[0][0];
  
  //fetch the data from the file and put it in something pointing at hm
  NcBool test;
  test = timer->get(phm,*edges,*(edges+1),0,0,0);

  //repack the data into vector form
  int i;
  int j;
  for(i=0;i<*edges;i++){
    for(j=0;j<*(edges+1);j++){
      tData[*(edges+1)*i+j]=hm[i][j];
    }
  }

  //cleanup
  delete [] edges;
  ncfid.close();
}
  return 1;
}


//----------------------------- o ---------------------------------------


///find the unique values of telUtc
bool TimePlace::findUniqueTelUtc()
{
  //telUtc should be increasing or constant, so let's test for this first
  for(int i=1;i<nSamples;i++){
    if(telUtc[i] < telUtc[i-1]){
      cerr << "TimePlace::findUniqueTelUtc(): ";
      cerr << "telUtc is not constant or increasing." << endl;
      cerr << "Measuring glitch ..." << endl;

      //measure the glitch.  If it's less than 32 samples we will
      //ignore it for now
      int bad=1;
      int nglitch=0;
      for(int j=0;j<32;j++) if(telUtc[i+j] > telUtc[i-1]){
	  bad=0;
	  nglitch=j+1;
	  break;
	}
      if(bad){
	cerr << "Warning, there may be large glitches in telUtc ";
	cerr << "in this data file." << endl;
      }else{
	cerr << "UTC glitch is " << nglitch << " samples long.  ";
	cerr << "This is less than 32 samples.  Ignoring it." << endl;
	cerr << "Consider throwing out this data file." << endl;
      }
    }
  }
  
  //count up the unique values
  int count=1;
  for(int i=1;i<nSamples;i++) if(telUtc[i] > telUtc[i-1]) count++;

  //build the new arrays
  telUtcUnique.resize(count);
  locUnique.resize(count);
  telUtcUnique[0] = telUtc[0];
  locUnique[0] = 0;
  int counter=1;
  for(int i=1;i<nSamples;i++) if(telUtc[i] > telUtc[i-1]){
      telUtcUnique[counter] = telUtc[i];
      locUnique[counter] = i;
      counter++;
    }
  nUnique = count;

  cerr << "TimePlace::findUniqueTelUtc(): ";
  cerr << nUnique << " unique samples." << endl;

  return 1;
}


//----------------------------- o ---------------------------------------


///interpolate the timing signals
bool TimePlace::alignWithDetectors()
{
  //timing signals start out referenced to telUtcUnique
  //here we interpolate them to detUtc to align them
  //with the detector signals
  //we also have to add an offset to the time to best align them

  //here's the model for this
  VecDoub xx(nUnique),yy(nUnique);       //original signals

  cerr << "TimePlace::alignWithDetectors() with " << nSamples;
  cerr << " samples." << endl;

  //lst
  for(int i=0;i<nUnique;i++){
    xx[i] = telUtcUnique[i]+timeOffset;
    yy[i] = lst[locUnique[i]];
  }
  interpolateLinear(xx, yy, nUnique,
		    &detUtc[0], &lst[0], nSamples);

  cerr << "TimePlace::alignWithDetectors(): done." << endl;

  return 1;
}


//----------------------------- o ---------------------------------------


TimePlace::~TimePlace()
{

}
