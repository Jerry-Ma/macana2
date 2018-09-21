#include <netcdfcpp.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <omp.h>
using namespace std;

#include "nr3.h"
#include "tinyxml2.h"
#include "AnalParams.h"
#include "Source.h"
#include "Telescope.h"
#include "TimePlace.h"
#include <gsl/gsl_statistics.h>
#include "astron_utilities.h"
#include "vector_utilities.h"

///Telescope constructor
/** Here we determine the observatory, collect the pointing data, 
    signal condition it, and generate some delta-source pointings.
**/
Telescope::Telescope(AnalParams* anal, TimePlace* tP, Source* source)
{
  ap = anal;
  dataFile = ap->getDataFile();
  string timeVarName = ap->getTimeVarName();
  string observatory = ap->getObservatory();
  timePlace = tP;

  //collect the UTDate
  utDate = timePlace->utDate;



  long int * edges = NULL;
  //set the number of available samples
#pragma omp critical (dataio)
{
  NcFile ncfid(dataFile, NcFile::ReadOnly);
  NcVar* timeVar = ncfid.get_var(timeVarName.c_str());
  edges = timeVar->edges();
  nSamples = edges[0]*edges[1];
  ncfid.close();
}
  delete [] edges;

  //get the pointing signals
  if(observatory.compare("LMT") == 0){
    cerr << "Telescope:Telescope(): This is an LMT data file." << endl;
    getLMTPointing();
  } else if(observatory.compare("ASTE") == 0){
    cerr << "Telescope:Telescope(): This is an ASTE data file." << endl;
    getASTEPointing();
  } else {
    cerr << "JCMT files not yet supported." << endl;
    exit(1);
  }

  //calculate az/el physical pointing
  cerr << "Telescope(): generating az/el delta-source pointings." << endl;
  absToPhysHorPointing(source);

  //define the scans
  cerr << "Telescope(): defining scans." << endl;
  double timeChunk=ap->getTimeChunk();
  defineScans(ap, timeChunk);

  //collect the M2 offsets
  getM2Offsets();

  //collect the M1 Zernike coefficients
  getM1Zernike();

  //collect the obsNum
  getObsNumFromFile();

  //get the Project ID
  getProjectIDFromFile();
}


//----------------------------- o ---------------------------------------

///collect and signal condition the LMT pointing signals
/** The LMT pointing is easy since it's all right there in the data
    file.  Collect it, do some signal conditioning to remove dropouts,
    NaNs, and spikes, correct the rollovers if needed, and align with
    the detector signals.
**/
bool Telescope::getLMTPointing()
{
  //in this case we just grab the pointing data directly from the file
  paraAngle.resize(nSamples);
  //getTelescopeValues("Data.AztecBackend.ParAng", &paraAngle[0]);
  getTelescopeValues("Data.AztecBackend.ParAng", &paraAngle[0]);
  for(int i=0;i<nSamples;i++) paraAngle[i] = TWO_PI/2.-paraAngle[i];
  hTelRa.resize(nSamples);
  getTelescopeValues("Data.AztecBackend.SourceRaAct", &hTelRa[0]);
  hTelDec.resize(nSamples);
  getTelescopeValues("Data.AztecBackend.SourceDecAct", &hTelDec[0]);
  hTelAzAct.resize(nSamples);
  getTelescopeValues("Data.AztecBackend.TelAzAct", &hTelAzAct[0]);
  hTelElAct.resize(nSamples);
  getTelescopeValues("Data.AztecBackend.TelElAct", &hTelElAct[0]);
  hTelAzDes.resize(nSamples);
  getTelescopeValues("Data.AztecBackend.TelAzDes", &hTelAzDes[0]);
  hTelElDes.resize(nSamples);
  getTelescopeValues("Data.AztecBackend.TelElDes", &hTelElDes[0]);
  hTelAzCor.resize(nSamples);
  getTelescopeValues("Data.AztecBackend.TelAzCor", &hTelAzCor[0]);
  hTelElCor.resize(nSamples);
  getTelescopeValues("Data.AztecBackend.TelElCor", &hTelElCor[0]);

  //calcParallacticAngle();
  cerr<<"Got Values from file"<<endl;
  //apply levels 0 signal conditioning, this removes NaNs and dropouts
  signalCondition0();

  //now take out spikes that are not NaNs
  signalCondition1();

  //apply the rollover correction
  fixRollover();

  //do the alignment with detector signals
  alignWithDetectors();

  //fetch the az and el user offsets
#pragma omp critical (dataio)
{
  NcFile ncfid(dataFile, NcFile::ReadOnly);
  NcVar* azo=ncfid.get_var("Header.PointModel.AzUserOff");
  NcVar* elo=ncfid.get_var("Header.PointModel.ElUserOff");
  azUserOff = azo->as_double(0);
  elUserOff = elo->as_double(0);
  ncfid.close();
}

  //calculate the parallactic angle used to make bolo pointings
//  calcParallacticAngle();
//  for(int i=0;i<nSamples;i++) paraAngle[i] = TWO_PI/2.-paraAngle[i];

  return 1;

}


//----------------------------- o ---------------------------------------

///collect and signal condition the ASTE pointing signals
/**For ASTE we have some work to do as less is precalculated.
   We calculate them here along with doing the same signal
   conditioning.
 **/
bool Telescope::getASTEPointing()
{
  //get Az and El from file, they are called Act but are actually Des
  //swap the sign on the error signals so that the sign convention 
  //matches the LMT's Cor signals.
  hTelAzDes.resize(nSamples);
  getTelescopeValues("aste_azact", &hTelAzDes[0]);
  for(int i=0;i<nSamples;i++) hTelAzDes[i] /= DEG_RAD;
  hTelElDes.resize(nSamples);
  getTelescopeValues("aste_elact", &hTelElDes[0]);
  for(int i=0;i<nSamples;i++) hTelElDes[i] /= DEG_RAD;
  hTelAzCor.resize(nSamples);
  getTelescopeValues("aste_azerr", &hTelAzCor[0]);
  hTelElCor.resize(nSamples);
  getTelescopeValues("aste_elerr", &hTelElCor[0]);  

  //the other signals need storage space with zeros
  paraAngle.assign(nSamples,1.);
  hTelAzAct.assign(nSamples,1.);
  hTelElAct.assign(nSamples,1.);
  hTelRa.assign(nSamples,1.);
  hTelDec.assign(nSamples,1.);

  //signal conditioning
  signalCondition0();
  signalCondition1();
  fixRollover();
  alignWithDetectors();

  //the Act signals will be the error corrected signals
  for(int i=0;i<nSamples;i++){
    hTelAzCor[i] *= (-1.);
    hTelElCor[i] *= (-1.);
    hTelAzAct[i] = hTelAzDes[i] - hTelAzCor[i]/hTelElDes[i];
    hTelElAct[i] = hTelElDes[i] - hTelElCor[i];
  }

  //generate the Ra/Dec signals from the Az/El Signals
  azElToRaDec2000(timePlace, &hTelAzAct[0], &hTelElAct[0],
		  &hTelRa[0], &hTelDec[0], nSamples);


  //calculate the parallactic angle used to make bolo pointings
  calcParallacticAngle();

  //set the user offsets to 0 - these will never be used
  azUserOff = 0.;
  elUserOff = 0.;

  return 1;
}


//----------------------------- o ---------------------------------------

///pulls the M2 offsets from the data file
bool Telescope::getM2Offsets()
{
  //open the netcdf file
#pragma omp critical (dataio)
{
  if (ap->getObservatory().compare("LMT")==0){
	  NcFile ncfid(dataFile, NcFile::ReadOnly);
	  NcVar* var=ncfid.get_var("Header.M2.XReq");

	  if (!var){
		  cerr<<"Telescope(): Warning no Header data in NetCDF file."<<endl;
		  exit(-1);
	  }

	  M2XReq = var->as_double(0);
	  var=ncfid.get_var("Header.M2.YReq");
	  M2YReq = var->as_double(0);
	  var=ncfid.get_var("Header.M2.ZReq");
	  M2ZReq = var->as_double(0);
	  ncfid.close();
  }else{ ///TODO: Set apropiate ASTE/JCMT variables if necessary
	  M2XReq = 0.0;
	  M2YReq = 0.0;
	  M2ZReq = 0.0;
  }
}
  return 1;
}


//----------------------------- o ---------------------------------------

///pulls the M2 offsets from the data file
bool Telescope::getM1Zernike()
{
  //open the netcdf file
#pragma omp critical (dataio)
{
  if (ap->getObservatory().compare("LMT")==0 && utDate>=2014.82235488798){
	  NcFile ncfid(dataFile, NcFile::ReadOnly);
	  NcVar* var=ncfid.get_var("Header.M1.ZernikeC");
	  if (!var){
		  cerr<<"Telescope(): Warning no Zernike Header data in NetCDF file."<<endl;
		  exit(-1);
	  }
	  M1ZernikeC0 = var->as_double(0);
	  M1ZernikeC1 = var->as_double(1);
	  M1ZernikeC2 = var->as_double(2);
	  ncfid.close();
  }else{ ///TODO: Set apropiate ASTE/JCMT variables if necessary
	  M1ZernikeC0 = 0.;
	  M1ZernikeC1 = 0.;
	  M1ZernikeC2 = 0.;
  }
}
  return 1;

}

//----------------------------- o ---------------------------------------

///pulls the obsnum from the data file
bool Telescope::getObsNumFromFile()
{
  //open the netcdf file
#pragma omp critical (dataio)
{
  if (ap->getObservatory().compare("LMT")==0){
	  NcFile ncfid(dataFile, NcFile::ReadOnly);
	  NcVar* var=ncfid.get_var("Header.Dcs.ObsNum");

	  if (!var){
		  cerr<<"Telescope(): Warning no Header data in NetCDF file."<<endl;
		  exit(-1);
	  }
	  obsNum = var->as_int(0);
	  ncfid.close();
  }else{ ///TODO: Set apropiate ASTE/JCMT variables if necessary
          obsNum = 0;
  }
}
  return 1;

}

//----------------------------- o ---------------------------------------

///pulls the project ID from the data file
bool Telescope::getProjectIDFromFile()
{
  //open the netcdf file
#pragma omp critical (dataio)
{
  if (ap->getObservatory().compare("LMT")==0 && utDate>=2014.75){
    NcFile ncfid(dataFile, NcFile::ReadOnly);
    NcVar* var=ncfid.get_var("Header.Dcs.ProjectId");
    
    if (!var){
      cerr<<"Telescope(): Warning no Header data in NetCDF file."<<endl;
      exit(-1);
    }
    projectID.assign(var->as_string(0));
    projectID.erase(projectID.find_last_not_of(" \n\r\t")+1);
    ncfid.close();
  }else{ ///TODO: Set apropiate ASTE/JCMT variables if necessary
    projectID.assign("");
  }
}

  return 1;
}


//------------------------------o----------------------------------------


///returns vector of Telescope values
bool Telescope::getTelescopeValues(const char* sigName, double *tData)
{

  //open the netcdf file
#pragma omp critical (dataio)
{
  NcFile ncfid(dataFile, NcFile::ReadOnly);
  NcVar* tele=ncfid.get_var(sigName);
  
  //need size of detector array
  long int* edges = tele->edges();
  
  //here's a vector for the data
  MatDoub hm(edges[0],edges[1]);
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
  ncfid.close();
}
  return 1;
}


//----------------------------- o ---------------------------------------


///take out NaNs and dropouts from dataset
bool Telescope::signalCondition0()
{
  //need to remove bad values a la aztec_makeLMTPointing.pro
  removeDropouts(&paraAngle[0],nSamples);
  removeDropouts(&hTelRa[0],nSamples);
  removeDropouts(&hTelDec[0],nSamples);
  removeDropouts(&hTelAzAct[0],nSamples);
  removeDropouts(&hTelElAct[0],nSamples);
  removeDropouts(&hTelAzDes[0],nSamples);
  removeDropouts(&hTelElDes[0],nSamples);
  removeDropouts(&hTelAzCor[0],nSamples);
  removeDropouts(&hTelElCor[0],nSamples);
  deNanSignal(&paraAngle[0],nSamples);
  deNanSignal(&hTelDec[0],nSamples);
  deNanSignal(&hTelAzAct[0],nSamples);
  deNanSignal(&hTelElAct[0],nSamples);
  deNanSignal(&hTelAzDes[0],nSamples);
  deNanSignal(&hTelElDes[0],nSamples);
  deNanSignal(&hTelAzCor[0],nSamples);
  deNanSignal(&hTelElCor[0],nSamples);
  return 1;
}


//----------------------------- o ---------------------------------------


///reset any signals that are out of bounds to mean of adjacent samples
bool Telescope::signalCondition1()
{
  double low=1.e-9;
  double high=5.*PI;
  correctOutOfBounds(&paraAngle[0], low, high, nSamples);
  correctOutOfBounds(&hTelRa[0], low, high, nSamples);
  correctOutOfBounds(&hTelDec[0], low, high, nSamples);
  correctOutOfBounds(&hTelAzAct[0], low, high, nSamples);
  correctOutOfBounds(&hTelElAct[0], low, high, nSamples);
  correctOutOfBounds(&hTelAzDes[0], low, high, nSamples);
  correctOutOfBounds(&hTelElDes[0], low, high, nSamples);
  correctOutOfBounds(&hTelAzCor[0], low, high, nSamples);
  correctOutOfBounds(&hTelElCor[0], low, high, nSamples);
  return 1;
}


//----------------------------- o ---------------------------------------

///apply a rollover correction to signals with periodic boundary conditions
bool Telescope::fixRollover()
{
  correctRollover(&hTelRa[0], PI, 1.99*PI, 2.0*PI, nSamples);
  correctRollover(&hTelDec[0], PI, 1.99*PI, 2.0*PI, nSamples);
  correctRollover(&hTelAzAct[0], PI, 1.99*PI, 2.0*PI, nSamples);
  correctRollover(&hTelElAct[0], PI, 1.99*PI, 2.0*PI, nSamples);
  correctRollover(&hTelAzDes[0], PI, 1.99*PI, 2.0*PI, nSamples);
  correctRollover(&hTelElDes[0], PI, 1.99*PI, 2.0*PI, nSamples);
  correctRollover(&hTelAzCor[0], PI, 1.99*PI, 2.0*PI, nSamples);
  correctRollover(&hTelElCor[0], PI, 1.99*PI, 2.0*PI, nSamples);
  return 1;
}

//----------------------------- o ---------------------------------------


///interpolate the pointing signals
bool Telescope::alignWithDetectors()
{
  //pointing signals start out referenced to telUtc
  //here we interpolate them to detUtc to align them
  //with the detector signals

  //here's the model for this
  VecDoub xx(timePlace->nUnique),yy(timePlace->nUnique);   //original signals

  cerr << "Telescope::alignWithDetectors() with " << nSamples;
  cerr << " samples." << endl;

  //deal with the xx vector first
  for(int i=0;i<timePlace->nUnique;i++)
    xx[i] = timePlace->telUtcUnique[i]+timePlace->timeOffset;

  //interpolate paraAngle
  for(int i=0;i<timePlace->nUnique;i++) 
    yy[i]=paraAngle[timePlace->locUnique[i]];
  interpolateLinear(xx, yy, timePlace->nUnique,
		    &timePlace->detUtc[0], &paraAngle[0], nSamples);

  //interpolate hTelRa
  for(int i=0;i<timePlace->nUnique;i++) 
    yy[i]=hTelRa[timePlace->locUnique[i]];
  interpolateLinear(xx, yy, timePlace->nUnique,
		    &timePlace->detUtc[0], &hTelRa[0], nSamples);

  //interpolate hTelDec
  for(int i=0;i<timePlace->nUnique;i++) 
    yy[i]=hTelDec[timePlace->locUnique[i]];
  interpolateLinear(xx, yy, timePlace->nUnique,
		    &timePlace->detUtc[0], &hTelDec[0], nSamples);

  //interpolate hTelAzAct
  for(int i=0;i<timePlace->nUnique;i++) 
    yy[i]=hTelAzAct[timePlace->locUnique[i]];
  interpolateLinear(xx, yy, timePlace->nUnique,
		    &timePlace->detUtc[0], &hTelAzAct[0], nSamples);

  //interpolate hTelElAct
  for(int i=0;i<timePlace->nUnique;i++) 
    yy[i]=hTelElAct[timePlace->locUnique[i]];
  interpolateLinear(xx, yy, timePlace->nUnique,
		    &timePlace->detUtc[0], &hTelElAct[0], nSamples);

  //interpolate hTelAzDes
  for(int i=0;i<timePlace->nUnique;i++) 
    yy[i]=hTelAzDes[timePlace->locUnique[i]];
  interpolateLinear(xx, yy, timePlace->nUnique,
		    &timePlace->detUtc[0], &hTelAzDes[0], nSamples);

  //interpolate hTelElDes
  for(int i=0;i<timePlace->nUnique;i++) 
    yy[i]=hTelElDes[timePlace->locUnique[i]];
  interpolateLinear(xx, yy, timePlace->nUnique,
		    &timePlace->detUtc[0], &hTelElDes[0], nSamples);

  //interpolate hTelAzCor
  for(int i=0;i<timePlace->nUnique;i++) 
    yy[i]=hTelAzCor[timePlace->locUnique[i]];
  interpolateLinear(xx, yy, timePlace->nUnique,
		    &timePlace->detUtc[0], &hTelAzCor[0], nSamples);

  //interpolate hTelElCor
  for(int i=0;i<timePlace->nUnique;i++) 
    yy[i]=hTelElCor[timePlace->locUnique[i]];
  interpolateLinear(xx, yy, timePlace->nUnique,
		    &timePlace->detUtc[0], &hTelElCor[0], nSamples);


  cerr << "Telescope::alignWithDetectors(): done." << endl;
  return 1;
}


//----------------------------- o ---------------------------------------


///define the scans in the dataFile
/** Here we use the same techniques found in the idl utilities.  I've
    checked a handful of cases that the scan boundaries agree but this
    could always be checked further.
**/
bool Telescope::defineScans(AnalParams* ap, double timeChunk)
{
  string timeVarName = ap->getTimeVarName();
  string observatory = ap->getObservatory();

  bool cutTurnArounds = ap->getCutTurnArounds();
  int cutSamplesAtEndOfScans = ap->getCutSamplesAtEndOfScans();

  //some of this is observatory-dependent
  bool LMT=0;
  bool ASTE=0;
  if(observatory.compare("LMT") == 0){
    LMT=1;
  } else if(observatory.compare("ASTE") == 0){
    ASTE=1;
  }

  //need samplerate from dataFile
  int samplerate;
  string dimname;
#pragma omp critical (dataio)
{
  NcFile ncfid(dataFile, NcFile::ReadOnly);
  if(LMT) dimname.assign("Data.AztecBackend.time_xlen");
  if(ASTE) dimname.assign("mfpersf");
  NcDim* mfPerSf=ncfid.get_dim(dimname.c_str());
  samplerate = mfPerSf->size();
  ncfid.close();
}
  //the goal is to pack up the turning array
  VecBool turning(nSamples);
  VecBool flagT2(nSamples,1);

  //set default turning state (telescope-specific)
  //turning: LMT
  if(LMT){
    //for LMT, scans are defined using the Hold signal from the
    //dataFile unless timeChunk is nonzero.  Get it anyway.  
    VecDoub holdDouble(nSamples);
    getTelescopeValues("Data.AztecBackend.Hold", &holdDouble[0]);  
    
    //recast hold from double to int
    VecInt hold(nSamples);
    for(int i=0;i<nSamples;i++) hold[i] = int(holdDouble[i]);

    //the telescope is turning at end of scan when hold&8=1
    for(int i=0;i<nSamples;i++) turning[i] = (hold[i]&8);
  }
    
  //turning: ASTE
  if(ASTE){
    //for ASTE we use the scantype signal
    VecDoub scantype(nSamples);
    getTelescopeValues("aste_scantype", &scantype[0]); 

    //now mimic turning signal
    for(int i=0;i<nSamples;i++) turning[i] = (scantype[i] != 3) ? 1 : 0;
  }  


  // Raster scan mode for timeChunk=0.
  if (timeChunk == 0){

    cerr << "Telescope::defineScans(): ";
    cerr << "timeChunk is zero. Raster scan mode enabled." << endl;

    // cutTurnArounds forces a search for stops direction changes that
    // the hold signal is not catching in LMT raster scans
    if (cutTurnArounds){
      VecDoub velMag(nSamples);
      VecDoub azVel(nSamples);
      VecDoub elVel(nSamples);
      VecDoub angles(nSamples);
      VecDoub absAngles(nSamples);

      //Use velocities to solve hold signal anomalies
      azVel = *derivate(timePlace->detUtc,hTelAzPhys);
      elVel = *derivate(timePlace->detUtc,hTelElPhys);
     
      for(int i=0;i<nSamples;i++) 
        velMag[i] = sqrt(pow(azVel[i],2)+pow(elVel[i],2));

      angles =  *getAngle(hTelAzPhys,hTelElPhys);
      for(int i=0;i<nSamples;i++) absAngles[i] = abs(angles[i]);

      double medVelMag = median(velMag);
      double medAbsAng = median(absAngles);

      double madVelMag = medabsdev(velMag);
      double madAbsAng = medabsdev(absAngles);

      for(int i=0;i<nSamples;i++) turning[i] = ((velMag[i] < (medVelMag-3.0*madVelMag)) || (velMag[i] < (0.2*medVelMag)) || (abs(absAngles[i]-medAbsAng) > (10*madAbsAng)) );

      double sum = 0; 
      for (long icta=0; icta<nSamples; icta++){
        if (turning[icta]==1){
          long begin = icta - cutSamplesAtEndOfScans;
          long end = icta + cutSamplesAtEndOfScans;
          if (begin < 0 ) begin = 0;
          if (end > (nSamples)) end = nSamples;
          sum = 1; // reset, and ensure we keep edge points of big turn blocks
          for (long jcta=begin; jcta <end; jcta++) sum += turning[jcta];
          if(sum > 0.5*(end-begin)) for (long jcta=begin; jcta <end; jcta++) flagT2[jcta]=0;
        }
      }

      for (long icta = 0; icta < nSamples; icta++) turning[icta] = !flagT2[icta];
  
    }

    //count up the number of scans and keep track of the scan number
    nScans=0;
    for(int i=1;i<nSamples;i++){
      if(turning[i]-turning[i-1]>0) nScans++;
    }

    if (turning[nSamples-1] == 0) nScans++;

    cerr << "Telescope::defineScans(): Found Scans: " << nScans << endl;


    //now build the scanIndex array
    //1st row is index of start of scan
    //2nd row is index of end of scan
    scanIndex.resize(2,nScans);
    int counter=-1;
    if(!turning[0]){
      scanIndex[0][0] = 1;
      counter++;
    }

    for(int i=1;i<nSamples;i++){
      if(turning[i]-turning[i-1] < 0){
        //this is the beginning of the next scan
        //intentionally sacrificing a sample to keep
        //things identical to aztec_defineScans.pro
        counter++;
        scanIndex[0][counter] = i+1;
      }
      if(turning[i]-turning[i-1] > 0){
        //one sample ago was the end of the scan
        scanIndex[1][counter] = i-1;
      }
    }
    //the last sample is the end of the last scan
    scanIndex[1][nScans-1] = nSamples-1;


  // Lissajous mode for timeChunk>0.
  } else {

    cerr << "Telescope::defineScans(): ";
    cerr << "timeChunk is nonzero. Lissajous mode enabled." << endl;

    //in this case we've got to take the bounds and redefine the scans
    //start with the first and last scan index
    //reduce by one scan since we seem to have a problem at the moment

    cerr << "Telescope::defineScans(): " << 
      "Redefining scans to use all data." << endl;
    int firstScanI=0;//scanIndex[0][0];
    int lastScanI =nSamples-1;//scanIndex[1][nScans-1];
    double period = floor(timeChunk*samplerate);
    int newNScans = floor((lastScanI-firstScanI+1)*1./period);
    scanIndex.resize(2,newNScans);
    for(int i=0;i<newNScans;i++){
      scanIndex[0][i] = i*period + firstScanI;
      scanIndex[1][i] = scanIndex[0][i] + period - 1;
    }
    nScans = newNScans;
    cerr << "Telescope::defineScans(): " << 
      "Now we have " << nScans << " scans." << endl;
  }

  //do a final check of scan length.  If a scan is
  //less than 2s of data then delete it
  VecBool obsFlags(nSamples,1);
  checkScanLengths(obsFlags,samplerate);    

  cerr << "Telescope::defineScans(): Found " << nScans << " scans." << endl;
  return 1;
}
  
void Telescope::checkScanLengths(VecBool obsFlags, int samplerate)
{
    MatInt tmpSI(2,nScans); 
    tmpSI = scanIndex; 

    //do a final check of scan length.  If a scan is
    //less than 2s of data then delete it
    int nBadScans=0;
    int sum=0;
    VecBool isBadScan(nScans);
    for(int i=0;i<nScans;i++)
    {
      sum=0;
      for(int j=tmpSI[0][i];j<(tmpSI[1][i]+1);j++) sum+=obsFlags[j];
      if(sum < 2.*samplerate){
        nBadScans++;
        isBadScan[i]=1;
      } else isBadScan[i]=0;
    }

    if(nBadScans > 0){
      cerr << "Telescope::checkScanLengths():";
      cerr << nBadScans<<" scans with duration less than 2 seconds detected...";
      cerr << "Ignoring." << endl;
    }

    //pack up the scan indices
    int c=0;
    scanIndex.resize(2,nScans-nBadScans);
    for(int i=0;i<nScans;i++)
    {
      if(!isBadScan[i]){
        scanIndex[0][c] = tmpSI[0][i];
        scanIndex[1][c] = tmpSI[1][i];
        c++;
      }
    }
    nScans = nScans-nBadScans;
 
}

//----------------------------- o ---------------------------------------

///transform Ra/Dec pointing signals from absolute to physical coordinates
/**I have checked this against aztec_abs2phys.pro and it looks good.
 **/
bool Telescope::absToPhysEqPointing()
{
  //create the phys pointing arrays
  hTelRaPhys.resize(nSamples);
  hTelDecPhys.resize(nSamples);

  //get the centerRa and centerDec of the map from ap
  double* pmg;
  pmg = ap->getMasterGridJ2000();
  double centerRa = pmg[0];
  double centerDec = pmg[1];

  //call the tangential projection 
  absToPhys(&hTelRa[0], &hTelDec[0], centerRa, centerDec, 
	    &hTelRaPhys[0], &hTelDecPhys[0], nSamples);

  return 1;
}


//----------------------------- o ---------------------------------------

///transform the AzEl pointing signals from absolute to physical coordinates
bool Telescope::absToPhysHorPointing(Source* source)
{
  //create the phys pointing arrays
  hTelAzPhys.resize(nSamples);
  hTelElPhys.resize(nSamples);

  //get the source az and el
  double* saz = &source->hSourceAz[0];
  double* sel = &source->hSourceEl[0];

  //let's try this though it may not be robust for all cases
  //The problem here is that the source az is limited to [0,2pi]
  //while the telescope az ranges from [0,2.5pi].
  for(int i=0;i<nSamples;i++){
    if((hTelAzAct[i]-saz[i]) > 0.9*TWO_PI) hTelAzAct[i] -= TWO_PI;
  }

  for(int i=0;i<nSamples;i++){
    hTelAzPhys[i] = (hTelAzAct[i] - saz[i])*cos(hTelElDes[i]) - hTelAzCor[i];
    hTelElPhys[i] = hTelElAct[i] - sel[i] - hTelElCor[i];
  }

  return 1;
}


//----------------------------- o ---------------------------------------

///Calculates parallactic angle with respect to telescope boresight
/**Note that this isn't really the parallactic angle since it also
   rotates through the change in epoch to J2000.
**/
bool Telescope::calcParallacticAngle()
{
  VecDoub daz(nSamples,100./3600.*TWO_PI/360.);
  VecDoub del(nSamples,0.);
  VecDoub sratmp(nSamples);
  VecDoub sdectmp(nSamples);
  VecDoub naz(nSamples);
  VecDoub nel(nSamples);
  VecDoub nra(nSamples);
  VecDoub ndec(nSamples);
  VecDoub dra(nSamples);
  VecDoub ddec(nSamples);
  double* saz = &hTelAzAct[0];
  double* sel = &hTelElAct[0];
  physToAbs(&daz[0],&del[0],&saz[0],&sel[0],&naz[0],&nel[0],nSamples);
  azElToRaDec2000(timePlace, &naz[0], &nel[0], &nra[0], &ndec[0], nSamples);
  for(int i=0;i<nSamples;i++)
    absToPhys(&nra[i], &ndec[i], hTelRa[i], hTelDec[i], &dra[i], &ddec[i], 1);
  for(int i=0;i<nSamples;i++) paraAngle[i] = atan2(ddec[i],dra[i]);
  return 1;
}

//----------------------------- o ---------------------------------------

int Telescope::getObsNum()
{
  return obsNum;
}

//----------------------------- o ---------------------------------------

double Telescope::getUTDate()
{
  return utDate;
}

//----------------------------- o ---------------------------------------

double Telescope::getM2XReq()
{
  return M2XReq;
}

//----------------------------- o ---------------------------------------

double Telescope::getM2YReq()
{
  return M2YReq;
}

//----------------------------- o ---------------------------------------

double Telescope::getM2ZReq()
{
  return M2ZReq;
}

//----------------------------- o ---------------------------------------

double Telescope::getAzUserOff()
{
  return azUserOff;
}

//----------------------------- o ---------------------------------------

double Telescope::getElUserOff()
{
  return elUserOff;
}

//----------------------------- o ---------------------------------------

double Telescope::getM1ZernikeC0()
{
  return M1ZernikeC0;
}

//----------------------------- o ---------------------------------------

double Telescope::getM1ZernikeC1()
{
  return M1ZernikeC1;
}

//----------------------------- o ---------------------------------------

double Telescope::getM1ZernikeC2()
{
  return M1ZernikeC2;
}


//----------------------------- o ---------------------------------------

string Telescope::getProjectID()
{
  return projectID;
}

//----------------------------- o ---------------------------------------

Telescope::~Telescope()
{

}
