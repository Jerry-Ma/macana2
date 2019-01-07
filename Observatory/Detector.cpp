#include <netcdfcpp.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <fftw3.h>
using namespace std;

#include "nr3.h"
#include "tinyxml2.h"
#include "Detector.h"
#include "astron_utilities.h"
#include "vector_utilities.h"

#include "../port/include/port_timestream_sensitivity.h"

///default constructor
/** All that happens here is that all pointers are initilized.
    All the real action happens in initialize() below.
**/
Detector::Detector()
{
  //initialize the unset pointers to NULL
  isDespiked=0;  
  isCleaned=0;  
  isLowpassed=0;  
  isCalibrated=0;  
  isDownsampled=0;  
  isPointingGenerated=0;
}


//----------------------------- o ---------------------------------------


///constructor that fills detector attributes
/** Here is where all the detector data is extracted from both
    the netcdf data file as well as the bolostats.xml file.
    Warning: the ncdfLocation variable may be wrong for some
    data files and so is not used.
**/
void Detector::initialize(AnalParams* analParams, const char* dId)
{
  //get our own pointer to analysis parameters
  ap = analParams;
  dataFile = ap->getDataFile();
  bolostatsFile = ap->getBolostatsFile();

  //get ancillary data from bolostats xml file keying off dId
  tinyxml2::XMLDocument bolostats;
  bolostats.LoadFile(bolostatsFile);

  //xml elements
  tinyxml2::XMLElement* xdet = bolostats.FirstChildElement(dId);
  tinyxml2::XMLElement* xtmp;

  //name, goodFlag, psf info, etc.
  string param;
  string n;

  param = "name";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  n.assign(xtmp->GetText());

  param = "ncdf_location";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  ncdfLocation = atol(xtmp->GetText());

  param = "az_fwhm";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  beamSigAz = atof(xtmp->GetText())/2.3548;

  param = "el_fwhm";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  beamSigEl = atof(xtmp->GetText())/2.3548;

  param = "az_offset";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  azOffset = atof(xtmp->GetText());

  param = "el_offset";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  elOffset = atof(xtmp->GetText());

  param = "bologain";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  fcf = atof(xtmp->GetText());

  param = "bolosens";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  sensitivity = atof(xtmp->GetText());

  param = "offset_dc2tau";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  dc2tau[0] = atof(xtmp->GetText());

  param = "offset_err_dc2tau";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  dc2tauErr[0] = atof(xtmp->GetText());

  param = "slope_dc2tau";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  dc2tau[1] = atof(xtmp->GetText());

  param = "slope_err_dc2tau";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  dc2tauErr[1] = atof(xtmp->GetText());

  param = "quad_dc2tau";
  xtmp = xdet->FirstChildElement(param.c_str());
  if(xtmp==0) dc2tau[2]=0.; else dc2tau[2] = atof(xtmp->FirstChildElement("value")->GetText());

  param = "quad_err_dc2tau";
  xtmp = xdet->FirstChildElement(param.c_str());
  if(xtmp==0) dc2tauErr[2]=0.; else dc2tauErr[2] = atof(xtmp->FirstChildElement("value")->GetText());

  param = "offset_dc2responsivity";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  dc2responsivity[0] = atof(xtmp->GetText());

  param = "offset_err_dc2responsivity";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  dc2responsivityErr[0] = atof(xtmp->GetText());

  param = "slope_dc2responsivity";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  dc2responsivity[1] = atof(xtmp->GetText());

  param = "slope_err_dc2responsivity";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  dc2responsivityErr[1] = atof(xtmp->GetText());

  param = "goodflag";
  xtmp = xdet->FirstChildElement(param.c_str())->FirstChildElement("value");
  if(xtmp==0) throwXmlError(param);
  goodFlag = (strcmp(xtmp->GetText(),"1") == 0) ? 1 : 0;


  //strip the 'd' off of the detector ID and make it d.id as int
  int ntmp;
  ntmp = strlen(dId);
  char stmp[10];
  id = atol(strcpy(stmp,dId+1));
  int test;
//#pragma omp critical (dataio)
//{
  //open the file 
  NcFile ncfid(dataFile, NcFile::ReadOnly);
  NcVar* bolo=ncfid.get_var(n.c_str());
  if(!bolo->is_valid()){
    cerr << "Detector:: bolometer variable not fetched from ncfile." << endl;
    exit(1);
  }

  //assign the detector name
  //strcpy(name,bolo->name());
  name.assign(bolo->name());

  //number of samples and variable edges
  long int* edges = bolo->edges();
  nSamples = edges[0]*edges[1];
  rawSamplerate = edges[1];
  samplerate=rawSamplerate;
  delete [] edges;

  //create the detector value array and pack it
  hValues.resize(nSamples);
  getBoloValues(&hValues[0]);
  hSampleFlags.resize(nSamples);
  atmTemplate.resize(0);
  for(int i=0;i<nSamples;i++) hSampleFlags[i]=1;

  //get the detector dc value
  Doub dctmp=0;
  for(int i=0;i<nSamples;i++) dctmp+=hValues[i];
  dc = dctmp/nSamples;

  //estimate opacity and responsivity, initialize extinction to 0 to remind
  //us we've got to do this in main
  estimateTau();
  extinction=0;
  estimateResponsivity();

  //get the electronics gain
  NcError ncerror(NcError::silent_nonfatal);
  NcVar* fec1cmd=ncfid.get_var("Data.AztecBackend.fec1_cntl");
  if(!fec1cmd){
    fec1cmd=ncfid.get_var("fec1_cntl");
    if(!fec1cmd){
      cerr << "Can't find fec1_cntl in data file. Aborting." << endl;
      exit(1);
    }
  }
  long int* cmdEdge = fec1cmd->edges();
  VecInt fec1cmdVec(*cmdEdge);
  test = fec1cmd->get(&fec1cmdVec[0],cmdEdge[0],0,0,0,0);
  if(!test){
    cerr << "Detector:: Failed to get fec1_cntl: " << name << 
      " data from datafile" << endl;
    exit(1);
  }
  //take the value half-way through the data file
  int cmd = fec1cmdVec[*cmdEdge/2.];
  fecGain = cmdToGain(cmd);
  
  //cleanup
  delete [] cmdEdge;

  //close the ncfile just to be sure
  test = ncfid.close();
//}
  if(!test){
    cerr << "Detector:: Failed to close the netcdf datafile." << endl;
    exit(1);
  }
}


//----------------------------- o ---------------------------------------


///returns vector of detector values
bool Detector::getBoloValues(double *bData)
{
  //open the netcdf file
//#pragma omp critical (dataio)
//{
  NcFile ncfid(dataFile, NcFile::ReadOnly);
  NcVar* bolo=ncfid.get_var(name.c_str());
  
  //need size of detector array
  long int* edges = bolo->edges();
  
  //here's a vector for the data
  MatDoub hm(*edges,*(edges+1));
  Doub *phm = &hm[0][0];
  
  //fetch the data from the file and put it in something pointing at hm
  NcBool test;
  test = bolo->get(phm,*edges,*(edges+1),0,0,0);
  
  //repack the data into vector form
  int i;
  int j;
  for(i=0;i<*edges;i++){
    for(j=0;j<*(edges+1);j++){
      bData[*(edges+1)*i+j]=hm[i][j];
    }
  }

  //cleanup
  delete [] edges;
  ncfid.close();
//}

  return 1;
}


//----------------------------- o ---------------------------------------


///Despiking - detector-only steps of despiking.
/** Despiking.  Here we are executing the despiking steps that do
    not depend on the values of other bolometers: finding hte spikes
    and flagging the timestream values around the spikes.
    The despiking algorithm here is taken from the AzTEC despiking
    routines written by Jay Austermann but with my own modifications.
    Of course, I am attempting to preserve critical functionality and
    I have tried to make comments where I deviate substantially.
    The biggest difference now is that I despike the entire timestream
    at once rather than breaking it up by scans.
    This has been fairly thoroughly tested but could always used
    some additional comparisons to the IDL utilities.
    /todo Implement neighboring spike detection.
**/
bool Detector::despike(double nSigmaSpikes)
{
  //here we do all the despiking we can do without the 
  //use of other detectors.  
  //1) find the spikes
  //2) flag the timestream values around the spikes
  
  
  //FIND THE SPIKES
  int nSpikes=0;
  
  //calculate vector of adjacent signals and its moments
  Doub deltaMean=0;
  Doub deltaStddev=0;
  VecDoub delta(nSamples-1);
  Doub *pDelta = &delta[0];
  int nDelta = nSamples-1;
  for(int i=0;i<nDelta;i++) delta[i]=hValues[i+1]-hValues[i];
  deltaMean = mean(pDelta,nDelta);
  deltaStddev = stddev(pDelta,nDelta,deltaMean);

  
  //search the delta array for spikes that are nSigmaSpikes
  //times larger than the stddev
  for(int i=1;i<nDelta;i++){
    if(abs(delta[i]-deltaMean) > nSigmaSpikes*deltaStddev){
      hSampleFlags[i]=0;
      delta[i]=0.;
      nSpikes++;
    }
  }

  
  //do this recursively since big spikes skew the mean and stddev
  int newfound=1;
  bool newSpikes;
  while(newfound == 1){
    newSpikes=0;
    deltaMean = mean(pDelta,nDelta);
    deltaStddev = stddev(pDelta,nDelta,deltaMean);      
    for(int i=1;i<nDelta;i++){
      if(abs(delta[i]-deltaMean) > nSigmaSpikes*deltaStddev){
	hSampleFlags[i]=0;
	delta[i]=0.;
	nSpikes++;
	newSpikes=1;
      }
    } 
    if(!newSpikes) newfound=0;
  }
  
  /*
  //if there are more than 100 spikes for a detector then something
  //is very wrong
  if(nSpikes > 100){
    cerr << "Detector::despike(): More than 100 spikes found for";
    cerr << " detector: " << name << endl;
    cerr << "Setting this detector's goodflag to 0" << endl;
    goodFlag=0;
  }
  */

  //if there are multiple spikes in a window of 10 samples then
  //call it a single spike at the location of the middle of the window
  for(int i=0;i<nSamples;i++){
    if(hSampleFlags[i] == 0){
      int c=0;
      int indx=0;
      for(int j=1;j<=10 && i+j<nSamples;j++){
	indx = i+j;
	if(hSampleFlags[indx] == 0) c++;
      }
      if(c>0){
	//reduce the number nSpikes c times
	nSpikes = nSpikes - c;
	//fix hSampleFlags in the window 
	//(11/11/15 - GW - this next line incorrectly
	//used to start j with 0 instead of 1)
	for(int j=1;j<=10 && i+j<nSamples;j++){
	  indx = i+j;
	  hSampleFlags[indx] = 1;
	}
	hSampleFlags[i+5] = 0;
      }
      i = i+9;
    }
  }

  //FLAG THE SAMPLES AROUND ANY SPIKES
  if(nSpikes > 0){
    //flag the samples around the spikes
    //note that I'm assuming a fixed time constant for the detectors
    //here.  I am also not implementing the decay_samples feature
    //of aztec_despike.pro
    double timeConstant = 0.015;   //assumed time constant [seconds]

    //recount the spikes to avoid pathologies
    nSpikes = 0;
    for(int i=0;i<nSamples-1;i++) if(hSampleFlags[i] == 0) nSpikes++;

    //get a list of spike locations and det values
    VecInt spikeLoc(nSpikes);
    VecDoub spikeVals(nSpikes);
    int count=0;
    for(int i=0;i<nSamples-1;i++){
      if(hSampleFlags[i] == 0){
	spikeLoc[count]=i+1;
	spikeVals[count]=hValues[i+1]-hValues[i];
	count++;
      }
    }

    //find biggest spike-free window
    //first deal with the harder case of multiple spikes
    //remember that there are nSpikes+1 possible windows
    int winIndex0;
    int winIndex1;
    int winSize;
    if(nSpikes > 1){
      VecInt deltaSpikeLoc(nSpikes+1);
      deltaSpikeLoc[0] = spikeLoc[0]-0;  //distance to first spike
      deltaSpikeLoc[nSpikes] = nSamples-spikeLoc[nSpikes-1];
      for(int i=1;i<nSpikes;i++) deltaSpikeLoc[i] = spikeLoc[i]-spikeLoc[i-1];
      int mxWindow=deltaSpikeLoc[0];
      int mxWindowIndex=0;
      for(int i=1;i<=nSpikes;i++){
	if(deltaSpikeLoc[i] > mxWindow){
	  mxWindow = deltaSpikeLoc[i];
	  mxWindowIndex = i;
	}
      }
      //set the starting and ending indices for the window
      //leave a 2 second pad after the spike beginning the
      //window and a 1 second pad before the spike ending the window
      if(mxWindowIndex == 0){
    	  winIndex0=0;
    	  winIndex1=spikeLoc[1] - samplerate;
      }
      else{
    	  if(mxWindowIndex == nSpikes){
    		  winIndex0 = spikeLoc[nSpikes-1];
    		  winIndex1 = nSamples;
    	  } else {
    		  winIndex0=spikeLoc[mxWindowIndex-1] + 2*samplerate;
    		  winIndex1=spikeLoc[mxWindowIndex] - 1*samplerate;
    	  }
      }
      //do some error checking
      if(winIndex0 > winIndex1 || winIndex0 < 0 || winIndex0 > nSamples ||
	 winIndex1 < 0 ||  winIndex1 > nSamples){
	cerr << "Detector::despike(): detector:" << name << endl;
	cerr << "Either spikes everywhere or something else is"; 
	cerr << " terribly wrong." << endl;
	cerr << "Data file: " << dataFile << endl;
	cerr << "winIndex0 = " << winIndex0 << endl;
	cerr << "winIndex1 = " << winIndex1 << endl;
	cerr << "nSamples = " << nSamples << endl;
	for(int i=0;i<nSpikes;i++) 
	  cerr << "spikeLoc["<<i<<"] = " << spikeLoc[i] << endl;
	cerr << "Overall there are " << nSpikes << " spikes." << endl;
	exit(1);
      }
    } 
    else {
      //in this case there is only one spike
      if(nSamples-spikeLoc[0] > spikeLoc[0]){
	winIndex0 = spikeLoc[0] + 2*samplerate;
	winIndex1 = nSamples-1;
      }
      else{
	winIndex0 = 0;
	winIndex1 = spikeLoc[0] - 1*samplerate;
      }
    }
    //and after all that, let's limit the maximum window size to 10s
    if((winIndex1-winIndex0-1)/samplerate > 10.){
      winIndex1 = winIndex0+10*samplerate+1;
    }
    winSize = winIndex1-winIndex0-1;

    //make a sub-array with values from biggest spike-free window
    VecDoub subVals(winSize);
    for(int i=winIndex0;i<winIndex1-1;i++){
      //cerr << i-winIndex0 << "   " << i << endl;
      subVals[i-winIndex0]=hValues[i];
    }

    //make a boxcar smoothed copy of the subarray
    VecDoub smoothedSubVals(winSize);
    smooth(subVals,smoothedSubVals,10);

    //here is the estimate of the stardard deviation
    double sigest;
    for(int i=0;i<winSize;i++) subVals[i] -= smoothedSubVals[i];
    sigest = stddev(subVals);
    if(sigest < 1.e-8) sigest=1.e-8;


    //calculate for each spike the time it takes to decay to sigest
    VecDoub decayLength(nSpikes);
    for(int i=0;i<nSpikes;i++){
      decayLength[i] = -samplerate*timeConstant*log(abs(sigest/spikeVals[i]));
      if(decayLength[i] < 6) decayLength[i]=6;
      if(decayLength[i] > samplerate*10.){
	cerr << "Detector::despike(): Decay length is too int.";
	cerr << "  Something's broken in the despiker." << endl;
	exit(1);
      }
    }

    //now flag samples, 1 decayLength before and 2 decayLengths after
    //if not lowpassed, otherwise extend window by length of lowpass
    //filter/2 on either side of spike
    for(int i=0;i<nSpikes;i++){
      if(!isLowpassed){
	for(int j=-decayLength[i];j<2*decayLength[i];j++){
	  if(spikeLoc[i]+j >= 0 && spikeLoc[i]+j < nSamples)
	    hSampleFlags[spikeLoc[i]+j] = 0;
	}
      } else {
	for(int j=-(despikeWindow-1)/2;j<(despikeWindow-1)/2;j++){
	  if(spikeLoc[i]+j >= 0 && spikeLoc[i]+j < nSamples)
	    hSampleFlags[spikeLoc[i]+j] = 0;
	}
      }
    }
  }//if nSpikes gt 0 then flag samples around spikes
    
  isDespiked=1;
  return 1;
}


//----------------------------- o ---------------------------------------


//despiking self (gpu version)
bool Detector::despike(int nSigmaSpikes, int gpuId)
{
  //here we do all the despiking we can do without the 
  //use of other detectors.  
  //1) find the spikes
  //2) flag the timestream values around the spikes
  // TODO
  (void) nSigmaSpikes;
  (void) gpuId;
  
  isDespiked=1;
  return 1;
}


//----------------------------- o ---------------------------------------


///Lowpass - lowpass the detector timestreams
/** This is a simple lowpass filter being applied to the
    timestream samples.  The digital filter coefficients
    are calculated in the Array class since they are common
    across all detectors.  This routine is one of the slower
    ones in the group.
**/
bool Detector::lowpass(double* digFiltTerms, int nTerms)
{
  //now just do the convolution
  //I'm not going to lowpass the first or last nTerms of the detector
  //signals.  This only corresponds to 2s of data on either end.  I
  //will set those sampleFlags to 0.
  VecDoub storage(nSamples);
  double tmp;
  int nCoef=2*nTerms+1;
  int tt;
  int th=nSamples-nTerms;
  int t;
  int i;
 
  for(t=nTerms;t<th;t++){
    tmp=0.;
    tt = t-nTerms;
    for(i=0;i<nCoef;i++){
      tmp += hValues[tt+i]*digFiltTerms[i];
    }
    storage[t] = tmp;
  }

  for(i=0;i<nTerms;i++) hSampleFlags[i]=0;
  for(i=th;i<nSamples;i++) hSampleFlags[i]=0;
  
  //now replace hValues with the filtered values
  for(i=nTerms;i<th;i++) hValues[i] = storage[i];


  //the despiking window size should be set by the extent of
  //the filter
  despikeWindow = nCoef;

  isLowpassed=1;
  return 1;
}


//----------------------------- o ---------------------------------------


//lowpass the timestream signals (gpu version)
bool Detector::lowpass(int gpuId)
{
  //this is a digital FIR filter
  // TODO
  (void) gpuId;
  isLowpassed=1;
  return 1;
}


//----------------------------- o ---------------------------------------


//downsample the data
bool Detector::downsample(double desiredSamplerate)
{
  //downsample the data - this is trivial so can
  //be done on the cpu.  make sure to update
  //the data values and sampleflags pointers if needed
  
  //the desired samplerate may not be practical so use
  //the nearest best value
  //samplerate = best value;

  // TODO
  (void) desiredSamplerate;
  
  isDownsampled=1;
  return 1;
}


//----------------------------- o ---------------------------------------


//estimate opacity
bool Detector::estimateTau()
{
  if(dc2tau[0] < 9998){
    estimatedTau = dc2tau[2]*pow(dc,2) + dc2tau[1]*dc + dc2tau[0]; 
    return 1;
  }else{
    estimatedTau = -9999.;
    return 1;
  }
}


//----------------------------- o ---------------------------------------


//estimate extinction
//inputting tau here since array averaged tau is more accurate
//than single detector's estimatedTau.
//the following is copied from aztec_correct_extinction.pro
bool Detector::estimateExtinction(double tau)
{
  double tau265 = tau*1.67;
  extinction = exp(tau265);
  return 1;
}


//----------------------------- o ---------------------------------------


//estimate responsivity
bool Detector::estimateResponsivity()
{
  if(dc2responsivity[0] < 9998){
    responsivity = dc2responsivity[1]*dc + dc2responsivity[0];
  }else{
    cerr << "Detector::estimateResponsivity(): " << name;
    cerr << " does not have valid dc2responsivity conversion." << endl;
    exit(1);
  }
  return 1;
}


//----------------------------- o ---------------------------------------


//convert from input command to electronics gain
double Detector::cmdToGain(int cmd)
{
  if(cmd < 128 || cmd > 255){
    cerr << "Detector::cmdToGain(): cmd must be between 128 and 255: ";
    cerr << cmd << endl;
    //    exit(1);
  }
  return (54.0244 + 30.6951*(cmd-128))*2.4662;
}


//----------------------------- o ---------------------------------------


//calibrate the timestream signals
bool Detector::calibrate()
{
  //this should be straightforward now that the constructor
  //has calculated everything we need
  //can do this on the cpu since it's fast
  if (isCalibrated) 
    return isCalibrated;

  //correct for extinction
//  if(extinction == 0)
//  for(int i=0;i<nSamples;i++) hValues[i] *= extinction;
//
//  //correct nonlinearity
//  for(int i=0;i<nSamples;i++) hValues[i] /= responsivity;
//
////  //correct electronics gain
//  for(int i=0;i<nSamples;i++) hValues[i] /= fecGain;
//
//
//  for(int i=0;i<nSamples;i++) hValues[i] *= abs(fcf)*1.e-3;
//
//  //multiply by absolute value of bolometer gain measured in beammaps
  double calFactor = getCalibrationFactor();
  for(int i=0;i<nSamples;i++) hValues[i]*=calFactor;
  isCalibrated=1;
  return 1;
}


double Detector::getCalibrationFactor(){
	if (extinction ==0){
	    cerr << "Detector::calibrate(): Must call estimateExtinction(tau)";
	    cerr << " from main" << endl;
	    exit(1);
	  }
	return (extinction*abs(fcf)*1e-3/(responsivity*fecGain));

}

//----------------------------- o ---------------------------------------

bool Detector::calculateScanWeight(Telescope* tel)
{
  //calculate the inverse variance of the unflagged samples
  //for each scan
  int nScans = tel->scanIndex.ncols();
  scanWeight.resize(nScans);

  //loop through the scans
  for(int i=0;i<nScans;i++){
    int si = tel->scanIndex[0][i];
    int ei = tel->scanIndex[1][i];
    double tmp=0.;
    double count=0.;
    tmp = stddev(hValues, hSampleFlags, si, ei, &count);
 
    //insist on at least 1s of good samples
    if (tmp!=tmp || count < samplerate)
    	scanWeight[i]=0.0;
    else
    	scanWeight[i]= 1.0/pow(tmp,2.0);
  }

  return 1;
}


//----------------------------- o ---------------------------------------

void Detector::setFcf(double f)
{
  fcf = f;
}

//----------------------------- o ---------------------------------------
///calculates the sensitivity of the Detector comparable to the IDL utilities. 
//This is still a work in progress and lacks comments because I don't know what anything does
//-Vishnu
double Detector::calculateSensitivity(Telescope* tel)
{
  using namespace Eigen;
  // prepare eigen containers
  VectorXd scans = Map<VectorXd>(&hValues[0], hValues.size());
  MatrixXI scanindex = Map<Matrix<int, Dynamic, Dynamic, RowMajor>>(tel->scanIndex[0], tel->scanIndex.nrows(), tel->scanIndex.ncols()).cast<Index>();
  auto logger = logging::createLogger("detector.sensitivity", this);
  SPDLOG_LOGGER_TRACE(logger, "scans{}", logging::pprint(scans));
  SPDLOG_LOGGER_TRACE(logger, "scanindex{}", logging::pprint(scanindex));
  MatrixXd sensitivities;
  MatrixXd noisefluxes;
  timestream::sensitivity(scans, scanindex, sensitivities, noisefluxes, samplerate, {3., 5.});
  sensitivity = sensitivities.mean();
  return sensitivity;
}

//----------------------------- o ---------------------------------------

void Detector::calibrateSensitivity()
{
  double calibrationFactor = getCalibrationFactor();
  sensitivity *= calibrationFactor;
}

//----------------------------- o ---------------------------------------
///adds a gaussian to the Detector timestream given a MatDoub of parameters for a list of detectors
void Detector::addGaussian(MatDoub &params, int detNumber)
{
  //extract parameters
  double amplitude = params[detNumber][1];
  double sigma_x = params[detNumber][2];
  double sigma_y = params[detNumber][3];
  double offset_x = params[detNumber][4];
  double offset_y = params[detNumber][5];

  //calculate the value to be added at each position
  for(int i=0;i<nSamples;i++){
    double toAdd = amplitude*exp(-1.*(pow(hRa[i] - offset_x, 2) / (2.*pow(sigma_x,2))
                                  + pow(hDec[i] - offset_y, 2) / (2.*pow(sigma_y,2))));
    //use subtraction sign because hValues are negative
		hValues[i] -= toAdd;
  }
}

//----------------------------- o ---------------------------------------



//generate the pointing data for this detector
//this version calculates it from scratch and
//requires a telescope object
bool Detector::getPointing(Telescope* tel, TimePlace* tp, Source* source)
{
  //some of this is observatory-dependent
  bool LMT=0, ASTE=0;
  string observatory = ap->getObservatory();
  if(observatory.compare("LMT") == 0){
    LMT=1;
  } else if(observatory.compare("ASTE") == 0){
    ASTE=1;
  }

  //is this an azel map?
  int azelMap = ap->getAzelMap();

  /*
  Here is the strategy:
  1) get map center ra/dec in J2000  (this held in AnalParams but possibly
     pulled from source coordinates in Source constructor.
  2) get telescope bsight ra/dec J2000 and convert from abs to phys
     (this done in Telescope since it only needs to happen once)
  3) get parallactic angle for each az/el sample
     (this done in Telescope since it only needs to happen once)
     It's unclear if this should be calculated for J2000 or current epoch.
  4) get detector offsets in delta az/el phys, rotate according 
     to elevation angle and rotate to delta ra/dec phys
  5) compute detector ra/dec
  */

  hRa.resize(nSamples);
  hDec.resize(nSamples);

  //(2) Telescope boresight ra/dec
  //These are in physical coordinates.
  double* raPhys=NULL;
  double* decPhys=NULL;
  if(azelMap == 0){
    raPhys = &tel->hTelRaPhys[0];
    decPhys = &tel->hTelDecPhys[0];
  } else {
    raPhys = &tel->hTelAzPhys[0];
    decPhys = &tel->hTelElPhys[0];
  }

  //(3) Parallactic Angle
  double* pA = &tel->paraAngle[0];

  //(4) Detector Offsets + ra/dec computation
  //start by rotating offsets by elevation angle to
  //counteract field rotation

 // if(1){
    VecDoub raOff(nSamples);
    VecDoub decOff(nSamples);
    double* bsOffset = ap->getBsOffset();
    double ratmp;
    double dectmp;
    double azOfftmp;
    double elOfftmp;
    double* hElDes = &tel->hTelElDes[0];
    double pa2;
    for(int i=0;i<nSamples;i++){
      if(LMT){
    	  azOfftmp = cos(hElDes[i])*azOffset - sin(hElDes[i])*elOffset;
    	  elOfftmp = cos(hElDes[i])*elOffset + sin(hElDes[i])*azOffset;
      } else {
    	  azOfftmp = azOffset;
    	  elOfftmp = elOffset;
      }

      //apply the bs offset assuming it is in arcseconds like azOffset
      //and elOffset
      azOfftmp += bsOffset[0];
      elOfftmp += bsOffset[1];

      //don't apply the pa transformation if azelMap is selected
      if(!azelMap){
	pa2 = pA[i]-PI;
	azOfftmp *= -1.;
	ratmp = azOfftmp*cos(pa2) - elOfftmp*sin(pa2);
	dectmp= azOfftmp*sin(pa2) + elOfftmp*cos(pa2);
	hRa[i] = ratmp*RAD_ASEC + raPhys[i];
	hDec[i] = dectmp*RAD_ASEC + decPhys[i];
      } else {
	hRa[i] = azOfftmp*RAD_ASEC + raPhys[i];
	hDec[i] = elOfftmp*RAD_ASEC + decPhys[i];
      }
    }
//  }//elegant approach

/*
  if(0){
    //here's the brute force approach
    double* saz = &source->hSourceAz[0];
    double* sel = &source->hSourceEl[0];
    VecDoub azTmp(nSamples);
    VecDoub elTmp(nSamples);
    VecDoub azTmp2(nSamples);
    VecDoub elTmp2(nSamples);
    double cdes;
    double sdes;
    azOffset = azOffset/3600.*TWO_PI/360.;
    elOffset = elOffset/3600.*TWO_PI/360.;
    for(int i=0;i<nSamples;i++){
      cdes = cos(tel->hTelElDes[i]);
      sdes = sin(tel->hTelElDes[i]);
      azTmp[i] = (tel->hTelAzAct[i]-saz[i])*cdes-tel->hTelAzCor[i];
      elTmp[i] = tel->hTelElAct[i] -sel[i] -tel->hTelElCor[i];
      azTmp[i] += cdes*azOffset - sdes*elOffset; 
      elTmp[i] += cdes*elOffset + sdes*azOffset;
    }
    physToAbs(&azTmp[0],&elTmp[0],saz,sel,&azTmp2[0],&elTmp2[0],nSamples);
    azElToRaDec2000(tp, &azTmp2[0], &elTmp2[0], &hRa[0], &hDec[0], nSamples);
    double* pmg;
    pmg = ap->getMasterGridJ2000();
    double centerRa = pmg[0];
    double centerDec = pmg[1];
    absToPhys(&hRa[0], &hDec[0], centerRa, centerDec, 
	      &hRa[0], &hDec[0], nSamples);
  }//brute force approach
*/
  //get max and min values of each coordinate during valid scans
  int nScans = tel->scanIndex.ncols();
  minX = 8*PI;
  maxX = -8*PI;
  minY = 8*PI;
  maxY = -8*PI;

  //loop through the scans
  for(int i=0;i<nScans;i++){
    int si = tel->scanIndex[0][i];
    int ei = tel->scanIndex[1][i]+1;
    for(int j=si;j<ei;j++){
      if(hSampleFlags[j]){
    	  if(hRa[j] < minX) minX = hRa[j];
    	  if(hRa[j] > maxX) maxX = hRa[j];
    	  if(hDec[j] < minY) minY = hDec[j];
    	  if(hDec[j] > maxY) maxY = hDec[j];
      }
    }
  }

  isPointingGenerated=1;
  return isPointingGenerated;
}


//---------------------------- o ----------------------------------------
bool Detector::getAzElPointing (Telescope *tel){
  azElRa.resize(nSamples);
  azElDec.resize(nSamples);
  azElRaPhys.resize(nSamples);
  azElDecPhys.resize(nSamples);
  
  bool LMT = (ap->getObservatory().compare("LMT") == 0);

  double azOff, elOff;

  for (int i=0; i<nSamples; i++){

	  if(LMT){
	      	  azOff = cos(tel->hTelElAct[i])*azOffset - sin(tel->hTelElAct[i])*elOffset;
	      	  elOff = cos(tel->hTelElAct[i])*elOffset + sin(tel->hTelElAct[i])*azOffset;
      } else {
	      	  azOff = azOffset;
	      	  elOff = elOffset;
	  }


	  azElRa[i] = tel->hTelAzAct[i]*360.0/TWO_PI + azOff/3600.0;
	  azElDec[i] = tel->hTelElAct[i]*360.0/TWO_PI + elOff/3600.0;
	  azElRaPhys[i] = tel->hTelAzPhys[i]*360.0/TWO_PI + azOff/3600.0;
	  azElDecPhys[i]= tel->hTelElPhys[i]*360.0/TWO_PI + elOff/3600.0;
  }
  return 1;
}

//----------------------------- o ---------------------------------------


bool Detector::makeKernelTimestream(Telescope* tel)
{
  //The idea here is to make a fake timestream representing the
  //detector's response to a 1Jy point source located at the 
  //pointing center.
  //Do this entirely in physical coordinates to make it easy.

  //There is a choice here.  We could do this in az/el coordinates and
  //then transform back to the timestreams, or use Physical Ra/Dec
  //coordinates and put the kernel at the mastergrid location.  I'm
  //going to make these selectable but default to the Ra/Dec approach
  //as it seems to match the IDL outputs a little better.

  //The goal is to calculate the distance for use in the gaussian kernel
  VecDoub dist(nSamples);
  
  if(0){
    //This is the Az/El approach
    //some of what follows is observatory-dependent
    bool LMT=0, ASTE=0;
    string observatory = ap->getObservatory();
    if(observatory.compare("LMT") == 0){
      LMT=1;
    } else if(observatory.compare("ASTE") == 0){
      ASTE=1;
    }
    
    //rotated detector offsets
    VecDoub azo(nSamples);
    VecDoub elo(nSamples);
    for(int i=0;i<nSamples;i++){
      if(LMT){
	azo[i] = cos(tel->hTelElDes[i])*azOffset-
	  sin(tel->hTelElDes[i])*elOffset;
	elo[i] = cos(tel->hTelElDes[i])*elOffset+
	  sin(tel->hTelElDes[i])*azOffset;
      } else {
	azo[i] = azOffset;
	elo[i] = elOffset;
      }
    }

    //distance from detector to source in physical coordinates
    //make sure to use bsOffsets
    double* bsOffset = ap->getBsOffset();
    for(int i=0;i<nSamples;i++){
      azo[i] = azo[i]/3600./360.*TWO_PI;
      elo[i] = elo[i]/3600./360.*TWO_PI;
      azo[i] += tel->hTelAzPhys[i] + bsOffset[0]*RAD_ASEC;
      elo[i] += tel->hTelElPhys[i] + bsOffset[1]*RAD_ASEC;
      dist[i] = sqrt(pow(azo[i],2)+pow(elo[i],2));
    }
  } else {
    //this makes use of the previously calculate physical Ra/Dec values
    for(int i=0;i<nSamples;i++){
      dist[i] = sqrt(pow(hRa[i],2)+pow(hDec[i],2));
    }  
  }

  //now run through the distances and generate a signal.  If the 
  //source is more than 3 beam sigmas away, call it 0.
  hKernel.resize(nSamples);
  double sigma = (beamSigAz+beamSigEl)/2./3600./360.*TWO_PI;
  for(int i=0;i<nSamples;i++){
    hKernel[i] = 0.;
    if(dist[i] <= 3.*sigma){
      hKernel[i] = exp(-0.5*pow(dist[i]/sigma,2));
    }
  }

  //check if altKernel is set.  If so, redefine the kernel
  if(ap->getAltKernel()){
    if(ap->getKernelName() == "core"){
      //This is the starless core kernel defined by Rob Gutermuth
      double r0 = ap->getCoreR0();
      double p = ap->getCoreP();
      double ar = ap->getCoreAxisRatio();
      for(int i=0;i<nSamples;i++){
	hKernel[i] = 0.;
	if(dist[i] <= 6.*sigma){
	  double thetas = atan(hDec[i]/hRa[i]);
	  double rads = dist[i];
	  double reff = sqrt(pow(rads*cos(thetas)*sqrt(ar),2) +
			     pow(rads*sin(thetas)*sqrt(ar),2));
	  hKernel[i] = pow(1.+(reff/r0),-p);
	}
      }     
    }
  }


  return 1;
}


//----------------------------- o ---------------------------------------


void Detector::throwXmlError(string p)
{
  cerr << "XML error in " << bolostatsFile << ": " << p << endl;
  exit(-1);
}


//----------------------------- o ---------------------------------------


int Detector::getNSamples()
{
  return nSamples;
}

double Detector::getSamplerate()
{
  return samplerate;
}

double Detector::getMaxX()
{
  return maxX;
}

double Detector::getMinX()
{
  return minX;
}

double Detector::getMaxY()
{
  return maxY;
}

double Detector::getMinY()
{
  return minY;
}

double Detector::getFcf()
{
  return fcf;
}

double Detector::getSensitivity()
{
  return sensitivity;
}

string Detector::getName()
{
  return name;
}

double Detector::getExtinction()
{
  return extinction;
}

double Detector::getFecGain()
{
  return fecGain;
}

int Detector::getId()
{
  return id;
}



bool Detector::setAtmTemplate(VecDoub temp) {
	if (temp.size() != size_t(nSamples)){
		cerr<<"SetAtmTemplate(). Wrong data size. Imploding"<<endl;
		exit(-1);
	}
	atmTemplate = temp;
	return true;
}

//----------------------------- o ---------------------------------------


//destructor
Detector::~Detector()
{

}

