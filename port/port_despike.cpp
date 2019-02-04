#include "port_despike.h"

namespace despike {

/*----------------------------------------------------------------------------------------------*/

//calculates mean of an array of length nsamp
double mean(double* arr, int nSamp)
{
  double mn=0.;
  for(int i=0;i<nSamp;i++) mn+=arr[i];
  return mn/nSamp;
}

/*----------------------------------------------------------------------------------------------*/

//calculates mean of a VecDoub array
double mean(VecDoub &arr)
{
  double mn=0.;
  int nSamp=arr.size();
  int nValid = 0;
  for(int i=0;i<nSamp;i++){
    if (!isnan(arr[i]) && !isinf(arr[i])){
        mn+=arr[i];
        nValid++;
    }

  }
  return mn/nValid;
}

/*----------------------------------------------------------------------------------------------*/

//calculates standard deviation of an array of length nsamp
//and mean of mn
double stddev(double* arr, int nSamp, double mn)
{
  double std=0.;
  for(int i=0;i<nSamp;i++)
    std+=(arr[i]-mn)*(arr[i]-mn);
  std = std/(nSamp-1.);
  return sqrt(std);
}

/*----------------------------------------------------------------------------------------------*/

//calculates standard deviation of a VecDoub array
double stddev(VecDoub &arr)
{
  double mn=mean(arr);
  int nSamp = arr.size();
  int nValid = 0;
  double std=0.;
  for(int i=0;i<nSamp;i++){
     if (!isnan(arr[i]) && !isinf(arr[i])){
        std+=(arr[i]-mn)*(arr[i]-mn);
        nValid++;
     }
  }
  std = std/(nValid-1.);
  return sqrt(std);
}

/*----------------------------------------------------------------------------------------------*/

//implements the boxcar smoothing algorithm used by IDL
//inArr is the original array to be smoothed
//outArr is the smoothed version
//w is the boxcar width in samples
void smooth(VecDoub &inArr, VecDoub &outArr, int w)
{
  int nIn = inArr.size();
  int nOut = outArr.size();
  if(nIn != nOut){
    cerr << "vector_utilities::smooth() the input array must be";
    cerr << " the same size as output array" << endl;
    exit(1);
  }

  //as with idl, if w is even then add 1
  if(w % 2 == 0) w++;

  //first deal with the end cases where the output is the input
  for(int i=0;i<(w-1)/2.;i++) outArr[i] = inArr[i];
  for(int i=nIn-(w+1)/2.+1;i<nIn;i++) outArr[i] = inArr[i];

  //here is the smoothed part
  double winv = 1./w;
  int wm1d2 = (w-1)/2.;
  int wp1d2 = (w+1)/2.;
  double tmpsum;
  for(int i=wm1d2;i<=nIn-wp1d2;i++){
    tmpsum=0;
    for(int j=0;j<w;j++) tmpsum += inArr[i+j-wm1d2];
    outArr[i] = winv*tmpsum;
  }
}

/*----------------------------------------------------------------------------------------------*/

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
bool Macanadespike(VecDoub &hValues, VecDoub &hSampleFlags, double nSigmaSpikes, const int nSamples,
                   const double samplerate, const int despikeWindow, bool isLowpassed, bool isDespiked)
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
        cerr << "Found a spike at " << i << endl;
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
        //cerr << "Detector::despike(): detector:" << name << endl;
        cerr << "Either spikes everywhere or something else is";
        cerr << " terribly wrong." << endl;
        //cerr << "Data file: " << dataFile << endl;
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
} //namespace
