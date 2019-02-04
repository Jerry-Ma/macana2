#ifndef PORT_DESPIKE_H
#define PORT_DESPIKE_H

//Eigen Includes
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/CXX11/Tensor>

//Macana Includes
#include "nr3.h"
#include <fftw3.h>
#include <gsl/gsl_sf_bessel.h> //Needed for filter

namespace despike{

using namespace Eigen;

//Funciton to write vectors to file
template<typename Scalar1>
void write_to_file(Scalar1 &vec, string filepath){
    ofstream filename;
       filename.open(filepath);
       for(int i = 0; i < vec.size(); i++){
           stringstream sstm;
           filename << vec[i];
           filename << endl;
     }
       filename.close();
}

/*----------------------------------------------------------------------------------------------*/

//Functions for Macana
bool Macanadespike(VecDoub &hValues, VecDoub &hSampleFlags, double nSigmaSpikes, const int nSamples,
                   const double samplerate, const int despikeWindow, bool isLowpassed, bool isDespiked);
double mean(double* arr, int nSamp);
double mean(VecDoub &arr);
double stddev(double* arr, int nSamp, double mn);
double stddev(VecDoub &arr);
void smooth(VecDoub &inArr, VecDoub &outArr, int w);


/*----------------------------------------------------------------------------------------------*/

//Function to generate a time series of random data points.
template <typename Scalar>
void MakeData(Scalar &hValues, const int nSamples, const double distmean, const double diststddev){

    std::mt19937 generator;
    std::normal_distribution<double> dist(distmean, diststddev);

    // Add Gaussian noise
    for(int i=0; i<nSamples; i++){
        hValues[i] = hValues[i] + dist(generator);
    }
}

/*----------------------------------------------------------------------------------------------*/

//Function to generate a time series of random data points.
template <typename Scalar>
void AddSpikes(Scalar &hValues, const int nSamples, const int nSpikes, const int SpikeSigma, const double diststddev){

    std::mt19937 generator;

    std::uniform_int_distribution<int> gen1(0, nSamples-1);
    std::uniform_int_distribution<int> gen2(0, SpikeSigma);

    int SpikeIndex;

    // Add Gaussian noise
    for(int i=0; i<nSpikes; i++){
        SpikeIndex = gen1(generator);
        hValues[SpikeIndex] = gen2(generator)*diststddev;
    }
}

/*----------------------------------------------------------------------------------------------*/

//implements the boxcar smoothing algorithm used by IDL
//inArr is the original array to be smoothed
//outArr is the smoothed version
//w is the boxcar width in samples
template <typename Scalar1, typename Scalar2>
void Eigensmooth(Eigen::DenseBase<Scalar1> &inArr, Eigen::DenseBase<Scalar2> &outArr, int w)
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
  outArr.head((w-1)/2) = inArr.head((w-1)/2);
  outArr.segment(nIn-(w+1)/2.+1,(w+1)/2.-1) = inArr.segment(nIn-(w+1)/2.+1,(w+1)/2.-1);

  //here is the smoothed part
  double winv = 1./w;
  int wm1d2 = (w-1)/2.;
  int wp1d2 = (w+1)/2.;
  for(Index i=wm1d2;i<=nIn-wp1d2;i++){
    outArr[i]=winv*inArr.segment(i-wm1d2,w).sum();
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

template <typename Scalar1, typename Scalar2>
bool Eigendespike(Eigen::DenseBase<Scalar1> &hValues, Eigen::DenseBase<Scalar2> &hSampleFlags,
                  double nSigmaSpikes, const int nSamples, const double samplerate,
                  const int despikeWindow, const double timeConstant, bool isLowpassed, bool isDespiked)
{
  //here we do all the despiking we can do without the
  //use of other detectors.
  //1) find the spikes
  //2) flag the timestream values around the spikes

  //FIND THE SPIKES
  int nSpikes=0;

  //calculate vector of adjacent signals and its moments
  double deltaMean=0;
  double deltaStddev=0;
  Eigen::VectorXf delta(nSamples-1);
  //Doub *pDelta = &delta[0];
  int nDelta = nSamples-1;

  Eigen::VectorXf deltaLower(nSamples -1);
  Eigen::VectorXf deltaUpper(nSamples -1);

  deltaLower = hValues.head(nSamples-1);
  deltaUpper = hValues.tail(nSamples-1);

  delta = deltaUpper - deltaLower;

  deltaMean = delta.mean();
  deltaStddev = std::sqrt((delta.array() - deltaMean).square().sum()/(nDelta-1));

  //search the delta array for spikes that are nSigmaSpikes
  //times larger than the stddev
  for(Index i=1;i<nDelta;i++){
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
    deltaMean = delta.mean();
    deltaStddev = std::sqrt((delta.array() - deltaMean).square().sum()/(nDelta-1));
    for(Index i=1;i<nDelta;i++){
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
  int c;
  int indx;
  for(Index i=0;i<nSamples;i++){
    if(hSampleFlags[i] == 0){
      c=0;
      indx=0;
      for(Index j=1;j<=10 && i+j<nSamples;j++){
        indx = i+j;
        if(hSampleFlags[indx] == 0) c++;
      }
      if(c>0){
        //reduce the number nSpikes c times
        nSpikes = nSpikes - c;
        //fix hSampleFlags in the window
        //(11/11/15 - GW - this next line incorrectly
        //used to start j with 0 instead of 1)
        for(Index j=1;j<=10 && i+j<nSamples;j++){
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

    //recount the spikes to avoid pathologies
    nSpikes = 0;
    for(int i=0;i<nSamples-1;i++) if(hSampleFlags[i] == 0) nSpikes++;

    //get a list of spike locations and det values
    Eigen::VectorXi spikeLoc(nSpikes);
    Eigen::VectorXf spikeVals(nSpikes);
    int count=0;

    for(Index i=0;i<nSamples-1;i++){
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
      Eigen::VectorXi deltaSpikeLoc(nSpikes+1);
      deltaSpikeLoc[0] = spikeLoc[0]-0;  //distance to first spike
      deltaSpikeLoc[nSpikes] = nSamples-spikeLoc[nSpikes-1];
      deltaSpikeLoc.tail(nSpikes) = spikeLoc.tail(nSpikes) - spikeLoc.head(nSpikes);

      int mxWindow=deltaSpikeLoc[0];
      int mxWindowIndex=0;
      for(Index i=1;i<=nSpikes;i++){
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
    VectorXf subVals(winSize);
    //for(int i=winIndex0;i<winIndex1-1;i++){
      //cerr << i-winIndex0 << "   " << i << endl;
      //subVals[i-winIndex0]=hValues[i];
      subVals.head(winIndex1-1 - winIndex0) = hValues.head(winIndex1-1 - winIndex0);
    //}

    //make a boxcar smoothed copy of the subarray
    Eigen::VectorXf smoothedSubVals(winSize);
    Eigensmooth<Eigen::VectorXf, Eigen::VectorXf>(subVals,smoothedSubVals,10);

    //here is the estimate of the stardard deviation
    double sigest;
    subVals.head(winSize) =  subVals.head(winSize) - smoothedSubVals.head(winSize);
    sigest = std::sqrt((subVals.array() - subVals.mean()).square().sum()/(subVals.size()-1));
    if(sigest < 1.e-8) sigest=1.e-8;

    //calculate for each spike the time it takes to decay to sigest
    Eigen::VectorXf decayLength(nSpikes);
    decayLength = -samplerate*timeConstant*log(abs(sigest/spikeVals.array()));
    for(Index i=0;i<nSpikes;i++){
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
    for(Index i=0;i<nSpikes;i++){
      if(!isLowpassed){
        for(Index j=-decayLength[i];j<2*decayLength[i];j++){
          if(spikeLoc[i]+j >= 0 && spikeLoc[i]+j < nSamples)
            hSampleFlags[spikeLoc[i]+j] = 0;
        }
      } else {
        for(Index j=-(despikeWindow-1)/2;j<(despikeWindow-1)/2;j++){
          if(spikeLoc[i]+j >= 0 && spikeLoc[i]+j < nSamples)
            hSampleFlags[spikeLoc[i]+j] = 0;
        }
      }
    }
  }//if nSpikes gt 0 then flag samples around spikes

  isDespiked=1;
  return 1;
}
} // namespace

#endif // PORT_DESPIKE_H
