#include <netcdfcpp.h>
#include <cmath>
#include <algorithm>
#include <stdio.h>
using namespace std;
#include "nr3.h"
#include "astron_utilities.h"
#include "vector_utilities.h"
#include <gsl/gsl_vector.h>
//#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>


#include "mpfit.h"

double select(vector<double> input, int index){
 //Partition input based on a selected pivot
    //More on selecting pivot later
  //you now know the index of the pivot
  //if that is the index you want, return
  //else, if it's larger, recurse on the larger partition
  //if it's smaller, recurse on the smaller partition
  unsigned int pivotIndex = rand() % input.size();
  double pivotValue = input[pivotIndex];
  vector<double> left;
  vector<double> right;
  for(unsigned int x = 0; x < input.size(); x++){
    if(x != pivotIndex){
      if(input[x] > pivotValue){
        right.push_back(input[x]);
      }
      else{
        left.push_back(input[x]);
      }
    }
  }
  if((int) left.size() == index){
    return pivotValue;
  }
  else if((int) left.size() < index){
    return select(right, index - left.size() - 1);
  }
  else{
    return select(left, index);
  }
}


//calculates mean of an array of length nsamp
double mean(double* arr, int nSamp)
{
  double mn=0.;
  for(int i=0;i<nSamp;i++) mn+=arr[i];
  return mn/nSamp;
}


double mean(double* arr, bool *flags, int nSamp)
{
  double mn=0.;
  size_t valid = 0;

  for(int i=0;i<nSamp;i++){
	  if (!isnan(arr[i]) && !isinf(arr[i]) && flags[i]){
		  mn+=arr[i];
		  valid++;
	  }

  }
  if (valid == 0){
	  cerr<<"mean(). No valid points specified"<<endl;
	  exit(-1);
  }
  return mn/valid;
}


//----------------------------- o ---------------------------------------


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

//calculates mean of a VecDoub array
double mean(VecDoub &arr, VecBool &flags)
{
  double mn=0.;
  int nSamp=arr.size();
  int nValid = 0;
  for(int i=0;i<nSamp;i++){
    if (!isnan(arr[i]) && !isinf(arr[i]) &&flags[i]){
	mn+=arr[i];
	nValid++;
    }
    
  }
  return mn/nValid;
}


//----------------------------- o ---------------------------------------


//calculates mean of a portion of a VecDoub array
//samples start at start and end at end
double mean(VecDoub &arr, int start, int end)
{
  double mn=0;
  int nSamp=arr.size();
  if(start < 0 || end >= nSamp){
    cerr << "vector_utilities::mean(): Out of bounds indices." << endl;
    exit(1);
  }
  for(int i=start;i<=end;i++) mn+=arr[i];
  return mn/(end-start);
}

//calculates mean of a portion of a VecDoub array
//samples start at start and end at end
double mean(VecDoub &arr, VecBool &flags, int start, int end, double *ncount )
{
  double mn=0;
  int nSamp=arr.size();
  if(start < 0 || end >= nSamp){
    cerr << "vector_utilities::mean(): Out of bounds indices." << endl;
    exit(1);
  }
  double ngood = 0;
  for(int i=start;i<=end;i++)
	  if (flags[i]){
		  mn+=arr[i];
		  ngood++;
	  }

  if (ncount != NULL)
	  *ncount = ngood;

  if (ngood ==0)
  {
      cerr << "Warning. Not valid points" << endl;
      return std::numeric_limits<double>::quiet_NaN();
  }

  return mn/(end-start);
}

//----------------------------- o ---------------------------------------


//calculates mean of a portion of a VecDoub array
//samples start at start and end at end
double mean(double *arr, int start, int end)
{
  double mn=0;
  for(int i=start;i<=end;i++) mn+=arr[i];
  return mn/(end-start);
}

//----------------------------- o ---------------------------------------


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


//----------------------------- o ---------------------------------------


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

//----------------------------- o ---------------------------------------


//calculates standard deviation of a VecDoub array
double stddev(VecDoub &arr, VecBool &flags)
{
  double mn=mean(arr);
  int nSamp = arr.size();
  int nValid = 0;
  double std=0.;
  for(int i=0;i<nSamp;i++){
     if (!isnan(arr[i]) && !isinf(arr[i]) && flags[i]){
	std+=(arr[i]-mn)*(arr[i]-mn);
	nValid++;
     }
  }
  std = std/(nValid-1.);
  return sqrt(std);
}


//----------------------------- o ---------------------------------------


//calculates standard deviation of a portion of a VecDoub array
//samples start at start and end at end
double stddev(VecDoub &arr, int start, int end)
{
  int nSamp=arr.size();
  if(start < 0 || end >= nSamp){
    cerr << "vector_utilities::stddev(): Out of bounds indices." << endl;
    exit(1);
  }
  double mn=mean(arr,start,end);
  double std=0.;
  for(int i=start;i<=end;i++) std+=(arr[i]-mn)*(arr[i]-mn);
  std = std/(end-start);
  return sqrt(std);
}

double stddev(VecDoub &arr, VecBool &flags, int start, int end, double *ncount)
{
  int nSamp=arr.size();
  if(start < 0 || end >= nSamp){
    cerr << "vector_utilities::stddev(): Out of bounds indices." << endl;
    exit(1);
  }
  double mn=mean(arr,flags,start,end, ncount);

  if (mn!=mn){
	  cerr<<"stddev():: Invalid mean value"<<endl;
	  return mn;
  }

  double std=0.;
  double ngood=0;
  for(int i=start;i<=end;i++)
	  if (flags[i]){
		  std+=(arr[i]-mn)*(arr[i]-mn);
		  ngood++;
	  }
  std = std/(ngood+1.0);
  return sqrt(std);
}

//----------------------------- o ---------------------------------------


//calculates MAD, median absolute deviation (from the median), of a VecDoub
//array. MAD is a relatively robust outlier-resistant replacement for stddev.

double medabsdev(VecDoub &arr)
{
  int nSamp=arr.size();
  double med = median(arr);
  VecDoub absDelt(nSamp);
  for(int i=0;i<nSamp;i++)
    absDelt[i] = abs(arr[i]-med);
  return median(absDelt);
}




//----------------------------- o ---------------------------------------


//finds maximum and minimum elements of a VecDoub array
void maxmin(VecDoub &arr, double* max, double* min)
{
  double mxtmp = arr[0];
  double mntmp = arr[0];
  double npts = arr.size();
  for(int i=1;i<npts;i++){
    mxtmp = (arr[i] > mxtmp) ? arr[i] : mxtmp;
    mntmp = (arr[i] < mntmp) ? arr[i] : mntmp;
  }
  *max = mxtmp;
  *min = mntmp;
}

//----------------------------- o ---------------------------------------


//finds maximum and minimum elements of a VecInt array
void maxmin(VecInt &arr, int* max, int* min)
{
  int mxtmp = arr[0];
  int mntmp = arr[0];
  int npts = arr.size();
  for(int i=1;i<npts;i++){
    mxtmp = (arr[i] > mxtmp) ? arr[i] : mxtmp;
    mntmp = (arr[i] < mntmp) ? arr[i] : mntmp;
  }
  *max = mxtmp;
  *min = mntmp;
}


//----------------------------- o ---------------------------------------


//finds maximum and minimum elements of a double* array
void maxmin(double* arr, double* max, double* min, int npts)
{
  double mxtmp = arr[0];
  double mntmp = arr[0];
  for(int i=1;i<npts;i++){
    mxtmp = (arr[i] > mxtmp) ? arr[i] : mxtmp;
    mntmp = (arr[i] < mntmp) ? arr[i] : mntmp;
  }
  *max = mxtmp;
  *min = mntmp;
}

//----------------------------- o ---------------------------------------

double median (double *data, size_t nSamp){

	double *tmpData = new double[nSamp];
	double median = 0.0;

	for (register size_t i=0; i<nSamp; i++){
		tmpData[i]=data[i];
	}
	gsl_sort(tmpData,1,nSamp);
	median = gsl_stats_median_from_sorted_data(tmpData,1,nSamp);
	delete [] tmpData;
	return median;
}

double median (double *data, bool *flags, size_t nSamp){

	size_t ngood = 0;

	for (register size_t i=0; i<nSamp; i++)
		if (flags[i])
			ngood++;
	double *tmpData = new double[ngood];
	double median = 0.0;
	size_t ix=0;
	for (register size_t i=0; i<nSamp; i++){
		if (flags[i])
			tmpData[ix++]=data[i];
	}
	gsl_sort(tmpData,1,ngood);
	median = gsl_stats_median_from_sorted_data(tmpData,1,ngood);
	delete [] tmpData;
	return median;
}


//----------------------------- o ---------------------------------------

double percentile (double *data, size_t nSamp, double percentile){

	double *tmpData = new double[nSamp];
	double p = 0.0;

	for (register size_t i=0; i<nSamp; i++){
		tmpData[i]=data[i];
	}
	gsl_sort(tmpData,1,nSamp);
	p = gsl_stats_quantile_from_sorted_data(tmpData,1,nSamp, percentile);
	delete [] tmpData;
	return p;
}

//----------------------------- o ---------------------------------------


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



//----------------------------- o ---------------------------------------


//implements the boxcar smoothing algorithm used by IDL
//using the edge_truncate keyword activated
//inArr is the original array to be smoothed
//outArr is the smoothed version
//w is the boxcar width in samples
void smooth_edge_truncate(VecDoub &inArr, VecDoub &outArr, int w)
{
  int nIn = inArr.size();
  int nOut = outArr.size();
  if(nIn != nOut){
    cerr << "vector_utilities::smooth_edge_truncate() ";
    cerr << "the input array must be ";
    cerr << "the same size as output array" << endl;
    exit(1);
  }

  //as with idl, if w is even then add 1
  if(w % 2 == 0) w++;

  //do this all at once
  double winv = 1./w;
  int wm1d2 = (w-1)/2;
  double tmpsum;
  for(int i=0;i<nIn;i++){
    tmpsum=0;
    for(int j=0;j<w;j++){
      int addindex = i+j-wm1d2;
      if(addindex < 0) addindex=0;
      if(addindex > nIn-1) addindex=nIn-1;
      tmpsum += inArr[addindex];
    }    
    outArr[i] = winv*tmpsum;
  }
}

void convolve (double *data, size_t nData, double *kernel, bool first)
{
	gsl_fft_real_wavetable * real = gsl_fft_real_wavetable_alloc (nData);
	gsl_fft_halfcomplex_wavetable * hc;
	gsl_fft_real_workspace * work =gsl_fft_real_workspace_alloc (nData);


	gsl_fft_real_transform (data, 1, nData, real, work);
	if (first){
		double totalKernel =0.0;
		for (size_t i=0; i<nData; i++)
			totalKernel+=kernel[i];
		for (size_t i=0; i<nData; i++)
			kernel[i]/=totalKernel;
		gsl_fft_real_transform (kernel,1,nData, real, work);
	}
	gsl_fft_real_wavetable_free (real);

	for (size_t i=0; i<nData; i++)
		data[i]*=kernel[i];
	hc = gsl_fft_halfcomplex_wavetable_alloc (nData);
	gsl_fft_halfcomplex_inverse (data, 1, nData, hc, work);

	gsl_fft_halfcomplex_wavetable_free (hc);
	gsl_fft_real_workspace_free (work);

}


//----------------------------- o ---------------------------------------


//debug utility - write out vector data as txt file
//template <class T> bool writeVecOut(const char* outFile,  T* outData, size_t nOut)
//{
//  ofstream out(outFile);
//  if(out.bad()){
//    cerr << "vector_utilities::writeVecOut():";
//    cerr << " Error opening " << outFile << endl;
//    return 0;
//  }
//
//  //set precision to something reasonable for det values
//  out.precision(14);
//
//  //and the output
//  for(size_t i=0;i<nOut;i++) out << outData[i] << "\n";
//  if(out.bad()){
//    cerr << "vector_utilities::writeVecOut():";
//    cerr << " Error writing to " << outFile << endl;
//    return 0;
//  }
//
//  //flush the buffer just to be sure
//  out.flush();
//
//  //announce what you just did
//  cout << "Wrote out " << outFile << endl;
//
//  return 1;
//}

//----------------------------- o ---------------------------------------


//debug utility - write out matrix data as txt file
bool writeMatOut(const char* outFile, MatDoub &outMat)
{
  ofstream out(outFile);
  if(out.bad()){
    cerr << "vector_utilities::writeMatOut():";
    cerr << " Error opening " << outFile << endl;
    return 0;
  }

  //set precision to something reasonable for det values
  out.precision(14);

  //and the output
  for(int i=0;i<outMat.nrows();i++)
  for(int j=0;j<outMat.ncols();j++)
    {
      out << outMat[i][j] << "\n";
      if(out.bad()){
	cerr << "vector_utilities::writeMatOut():";
	cerr << " Error writing to " << outFile << endl;
	return 0;
      }
    }

  //flush the buffer just to be sure
  out.flush();

  //announce what you just did
  cout << "Wrote out " << outFile <<  " with nrows=" << outMat.nrows();
  cout << " and ncols=" << outMat.ncols() << endl;

  return 1;
}


//----------------------------- o ---------------------------------------


//This is a direct translation of digital_filter.pro from the idl
//distribution.  Note that coefOut must be allocated with
//2*nTerms+1 samples;
//fLow and fHigh are fractions of the Nyquist frequency
//set aGibbs to 50 (according to IDL help)
//the returned filter will be centered (this is unlike the idl version)
bool digitalFilter(double fLow, double fHigh, double aGibbs, 
		   int nTerms, double* coefOut)
{
  // computes Kaiser weights W(N,K) for digital filters.
  // W = COEF = returned array of Kaiser weights
  // N = value of N in W(N,K), i.e. number of terms
  // A = Size of gibbs phenomenon wiggles in -DB.
  
  double alpha;
  if(aGibbs <= 21.){
    alpha = 0.;
  }
  if(aGibbs >= 50){
    alpha = 0.1102*(aGibbs-8.7);
  }
  if(aGibbs > 21. && aGibbs < 50){
    alpha = 0.5842*pow(aGibbs-21.,0.4) + 0.07886*(aGibbs-21.);
  }

  //Band stop?
  double fStop = (fHigh < fLow) ? 1. : 0.;

  //Arg
  VecDoub arg(nTerms);
  for(double i=0;i<nTerms;i++) arg[i]=(i+1.)/nTerms;

  //Coef 
  VecDoub coef(nTerms);
  for(int i=0;i<nTerms;i++){
    coef[i] = gsl_sf_bessel_I0(alpha*sqrt(1.-pow(arg[i],2))) /
      gsl_sf_bessel_I0(alpha);
  }

  //t
  VecDoub t(nTerms);
  for(double i=0;i<nTerms;i++) t[i]=(i+1.)*PI;

  //here we go
  for(int i=0;i<nTerms;i++){
    coef[i] *= (sin(t[i]*fHigh)-sin(t[i]*fLow))/t[i];
  }

  //build the uncentered version
  for(int i=0;i<nTerms;i++) coefOut[i] = coef[nTerms-i-1];
  coefOut[nTerms] = fHigh-fLow-fStop;
  for(int i=0;i<nTerms;i++) coefOut[i+nTerms+1] = coef[i];

  //renormalize
  double area=0;
  for(int i=0;i<2*nTerms+1;i++) area += coefOut[i];
  for(int i=0;i<2*nTerms+1;i++) coefOut[i] /= area;

  return 1;

}


//----------------------------- o ---------------------------------------


Int Base_interp::locate(const Doub x)
{
  Int ju, jm, jl;
  if(n<2 || mm<2 || mm>n){
    cerr << "vector_utilities::Base_interp::locate(): locate size error";
    exit(1);
  }
  bool ascnd=(xx[n-1] >= xx[0]);
  jl=0;
  ju=n-1;
  while(ju-jl>1){
    jm=(ju+jl) >> 1;
    if((x >= xx[jm]) == ascnd)
      jl=jm;
    else
      ju=jm;
  }
  cor=abs(jl-jsav)>dj ? 0:1;
  jsav=jl;
  return max(0,MIN(n-mm,jl-((mm-2)>>1)));
}


//----------------------------- o ---------------------------------------


Int Base_interp::hunt(const Doub x)
/*Given a value x, return a value j such that x is (insofar as
possible) centered in the subrange xx[j .. j+mm-1] , where xx is the
stored pointer. The values in xx must be monotonic, either increasing
or decreasing. The returned value is not less than 0, nor greater than
n-1.  */
{
  Int jl=jsav, jm, ju, inc=1;
  if (n<2 || mm<2 || mm>n){
    cerr << "vector_utilities::Base_interp(): hunt size error.";
    exit(1);
  }
  ju=n-1;           //this is my line since otherwise it is uninitialized
  bool ascnd=(xx[n-1] >= xx[0]);
  if(jl<0 || jl>n-1){
    jl=0;
    ju=n-1;
  }else{
    if((x>=xx[jl]) == ascnd){
      for(;;){
	if(ju>=n-1){ju=n-1; break;}
	else if((x<xx[ju]) == ascnd) break;
	else{
	  jl=ju;
	  inc+=inc;
	}
      }
    }else{
      ju=jl;
      for(;;){
	jl=jl-inc;
	if(jl<=0){jl=0;break;}
	else if ((x>=xx[jl]) == ascnd) break;
	else{
	  ju=jl;
	  inc+=inc;
	}
      }
    }
  }
  while (ju-jl > 1){
    jm=(ju+jl) >> 1;
    if((x>xx[jm]) == ascnd)
      jl=jm;
    else
      ju=jm;
  }
  cor=abs(jl-jsav)>dj ? 0 : 1;
  jsav=jl;
  return max(0,MIN(n-mm,jl-((mm-2)>>1)));
}


//----------------------------- o ---------------------------------------


//linear interpolation using the routines above
bool interpolateLinear(VecDoub xx, VecDoub yy, int nUnique,
		       double* x, double *result, int nSamples)
{
  //preliminaries
  VecDoub tmp(nSamples,0.);   //storage
  double mx=xx[nUnique-1];
  double mn=xx[0];
  int locMax=nUnique-1;
  int locMin=0;        

  //build interpolation function
  Linear_interp interpme(xx,yy);

  //here's the interpolation
  for(int i=0;i<nSamples;i++){
    if(x[i] > mn && x[i] < mx)
      tmp[i]=interpme.interp(x[i]);
    if(x[i] <= mn)
      tmp[i]=yy[locMin];
    if(x[i] >= mx)
      tmp[i]=yy[locMax];
  }

  //edges are not guaranteed to come out right so set them to adjacent
  //values
  tmp[0] = tmp[1];
  tmp[nSamples-1]=tmp[nSamples-2];
  for(int i=0;i<nSamples;i++) result[i]=tmp[i];
  return 1;
}


//----------------------------- o ---------------------------------------


bool deNanSignal(double* sig, int nSamples)
{
  //need to find single point NaNs and replace with average of
  //bracketing values
  for(int i=1;i<nSamples-1;i++){
    if(sig[i] != sig[i]){
      //this is a NaN
      sig[i] = (sig[i-1]+sig[i+1])/2.;
    }
  }
  if(sig[0] != sig[0]){
    cerr << "vector_utilities::deNanSignal():";
    cerr << " First data point is corrupted as a NaN." << endl;
    cerr << "Setting equal to adjacent data point." << endl;
    sig[0] = sig[1];
  }
  if(sig[nSamples-1] != sig[nSamples-1]){
    cerr << "vector_utilities::deNanSignal():";
    cerr << " Last data point is corrupted as a NaN." << endl;
    cerr << "Setting equal to adjacent data point." << endl;
    sig[nSamples-1] = sig[nSamples-2];
  }
  return 1;
}


//----------------------------- o ---------------------------------------


bool removeDropouts(double* sig, int nSamples)
{
  //sometimes the LMT signals drop out.  Replace these with average of
  //adjacent signals
  for(int i=1;i<nSamples-1;i++){
    if(sig[i] <= 1.e-9){
      //this is a dropout
      sig[i] = (sig[i-1]+sig[i+1])/2.;
    }
  }
  if(sig[0] != sig[0]){
    cerr << "vector_utilities::removeDropouts():";
    cerr << " First data point is corrupted as a dropout." << endl;
    cerr << "Setting equal to adjacent data point." << endl;
    sig[0] = sig[1];
  }
  if(sig[nSamples-1] != sig[nSamples-1]){
    cerr << "vector_utilities::removeDropouts():";
    cerr << " Last data point is corrupted as a dropout." << endl;
    cerr << "Setting equal to adjacent data point." << endl;
    sig[nSamples-1] = sig[nSamples-2];
  }
  return 1;
}


//----------------------------- o ---------------------------------------


//turn wrapped signals to monotonically increasing
bool correctRollover(double* sig, double low, double high, 
		     double ulim, int nSamples)
{
  double mx=sig[0];
  double mn=sig[0];
  for(int i=1;i<nSamples;i++){
    mx = (mx > sig[i]) ? mx : sig[i];
    mn = (mn < sig[i]) ? mn : sig[i];
  }
  if(mx > high && mn < low){
    for(int i=0;i<nSamples;i++){
      if(sig[i]<low) sig[i]+=ulim;
    }
  }
  return 1;
}


//----------------------------- o ---------------------------------------


//reset any signals that our out of bounds to the mean of the adjacent samples
bool correctOutOfBounds(double* sig, double low, double high, int nSamples)
{
  for(int i=1;i<nSamples-1;i++)
    if(sig[i]<low || sig[i]>high) sig[i] = (sig[i-1]+sig[i+1])/2.;
  if(sig[0]<low || sig[0]>high){
    cerr << "correctOutOfBounds():";
    cerr << " First data point is out of bounds." << endl;
    cerr << "Setting equal to adjacent data point." << endl;
    sig[0] = sig[1];
  }
  if(sig[nSamples-1]<low || sig[nSamples-1]>high){
    cerr << "correctOutOfBounds():";
    cerr << " First data point is out of bounds." << endl;
    cerr << "Setting equal to adjacent data point." << endl;
    sig[nSamples-1] = sig[nSamples-2];
  }
  return 1;
}


//----------------------------- o ---------------------------------------


//a rewrite of IDL's hanning function
//but with alpha fixed to 0.5
MatDoub hanning(int n1in, int n2in)
{
  double n1 = (double) n1in;
  double n2 = (double) n2in;
  double a = 2.*PI/n1;
  VecDoub index(n1in);
  for(int i=0;i<n1in;i++) index[i] = (double) i;
  double b = 2.*PI/n2;
  VecDoub row(n1in);
  for(int i=0;i<n1in;i++) row[i] = -0.5 * cos(index[i]*a) + 0.5;
  index.resize(n2in);
  for(int i=0;i<n2in;i++) index[i] = (double) i;
  VecDoub col(n2in);
  for(int i=0;i<n2in;i++) col[i] = -0.5 * cos(index[i]*b) + 0.5;

  MatDoub han(n1in,n2in);
  for(int i=0;i<n1in;i++) 
    for(int j=0;j<n2in;j++)
      han[i][j] = row[i]*col[j];

  return han;
}

VecDoub hanning(int n1in){
  double n1 = (double) n1in;
  double a = 2.*PI/n1;
  
  VecDoub han(n1in);
  VecDoub index(n1in);
  
  for(int i=0;i<n1in;i++)
     index[i]=(double)i;
  for (int i=0;i<n1in; i++)
    han[i] = -0.5 * cos (index[i]*a) +0.5;
  
  return han;
}


//----------------------------- o ---------------------------------------


//a generic tool to histogram a map
//this uses the gsl histogramming package
//image - the image to be histogrammed, apply the coverage cut first
//nbins - the number of bins in the output histogram
//binloc - the lower range of each bin (must have size nbins)
//hist - the corresponding histogram values (must have size nbins)
bool histogramImage(MatDoub &image, int nbins, VecDoub &binloc, VecDoub &hist)
{
  //allocate memory for the histogram, binloc, and hist
  binloc.resize(nbins);
  hist.resize(nbins);
  gsl_histogram *h = gsl_histogram_alloc(nbins);

  //find the min and max of image and set the ranges
  double min=image[0][0];
  double max=image[0][0];
  for(int i=0;i<image.nrows();i++)
    for(int j=0;j<image.ncols();j++){
      if(image[i][j] < min) min = image[i][j];
      if(image[i][j] > max) max = image[i][j];
    }

  //force the histogram to be symmetric about 0.
  double rg = (abs(min) > abs(max)) ? abs(min) : abs(max);
  gsl_histogram_set_ranges_uniform(h, -rg, rg);

  //fill up the histogram
  for(int i=0;i<image.nrows();i++)
    for(int j=0;j<image.ncols();j++){
      gsl_histogram_increment(h, image[i][j]);
    }
  for(int i=0;i<nbins;i++){
    binloc[i] = h->range[i];
    hist[i] = h->bin[i];
  }

  //free resources
  gsl_histogram_free(h);
  return 1;
}

//----------------------------- o ---------------------------------------

//this is a direct knockoff of IDL's shift function
//n is the index shift value
bool shift(VecDoub &vec, int n){
  int nx = vec.size();
  VecDoub vec2(nx);
  for(int i=0;i<nx;i++){
  	int ti = (i+n)%nx;
	int shifti = (ti < 0) ? nx+ti : ti;
	vec2[shifti] = vec[i];
  }
  for(int i=0;i<nx;i++) vec[i] = vec2[i];
  return 1;
}

//same but for a matrix
//n1 and n2 are the index shifting values in each dim
bool shift(MatDoub &mat, int n1, int n2){
  int nx = mat.nrows();
  int ny = mat.ncols();
  MatDoub mat2(nx, ny);
  for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++){
      int ti = (i+n1)%nx;
      int tj = (j+n2)%ny;
      int shifti = (ti < 0) ? nx+ti : ti;
      int shiftj = (tj < 0) ? ny+tj : tj;
      mat2[shifti][shiftj] = mat[i][j];
    }
  for(int i=0;i<nx;i++) for(int j=0;j<ny;j++) mat[i][j] = mat2[i][j];
  return 1;
}


//same yet again but returns just a single value of the shifted
//matrix at location (i,j)
double shift(MatDoub &mat, int n1, int n2, int i, int j){
  int nx = mat.nrows();
  int ny = mat.ncols();
  int ti = (i+n1)%nx;
  int tj = (j+n2)%ny;
  int shifti = (ti < 0) ? nx+ti : ti;
  int shiftj = (tj < 0) ? ny+tj : tj;
  return mat[shifti][shiftj];
}



//----------------------------- o ---------------------------------------

// Simple derivation algorithm
// Uses Ridders's method to estimate the derivate of a sequence of non-uniform tabulated data
// First element  is calculated by the forward approximation
// Last element is calculated using the backwards approximation
// x - Independent variable
// y - dependent variable
// Returns VecDoub pointer to the derivate values 
VecDoub* derivate (VecDoub x, VecDoub y)
{
  if (x.size() != y.size ())
    return NULL;
  long len = x.size();
  if (len < 3){
    cerr<<"Derivate Error. Input Array must have at least three elements to compute the derivates"<<endl;
    return NULL;
  }
  VecDoub *deriv = new VecDoub(len);
  for (long i =1; i< len-1; i++)
    (*deriv)[i] = (y[i+1]-y[i-1])/(x[i+1]-x[i-1]);
  
  (*deriv)[0] = (y[1]-y[0])/(x[1]-y[0]);
  (*deriv)[len-1] = (y[len-1]-y[len-2])/(x[len-1]-y[len-2]);
  
  return deriv;
}

// Simple derivation algorithm
// Uses Element to element differences to estimate the local angle of a 2D (az-el) trajectory
// Last element is calculated using the backwards approximation
// x - az variable
// y - el variable
// Returns VecDoub pointer to the derivate values 
VecDoub* getAngle (VecDoub x, VecDoub y)
{
  if (x.size() != y.size ())
    return NULL;
  long len = x.size();
  if (len < 3){
    cerr<<"Angle Estimation Error. Input Array must have at least three elements to compute angle"<<endl;
    return NULL;
  }
  VecDoub *angle = new VecDoub(len);
  for (long i =0; i< len-1; i++)
    (*angle)[i] = atan2(y[i+1]-y[i],x[i+1]-x[i]);
  
  (*angle)[len-1] = atan2(y[len-1]-y[len-2],x[len-1]-y[len-2]);
  
  return angle;
}


//----------------------------- o ---------------------------------------


///writes a MatDoub object to a netcdf file 
bool writeMatDoubToNcdf(MatDoub &mat, string ncdfFilename)
{
  //create the file
  NcFile ncfid = NcFile(ncdfFilename.c_str(), NcFile::Replace);
  if (!ncfid.is_valid()){
    cerr << "writeMatDoubToNcdf(): Couldn't open netcdf file for writing!\n";
    return 0;
  }

  //create dimensions
  NcDim* rowDim = ncfid.add_dim("nrows", mat.nrows());
  NcDim* colDim = ncfid.add_dim("ncols", mat.ncols());

  //define variables for maps
  NcVar *matVar = ncfid.add_var("mat", ncDouble, rowDim, colDim);

  //and write the maps
  matVar->put(&mat[0][0], mat.nrows(), mat.ncols());
  cerr << "writeMatDoubToNcdf(): Matrix written to ";
  cerr << ncdfFilename << endl;

  return 1;
}

size_t robustMedian(double* arr, size_t nSamp, double cutStd, double* outMedian,
		double* outDev) {

	*outMedian = median (arr,nSamp);
	*outDev = stddev(arr, nSamp, *outMedian);

	bool flags [nSamp];
	size_t nRemoved;
	for (size_t i=0; i<nSamp; i++)
		flags[i]= true;
	size_t ngood = nSamp;
	do{
		nRemoved = 0;
		for (size_t i=0; i< nSamp; i++)
			if (abs(arr[i]-*outMedian)/ *outDev > cutStd && flags[i]){
				nRemoved++;
				ngood--;
				flags[i]= false;
			}
		if (nRemoved == 0)
			break;
		double goodData [ngood];

		size_t ii=0;
		for (size_t i=0; i< nSamp; i++)
					if (flags[i])
						goodData[ii++]= arr[i];
		*outMedian = median (goodData,ngood);
		*outDev = stddev (goodData, ngood, *outMedian);

	}while (1);

	return ngood;
}

bool writeGslMatrix(const char* outFile, void *m){
	  gsl_matrix * matrix = (gsl_matrix*)m;
	  ofstream out(outFile);
	  if(out.bad()){
	    cerr << "vector_utilities::writeVecGslMatrix():";
	    cerr << " Error opening " << outFile << endl;
	    return 0;
	  }

	  //set precision to something reasonable for det values
	  out.precision(14);

	  //and the output
	  for(size_t i=0;i<matrix->size1;i++)
		  for(size_t j=0;j<matrix->size2;j++){
			  out << gsl_matrix_get(matrix,i,j);
			  if (j < matrix->size2-1)
				  out<<" ";
			  else
				  out<<endl;

		  }
	  out.close();
	  cerr<<"Written file:"<<outFile<<endl;
	  return 1;
}

//----------------------------- o ---------------------------------------


///writes a MatDoub object to a netcdf file 
bool writeVecDoubToNcdf(VecDoub &vec, string ncdfFilename)
{
  //create the file
  NcFile ncfid = NcFile(ncdfFilename.c_str(), NcFile::Replace);
  if (!ncfid.is_valid()){
    cerr << "writeVecDoubToNcdf(): Couldn't open netcdf file for writing!\n";
    return 0;
  }

  //create dimensions
  NcDim* sizeDim = ncfid.add_dim("size", vec.size());

  //define variable for vector
  NcVar *vecVar = ncfid.add_var("vec", ncDouble, sizeDim);

  //and write the maps
  vecVar->put(&vec[0], vec.size());
  cerr << "writeVecDoubToNcdf(): Vector written to ";
  cerr << ncdfFilename << endl;

  return 1;
}



//----------------------------- o ---------------------------------------


///finds a weight threshold for a given coverage cut value
///this is a translation of the IDL technique
double findWeightThreshold(MatDoub &myweight, double cov)
{
  int nr = myweight.nrows();
  int nc = myweight.ncols();
  vector<double> og;
  for(int x = 0; x < nr; x++){
    for(int y = 0; y < nc; y++){
      if(myweight[x][y] > 0.){
	og.push_back(myweight[x][y]);
      }
    }
  }
  //using gsl vector sort routines to do the coverage cut
  //so we need to repack the maps
  //start with the weight map (keep idl utils nomenclature)
  
  //find the point where 25% of nonzero elements have greater weight
  double covlim;
  int covlimi;
  covlimi = 0.75*og.size();
  covlim = select(og, covlimi);
  double mval;
  double mvali;
  mvali = floor((covlimi+og.size())/2.);
  mval = select(og, mvali);
  
  //return the weight cut value
  return cov*mval;
}


int linearDeviates (int m, int n,  double *p, double *deviates, double **derivs, void * private_data){
	mpfit_basic_data *data = (mpfit_basic_data *) private_data;
	for (int i =0; i< m; i++)
		deviates[i] = abs((data->y[i]-p[0])/p[1]- data->x[i]);

	return 0;
}

//double linfit_flags (const double *x,const bool *flagsx, const double *y, const bool *flagsy, size_t nSamples, double *c0, double *c1, bool useMpfit){
	//double tol = 1e-6;
	//double chisq = 0;
	//size_t rank = 1;
	//gsl_vector  *weights = gsl_vector_alloc(nSamples);
	//gsl_matrix *X = gsl_matrix_alloc(nSamples,2);
	//gsl_vector *vy = gsl_vector_alloc(nSamples);
	//gsl_matrix *cov = gsl_matrix_alloc (2,2);
	//gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(nSamples,2);

	//gsl_vector *c =gsl_vector_alloc(2);
	//gsl_vector_set (c,0,0.0);
	//gsl_vector_set (c,1,1.0);

	//gsl_matrix_set_all (X,1.0);
	//gsl_vector_set_all(weights,1.0);
	//for (size_t i = 0; i<nSamples; i++){
		//if (!flagsx[i] || !flagsy[i])
			//gsl_vector_set (weights,i,1e-3);
		//gsl_matrix_set (X,i,1,x[i]);
		//gsl_vector_set (vy,i,y[i]);
	//}

	//gsl_multifit_wlinear_svd(X,weights,vy,tol,&rank,c,cov,&chisq,work);

	//*c0 = gsl_vector_get (c,0);
	//*c1 = gsl_vector_get (c,1);


	//gsl_multifit_linear_free(work);
	//gsl_matrix_free (X); gsl_matrix_free(cov);
	//gsl_vector_free (weights); gsl_vector_free (vy); gsl_vector_free (c);
	//return chisq;

//}


//Uses mpfit to provide a linear fit for data x,y ignoring flagged data
double linfit_flags (const double *x,const bool *flagsx, const double *y, const bool *flagsy, size_t nSamples, double *c0, double *c1, bool useMpfit){

	size_t ngood=0;


		for (size_t i =0; i<nSamples; i++){
			if (flagsx[i] && flagsy[i])
				ngood++;
		}

	if (ngood == 0){
		cerr<<"linfit_flags()::No good valid points"<<endl;
		*c0=0.0;
		*c1=1.0;
		return -1.0;
	}
	//cerr <<"linfit_flags():: Using a factor points of: "<<double(ngood)/double(nSamples)<<endl;
	//Copy data
	double *newx = new double [ngood];
	double *newy = new double [ngood];

	size_t iSamples = 0;


	for (size_t i =0; i<nSamples; i++){
		if (flagsx[i] && flagsy[i]){
			newx[iSamples] = x[i];
			newy[iSamples++] = y[i];
		}
	}
	double tResult;
	if (!useMpfit){
		double c00, c01, c11, sumsq;

		//gsl_matrix *X = gsl_matrix_alloc (nSamples,2);
		//gsl_matrix_set_all (X,1.0);

		//gsl_multifit_linear_svd()
		gsl_fit_linear(newx,1,newy,1,ngood,c0,c1,&c00,&c01,&c11,&sumsq);
		tResult = sumsq;
	}else{
		double pars [2] = {0.0,1.0};
		mpfit_basic_data bd;
		bd.x = newx;
		bd.y = newy;
		mp_result result;
		memset(&result,0,sizeof(result));
		int status = mpfit (&linearDeviates,ngood, 2,pars,0,0,&bd,&result);
		*c0 = pars[0];
		*c1 = pars[1];
		if (status < 0){
			*c0=0.0;
			*c1=1.0;
			cerr<<"Bad fit in data"<<endl;
		}
		tResult = status;
	}
	delete [] newx;
	delete [] newy;


	return tResult;
	//return sumsq;
}

double linfit_flags (const double *x,const double *flagsx, const double *y, const double *flagsy, size_t nSamples, double *c0, double *c1, bool useMpfit){
	bool *fx = new bool [nSamples];
	bool *fy = new bool [nSamples];

	for (size_t i=0; i<nSamples; i++){
		fx[i]=flagsx[i]!=0?true:false;
		fy[i]=flagsy[i]!=0?true:false;
	}

	double retval = linfit_flags (x,fx,y,fy,nSamples,c0,c1,useMpfit);

	delete []fx;
	delete []fy;

	return retval;

}



double flagCorrelation(const double* x, const double* fx, const double* y,const double* fy, size_t nSamples) {
	size_t ngood=0;
	//double *xx, *yy;

	for (size_t i=0; i< nSamples; i++)
		if (fx[i] != 0 && fy[i] != 0)
			ngood++;

	if (ngood ==0){
		cerr<<"Not valid scan"<<endl;
		return 0.0;
	}

	double xx[ngood];
	double yy[ngood];
	size_t ii = 0;


	for (size_t i=0; i< nSamples; i++)
		if (fx[i] !=0 && fy[i] != 0){
			xx[ii]= x[i];
			yy[ii++]= y[i];
		}

	double correlation = gsl_stats_correlation(xx,1,yy,1 ,ngood);

	//delete [] xx;
	//delete [] yy;

	return correlation;
}

double flagCorrelation(const double* x, const bool* fx, const double* y,const bool* fy, size_t nSamples) {
	size_t ngood=0;
	//double *xx, *yy;

	for (size_t i=0; i< nSamples; i++)
		if (fx[i] && fy[i] )
			ngood++;

	if (ngood ==0){
		cerr<<"Not valid scan"<<endl;
		return 0.0;
	}

	double xx[ngood];
	double yy[ngood];
	size_t ii = 0;


	for (size_t i=0; i< nSamples; i++)
		if (fx[i] !=0 && fy[i] != 0){
			xx[ii]= x[i];
			yy[ii++]= y[i];
		}

	double correlation = gsl_stats_correlation(xx,1,yy,1 ,ngood);

	//delete [] xx;
	//delete [] yy;

	return correlation;
}


double flagCovariance(const double* x, const double* fx, const double* y,const double* fy, size_t nSamples) {

	double xx[nSamples];
	double yy[nSamples];

	for (size_t i=0; i< nSamples; i++)
		if (fx[i] !=0 && fy[i] != 0){
			xx[i]= x[i];
			yy[i]= y[i];

		}else
			xx[i]=yy[i]=0.0;

	double correlation = gsl_stats_covariance(xx,1,yy,1 ,nSamples);


	return correlation;
}


double median (VecDoub &arr){
	return median (&arr[0], arr.size());
}

