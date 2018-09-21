#ifndef _vector_utilities_h_
#define _vector_utilities_h_

//#include <gsl/gsl_matrix.h>
#include "nr3.h"
#include <iostream>



double select(vector<double> input, int index);
double mean(double* arr, int nSamp);
double mean(double* arr, bool *flags, int nSamp);
double mean(VecDoub &arr);
double mean(VecDoub &arr, VecBool &flags);
double mean(VecDoub &arr, int start, int end);
double mean(VecDoub &arr, VecBool &flags, int start, int end, double *ncount = NULL);
double mean(double *arr, int start, int end);
double stddev(double* arr, int nSamp, double mn);
double stddev(VecDoub &arr);
double sttdev(VecDoub &arr, VecBool & flags);
double stddev(VecDoub &arr, int start, int end);
double stddev(VecDoub &arr, VecBool &flags, int start, int end, double *ncount = NULL);
double medabsdev(VecDoub &arr);
double median (double*arr, size_t nSamp);
double median (double *data, bool *flags, size_t nSamp);
double median (VecDoub &arr);
size_t robustMedian (double *arr, size_t nSamp, double cutStd, double *outMedian, double *outDev);
double percentile (double *data, size_t nSamp, double percentile);
void maxmin(VecDoub &arr, double* max, double* min);
void maxmin(VecInt &arr, int* max, int* min);
void maxmin(double* arr, double* max, double* min, int npts);
void smooth(VecDoub &inArr, VecDoub &outArr, int w);
void smooth_edge_truncate(VecDoub &inArr, VecDoub &outArr, int w);
//bool writeVecOut(const char* outFile, double* outData, int nOut);
//template <class T> bool writeVecOut(const char* outFile,  T* outData, size_t nOut);
bool writeMatOut(const char* outFile, MatDoub &outMat);
bool digitalFilter(double fLow, double fHigh, double aGibbs, 
		   int nTerms, double* coefOut);
bool interpolateLinear(VecDoub xx, VecDoub yy, int nUnique,
		       double* x, double *result, int nSamples);
bool deNanSignal(double* sig, int nSamples);
bool removeDropouts(double* sig, int nSamples);
bool correctRollover(double* sig, double low, double high, 
		     double ulim, int nSamples);
bool correctOutOfBounds(double* sig, double low, double high, int nSamples);
MatDoub hanning(int n1in, int n2in);
VecDoub hanning(int n1in);
bool histogramImage(MatDoub &image, int nbins, VecDoub &binloc, VecDoub &hist);
bool shift(VecDoub &vec, int n);
bool shift(MatDoub &mat, int n1, int n2);
double shift(MatDoub &mat, int n1, int n2, int i, int j);
VecDoub* derivate (VecDoub x, VecDoub y);
VecDoub* getAngle (VecDoub x, VecDoub y);
bool writeMatDoubToNcdf(MatDoub &mat, string ncdfFilename);
bool writeVecDoubToNcdf(VecDoub &vec, string ncdfFilename);
void convolve (double *data, size_t nData, double *kernel, bool first);
double findWeightThreshold(MatDoub &myweight, double cov);
double linfit_flags (const double *x, const bool *flagsx, const double *y, const bool *flagsy, size_t nSamples, double *c0, double *c1, bool useMpfit);
double linfit_flags (const double *x, const double *flagsx, const double *y, const double *flagsy, size_t nSamples, double *c0, double *c1, bool useMpfit);
double flagCorrelation (const double *x, const double *fx, const double *y, const double *fy, size_t nSamples);
double flagCorrelation(const double* x, const bool* fx, const double* y,const bool* fy, size_t nSamples);
double flagCovariance(const double* x, const double* fx, const double* y,const double* fy, size_t nSamples);
bool writeGslMatrix(const char * outFile, void *matrix);
//gsl_matrix *castMatDoub (MatDoub matrix);
struct Base_interp
{
  Int n, mm, jsav, cor, dj;
  const Doub *xx, *yy;
  Base_interp(VecDoub_I &x, const Doub *y, Int m)
    :n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y){
    dj=max(1,(int)pow((Doub)n,0.125));
  }

  Doub interp(Doub x){
    Int jlo=cor ? hunt(x) : locate(x);
    return rawinterp(jlo,x);
  }

  Int locate(const Doub x);
  Int hunt(const Doub x);
  Doub virtual rawinterp(Int jlo, Doub x) = 0;
};


struct Linear_interp : Base_interp
{
  Linear_interp(VecDoub_I &xv, VecDoub_I &yv)
    : Base_interp(xv,&yv[0],2){}
  Doub rawinterp(Int j, Doub x){
    if(xx[j]==xx[j+1]) return yy[j];
    else return yy[j]+((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
  }
};


struct Bilin_interp
{
  int m,n;
  const MatDoub &y;
  Linear_interp x1terp, x2terp;

  Bilin_interp(VecDoub_I &x1v, VecDoub_I &x2v, MatDoub_I &ym)
  : m(x1v.size()), n(x2v.size()), y(ym),
    x1terp(x1v,x1v), x2terp(x2v,x2v) {}

  Doub interp(Doub x1p, Doub x2p)
  {
    int i,j;
    Doub yy, t, u;
    i = x1terp.cor ? x1terp.hunt(x1p) : x1terp.locate(x1p);
    j = x2terp.cor ? x2terp.hunt(x2p) : x2terp.locate(x2p);

    t = (x1p-x1terp.xx[i])/(x1terp.xx[i+1]-x1terp.xx[i]);
    u = (x2p-x2terp.xx[j])/(x2terp.xx[j+1]-x2terp.xx[j]);
    yy = (1.-t)*(1.-u)*y[i][j] + t*(1.-u)*y[i+1][j] + 
         (1.-t)*u*y[i][j+1] + t*u*y[i+1][j+1];

    return yy;
  }
};

template <class T> bool writeVecOut(const char* outFile,  T* outData, size_t nOut)
{
  ofstream out(outFile);
  if(out.bad()){
    cerr << "vector_utilities::writeVecOut():";
    cerr << " Error opening " << outFile << endl;
    return 0;
  }

  //set precision to something reasonable for det values
  out.precision(14);

  //and the output
  for(size_t i=0;i<nOut;i++) out << outData[i] << "\n";
  if(out.bad()){
    cerr << "vector_utilities::writeVecOut():";
    cerr << " Error writing to " << outFile << endl;
    return 0;
  }

  //flush the buffer just to be sure
  out.flush();

  //announce what you just did
  cout << "Wrote out " << outFile << endl;

  return 1;
}

typedef struct{
	double *x;
	double *y;
}mpfit_basic_data;


#endif
