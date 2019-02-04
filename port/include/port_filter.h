#ifndef PORT_FILTER_H
#define PORT_FILTER_H

//Eigen Includes
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/CXX11/Tensor>

//Macana Includes
#include "nr3.h"
#include <fftw3.h>
#include <gsl/gsl_sf_bessel.h> //Needed for filter

namespace filter{

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

//Declare function for Macana Filter
bool digiFilter(const double fLow, const double fHigh, const double aGibbs,
                   const int nTerms, VecDoub &coefOut);

//Declare function for Macana Convolution
void MacanaFilter(VecDoub &hValues, const VecDoub &digFiltTerms, VecDoub &storage,
                     const int nSamples, const int nTerms, const int nCoef);

/*----------------------------------------------------------------------------------------------*/

//Function to convert from Eigen::VectorXf to Eigen::Tensor
template<typename Scalar1>
auto Vector_to_Tensor(Scalar1 vec, int size){
   auto mapped_t = Eigen::TensorMap<Eigen::Tensor<float, 1, Eigen::RowMajor>>(&vec.data()[0], size);
   Eigen::Tensor<float, 1> t = Eigen::TensorLayoutSwapOp<Eigen::Tensor<float, 1, Eigen::RowMajor>>(mapped_t);
   return t;
}

/*----------------------------------------------------------------------------------------------*/

//Function to generate a time series of random data points.
template <typename Scalar>
void MakeData(Scalar &hValues, const int nSamples, const double mean, const double stddev){

    std::mt19937 generator;
    std::normal_distribution<double> dist(mean, stddev);

    // Add Gaussian noise
    for(int i=0; i<nSamples; i++){
        hValues[i] = hValues[i] + dist(generator);
    }
}
/*----------------------------------------------------------------------------------------------*/

//Filter function
template<typename Scalar1, typename Scalar2,typename Scalar3,
         typename Scalar4, typename Scalar5>
bool EigenDigitalFilter(const Scalar1 fLow, const Scalar2 fHigh, const Scalar3 aGibbs,
                        const Scalar4 nTerms, Eigen::DenseBase<Scalar5> &coefOut)
{
    //Band stop?
    double fStop = (fHigh < fLow) ? 1. : 0.;

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

    //Arg
    VectorXf arg(nTerms);
    arg.setLinSpaced(nTerms,1,nTerms);
    arg = arg/nTerms;

    //Coef
    VectorXf coef(nTerms);
    for(int i=0;i<nTerms;i++){
     coef[i] = gsl_sf_bessel_I0(alpha*sqrt(1.-pow(arg[i],2))) /
     gsl_sf_bessel_I0(alpha);
    }

    //t
    VectorXf t(nTerms);
    t.setLinSpaced(nTerms,1,nTerms);
    t = t*M_PI;

    //here we go
    coef = coef.array()*(sin(t.array()*fHigh)-sin(t.array()*fLow))/t.array();

    //build the uncentered version
    coefOut.head(nTerms) = coef.reverse();
    coefOut[nTerms] = fHigh-fLow-fStop;
    coefOut.tail(nTerms) = coef;

    //renormalize
    double area = 0;
    area = coefOut.sum();
    //coefOut = coefOut/area;

    return 1;
}

/*----------------------------------------------------------------------------------------------*/

//Convoltuion code using Eigen FFT
template <typename Scalar1, typename Scalar2, typename Scalar3>
void EigenFFTFilter(const Eigen::MatrixBase<Scalar1>& hValues,
                    const Eigen::MatrixBase<Scalar2>& digitalFilterTerms,
                    Eigen::MatrixBase<Scalar3>& convolution)
{
   //do forward fft
   Eigen::FFT<float> fft;

   Eigen::VectorXcf fft_hValues(hValues.size());
   Eigen::VectorXcf fft_digitalFilterTerms(digitalFilterTerms.size());
   Eigen::VectorXcf fft_convolution(convolution.size());

   fft.fwd(fft_hValues,hValues.eval());
   fft.fwd(fft_digitalFilterTerms,digitalFilterTerms.eval());

   //Multiply the FFT conjugates
   fft_convolution = fft_hValues.array() * fft_digitalFilterTerms.array();

   //Take the inverse FFT of the product of FFT conjugates
   fft.inv(convolution, fft_convolution);
}

/*----------------------------------------------------------------------------------------------*/

//Convoltuion code using Eigen FFT
template <typename Scalar1, typename Scalar2>
void EigenTensorFilter(const Scalar1 &hValues,
                       const Scalar1 &digitalFilterTerms,
                       const Scalar2 &dims,
                       Scalar1 &convolution){

    convolution = hValues.convolve(digitalFilterTerms,dims);
}

} //namespace
#endif // PORT_FILTER_H
