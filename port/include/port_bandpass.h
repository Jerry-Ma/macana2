#ifndef PORT_BANDPASS_H
#define PORT_BANDPASS_H

//Eigen
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/CXX11/Tensor>

#define EIGEN_USE_MKL_ALL

//MKL
#include <mkl_dfti.h>

//macana
#include "nr3.h"
#include <fftw3.h>

namespace bandpass{

using namespace Eigen;

template<typename Scalar1>
auto Vector_to_Tensor(Scalar1 vec, int size){
   auto mapped_t = Eigen::TensorMap<Eigen::Tensor<float, 1, Eigen::RowMajor>>(&vec.data()[0], size);
   Eigen::Tensor<float, 1> t = Eigen::TensorLayoutSwapOp<Eigen::Tensor<float, 1, Eigen::RowMajor>>(mapped_t);
   return t;
}

template<typename Scalar1>
auto Vector_to_VecDoub(Scalar1 vec, int size){
   VecDoub t(vec.size());
   for(int i=0;i<vec.size();i++) t[i] = vec[i];
   return t;
}


template <typename Scalar1, typename Scalar2>
void EigenFFT(const MatrixBase<Scalar1>& dlt_func,
              const MatrixBase<Scalar1>& kernel,
              MatrixBase<Scalar2>& fft_dlt_func,
              MatrixBase<Scalar2>& fft_kernel,
              MatrixBase<Scalar2>& fft_dlt_kernel,
              MatrixBase<Scalar1>& convolution){
   
   //do forward fft
   Eigen::FFT<float> fft;

   fft.fwd(fft_dlt_func, dlt_func);
   fft.fwd(fft_kernel, kernel);

   //Multiply the FFT conjugates
   fft_dlt_kernel = fft_dlt_func.array() * fft_kernel.array();

   //Take the inverse FFT of the product of FFT conjugates
   fft.inv(convolution, fft_dlt_kernel);
}

template <typename Scalar1>
void EigenTensor(const Scalar1 &dlt_tensor,
                 const Scalar1 &kernel_tensor,
                 const Eigen::array<ptrdiff_t, 1> &dims,
                 Scalar1 &conv_tensor){

    conv_tensor = dlt_tensor.convolve(kernel_tensor,dims);
}

void MacanaBandpass(VecDoub &hValues, const VecDoub &digFiltTerms, VecDoub &storage,
                     const int &th, const size_t &nSamples, const size_t &nTerms,
                     const size_t &nCoef, int &tt, double &tmp);

} //namespace

#endif // PORT_BANDPASS_H
