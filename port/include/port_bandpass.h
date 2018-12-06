#ifndef BANDPASS_H
#define BANDPASS_H

//Eigen
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/CXX11/Tensor>

#define EIGEN_USE_MKL_ALL

//MKL
#include <mkl_dfti.h>

//Other
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <string>
#include <fftw3.h>

#include "nr3.h"
#include "tinyxml2.h"
#include "Detector.h"
#include "astron_utilities.h"
#include "vector_utilities.h"

#include <chrono>

namespace generic {

using namespace Eigen;

template<typename T>
using  MatrixType = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>;

//Test template to convert from Tensor to Matrix
template<typename Scalar,int rank, typename sizeType>
auto Tensor_to_Matrix(const Eigen::Tensor<Scalar,rank> &tensor,const sizeType rows,const sizeType cols)
{
    return Eigen::Map<const MatrixType<Scalar>> (tensor.data(), rows,cols);
}

//Test template to convert from Matrix to Tensor
template<typename Scalar, typename... Dims>
auto Matrix_to_Tensor(const MatrixType<Scalar> &matrix, Dims... dims)
{
    constexpr int rank = sizeof... (Dims);
    return Eigen::TensorMap<Eigen::Tensor<const Scalar, rank>>(matrix.data(), {dims...});
}

//Eigen FFT Template
template <typename Derived1, typename Derived2, typename Derived3,
          typename Derived4,typename Derived5, typename Derived6>
double eigen_fft_bandpass(const MatrixBase<Derived1>& dlt_func,
                          const MatrixBase<Derived2>& kernel,
                          MatrixBase<Derived3>& fft_dlt_func,
                          MatrixBase<Derived4>& fft_kernel,
                          MatrixBase<Derived5>& fft_dlt_kernel,
                          MatrixBase<Derived6>& convolution){

    FFT<float> fft;

    //Start the clocks!
    auto t1 = std::chrono::high_resolution_clock::now();

    //Carry out the forward FFTs
    fft.fwd(fft_dlt_func, dlt_func);
    fft.fwd(fft_kernel, kernel);

    //Multiply the FFT conjugates
    fft_dlt_kernel = fft_dlt_func.array() * fft_kernel.array();

    //Take the inverse FFT of the product of FFT conjugates
    fft.inv(convolution, fft_dlt_kernel);

    //Time out!
    auto t2 = std::chrono::high_resolution_clock::now();

    return std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
}


double eigen_convolve_bandpass(const Tensor<float, 1, ColMajor> &dlt_tensor,
                               const Tensor<float, 1, ColMajor> &kernel_tensor,
                               const Eigen::array<ptrdiff_t, 1> &dims,
                               Tensor<float, 1, ColMajor> &conv_tensor);

double macana_bandpass(VecDoub &hValues, const VecDoub &digFiltTerms, VecDoub &storage, const int &th,
                       const size_t &nSamples, const size_t &nTerms, const size_t &nCoef, int &tt,
                       double &tmp);

} // namespace generic

#endif // BANDPASS_H
