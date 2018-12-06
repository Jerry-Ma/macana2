//Google test
#include <gtest/gtest.h>
#include <gmock/gmock.h>

//Eigen
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/CXX11/Tensor>

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

#include <port_bandpass.h>

namespace {

using namespace Eigen;
using namespace std;
using namespace generic;

class kernel_test: public ::testing::Test{
protected:
    //Data length
    size_t nSamples=1000;

    //Gaussian kernel parameters
    size_t nTerms=32;
    size_t nCoef=2*nTerms+1;;
    double lower=0.0;
    double upper=10.0;
    double sigma=0.5;
    double center=(upper+lower)/2.0;

    //Macana bandpass params
    double tmp;
    int th=nSamples-nTerms;
    int tt;

    //Eigen Vectors
    VectorXf dlt_func{nSamples};
    VectorXcf fft_dlt_func{nSamples};

    VectorXf kernel{nSamples};
    VectorXcf fft_kernel{nSamples};

    VectorXcf fft_dlt_kernel{nSamples};
    VectorXf convolution{nSamples};

    ArrayXf x = ArrayXf::LinSpaced(nSamples, lower, upper);

    //Macana Vectors
    VecDoub hValues{nSamples};
    VecDoub digFiltTerms{nCoef};
    VecDoub xm{nCoef};
    VecDoub storage{nSamples};

    //Eigen Tensors
    Tensor<float, 1, ColMajor> dlt_tensor{static_cast<Index>(nSamples)};
    Tensor<float, 1, ColMajor> kernel_tensor{static_cast<Index>(nCoef)};
    Tensor<float, 1, ColMajor> conv_tensor{static_cast<Index>(nSamples - nCoef + 1)};


    Eigen::array<ptrdiff_t, 1> dims{0};

    kernel_test() {};
    ~kernel_test() override {}
    void TearDown() override {}
    void SetUp() override {

    //Eigen Vector setup
    /*-----------------------------------------------------------------------------------*/
    dlt_func.setZero(nSamples);
    convolution.setZero(nSamples);
    dlt_func[0] = 1.0;

    kernel = exp(-pow(x.array() - center,2.)/(2.*pow(sigma,2.0)));

    //Macana setup
    /*-----------------------------------------------------------------------------------*/
    xm[0] = lower;

    for (int i = 1; i<nCoef; i++){
        xm[i] = xm[i-1] + (upper - lower)/nCoef;
     }

    for(int i = 0; i<nSamples; i++){
        hValues[i] = 0.0;
     }

    hValues[nCoef] = 1.0;

     for(int i = 0; i<nCoef; i++){
         digFiltTerms[i] = exp(-pow(xm[i] - center,2.)/(2.*pow(sigma,2.0)));
     }

     //Eigen Tensor Setup
     /*-----------------------------------------------------------------------------------*/
     dlt_tensor.setConstant(0.0);
     dlt_tensor[nCoef] = 1.0;
     dims[0] = 0.0;

     for(int i=0;i<nCoef; i++){
         kernel_tensor[i] = exp(-pow(xm[i] - center,2.)/(2.*pow(sigma,2.0)));
     }
    }

};

MATCHER_P(NearWithPrecision, precision, "") {
    return abs(get<0>(arg) - get<1>(arg)) < precision;
}

using namespace testing;

//Testing the Eigen FFT routine.
TEST_F(kernel_test, eigen_fft_bandpass) {
    double err = 1e-5;
    double time = 0;

    time = eigen_fft_bandpass(dlt_func, kernel, fft_dlt_func, fft_kernel, fft_dlt_kernel,
                              convolution);

    std::vector<double> conv(nSamples);
    std::vector<double> ker(nSamples);

    for(int bb=0;bb<nSamples;bb++){
        conv[bb] = convolution[bb];
        ker[bb] = kernel[bb];
    }

    EXPECT_THAT(conv, Pointwise(NearWithPrecision(err), ker));
}

//Testing the Eigen convolve routine
TEST_F(kernel_test, eigen_convolve_bandpass) {
    double err = 1e-5;
    double time;

    time = eigen_convolve_bandpass(dlt_tensor, kernel_tensor, dims, conv_tensor);

    std::vector<double> conv(nCoef);
    std::vector<double> ker(nCoef);

    for(int bb=0;bb<nCoef;bb++) {
        conv[bb-0] = conv_tensor[bb];
        ker[bb] = kernel_tensor[bb];
    }

    EXPECT_THAT(conv, Pointwise(NearWithPrecision(err), ker));
}

//Testing the macana convolve routine
TEST_F(kernel_test, macana_bandpass) {
    double err = 1e-5;
    double time = 0;

    time = macana_bandpass(hValues, digFiltTerms, storage, th, nSamples, nTerms, nCoef, tt, tmp);

    std::vector<double> conv(nCoef);
    std::vector<double> ker(nCoef);

    for(int bb=nTerms;bb<nCoef + nTerms;bb++) {
        conv[bb-nTerms] = hValues[bb];
        ker[bb-nTerms] = digFiltTerms[bb-nTerms];
        //cerr << hValues[bb] << "\n";
    }

    EXPECT_THAT(conv, Pointwise(NearWithPrecision(err), ker));
}
} //namespace

