//Gtest
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "port_bandpass.h"

namespace{

using namespace Eigen;
using namespace bandpass;

class BandpassTest:public::testing::Test
{
   public:
      size_t nSamples=1000;
      size_t nTerms=32;
      size_t nCoef=2*nTerms+1;;
      double lower=0.0;
      double upper=10.0;
      double sigma=0.5;
      double center=(upper+lower)/2.0;
      double tmp;
      int th=nSamples-nTerms;
      int tt;

      double err = 1e-5;
      
      BandpassTest() {};
      ~BandpassTest() override {}
      void TearDown() override {}
      void SetUp() override {}

      template <typename VectorType>
         VectorType generate_ker(VectorType& dlt_func, int PeakLoc, int size){
            dlt_func.setZero(nSamples);
            dlt_func[PeakLoc] = 1.0;
            VectorType kernel(size);
            ArrayXf x = ArrayXf::LinSpaced(size, lower, upper);

            kernel = exp(-pow(x.array() - center,2.)/(2.*pow(sigma,2.0)));
            return kernel;
         }
};

using namespace testing;

MATCHER_P(NearWithPrecision, precision, "") {
    return abs(get<0>(arg) - get<1>(arg)) < precision;
}

TEST_F(BandpassTest,EigenFFTTest){
   VectorXf dlt_func, convolution;
   VectorXcf fft_dlt_func, fft_kernel, fft_dlt_kernel;

   VectorXf kernel = generate_ker<VectorXf>(dlt_func, 0, nSamples);
   EigenFFT(dlt_func, kernel, fft_dlt_func, fft_kernel, fft_dlt_kernel,convolution);

   std::vector<double> c(convolution.data(),convolution.data() + convolution.size());
   std::vector<double> k(kernel.data(),kernel.data() + kernel.size());

   EXPECT_THAT(c, Pointwise(NearWithPrecision(err), k));
}

TEST_F(BandpassTest,EigenTensorTest){
   VectorXf dlt_func;

   //Eigen Tensor
   Tensor<float, 1, ColMajor> conv_tensor(static_cast<Index>(nSamples - nCoef + 1));

   Eigen::array<ptrdiff_t, 1> dims{0};

   VectorXf kernel = generate_ker<VectorXf>(dlt_func, nCoef-1, nCoef);

   Eigen::Tensor<float, 1> dlt_tensor = Vector_to_Tensor<VectorXf>(dlt_func, nSamples);
   Eigen::Tensor<float, 1> kernel_tensor = Vector_to_Tensor<VectorXf>(kernel, nCoef);

   EigenTensor(dlt_tensor, kernel_tensor, dims, conv_tensor);

   //for(int i=0;i<nCoef;i++) cerr << conv_tensor[i] << "," << kernel_tensor[i] << "\n";

   std::vector<double> c(conv_tensor.data(),conv_tensor.data() + nCoef);
   std::vector<double> k(kernel_tensor.data(),kernel_tensor.data() + nCoef);

   EXPECT_THAT(c, Pointwise(NearWithPrecision(err), k));
}

TEST_F(BandpassTest,MacanaTest){
    VectorXf dlt_func, convolution;
    VectorXcf fft_dlt_func, fft_kernel, fft_dlt_kernel;

    VectorXf kernel = generate_ker<VectorXf>(dlt_func, nCoef - 1, nCoef);

    VecDoub hValues = Vector_to_VecDoub<VectorXf>(dlt_func, nSamples);
    VecDoub digFiltTerms = Vector_to_VecDoub<VectorXf>(kernel, nCoef);

    VecDoub storage(nSamples);

    MacanaBandpass(hValues, digFiltTerms, storage, th, nSamples, nTerms, nCoef, tt, tmp);

        std::vector<double> c(nCoef);
        std::vector<double> k(nCoef);

    for(int i=nTerms;i<nCoef + nTerms;i++) {
            c[i-nTerms] = hValues[i];
            k[i-nTerms] = digFiltTerms[i-nTerms];
    }

    EXPECT_THAT(c, Pointwise(NearWithPrecision(err), k));
}


} //namespace
