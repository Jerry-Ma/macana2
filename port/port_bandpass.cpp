#include <port_bandpass.h>
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;

namespace generic{


double eigen_convolve_bandpass(const Tensor<float, 1, ColMajor> &dlt_tensor,
                               const Tensor<float, 1, ColMajor> &kernel_tensor,
                               const Eigen::array<ptrdiff_t, 1> &dims,
                               Tensor<float, 1, ColMajor> &conv_tensor){
    //Start the clocks!
    auto t1 = std::chrono::high_resolution_clock::now();

    //Convolve it!
    conv_tensor = dlt_tensor.convolve(kernel_tensor,dims);

    //Time out!
    auto t2 = std::chrono::high_resolution_clock::now();

    return std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
}

double macana_bandpass(VecDoub &hValues, const VecDoub &digFiltTerms, VecDoub &storage, const int &th,
                       const size_t &nSamples, const size_t &nTerms, const size_t &nCoef, int &tt,
                       double &tmp){

    //Start the clocks!
    auto t1 = std::chrono::high_resolution_clock::now();

    //Convolve it!
    for(int t=nTerms;t<th;t++){
        tmp=0.;
        tt = t-nTerms;
        for(int i = 0;i<nCoef;i++){
            tmp += hValues[tt+i]*digFiltTerms[i];
        }
        storage[t] = tmp;
    }

    //Replace old values
    for(int i=nTerms;i<th;i++) {
      hValues[i] = storage[i];
    }

    //Time out!
    auto t2 = std::chrono::high_resolution_clock::now();

    return std::chrono::duration_cast<std::chrono::milliseconds>(t1-t2).count();
}

} //namespace
