#ifndef PORT_DOWNSAMPLE_H
#define PORT_DOWNSAMPLE_H

//Eigen Includes
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/CXX11/Tensor>

//Macana Includes
#include <fstream>

namespace downsample{

using namespace Eigen;
using namespace std;

//Write vectors to file
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

template<typename Scalar1>
Scalar1 downsampler(Eigen::DenseBase<Scalar1>& hValues, int dsf){
    return Eigen::Map<Eigen::VectorXf, 0, Eigen::InnerStride<Eigen::Dynamic> >(&hValues[0], hValues.size()/dsf, Eigen::InnerStride<Eigen::Dynamic>(dsf));
}

} //namespace
#endif // PORT_DOWNSAMPLE_H
