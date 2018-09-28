# Macana2

## Install

### MacOS

    brew tap brewsci/science
    brew install libopm suitesparse gsl fftw ccfits netcdf
    # this is necessary because suitesparse headers are not linked by default
    ln -s /usr/local/Cellar/suite-sparse/5.3.0/include /usr/local/include/suitesparse
    git clone https://github.com/Jerry-Ma/macana2.git
    cd macana2
    mkdir build
    cd build
    cmake ..

## Ideas

* Git
* CMake
* setup.py ?
