# Macana2

## Build

### MacOS

```plain
brew tap brewsci/science
brew install libopm suitesparse gsl fftw ccfits netcdf
# this is necessary because suitesparse headers are not linked by default
ln -s /usr/local/Cellar/suite-sparse/5.3.0/include /usr/local/include/suitesparse
build 
git clone https://github.com/Jerry-Ma/macana2.git
cd macana2
mkdir build
cd build
cmake ..
```

## Run

```
./{build_dir}/bin/beammap <apb file>
```

## Run test
    TBD

## Documentation
    TBD
