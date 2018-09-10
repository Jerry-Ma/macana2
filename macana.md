
# The AzTEC C++ Pipeline

The AzTEC millimeter wavelength camera is a 144 element bolometer array configurated to carry out sensitive astronomical observations in the 1.1mm atmospheric window. It is currently installed on the Large Millimeter Telescope [(LMT)](www.lmtgtm.org) in Mexico.  The **AzTEC C++ pipeline** (codename: *macana*) is the prefered tool to convert the raw observed timestreams into science products (maps). A description of the instrument can be found in the [Wilson et al. 2008 paper](https://academic.oup.com/mnras/article/386/2/807/1056754). An overview of the data processing can be found in the [Scott et al 2012 paper](https://academic.oup.com/mnras/article/405/4/2260/1044926). 

## Installation

### Prerequisites

Macana has been extensively tested on the Long Term (LTS) versions of Ubuntu and moderately tested on another few other linux distributions. Also it is possible to install it on MacOs using the **MacPorts** system (see [this page](macanaOSX.html)). In general, you require to install following software and libraries:

 * GNUC C++ compiler (GCC)
 * OpenMP
 * NetCDF C++
 * GNU Scientific Library (GSL)
 * BLAS
 * FFTW3
 * CCFits
 * SuiteSparse (CXSparse)
 * Subversion

You can install the packages for these libraries using the APT tool on Debian based systems:

```
 sudo apt-get install build-essential libfftw3-dev libgsl0-dev libsuitesparse-dev libnetcdf-dev libccfits-dev subversion
```

Ubuntu 16.04 LTS, Debian Stretch and more recent linux distributions, a new version of the C++ NetCDF library has been released. Macana will still work using the *legacy* package:

```
 sudo apt-get install build-essential libfftw3-dev libgsl0-dev libsuitesparse-dev libnetcdf-cxx-legacy-dev libccfits-dev subversion
```

On RPM based system (like Fedora, OpenSuse, CentOS,...) you can use the YUM tool to install the compiler and libraries:

```
sudo yum install make automake gcc gcc-c++ kernel-devel fftw-devel gsl-devel suitesparse-devel CCfits-devel netcdf-cxx-devel
```

### Getting Macana

You can obtain the code from the LMT proposal system site [(www.lmtobservatory.org)](www.lmtobservatory.org) using the subversion tool:

```
svn co svn://www.lmtobservatory.org/aztec_c++
```

Now compile the code by entering to the *aztec_c++* directory and running the make utility:

```
cd aztec_c++
make
```

#### Troubleshooting

##### 	Error Fatal: fitsio.h No such file or directory
This error is frequent on RPM-based distribution when the cfitsio include files are stored in /usr/include/cfitsio instead of the default path for include files (/usr/include). There are two ways to workaround this problem:

 1. Make asymbolic link to header files (you need root privileges):
 
 ```
 sudo cp --symbolic-link /usr/include/cfitsio/* /usr/include
 ```
 2. Modify the macana Makefile to look for the cfitsio include file in the appropriate path. Open the Makefile in your favorite text editor and change the **IFLAGS** variable at line 3 to:
 
 ```
 IFLAGS=-I include/ -I /usr/include -I Sky/Novas/  -I /usr/include/cfitsio
 ```


### Setting up environment variables
In order to work correctly macana requires to set the `AZTEC_MACANA_PATH` environment variable. This variable must contain the absolute path to the folder containing macana distribution. For bash-like terminal interpreter you can add a the following  line a the end you `$HOME/.bashrc` file:

 ```
 export AZTEC_MACANA_PATH="/home/myuser/some-path/aztec_c++"
 ```
 
just replace the `/home/myuser/some-path` with the actual path where you downloaded macana. For a csh-like terminal interpreter, edit the `$HOME/.cshrc` file and add the line:

 ```
 setenv AZTEC_MACANA_PATH "/home/myuser/some-path/aztec_c++"
 ```
 
It is strongly recommended (but not necessary) to add the binaries of the AzTEC C++ pipeline to your `PATH` environment variable:

 ```
 export PATH=$PATH:$AZTEC_MACANA_PATH/bin				#bash-like interpreter
 
 setenv PATH $PATH:$AZTEC_MACANA_PATH/bin				#csh-like interpreter
 ```
 
this will allow you to run macana as a global command. 

### Setting up the python utilities
Macana include a python package useful to generate the directory structure, bolometer calibration files (bstats) and  pointing correction values necessary to produce a scientific-class image from the raw observed timestreams. These utilities require the following packages to be installed:

* numpy (1.11.0 or greater)
* scipy (0.17.0 or greater)
* matplotlib (1.5.1 or greater)
* astropy (1.1.1 or greater)
* lxml (3.5.0 or greater)

you can install this packages using the package manager for your linux distribution. For Ubuntu and other Debian-based distributions: 
```
apt-get install python-numpy python-scipy python-maptplotlib python-lxml python-astropy
```
for Fedora and other RPM-based distribution:

```
yum install python2-numpy python2-scipy python2-maptplotlib python2-lxml python2-astropy
```

After these python packages are installed you need to add the location of the macana python utilities to your PYTHONPATH environment variable:

```
export PYTHONPATH=$AZTEC_MACANA_PATH/python:$PYTHONPATH:$AZTEC_MACANA_PATH/bin				#bash-like interpreter

setenv PYTHONPATH $AZTEC_MACANA_PATH/python:$PYTHONPATH:$AZTEC_MACANA_PATH/bin				#csh-like interpreter
```

### Updating the pipeline

The pipeline is under continous development. It is strongly recommended to update your working copy from the repository:

```
cd $AZTEC_MACANA_PATH
svn update
```

## Running Macana

The input to the macana file is a xml file containing the parameter and the location of the observation files. The macana distribution contain a sample observation under the `sample_files/lmt` directory. You can run macana over the test data using the following commands:

    cd $AZTEC_MACANA_PATH
    macana apDefault.xml
