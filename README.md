# Macana2


## Get `Macana2`

    git clone https://github.com/toltec-astro/macana2.git

This will create a folder of name `macana2` in your current directory.


## Build `Macana2`

`Macana2` consists of with two sets of programs:

* Production tools

    * include `macanap` and `beammap`
    * for actuall data processing
    * run the same way as the original `macana` (a. k. a., `aztec_c++`).

* Testing tools

    * include `beammap_gui` and `macana_test`
    * for development and diagnostic purpose

While both of these depend on some common external libraries to compile and
run, the testing tools require additional software packages and settings.

Build guides:

* [macOS](#macos)

### macOS

#### Install common dependencies

`Homebrew` is recommended to install the common dependencies.

* Install/Update Homebrew

    see <https://brew.sh>.

* Install the dependencies

        brew tap brewsci/science
        brew install cmake libomp suitesparse gsl fftw ccfits netcdf

    Note: pay attention to the output of the `brew install` command, as there
    may be further instructions needed to complete the installation. One
    example is that you may need to run `brew link --overwrite xxx` if there is
    an old copy of `xxx` installed previously which may prevent `brew` from
    finishing its link step.

* Notes

    Some tweaks are needed because some of the Homebrew recipes does not install
    things correctly.

        # suitesparse headers are not linked by default
        # replace "x.x.x" with the actual version numbers
        ln -s /usr/local/Cellar/suite-sparse/x.x.x/include /usr/local/include/suitesparse

#### Build production tools

With the common dependencies installed, we can build the production tools
using `cmake`:

    cd macana2
    mkdir build
    cd build
    cmake ..


#### Install dependencies for testing tools

The testing GUI is built using `Qt`, and the unittests uses
`Google Test` framework. `qmake` from `Qt` is needed to generate
the `Makefile` to build both tools.

* Install `Qt` and `qmake`

    see <https://www.qt.io/download> and more on <http://doc.qt.io/qt-5/gettingstarted.html>.

    Once installed, `qmake` could be accessed from commandline with its full path:

        <Qt install path>/<Qt version>/clang_64/bin/qmake

    e.g.,

        /Applications/Qt/5.11.2/clang_64/bin/qmake


* Install `Google Test Framework`

    `Google Test` is not available in Homebrew by default, but we could still
    install it with a custom recipe by
    [Kronuz](https://gist.githubusercontent.com/Kronuz/96ac10fbd8472eb1e7566d740c4034f8/raw/gtest.rb):

        brew install --HEAD https://gist.githubusercontent.com/Kronuz/96ac10fbd8472eb1e7566d740c4034f8/raw/gtest.rb

#### Build testing tools

Once the additional dependencies are installed, we can build the testing
tools using `qmake`:

    cd macana2
    mkdir qtbuild
    cd qtbuild
    /path/to/qmake ../macana2.pro

Note: There might be issues running qmake if `Qt` is not installed in a
standard location (such as `/Applications/Qt` for macOS). It is recommended to
use `Qt Creator` to configure the building environment, as well as build and
run the project.

Caveat: In `Qt Creator`, the default "Run" settings under "Project" tab may have "Add build library search path to ..." checked.
Uncheck if the program could not find the shared/dylib libraries at runtime.

## Run `Macana2`

### Production tools

The executables including `macanap` and `beammap` reside in `build/bin/`, e.g.,

    /path/to/build_dir/bin/beammap apb.xml

### Testing tools

The `beammap_gui` executable is in `qtbuild/beammap_gui/`, and `macana_test`
is in `qtbuild/test/`.


## Documentation

    TBD
