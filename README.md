# Macana2


## Get `Macana2`

    git clone https://github.com/toltec-astro/macana2.git

This will create a folder of name `macana2` in your current directory.


## Build `Macana2`

`Macana2` consists of multiple sub-programs:

* Production programs

    * include `macanap` and `beammap`
    * for actuall data processing
    * run the same way as the original `macana` (a. k. a., `aztec_c++`).

* Development programs

    * include `beammap_gui`, `macana2_test`, and `port_test`
    * for development and diagnostic purpose

While all of these depend on some common external libraries to compile and
run, the development programs require additional software packages and settings.

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

#### Build production programs

With the common dependencies installed, we can build the production programs
with `cmake`:

    cd macana2
    mkdir build
    cd build
    cmake ..
    make

The binaries will be in `build/bin/`.

#### Install dependencies for development programs

##### `beammap_gui`

The `beammap_gui` is built using `Qt`. `qmake` from `Qt` is needed to generate
the `Makefile` to build the executable.

* Install `Qt` and `qmake`

    see <https://www.qt.io/download> and more on <http://doc.qt.io/qt-5/gettingstarted.html>.

    Once installed, `qmake` could be accessed from commandline with its full path:

        <Qt install path>/<Qt version>/clang_64/bin/qmake

    e.g.,

        /Applications/Qt/5.11.2/clang_64/bin/qmake


##### `macana2_test` and `port_test`

Both of these include unit tests that are written using the Google Test
Framework. To compile these programs, one can use either `qmake` or `cmake`

To use `qmake`, one has to install Google Test Framework separately.
Google Test, however, is not available in Homebrew by default, but we could still
install it with a custom recipe by
[Kronuz](https://gist.githubusercontent.com/Kronuz/96ac10fbd8472eb1e7566d740c4034f8/raw/gtest.rb):

    brew install --HEAD https://gist.githubusercontent.com/Kronuz/96ac10fbd8472eb1e7566d740c4034f8/raw/gtest.rb

The `cmake` does not use Google Test as an external dependency, instead it will
pull the source code during the build time. Therefore, we don't need to install
Google Test separately.

#### Build development programs

Once the additional dependencies are installed, we can build the development programs using `qmake`:

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


Or using cmake to build the unit tests:

    cd macana2/build
    cmake ..
    make macana2_test
    make port_test

## Run `Macana2`

### Production programs

The executables including `macanap` and `beammap` reside in `build/bin/`, e.g.,

    /path/to/build_dir/bin/beammap apb.xml

### Development programs

The `beammap_gui` executable is in `qtbuild/beammap_gui/`.

The unit tests `macana2_test` is in `qtbuild/test/` or `build/test`, depending on
which make-tool is used.

The unit tests `port_test` is in `qtbuild/port` or `build/port`, depending on
which make-tool is used.


## Documentation

    TBD
