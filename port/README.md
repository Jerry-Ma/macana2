# Porting Algorithms to Citlali

This is to address issue #5 (https://github.com/toltec-astro/macana2/issues/5)

This document is to walk through the steps to package the code in a way that
minimum changes will be required during the actual porting phase.

## Directories and Files

    macana2/
        port/
            port_timestream_despike.cpp
            port_timestream_lowpass.cpp
            port_map_fit_gaussian.cpp
            port_<stage>_<func>.cpp
            ...
            test/
                test_port.cpp  // test main
                test_port_timestream_despike.cpp
                test_port_timestream_lowpass.cpp
                test_port_map_fit_gaussian.cpp
                test_port_<stage>_<func>.cpp
                ...
            data/
                test_<sensible name>.nc
                ...
            include/
                params.h
                port_timestream_despike.h
                ...

## port_*.cpp

The `port_*.cpp` file should contain no class but a (stateless) function that
implements the algorithm to be ported, enclosed by a named namespace matches
the data stage (timestream, map, etc.) to avoid name clash:

    namespace timestream {

        void despike(const InDataType& in, OutDataType& out, [ArgType arg [ArgType arg2] ...]) {
        }

    }

The first two arguments should be the input and output data. The input data
shall be passed as const reference and the output data shall be passed as reference
and will be populated/modified at the end of the function call.

The rest of the arguments are all other tunable parameters. We prefer use a flat
list of arguments rather than a struct to pass the parameters to make the functionality
more explicit and easy to read.

Note that for some of the functions it may be depends on tens of parameters, which
may be too verbose. In that case, a plain struct defined in `port/include/params.h`
that groups related parameters is desired:

    // port/include/params.h
    struct TelecopeParams
    {
        double arg = 0.;
        in arg2 = 1;
        ...
    }

# test/test_port_*.cpp

The `test/test_port_*.cpp` is the place where tests reside.

## Function level unit tests

Each `port_<stage>_<func>.cpp` should comes with a
`test_port_<stage>_<func>.cpp`, in which a test fixture is defined. A test
fixture is a class that provide the run-time of a test case, and the test body
have full access to the members and methods to the fixture:

    #include <gtest/gtest.h>
    #include "AnalParams.h"

    namespace {

    class AnalParamsTest : public ::testing::Test
    {
    protected:
        AnalParamsTest(): ap(new AnalParams(apXmlFile, 1)) {}
        ~AnalParamsTest() override {delete ap;}
        void SetUp() override {}
        void TearDown() override {}
        std::string apXmlFile = "../data/test_apb.xml";
        AnalParams* ap = nullptr;
    };  // this is the test fixture class
        // note that here we are refering to the data folder for
        // some static test data files
 
    TEST_F(AnalParamsTest, AnalParamsFromApXml) {
        // test body of test named "AnalParamsFramApXml" that
        // test the read-in of ap.xml file
        EXPECT_EQ(ap->beammapping, 1);  // test assertion. there could be multiple of these.
    }

    TEST_F(AnalParamsTest, AnalParamsGetNFiles) {
        // yet another test
        EXPECT_EQ(ap->getNFiles(), 1);
    }

The google test framework has great documentation and please read through
and learn the capability of it, and use. some links:

* https://github.com/google/googletest/tree/master/googletest/docs
* https://github.com/google/googletest/blob/master/googletest/docs/primer.md
* https://github.com/google/googletest/blob/master/googletest/docs/advanced.md#how-to-write-value-parameterized-tests
* https://github.com/google/googletest/blob/master/googletest/docs/samples.md


## State level test

More integrated test also goes in the `port/test` folder, e.g.

    test_timestream.cpp

The content of these tests can be highly flexible, as it should mimic how
the data is handled in the pipeline. It is up to the programmers to implement
in whichever way they prefer. It may or may not depends on the google test framework.


# Container classes and other stuff

The `port/` folder shall be indepedant to anyof the old macana2 classes, except
that the algorithm body code are copied over from the original macana2 classes.
Therefore we should focus on exploring and utilizing new data container classes
that will be used in Citlali, such as `Eigen3`.
