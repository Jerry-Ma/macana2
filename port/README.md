# Porting Algorithms to Citlali

This is to address issue #5 (https://github.com/toltec-astro/macana2/issues/5)

This document is to walk through the steps to package the code in a way that
minimum changes will be required during the actual porting phase.

## Namespaces

To minimize name clash we use the following namespace scheme:

    ::citlali   // the root namespace. All public api should be wrapped within

        /// namespaces of high level algorithms, grouped by the data stage
        generic  // applicable to non-specific time of data, e.g, model fitting
        timestream  // timestream algorithms
        map  // map space algorithms

        /// namespaces for utility functions and helpers
        logging  // logging related
        math  // low level math and algorithms
        phys  // physics and astronomy related
        telescope  // telescope related
        ...  // could have many more of these
        misc  // random stuff that does not fit to any sensible category

        /// special namespace
        // any functions and helpers that are not the public API should
        // be wrapped in
        internal  // This could be nested within any of the namespace above

Examples:

    ::citlali::generic::leastsq()
    ::citlali::timestream::lowpass()
    ::citlali::timestream::internal::psd()  // non-public, used by sensitivity()
    ::citlali::timestream::sensitivity()
    ::citlali::map::gaussianfit()
    ::citlali::map::wienerfilter()
    ::citlali::logging::pprint()  // pretty print large data array
    ::citlali::math::sigma_clipped_stats()
    ::citlali::math::Interval  // simple tuple-like class for a value interval
    ::citlali::misc::ei_nullptr()  // generic pointer to use for optional out-argument


## Directories and Files

    macana2/
        port/
            timestream_despike.cpp
            timestream_lowpass.cpp
            generic_modelfit.cpp
            <stage>_<func>.cpp
            ...
            logging.cpp
            eigen_utils.cpp
            math_utils.cpp
            phys_utils.cpp

            test/
                test_port.cpp
                test_timestream_despike.cpp
                test_timestream_lowpass.cpp
                test_generic_modelfit.cpp
                ...
                test_<stage>_<func>.cpp
            data/
                testdata_{name}.nc
                ...
            include/
                timestream_lowpass.h
                timestream_despike.h
                ...
                logging.h
                eigen_utils.h
                math_utils.h
                phys_utils.h

### `<stage>_<func>.{h,cpp}`

These files implement functionality `func` for data stage `stage`.

Types/functions that are supposed to be used/called from other modules
(or to be short, public API) should be wrapped in namespace `<stage>`,
while the implementation details should be guarded with namespace `internal`,
e.g.

    // timestream_lowpass.h
    #pragma once

    namespace citlali {
    namespace timestream {
    namespace internal {

    Out foo(In in);   // function called by lowpass()

    } // namespace internal

    void lowpass(In in, Out out, LowpassParams lp);   // public API

    struct LowpassParams {   // public API
        enum {
           // some useful switch
        }
        // some other parameter
    };

    } // namespace timestream
    } // namespace citlali


Arguments passed to functions usually fall into several categories

* Input data
* Output data
* Compile-time switch/parameters
* Run-time switch/parameters

We use `Eigen::Matrix` as input/output data container. To allow passing
general Eigen types (expressions, views, etc.), one has to use
template, see https://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html

Compile-time switches/parameters can be used to eliminate branching during run-time.

Run-time switches/parameters can be optionally grouped into struct
`<Func>Params`. Although we prefer use a flat list of arguments rather than a
struct to make the functionality more explicit and easy to read.


Example:

    namespace timestream {

    enum Impl {
        impl1 = 0,
        impl2 = 1
    };

    template <Impl impl=impl1, typename Derived1, typename Derived2>
    void despike(
        const DenseBase<Derived1>& in
        DenseBase<Derived2>& out
        DespikeParams dp) {
        ...  // implement
        if constexpr (impl == impl1) {
        ...  // compile-time branching, C++17 needed
        } else if constexpr (impl == impl2) {
        ...
        }
        ...
    }
    }  // namespace timestream


About the function body, it should be the actual data handling part, and it should not
contain code that create objects to be used elsewhere outside of the scope of
the function. Instead, these objects need to be created outside and passed as
input and output data.

Temporary data are OK to be created with the scope of the function while being
aware of the extra memory pressure.


## test/test_port_*.cpp

The `test/test_*.cpp` is the place where tests reside.

### Function level unit tests

Each `<stage>_<func>.{h,cpp}` should comes with a
`test_<stage>_<func>.cpp`, in which a test fixture is defined. A test
fixture is a class that provide the run-time of a test case, and the test body
have full access to the members and methods to the fixture.

Any code that is responsible to setup the input/output data, and is to be
used by multiple different implementations of the same algorithm should be
put into the fixture class.

For each function, there should be at least one unit test covering the function
to make sure it works as expected.

To setup the inputs, it is recommended for function level unit tests to use
generated/simulated/theoretical data rather than external data. Although for
some cases, the latter could still be used as a last resort or by necessity.

In the test body, ideally there should only be the function call and assertion
clause, although in some cases, some massage of the data is needed prior to the
function call.

Google test also support parameterized tests, which basically
allows crating multiple tests with different parameters the runs the same code.

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


### Stage level test

More integrated test also goes in the `port/test` folder, e.g.

    test_timestream.cpp

The content of these tests can be highly flexible, as it should mimic how
the data is handled in the pipeline. It is up to the programmers to implement
in whichever way they prefer. It may or may not depends on the google test framework.


# Container classes and other new stuff

The `port/` folder shall be independent to any of the old `macana` classes,
except that the algorithm body codes are copied over from the original macana2
classes. Therefore we should focus on exploring and utilizing new data
container classes as well as new libraries that will be used or tested in
Citlali, such as `Eigen3`.

## External Dependencies

The external dependencies shall be put into `port/CMakeLists.txt`, where rules
to generate the `test_port` executable will be created.

TODO: instructions to add an external dependency in `cmake`.

# Git-flow

We use git-flow to manage the whole development process.

To port an algorithm:

1. (optional) Open an issue "Port some func in some stage ...", or some other sensible name.
   This step can be omitted, if you prepare
   to contribute code via pull request. In that case any discussion will be done in the body
   of the pull request and the pull request is essentially the issue. This is to reduce cluttering
2. Create a feature branch "feature/x-port-xxx-xxx". x should be the issue number in (1) if done,
   or 5, which is the number of "master" issue
3. Open a pull request of the feature branch to the develop branch, include
   "resolves #x" in the pull request body to signal a reference to the issue
   related.
4. Write and commit the code to the feature branch
5. Notify the code reviewer to review upon finish
6. Work back-and-forth with reviewer to address any problems
7. Merge the pull request to the develop branch

# Port back to macana

At some point, we may want to port the code back to `macana`. It could be
easily done by swapping out (rename) the original methods in the `macana`
classes, creating a new method with the original name, and implement using the
new functions in the `port` folder. We will do this as a practice before we do
the `citlali`.
