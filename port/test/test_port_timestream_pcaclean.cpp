#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Eigen/Core>

#include <port_timestream_pcaclean.h>

namespace {

using namespace Eigen;

TEST(TimeStreamPcaCleanTest, clean) {
    MatrixXd scans(2000, 7000);
    scans.topRows(1000).setConstant(1.);
    scans.bottomRows(1000).setConstant(2.);

    MatrixXd kernelscans(scans);

    MatrixXb scanflags(2000, 7000);
    scanflags.setConstant(true);

    MatrixXI scanindex(2, 2);
    scanindex <<     0, 1100,
                  1100, 2000;

    MatrixXd cleanedscans, cleanedkernelscans;

    timestream::pcaclean(scans, kernelscans, scanflags, scanindex,
                         cleanedscans, cleanedkernelscans,
                         10, -1.);
}

} // namespace
