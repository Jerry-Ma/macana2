#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Eigen/Core>

#include <port_timestream_sensitivity.h>

namespace {

using namespace Eigen;

TEST(TimeStreamSensitivityTest, psd) {
    VectorXd scan, freqs, psd;
    scan.setOnes(101);
    timestream::psd(scan, freqs, psd, 64., true);
}

TEST(TimeStreamSensitivityTest, sensitivity) {
    VectorXd scans(100);
    scans.head(50).setConstant(1.);
    scans.tail(50).setConstant(2.);

    MatrixXI scanindex(2, 2);
    scanindex <<  0, 55,
                 40, 100;
    timestream::sensitivity(scans, scanindex, 64., 0.);
}

} // namespace
