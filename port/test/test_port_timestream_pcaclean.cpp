#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <benchmark/benchmark.h>

#include <Eigen/Core>

#include <port_timestream_pcaclean.h>
#include <logging.h>

namespace {

using namespace Eigen;

TEST(TimeStreamPcaCleanTest, clean) {
    Index ndetectors = 100;
    MatrixXd scans(2000, ndetectors);
    scans.topRows(1000).setConstant(100.);
    scans.bottomRows(1000).setConstant(200.);
    MatrixXd kernelscans(scans);
    MatrixXb scanflags(2000, ndetectors);
    scanflags.setConstant(true);

    MatrixXI scanindex(2, 2);
    scanindex <<     0, 1100,
                  1100, 2000;

    MatrixXd cleanedscans, cleanedkernelscans;

    timestream::pcaclean<timestream::EigenBackend>(scans, kernelscans, scanflags, scanindex,
                         cleanedscans, cleanedkernelscans,
                         10, -1.);
    MatrixXd cleanedscans_alt, cleanedkernelscans_alt;
    timestream::pcaclean<timestream::SpectraBackend>(scans, kernelscans, scanflags, scanindex,
                     cleanedscans_alt, cleanedkernelscans_alt,
                     10, -1.);
    cleanedscans -= cleanedscans_alt;
    cleanedkernelscans -= cleanedkernelscans_alt;
    EXPECT_DOUBLE_EQ(cleanedscans.mean(), 0.);
}

template <timestream::EigenSolverBackend backend>
void BM_pcaclean(benchmark::State& state) {
    Index ndetectors = 4000;
    MatrixXd scans(2000, ndetectors);
    scans.topRows(1000).setConstant(1.);
    scans.bottomRows(1000).setConstant(2.);
    MatrixXd kernelscans(scans);
    MatrixXb scanflags(2000, ndetectors);
    scanflags.setConstant(true);

    MatrixXI scanindex(2, 2);
    scanindex <<     0, 1100,
                  1100, 2000;

    MatrixXd cleanedscans, cleanedkernelscans;

    for (auto _ : state) {
         timestream::pcaclean<backend>(scans, kernelscans, scanflags, scanindex,
                             cleanedscans, cleanedkernelscans,
                             10, -1.);
    }
}
BENCHMARK_TEMPLATE(BM_pcaclean, timestream::EigenBackend)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_pcaclean, timestream::SpectraBackend)->Unit(benchmark::kMillisecond);

} // namespace
