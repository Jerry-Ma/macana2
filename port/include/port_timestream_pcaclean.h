#pragma once

#include <logging.h>
#include <eigen_utils.h>
#include <Eigen/Dense>

namespace timestream {

using Eigen::ArrayXd;
using Eigen::DenseBase;
using Eigen::Index;
using Eigen::MatrixXd;
using Eigen::MatrixXI;
using Eigen::VectorXcd;
using Eigen::VectorXd;
using Eigen::VectorXI;
using Eigen::VectorXb;
using Eigen::MatrixXb;

namespace internal {

template <typename Derived>
std::tuple<double, double> sigma_clipped_stats(const DenseBase<Derived>& values,
                                               double sigma_lower, double sigma_upper) {
    double mean, stddev;
    return {mean, stddev};
}

template <typename Derived>
Index cutoffindex(const DenseBase<Derived>& evals, Index neigToCut, double cutStd) {
    if (neigToCut <= 0 && cutStd < 2.) {
        throw std::runtime_error("insufficient cut in eigen modes");
    }
    if (cutStd <= 0) {
        return neigToCut;
    }
    VectorXd ev = evals.derived().array().abs().log10();
    auto [mev, std] = sigma_clipped_stats(ev, cutStd, cutStd);
    double cut = pow(10, mev + cutStd * std);
    // return which ever the smaller, cut values beyond this point
    return neigToCut < cut? neigToCut: cut;
}

} // namespace internal

template <typename DerivedA, typename DerivedB, typename DerivedC>
void pcaclean(
    const Eigen::DenseBase<DerivedA> &scans,
    const Eigen::DenseBase<DerivedA> &kernelscans,
    const Eigen::DenseBase<DerivedB> &scanflags,
    const Eigen::DenseBase<DerivedC> &scanindex,
    Eigen::DenseBase<DerivedA> &cleanedscans,
    Eigen::DenseBase<DerivedA> &cleanedkernelscans,
    Index neigToCut,
    double cutStd
    ) {

    auto logger = logging::createLogger("timestream.pcaclean", nullptr);
    // logger->set_level(spdlog::level::trace);

    // scans, kernalscans, and scanflags are [ndata, ndetectors]
    // timestream data matrix from all detectors
    // scanindex is [2, nscans] matrix that store indexes to locate a scan
    Index ndetectors = scans.cols();
    Index nscans = scanindex.cols();
    VectorXI scanlengths = scanindex.row(1) - scanindex.row(0);

    // containers of pca [npts, ndetctors]
    MatrixXd det, ker, efdet, efker;
    MatrixXb flg;
    MatrixXd denom(ndetectors, ndetectors);
    MatrixXd pcaCorr(ndetectors, ndetectors);
    // return data
    cleanedscans.derived().resize(scans.rows(), scans.cols());
    cleanedkernelscans.derived().resize(kernelscans.rows(), kernelscans.cols());

    SPDLOG_LOGGER_TRACE(logger, "input data {}", logging::pprint(scans));
    // loop over all scans
    for (Index k = 0; k < nscans; ++k) {
        // Vector scans.segment(scanindex(0, i), scanlengths(i))
        Index npts = scanlengths(k);
        SPDLOG_LOGGER_TRACE(logger, "process scan {} out of {}: length={}", k + 1, nscans, npts);
        // prepare containers
        det.resize(ndetectors, npts);
        ker.resize(ndetectors, npts);
        flg.resize(ndetectors, npts);
        // populate
        for (Index i = 0; i < ndetectors; ++i) {
            det.row(i) = scans.col(i).segment(scanindex(0, k), npts);
            ker.row(i) = kernelscans.col(i).segment(scanindex(0, k), npts);
            flg.row(i) = scanflags.col(i).segment(scanindex(0, k), npts);
            // a function should be implemented to median-subtract the values
            // of det and ker inplace
            // it should possibly only include data of the center region
            // the same as in the case of macana, although it is implemented
            // to be true all the time there
        }
        // calculate denom
        // TODO figure out what this really does
        // denom = (flg * flg.adjoint()).array() - 1.;
        // noalias to force eigen evaluate into pcaCorr
        pcaCorr.noalias() = (det * det.adjoint()); // .cwiseQuotient(denom);
        // possibly make use of use additional corrMatrix from atmTemplate
        // pcaCorr = pcaCorr.cwiseProduct(corrMatrix)
        // compute eigen values and eigen vectors
        // these are sorted in increasing order, which is different from IDL
        Eigen::SelfAdjointEigenSolver<MatrixXd> solution(pcaCorr);
        auto evals = solution.eigenvalues();
        auto evecs = solution.eigenvectors();
        efdet.resize(ndetectors, npts);
        efker.resize(ndetectors, npts);
        efdet.noalias() = evecs.adjoint() * det;
        efker.noalias() = evecs.adjoint() * ker;
        // find out the cutoff index for evals and do the cut
        Index cut = internal::cutoffindex(evals, neigToCut, cutStd);
        SPDLOG_LOGGER_TRACE(logger, "cut {} largest modes", cut);
        efdet.bottomRows(cut).setConstant(0.);
        efker.bottomRows(cut).setConstant(0.);
        /*
        if (neigToCut > 0.) {
            efdet.bottomRows(neigToCut) = 0.;
            efker.bottomRows(neigToCut) = 0.;
        } else if (cutStd < 1.) {
            throw std::runtime_error("cutStd is too small");
        }
        Index cutAfterIndex = internal::
        VectorXd ev = evals.array().abs().log10();
        auto [mev, std] = sigma_clipped_stats(ev, cutStd, cutStd);
        double cut = pow(10, mev + cutStd * std);
        Index cutIndex =
        */
        // create cleaned data
        det.noalias() = evecs * efdet;
        ker.noalias() = evecs * efker;
        // update cleaned data
        for (Index i = 0; i < ndetectors; ++i) {
            // SPDLOG_LOGGER_TRACE(logger, "write cleaned data {}", logging::pprint(det.row(i)));
            cleanedscans.col(i).segment(scanindex(0, k), npts) = det.row(i);
            cleanedkernelscans.col(i).segment(scanindex(0, k), npts) = ker.row(i);
        }
    }
}

} // namespace timestream
