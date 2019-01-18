#pragma once

#include <logging.h>
#include <eigen_utils.h>
#include <Eigen/Dense>
#include <Spectra/SymEigsSolver.h>

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

enum EigenSolverBackend {
        EigenBackend = 0,
        SpectraBackend = 1
    };

template <EigenSolverBackend backend, typename DerivedA, typename DerivedB, typename DerivedC>
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

    SPDLOG_TRACE("input data {}", logging::pprint(scans));
    // loop over all scans
    for (Index k = 0; k < nscans; ++k) {
        // Vector scans.segment(scanindex(0, i), scanlengths(i))
        Index npts = scanlengths(k);
        SPDLOG_TRACE("process scan {} out of {}: length={}", k + 1, nscans, npts);
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
        if constexpr (backend == EigenBackend) {
            // these are sorted in increasing order, which is different from IDL
            Eigen::SelfAdjointEigenSolver<MatrixXd> solution(pcaCorr);
            const auto& evals = solution.eigenvalues();
            const auto& evecs_ = solution.eigenvectors();
            Eigen::MatrixXd evecs(evecs_);
            // find out the cutoff index for evals and do the cut
            Index cut = internal::cutoffindex(evals, neigToCut, cutStd);
            SPDLOG_TRACE("cut {} largest modes", cut);
            SPDLOG_TRACE("evals: {}", evals.tail(cut));
            evecs.rightCols(cut).setConstant(0);
            evecs.rightCols(cut).setConstant(0);
            efdet.resize(ndetectors, npts);
            efker.resize(ndetectors, npts);
            efdet.noalias() = evecs.adjoint() * det;
            efker.noalias() = evecs.adjoint() * ker;
            // efdet.bottomRows(cut).setConstant(0.);
            // efker.bottomRows(cut).setConstant(0.);
            // create cleaned data
            det.noalias() = evecs * efdet;
            ker.noalias() = evecs * efker;
       } else if constexpr (backend == SpectraBackend) {
            // Construct matrix operation object using the wrapper class DenseSymMatProd
            // int nev = ndetectors / 3;
            int nev = ndetectors <= 100?ndetectors - 1:100;
            int ncv = nev * 2.5 < ndetectors?int(nev * 2.5):ndetectors;
            SPDLOG_TRACE("spectra eigen solver nev={} ncv={}", nev, ncv);
            Spectra::DenseSymMatProd<double> op(pcaCorr);
            Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, Spectra::DenseSymMatProd<double>> eigs(&op, nev, ncv);
            // Initialize and compute
            eigs.init();
            int nconv = eigs.compute();  // the results will be sorted largest first
            // Retrieve results
            Eigen::VectorXd evals = Eigen::VectorXd::Zero(ndetectors);
            Eigen::MatrixXd evecs = Eigen::MatrixXd::Zero(ndetectors, ndetectors);
            if(eigs.info() == Spectra::SUCCESSFUL) {
                evals.head(nev) = eigs.eigenvalues();
                evecs.leftCols(nev) = eigs.eigenvectors();
            } else {
                throw std::runtime_error("failed to compute eigen values");
            }
            efdet.resize(ndetectors, npts);
            efker.resize(ndetectors, npts);
            efdet.noalias() = evecs.adjoint() * det;
            efker.noalias() = evecs.adjoint() * ker;
            // find out the cutoff index for evals and do the cut
            Index cut = internal::cutoffindex(evals, neigToCut, cutStd);
            if (cut > nev) {
                throw std::runtime_error("too few eigen values computed");
            }
            SPDLOG_TRACE("cut {} largest modes", cut);
            SPDLOG_TRACE("evals: {}", evals.head(cut));
            // here since we are computing the larget vectors, we first
            // construct the data from the larget modes, and then subtract
            efdet.bottomRows(ndetectors - cut).setConstant(0.);
            efker.bottomRows(ndetectors - cut).setConstant(0.);
            // create data to be cleaned and substract
            det.noalias() -= evecs * efdet;
            ker.noalias() -= evecs * efker;
        } else {
            static_assert(backend == EigenBackend, "UNKNOWN EIGEN SOLVER BACKEND");
        }
        // update cleaned data
        for (Index i = 0; i < ndetectors; ++i) {
            // SPDLOG_LOGGER_TRACE(logger, "write cleaned data {}", logging::pprint(det.row(i)));
            cleanedscans.col(i).segment(scanindex(0, k), npts) = det.row(i);
            cleanedkernelscans.col(i).segment(scanindex(0, k), npts) = ker.row(i);
        }
    }
}

} // namespace timestream
