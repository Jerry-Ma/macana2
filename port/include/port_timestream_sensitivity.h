#ifndef PORT_TIMESTREAM_SENSITIVITY_H
#define PORT_TIMESTREAM_SENSITIVITY_H

#include <cmath>
#include <tuple>

#include <Eigen/Core>
#include <unsupported/Eigen/FFT>

#include <eigen_utils.h>
#include <mlinterp.hpp>

#include <logging.h>

namespace timestream {

using Eigen::ArrayXd;
using Eigen::DenseBase;
using Eigen::Index;
using Eigen::MatrixXd;
using Eigen::MatrixXI;
using Eigen::VectorXcd;
using Eigen::VectorXd;
using Eigen::VectorXI;
using utils::Interval;

namespace internal {

//                         npts   nfreqs, dfreq
using FreqStat = std::tuple<Index, Index, double>;

FreqStat stat(Index scanlength, double samplerate);
VectorXd freq(Index npts, Index nfreqs, double dfreq);

// psd of individual scan npts/2 + 1 values with frequency [0, f/2];
template <typename DerivedA, typename DerivedB, typename DerivedC>
FreqStat psd(const Eigen::DenseBase<DerivedA> &scan,
             Eigen::DenseBase<DerivedB> &psd, Eigen::DenseBase<DerivedC> *freqs,
             double samplerate, bool hann = true) {
    auto logger = logging::createLogger("timestream.psd", nullptr);
    logger->set_level(spdlog::level::trace);

    auto stat = internal::stat(scan.size(), samplerate);
    auto [npts, nfreqs, dfreq] = stat;

    // create a copy of the data
    VectorXd timedata(scan.head(npts));

    if (hann) {
        SPDLOG_LOGGER_TRACE(logger, "apply hann window");
        const double PI = std::atan(1.0) * 4;
        // apply hann window
        // 0.5 - 0.5 cos([0, 2 * pi / (N - 1), ... 2 * pi * i/ (N - 1)])
        timedata.array() *=
            0.5 - 0.5 * ArrayXd::LinSpaced(npts, 0, 2. * PI).cos();
    }

    // do fft
    Eigen::FFT<double> fft;
    fft.SetFlag(Eigen::FFT<double>::HalfSpectrum);
    // fft.SetFlag(FFT<double>::Unscaled);
    VectorXcd freqdata;
    fft.fwd(freqdata, timedata);
    SPDLOG_LOGGER_TRACE(logger, "fft.fwd freqdata{}",
                        logging::pprint(&freqdata));

    // calcualte psd
    psd = freqdata.cwiseAbs2() /
          (samplerate * npts); // ! this is differnece from idl code
    // accound for the negative frequencies by an extra factor of 2. note: first
    // and last are 0 and nquist freq, so they only appear once
    psd.segment(1, nfreqs - 2) *= 2.;
    SPDLOG_LOGGER_TRACE(logger, "psd{}", logging::pprint(&psd));

    // make the freqency array when requested
    if (freqs) {
        freqs->operator=(internal::freq(npts, nfreqs, dfreq));
        SPDLOG_LOGGER_TRACE(logger, "freqs{}", logging::pprint(freqs));
    }
    return stat;
}

// a set of psd for multiple scans with potentialy different length.
// all individual psds are interpolated to a common frequency array
template <typename DerivedA, typename DerivedB, typename DerivedC,
          typename DerivedD>
FreqStat psds(const Eigen::DenseBase<DerivedA> &scans,
              const Eigen::DenseBase<DerivedB> &scanindex,
              Eigen::DenseBase<DerivedC> &_psds,
              Eigen::DenseBase<DerivedD> *freqs, double samplerate,
              bool hann = true) {
    auto logger = logging::createLogger("timestream.psds", nullptr);
    logger->set_level(spdlog::level::trace);

    typename Eigen::internal::ref_selector<DerivedC>::non_const_type psds(
        _psds.derived());

    Index nscans = scanindex.cols();
    VectorXI scanlengths = scanindex.row(1) - scanindex.row(0);
    // use the median length for computation
    Index len = utils::median(scanlengths);
    SPDLOG_LOGGER_TRACE(
        logger, "use median={} of scan lengths (min={} max={} nscans={})", len,
        scanlengths.minCoeff(), scanlengths.maxCoeff(), nscans);

    // get the common freq stat and freq array
    auto stat = internal::stat(len, samplerate);
    auto [npts, nfreqs, dfreq] = stat;
    VectorXd _freqs = internal::freq(npts, nfreqs, dfreq);
    SPDLOG_LOGGER_TRACE(logger, "use freqs{}", logging::pprint(&_freqs));

    // compute psd for each scan and interpolate onto the freq array
    psds.resize(nfreqs, nscans);

    // some temporaries
    VectorXd tfreqs, tpsd;
    VectorXI td(1);
    for (Index i = 0; i < nscans; ++i) {
        SPDLOG_LOGGER_TRACE(logger, "process scan {} out of {}", i + 1, nscans);
        internal::psd(scans.segment(scanindex(0, i), scanlengths(i)), tpsd,
                      &tfreqs, samplerate, hann);
        // interpolate (tfreqs, tpsd) on to _freqs
        td(0) = tfreqs.size();
        SPDLOG_LOGGER_TRACE(logger, "interpolate tfreqs{} to freqs{}",
                            logging::pprint(&tfreqs), logging::pprint(&_freqs));
        SPDLOG_LOGGER_TRACE(logger, "interpolate sizes{}",
                            logging::pprint(&td));
        SPDLOG_LOGGER_TRACE(logger, "interpolate ydata{}",
                            logging::pprint(&tpsd));
        mlinterp::interp(td.data(), nfreqs, tpsd.data(),
                         psds.data() + i * nfreqs, tfreqs.data(),
                         _freqs.data());
        SPDLOG_LOGGER_TRACE(logger, "interpolated y{}", logging::pprint(&psds));
    }

    SPDLOG_LOGGER_TRACE(logger, "calulated psds{}", logging::pprint(&psds));

    // update the freqs array if requested
    if (freqs) {
        freqs->operator=(_freqs);
    }
    return stat;
}

// convert psd to noise equivalent power
//  this will consume psd
template <typename Derived> Derived neps(Eigen::DenseBase<Derived> &&psds) {
    // convert from V^2/Hz to V*sqrt(s)
    return (psds.derived() / 2.).cwiseSqrt();
}

} // namespace internal

// this one will be implemented to incorporate the convertion factors FCF, FEC,
// etc. to calibrat the signal
//
// this will modify signal inplace
template <typename Derived>
void calibrate(Eigen::DenseBase<Derived> &signal, double gain) {}

template <typename DerivedA, typename DerivedB, typename DerivedC,
          typename DerivedD>
void sensitivity(
    const Eigen::DenseBase<DerivedA> &scans,
    const Eigen::DenseBase<DerivedB> &scanindex,
    Eigen::DenseBase<DerivedC> &sensitivities, // NEFD
    Eigen::DenseBase<DerivedD> &noisefluxes,   // sqrt(integ (NeFD^2 df))
    double samplerate, double gain, Interval<double> freqrange = {3., 5.}) {
    auto logger = logging::createLogger("timestream.sensitivity", nullptr);
    logger->set_level(spdlog::level::trace);

    // get psds
    MatrixXd tpsds;
    auto [npts, nfreqs, dfreq] = internal::psds(
        scans, scanindex, tpsds, utils::ei_nullptr(), samplerate, true);
    // get nep from psds
    MatrixXd neps = internal::neps(std::move(tpsds));
    SPDLOG_LOGGER_TRACE(logger, "neps{}", logging::pprint(&neps));

    // calibrate
    calibrate(neps, gain); // mJy/sqrt(Hz)

    // compute sensitivity with given freqrange
    // make use the fact that freqs = i * df to find index i
    auto i1 = static_cast<Eigen::Index>(freqrange.left() / dfreq);
    auto i2 = static_cast<Eigen::Index>(freqrange.right() / dfreq);
    auto nf = i2 + 1 - i1;
    sensitivities =
        neps.block(i1, 0, nf, scanindex.cols()).colwise().sum() / nf;
    SPDLOG_LOGGER_TRACE(logger, "sensitivities{}",
                        logging::pprint(&sensitivities));
    // take a mean on the sens for representitive sens to return
    auto meansens = sensitivities.mean();
    SPDLOG_LOGGER_TRACE(logger, "meansens={}", meansens);

    // compute noises by integrate over all frequencies
    noisefluxes = (neps.array().square() * 2. * dfreq).colwise().sum().sqrt();
    SPDLOG_LOGGER_TRACE(logger, "nosefluxes{}", logging::pprint(&noisefluxes));
    auto meannoise = noisefluxes.mean();
    SPDLOG_LOGGER_TRACE(logger, "meannoise={}", meannoise);
}

} // namespace timestream

#endif // PORT_TIMESTREAM_SENSITIVITY_H
