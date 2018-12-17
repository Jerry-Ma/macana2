#ifndef PORT_TIMESTREAM_SENSITIVITY_H
#define PORT_TIMESTREAM_SENSITIVITY_H

#include <eigen_utils.h>
#include <unsupported/Eigen/FFT>

#include <tuple>

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
using eiu::Interval;

namespace internal {

//                         npts   nfreqs, dfreq
using FreqStat = std::tuple<Index, Index, double>;

FreqStat stat(Index scanlength, double samplerate);
VectorXd freq(Index npts, Index nfreqs, double dfreq);

enum Window {
    NoWindow = 0,
    Hanning = 1
};

inline auto hann(Index npts)
{
    // hann window
    // N numbers starting from 0 to (include) 2pi/N * (N-1)
    // 0.5 - 0.5 cos([0, 2pi/N,  2pi/N * 2, ... 2pi/N * (N - 1)])
    // NENBW = 1.5 * df therefore we devide by 1.5 here to get the
    // equivalent scale as if no window is used
    return (0.5 - 0.5 * ArrayXd::LinSpaced(npts, 0, 2 * M_PI / npts * (npts - 1)).cos()).matrix() / 1.5;
}

// psd of individual scan npts/2 + 1 values with frequency [0, f/2];
template <Window win=Hanning, typename DerivedA, typename DerivedB, typename DerivedC>
FreqStat psd(const Eigen::DenseBase<DerivedA> &_scan,
             Eigen::DenseBase<DerivedB> &psd, Eigen::DenseBase<DerivedC> *freqs,
             double samplerate) {
    auto logger = logging::createLogger("timestream.psd", nullptr);
    logger->set_level(spdlog::level::trace);
    // decltype(auto) scan = _scan.derived();
    // decltype(auto) forward the return type of derived() so it declares a refernce as expected
    // if scan has storage, this is equivalent to:
    typename eiu::const_ref<DerivedA> scan(_scan.derived());

    auto stat = internal::stat(scan.size(), samplerate);
    auto [npts, nfreqs, dfreq] = stat;

    // prepare fft
    Eigen::FFT<double> fft;
    fft.SetFlag(Eigen::FFT<double>::HalfSpectrum);
    fft.SetFlag(Eigen::FFT<double>::Unscaled);
    VectorXcd freqdata;

    // branch according to whether applying hann
    if constexpr (win == Hanning) {
        SPDLOG_LOGGER_TRACE(logger, "apply hann window");
        // we need the eval() here per the requirement of fft.fwd()
        fft.fwd(freqdata, scan.head(npts).cwiseProduct(
                    internal::hann(npts)).eval());
    } else if (win == NoWindow) {
        fft.fwd(freqdata, scan.head(npts));
    } // note: at this point the freqdata is not normalized to NEBW yet

    SPDLOG_LOGGER_TRACE(logger, "fft.fwd freqdata{}",
                        logging::pprint(&freqdata));

    // calcualte psd
    // normalize to frequency resolution
    psd = freqdata.cwiseAbs2() / dfreq;  // V/Hz^0.5
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
template <Window win=Hanning, typename DerivedA, typename DerivedB, typename DerivedC,
          typename DerivedD>
FreqStat psds(const Eigen::DenseBase<DerivedA> &scans,
              const Eigen::DenseBase<DerivedB> &scanindex,
              Eigen::DenseBase<DerivedC> &_psds,
              Eigen::DenseBase<DerivedD> *freqs, double samplerate
              ) {
    auto logger = logging::createLogger("timestream.psds", nullptr);
    logger->set_level(spdlog::level::trace);

    // decltype(auto) psds = _psds.derived();
    typename eiu::ref<DerivedC> psds(_psds.derived());

    // prepare common freq grid
    Index nscans = scanindex.cols();
    VectorXI scanlengths = scanindex.row(1) - scanindex.row(0);
    // use the median length for computation
    Index len = eiu::median(scanlengths);
    SPDLOG_LOGGER_TRACE(
        logger, "use median={} of scan lengths (min={} max={} nscans={})", len,
        scanlengths.minCoeff(), scanlengths.maxCoeff(), nscans);

    // get the common freq stat and freq array
    auto stat = internal::stat(len, samplerate);
    auto [npts, nfreqs, dfreq] = stat;
    VectorXd _freqs = internal::freq(npts, nfreqs, dfreq);
    SPDLOG_LOGGER_TRACE(logger, "use freqs{}", logging::pprint(&_freqs));

    // compute psd for each scan and interpolate onto the common freq grid
    psds.resize(nfreqs, nscans);

    // make some temporaries
    VectorXd tfreqs, tpsd;
    VectorXI td(1);  // need this to store array ize for interp
    // get the psds
    for (Index i = 0; i < nscans; ++i) {
        SPDLOG_LOGGER_TRACE(logger, "process scan {} out of {}", i + 1, nscans);
        internal::psd<win>(scans.segment(scanindex(0, i), scanlengths(i)), tpsd,
                      &tfreqs, samplerate);
        // interpolate (tfreqs, tpsd) on to _freqs
        td << tfreqs.size();
        SPDLOG_LOGGER_TRACE(logger, "interpolate tpsd{} from tfreqs{} to freqs{}",
                           logging::pprint(&tpsd), logging::pprint(&tfreqs), logging::pprint(&_freqs));
        SPDLOG_LOGGER_TRACE(logger, "interpolate sizes{}",
                            logging::pprint(&td));
        // interp (tfreq, tpsd) on to freq and store the result in column i of psds
        mlinterp::interp(td.data(), nfreqs, tpsd.data(),
                         psds.data() + i * nfreqs, tfreqs.data(),
                         _freqs.data());
        SPDLOG_LOGGER_TRACE(logger, "updated psds{}", logging::pprint(&psds));
    }

    SPDLOG_LOGGER_TRACE(logger, "calulated psds{}", logging::pprint(&psds));

    // update the freqs array if requested
    if (freqs) {
        freqs->operator=(_freqs);
    }
    return stat;
}

// The helper class eiu::MoveEnabledUnaryOp is used to wrap a lambda function and
// provide overloaded calling signitures to allow optionally
// "consume" the input parameter if passed in as rvalue reference
// such as temporary objects and std::move(var).
// Data held by parameters passed this way is transfered to the returning
// variable at call site, e.g.
// call "auto ret = neps(std::move(in));" will move data in "in" to
// "ret" and apply the computation inplace on ret.
// Also note the auto&& type of input parameter, this will allow it
// work for different Eigen expression types
inline auto psd2sen = eiu::MoveEnabledUnaryOp([](auto&& psd){
    return (psd / 2.).cwiseSqrt();  // V * s^(1/2)
});

/*
inline auto sen2psd = eiu::MoveEnabledUnaryOp([](auto&& sen){
    return  sen.cwiseSquare() * 2.; // V^2 / Hz^(1/2)
});
*/

} // namespace internal

// this one will be implemented to incorporate the convertion factors FCF, FEC,
// etc. to calibrat the signal
//
/*
template <typename Derived>
inline auto calibrate = eiu::MoveEnabledUnaryOp([](
        auto&& signal,
        double gain, double gain_error, double fec, double fec_error, double fcf, double fcf_error
                                                ) {
    return signal * gain / fec * fcf;
}
);
*/

template <typename DerivedA, typename DerivedB, typename DerivedC,
          typename DerivedD>
void sensitivity(
    const Eigen::DenseBase<DerivedA> &scans,
    const Eigen::DenseBase<DerivedB> &scanindex,
    Eigen::DenseBase<DerivedC> &sensitivities, // V * s^(1/2)
    Eigen::DenseBase<DerivedD> &noisefluxes,   // V, = sqrt(\int PSD df)
    double samplerate,
    // double gain,
    Interval<double> freqrange = {3., 5.}) {

    auto logger = logging::createLogger("timestream.sensitivity", nullptr);
    logger->set_level(spdlog::level::trace);

    // get psds
    MatrixXd tpsds;
    auto [npts, nfreqs, dfreq] = internal::psds<internal::Hanning>(
        scans, scanindex, tpsds, eiu::ei_nullptr(), samplerate);

    // compute noises by integrate over all frequencies
    noisefluxes = (tpsds * dfreq).colwise().sum().cwiseSqrt();
    SPDLOG_LOGGER_TRACE(logger, "nosefluxes{}", logging::pprint(&noisefluxes));
    auto meannoise = noisefluxes.mean();
    SPDLOG_LOGGER_TRACE(logger, "meannoise={}", meannoise);

    // get sensitivity in V * s^(1/2)
    // this semantic is to indicate the tpsds is to be consumed herer
    // i.e., the data held by tpsds will be moved to neps after the call
    auto sens = internal::psd2sen(std::move(tpsds));
    // to create a copy, just call the following instead
    // MatrixXd sens = internal::psd2sen(tpsds * 2.);
    // to defer the computation, call the following
    // auto sens = internal::psd2sen(tpsds);

    SPDLOG_LOGGER_TRACE(logger, "consumed psds{}", logging::pprint(&tpsds));
    SPDLOG_LOGGER_TRACE(logger, "sens{}", logging::pprint(&sens));

    // calibrate
    // neps = calibrate(neps, gain); // mJy/sqrt(Hz)

    // compute sensitivity with given freqrange
    // make use the fact that freqs = i * df to find index i
    auto i1 = static_cast<Eigen::Index>(freqrange.left() / dfreq);
    auto i2 = static_cast<Eigen::Index>(freqrange.right() / dfreq);
    auto nf = i2 + 1 - i1;
    sensitivities =
        sens.block(i1, 0, nf, scanindex.cols()).colwise().sum() / nf;
    SPDLOG_LOGGER_TRACE(logger, "sensitivities{}",
                        logging::pprint(&sensitivities));
    // take a mean on the sens for representitive sens to return
    auto meansens = sensitivities.mean();
    SPDLOG_LOGGER_TRACE(logger, "meansens={}", meansens);

}

} // namespace timestream

#endif // PORT_TIMESTREAM_SENSITIVITY_H
