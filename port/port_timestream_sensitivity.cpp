#include <cmath>
#include <tuple>

#include <logging.h>
#include <mlinterp.hpp>
#include <port_timestream_sensitivity.h>
#include <unsupported/Eigen/FFT>

namespace timestream {

using namespace Eigen;

namespace details {
//         npts   nfreqs, dfreq
std::tuple<Index, Index, double> stat(Index scanlength, double samplerate,
                                      VectorXd *freqs = nullptr) {
    auto logger = logging::createLogger("timestream._stat", nullptr);
    logger->set_level(spdlog::level::trace);

    // make an even number of data points by rounding-down
    Index npts = scanlength;
    if (npts % 2 == 1)
        npts--;
    // prepare containers in frequency domain
    Index nfreqs = npts / 2 + 1; // number of one sided freq bins
    double dfreq = samplerate / npts;
    SPDLOG_LOGGER_TRACE(logger, "using npts={} (out of n={})", npts,
                        scanlength);
    SPDLOG_LOGGER_TRACE(logger, "samplerate={} dfreq={} nfreqs={}", samplerate,
                        dfreq, nfreqs);
    if (freqs)
        freqs->array() = dfreq * ArrayXd::LinSpaced(nfreqs, 0, npts / 2);
    SPDLOG_LOGGER_TRACE(logger, "freqs{}", logging::pprint(freqs));
    return {npts, nfreqs, dfreq};
};

} // namespace details

void psd(const VectorXd &scan, VectorXd &freqs, VectorXd &psd,
         double samplerate, bool hann) {
    auto logger = logging::createLogger("timestream.psd", nullptr);
    logger->set_level(spdlog::level::trace);

    auto [npts, nfreqs, dfreq] = details::stat(scan.size(), samplerate, &freqs);

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
    FFT<double> fft;
    fft.SetFlag(FFT<double>::HalfSpectrum);
    // fft.SetFlag(FFT<double>::Unscaled);
    VectorXcd freqdata;
    fft.fwd(freqdata, timedata);
    SPDLOG_LOGGER_TRACE(logger, "fft.fwd freqdata{}",
                        logging::pprint(&freqdata));

    // calcualte psd
    psd = freqdata.cwiseAbs2() / dfreq /
          npts; // ! this is differnece from idl code
    // accound for the negative frequencies by an extra factor of 2. note: first
    // and last are 0 and nquist freq, so they only appear once
    psd.segment(1, nfreqs - 2) *= 2.;
    SPDLOG_LOGGER_TRACE(logger, "psd{}", logging::pprint(&psd));
}

// this one will be implemented to incorporate the convertion factors FCF, FEC,
// etc.
void calibrate([[maybe_unused]] const VectorXd &signal,
               [[maybe_unused]] double gain) {}

// compute sensitivities for the given scans of same bolometer with optional
// gain here we assuming the scans are packed in 1d vector and scanindex (2,
// nscan) stores the start and end of each scan
double sensitivity(const VectorXd &scans, const MatrixXI &scanindex,
                   double samplerate, double gain,
                   std::pair<double, double> freqrange) {
    auto logger = logging::createLogger("timestream.sensitivity", nullptr);
    logger->set_level(spdlog::level::trace);

    Index nscans = scanindex.cols();
    VectorXI scanlengths = scanindex.row(1) - scanindex.row(0);
    // use the median length for computation
    Index len = details::median(scanlengths);
    SPDLOG_LOGGER_TRACE(
        logger, "use median={} of scan lengths (min={} max={} nscans={})", len,
        scanlengths.minCoeff(), scanlengths.maxCoeff(), nscans);

    // get the freq stat and freq array
    VectorXd freqs;
    auto [npts, nfreqs, dfreq] = details::stat(len, samplerate, &freqs);

    // compute psd for each scan and interpolate onto the freq array
    MatrixXd psds(nfreqs, nscans);
    VectorXd tfreqs, tpsd;
    VectorXI td(1);
    for (Index i = 0; i < nscans; ++i) {
        SPDLOG_LOGGER_TRACE(logger, "process scan {} out of {}", i + 1, nscans);
        psd(scans.segment(scanindex.coeffRef(0, i), scanlengths[i]), tfreqs,
            tpsd, samplerate, true);
        // interpolate
        td << tfreqs.size();
        mlinterp::interp(td.data(), nfreqs, tpsd.data(),
                         psds.data() + i * nfreqs, tfreqs.data(), freqs.data());
    }
    SPDLOG_LOGGER_TRACE(logger, "calculated psds{}", logging::pprint(&psds));
    // convert from V^2/Hz to V/sqrt(Hz)
    psds = psds.cwiseSqrt() / std::sqrt(2.);
    // calibrate
    calibrate(psds, gain);

    // compute sensitivity with given freqrange
    // make use the fact that freqs = i * df to find index i
    auto i1 = static_cast<Index>(freqrange.first / dfreq);
    auto i2 = static_cast<Index>(freqrange.second / dfreq);
    auto nf = i2 + 1 - i1;
    VectorXd sens = psds.block(i1, 0, nf, nscans).colwise().sum() / nf;
    SPDLOG_LOGGER_TRACE(logger, "sensitivity{}", logging::pprint(&sens));
    // take a mean on the sens for representitive sens to return
    auto s = sens.mean();
    SPDLOG_LOGGER_TRACE(logger, "mean sens={}", s);

    // integrate over all frequencies
    VectorXd intsens =
        (psds.array().square() * 2. * dfreq).colwise().sum().sqrt();
    SPDLOG_LOGGER_TRACE(logger, "integrated{}", logging::pprint(&intsens));
    auto i = intsens.mean();
    SPDLOG_LOGGER_TRACE(logger, "mean intsens={}", i);

    return s;
}

} // namespace timestream