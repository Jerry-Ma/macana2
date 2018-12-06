
#include <logging.h>
#include <mlinterp.hpp>
#include <port_timestream_sensitivity.h>

namespace timestream {

namespace internal {

std::tuple<Index, Index, double> stat(Index scanlength, double samplerate) {
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
    return {npts, nfreqs, dfreq};
}

VectorXd freq(Index npts, Index nfreqs, double dfreq) {
    return dfreq * VectorXd::LinSpaced(nfreqs, 0, npts / 2);
}

} // namespace internal

} // namespace timestream