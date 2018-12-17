#pragma once

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#define SPDLOG_DISABLE_DEFAULT_LOGGER
#define SPDLOG_FMT_EXTERNAL
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
// remove the extra filename:lineno in trace
#undef SPDLOG_LOGGER_TRACE
#define SPDLOG_LOGGER_TRACE(logger, ...)                                       \
    logger->log(spdlog::source_loc{}, spdlog::level::trace, __VA_ARGS__)
#include <fmt/format.h>
#include <fmt/ostream.h>

#include <Eigen/Core>

namespace logging {

// shared sink
inline const static auto console =
    std::make_shared<spdlog::sinks::stdout_color_sink_mt>();

inline std::shared_ptr<spdlog::logger> createLogger(std::string_view name,
                                                    const void *const ptr) {
    auto logger = std::make_shared<spdlog::logger>(
        fmt::format("{}@{:x}", name, reinterpret_cast<std::uintptr_t>(ptr)),
        logging::console);
    logger->set_level(spdlog::level::trace);
    return logger;
}

template <typename M> struct pprint {
    pprint(M *m) : m(m) {}
    template <typename OStream>
    friend OStream &operator<<(OStream &os, const pprint &pp) {
        const M &m = *(pp.m);
        if (m.size() == 0)
            return os << "(empty)";
        os << "(" << m.rows() << "," << m.cols() << ")";

        auto f = pp.matrix_formatter();
        print_matrix(os, m, f, 5, 5, 10);
        return os;
    }

  private:
    const M *m;
    Eigen::IOFormat matrix_formatter() const {
        if (m->cols() == 1)
            return Eigen::IOFormat(Eigen::StreamPrecision, Eigen::DontAlignCols,
                                   ", ", ", ", "", "", "[", "]");
        if (m->cols() < 3)
            return Eigen::IOFormat(Eigen::StreamPrecision, Eigen::DontAlignCols,
                                   ", ", " ", "[", "]", "[", "]");
        return Eigen::IOFormat(Eigen::StreamPrecision, 0, ", ", "\n", "[", "]",
                               "[\n", "]\n");
    }
    template <typename OStream, typename Derived>
    static OStream &print_matrix(OStream &s,
                                 const Eigen::DenseBase<Derived> &_m,
                                 const Eigen::IOFormat &fmt, int max_rows,
                                 int max_cols, int max_size = -1) {
        if (_m.size() == 0) {
            s << fmt.matPrefix << fmt.matSuffix;
            return s;
        }

        typename Derived::Nested m(_m.derived());

        if (max_size < 0)
            max_size = max_rows * max_cols;
        if (m.cols() == 1 || m.rows() == 1) {
            max_rows = max_size;
            max_cols = max_size;
        }

        Eigen::Index width = 0;

        std::streamsize explicit_precision;
        if (fmt.precision == Eigen::StreamPrecision) {
            explicit_precision = 0;
        } else {
            explicit_precision = fmt.precision;
        }
        std::streamsize old_precision = 0;
        if (explicit_precision)
            old_precision = s.precision(explicit_precision);

        bool align_cols = !(fmt.flags & Eigen::DontAlignCols);
        if (align_cols) {
            // compute the largest width
            for (Eigen::Index j = 0; j < m.cols(); ++j)
                for (Eigen::Index i = 0; i < m.rows(); ++i) {
                    std::stringstream sstr;
                    sstr.copyfmt(s);
                    sstr << m.coeff(i, j);
                    width = std::max<Eigen::Index>(
                        width, Eigen::Index(sstr.str().length()));
                }
        }

        auto print_row = [fmt, width, m, max_cols](OStream &s, Eigen::Index i) {
            if (i)
                s << fmt.rowSpacer;
            s << fmt.rowPrefix;
            if (width)
                s.width(width);
            s << m.coeff(i, 0);
            if (m.cols() <= max_cols) {
                for (Eigen::Index j = 1; j < m.cols(); ++j) {
                    s << fmt.coeffSeparator;
                    if (width)
                        s.width(width);
                    s << m.coeff(i, j);
                }
            } else {
                for (Eigen::Index j = 1; j < max_cols / 2; ++j) {
                    s << fmt.coeffSeparator;
                    if (width)
                        s.width(width);
                    s << m.coeff(i, j);
                }
                // s << fmt.coeffSeparator << ".." << m.cols() - max_cols / 2 *
                // 2 << " items..";
                s << fmt.coeffSeparator << "...";
                for (Eigen::Index j = m.cols() - max_cols / 2; j < m.cols();
                     ++j) {
                    s << fmt.coeffSeparator;
                    if (width)
                        s.width(width);
                    s << m.coeff(i, j);
                }
            }
            s << fmt.rowSuffix;
            if (i < m.rows() - 1)
                s << fmt.rowSeparator;
        };

        s << fmt.matPrefix;
        if (m.rows() <= max_rows) {
            for (Eigen::Index i = 0; i < m.rows(); ++i)
                print_row(s, i);
        } else {
            for (Eigen::Index i = 0; i < max_rows / 2; ++i)
                print_row(s, i);
            if (width)
                s.width(width);
            // s << ".." << m.rows() - max_rows / 2 * 2 << " items.." <<
            // fmt.rowSeparator;
            s << "..." << fmt.rowSeparator;
            for (Eigen::Index i = m.rows() - max_rows / 2; i < m.rows(); ++i)
                print_row(s, i);
        }
        s << fmt.matSuffix;
        if (explicit_precision)
            s.precision(old_precision);
        return s;
    }
};

} // namespace logging
