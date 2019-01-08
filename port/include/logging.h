#pragma once

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#define SPDLOG_FMT_EXTERNAL
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
// remove the extra filename:lineno in trace but use function instead
#undef SPDLOG_LOGGER_TRACE
#define SPDLOG_LOGGER_TRACE(logger, ...)                                       \
    logger->log(spdlog::source_loc{__FUNCTION__, __LINE__}, spdlog::level::trace, __VA_ARGS__)

#include <fmt/format.h>
#include <fmt/ostream.h>
#include <Eigen/Core>

namespace logging {

// shared sink
inline const static auto console =
    std::make_shared<spdlog::sinks::stdout_color_sink_mt>();

// this will create and return a runtime logger object to console
inline auto createLogger(std::string_view name) {
    return std::make_shared<spdlog::logger>(std::string(name), console);
}

// optionally we can use a unique name
inline auto createLogger(std::string_view name, const void* const ptr) {
    return createLogger(fmt::format("{}&{:x}", name, reinterpret_cast<std::uintptr_t>(ptr)));
}

// macro to insert private logger to class
// #define SPDLOG_LOGGER(name) \
//    private: const std::shared_ptr<spdlog::logger> logger = logging::logger(name);

using Eigen::DenseBase;
using Eigen::DontAlignCols;
using Eigen::Dynamic;
using Eigen::Index;
using Eigen::IOFormat;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::StreamPrecision;

namespace internal {

// pretty print larget Eigen array by only showing part of it
template <typename OStream, typename Derived>
OStream &pprint_matrix(OStream &s, const DenseBase<Derived> &_m,
                       const IOFormat &fmt, int max_rows, int max_cols,
                       int max_size = -1) {
    const Derived& m(_m.derived());

    if (m.size() == 0) {
        s << fmt.matPrefix << fmt.matSuffix;
        return s;
    }

    if (max_size < 0)
        max_size = max_rows * max_cols;
    if (m.cols() == 1 || m.rows() == 1) {
        max_rows = max_size;
        max_cols = max_size;
    }

    Index width = 0;

    std::streamsize explicit_precision;
    if (fmt.precision == StreamPrecision) {
        explicit_precision = 0;
    } else {
        explicit_precision = fmt.precision;
    }
    std::streamsize old_precision = 0;
    if (explicit_precision)
        old_precision = s.precision(explicit_precision);

    bool align_cols = !(fmt.flags & DontAlignCols);
    // if (align_cols) {
    if (false) {
        // compute the largest width
        for (Index j = 0; j < m.cols(); ++j)
            for (Index i = 0; i < m.rows(); ++i) {
                std::stringstream sstr;
                sstr.copyfmt(s);
                sstr << m.coeff(i, j);
                width = std::max<Index>(width, Index(sstr.str().length()));
            }
    }

    auto print_row = [fmt, width, m, max_cols](OStream &s, Index i) {
        if (i)
            s << fmt.rowSpacer;
        s << fmt.rowPrefix;
        if (width)
            s.width(width);
        s << m.coeff(i, 0);
        if (m.cols() <= max_cols) {
            for (Index j = 1; j < m.cols(); ++j) {
                s << fmt.coeffSeparator;
                if (width)
                    s.width(width);
                s << m.coeff(i, j);
            }
        } else {
            for (Index j = 1; j < max_cols / 2; ++j) {
                s << fmt.coeffSeparator;
                if (width)
                    s.width(width);
                s << m.coeff(i, j);
            }
            // s << fmt.coeffSeparator << ".." << m.cols() - max_cols / 2 * 2 <<
            // " items..";
            s << fmt.coeffSeparator << "...";
            for (Index j = m.cols() - max_cols / 2; j < m.cols(); ++j) {
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
        for (Index i = 0; i < m.rows(); ++i)
            print_row(s, i);
    } else {
        for (Index i = 0; i < max_rows / 2; ++i)
            print_row(s, i);
        if (width)
            s.width(width);
        // s << ".." << m.rows() - max_rows / 2 * 2 << " items.." <<
        // fmt.rowSeparator;
        s << "..." << fmt.rowSeparator;
        for (Index i = m.rows() - max_rows / 2; i < m.rows(); ++i)
            print_row(s, i);
    }
    s << fmt.matSuffix;
    if (explicit_precision)
        s.precision(old_precision);
    return s;
}

} // namespace internal

// wrapper class to log some objects that is mappable to eigen matrix
template <typename _Derived> struct pprint {
    using Derived =
        typename Eigen::internal::ref_selector<typename std::conditional<
            std::is_arithmetic<typename std::decay<_Derived>::type>::value,
            Map<const MatrixXd>, _Derived>::type>::type;
    // using a = typename Derived::nothing;

    pprint(const DenseBase<_Derived> &m) : m(m.derived()) {}

    // template <typename = std::enable_if< (NRows ==Dynamic || NCols ==
    // Dynamic)>>
    // pprint(const DenseBase<_Derived> &m) : m(m.derived()) {}

    pprint(const _Derived *m, Index nrows, Index ncols)
        : m(Derived(m, nrows, ncols)) {}
    pprint(const _Derived *m, Index size) : pprint(m, size, 1) {}

    template <typename OStream>
    friend OStream &operator<<(OStream &os, const pprint &pp) {
        auto &&m = pp.m;
        if (m.size() == 0)
            return os << "(empty)";
        os << "(" << m.rows() << "," << m.cols() << ")";

        auto f = pp.matrix_formatter();
        return internal::pprint_matrix(os, m, f, 5, 5, 10);
    }

  private:
    Derived m;
    IOFormat matrix_formatter() const {
        if (m.cols() == 1)
            return IOFormat(StreamPrecision, DontAlignCols, ", ", ", ", "", "",
                            "[", "]");
        else if (m.cols() < 3)
            return IOFormat(StreamPrecision, DontAlignCols, ", ", " ", "[", "]",
                            "[", "]");
        else
            return IOFormat(StreamPrecision, 0, ", ", "\n", "[", "]", "[\n",
                            "]\n");
    }
};

} // namespace logging

namespace logging1 {

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
