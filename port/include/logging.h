#ifndef LOGGING_H
#define LOGGING_H

#include <Eigen/Core>

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/fmt/ostr.h>
#undef SPDLOG_LOGGER_TRACE
#define SPDLOG_LOGGER_TRACE(logger, ...) logger->log(spdlog::source_loc{}, spdlog::level::trace, __VA_ARGS__)

namespace logging {

using namespace Eigen;

template <typename M>
struct pprint
{
    pprint(M* m): m(m) {}
    template <typename OStream>
    friend OStream &operator<<(OStream &os, const pprint& pp)
    {
        const M& m = *(pp.m);
        if (m.size() == 0) return os << "(empty)";
        os << "(" << m.rows() << "," << m.cols() << ")";

        auto f = pp.matrix_formatter();
        pp.print_matrix(os, m.eval(), f, 5, 5, 10);
        return os;
    }

private:
    const M* m;
    IOFormat matrix_formatter() const
    {
        if (m->cols() == 1)
            return IOFormat(StreamPrecision, DontAlignCols, ", ", ", ", "", "", "[", "]");
        else if (m->cols() < 3)
            return IOFormat(StreamPrecision, DontAlignCols, ", ", " ", "[", "]", "[", "]");
        else
            return IOFormat(StreamPrecision, 0, ", ", "\n", "[", "]", "[\n", "]\n");
    }
    template <typename OStream, typename Derived>
    OStream & print_matrix(OStream & s, const Derived& _m, const IOFormat& fmt, int max_rows, int max_cols, int max_size=-1) const
    {
      if(_m.size() == 0)
      {
        s << fmt.matPrefix << fmt.matSuffix;
        return s;
      }

      typename Derived::Nested m = _m;

      if (max_size < 0) max_size = max_rows * max_cols;
      if (m.cols() == 1 || m.rows() == 1)
      {
          max_rows = max_size;
          max_cols = max_size;
      }

      Index width = 0;

      std::streamsize explicit_precision;
      if(fmt.precision == StreamPrecision)
      {
        explicit_precision = 0;
      }
      else
      {
        explicit_precision = fmt.precision;
      }
      std::streamsize old_precision = 0;
      if(explicit_precision) old_precision = s.precision(explicit_precision);

      bool align_cols = !(fmt.flags & DontAlignCols);
      if(align_cols)
      {
        // compute the largest width
        for(Index j = 0; j < m.cols(); ++j)
          for(Index i = 0; i < m.rows(); ++i)
          {
            std::stringstream sstr;
            sstr.copyfmt(s);
            sstr << m.coeff(i,j);
            width = std::max<Index>(width, Index(sstr.str().length()));
          }
      }

      auto print_row = [fmt, width, m, max_cols] (OStream& s, Index i) {
        if (i)
          s << fmt.rowSpacer;
        s << fmt.rowPrefix;
        if(width) s.width(width);
        s << m.coeff(i, 0);
        if (m.cols() <= max_cols) {
            for(Index j = 1; j < m.cols(); ++j)
            {
              s << fmt.coeffSeparator;
              if (width) s.width(width);
              s << m.coeff(i, j);
            }
        } else {
            for (Index j = 1; j < max_cols / 2; ++j)
            {
              s << fmt.coeffSeparator;
              if (width) s.width(width);
              s << m.coeff(i, j);
            }
            // s << fmt.coeffSeparator << ".." << m.cols() - max_cols / 2 * 2 << " items..";
            s << fmt.coeffSeparator << "...";
            for (Index j = m.cols() -  max_cols / 2; j < m.cols(); ++j)
            {
              s << fmt.coeffSeparator;
              if (width) s.width(width);
              s << m.coeff(i, j);
            }
        }
        s << fmt.rowSuffix;
        if( i < m.rows() - 1)
          s << fmt.rowSeparator;
      };

      s << fmt.matPrefix;
      if (m.rows() <= max_rows) {
          for(Index i = 0; i < m.rows(); ++i)
              print_row(s, i);
      } else {
          for (Index i = 0; i < max_rows / 2; ++i)
              print_row(s, i);
          if(width) s.width(width);
          // s << ".." << m.rows() - max_rows / 2 * 2 << " items.." << fmt.rowSeparator;
          s << "..." << fmt.rowSeparator;
          for (Index i = m.rows() - max_rows / 2; i < m.rows(); ++i)
              print_row(s, i);
      }
      s << fmt.matSuffix;
      if(explicit_precision) s.precision(old_precision);
      return s;
    }
};

inline const static auto console = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();

} // namespace
#endif /* !LOGGING_H */
