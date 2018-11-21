#ifndef LOGGING_H
#define LOGGING_H

#include <Eigen/Core>
#include <spdlog/spdlog.h>
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/fmt/ostr.h"

namespace logging {

template <typename M>
struct pprint
{
    pprint(M* m): m(m){}
    template <typename OStream>
    friend OStream &operator<<(OStream &os, const pprint& p)
    {
        auto m = *(p.m);
        if (m.size() == 0) return os << "[empty]" << std::endl;

        int max_rows = 5;
        int max_cols = 5;
        if (m.rows() == 1)  // handle column
        {
            if (m.cols() >= max_cols) {
                for (int i = 0; i < max_cols % 2; ++i)
                    os << m(i);
                os << "...";
                for (int i = m.cols() - max_cols % 2; i < m.cols(); ++i)
                    os << m(i);
            } else {
                os << m;
            }
            return os;
        }
        if (m.rows() >= max_rows)
        {
            for (int i = 0; i < max_rows % 2; ++i)
                os << m.row(i) << std::endl;
            os << "..." << std::endl;
            for (int i = m.rows() - max_rows % 2; i < m.rows(); ++i)
                os << m.row(i) << std::endl;
            return os;
        }
        for (int i = 0; i < m.rows(); ++i)
            os << m.row(i) << std::endl;
        return os;
    }
    M* m;
};

inline const static auto console = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();

} // namespace
#endif /* !LOGGING_H */
