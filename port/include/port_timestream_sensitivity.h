#ifndef PORT_TIMESTREAM_SENSITIVITY_H
#define PORT_TIMESTREAM_SENSITIVITY_H

#include <Eigen/Core>

namespace Eigen {

using VectorXI = Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>;
using MatrixXI = Eigen::Matrix<Eigen::Index, Eigen::Dynamic, Eigen::Dynamic>;

}

namespace timestream {

namespace details {

template <typename Derived>
typename Derived::Scalar median(const Eigen::MatrixBase<Derived> &v_) {
    if (v_.size() == 0)
        throw std::runtime_error("cannot take median of empty vector");
    // make a copy of the original vector
    Derived v(v_);

    // convert to std vector
    std::vector<typename Derived::Scalar> vv;
    vv.resize(v.size());
    Derived::Map(&vv[0], v.size()) = v; // assign v to vv

    size_t n = vv.size() / 2;
    std::nth_element(vv.begin(), vv.begin() + n, vv.end());

    if (vv.size() % 2) {
        return vv[n];
    } else {
        // even sized vector -> average the two middle values
        auto max_it = std::max_element(vv.begin(), vv.begin() + n);
        return (*max_it + vv[n]) / 2.0;
    }
}

} // namespace details

void psd(const Eigen::VectorXd &scan, Eigen::VectorXd &freqs,
         Eigen::VectorXd &psd, double samplerate, bool hann = true);

double sensitivity(const Eigen::VectorXd &scans,
                   const Eigen::MatrixXI &scanindex, double samplerate,
                   double gain = 0.,
                   std::pair<double, double> freqrange = {3., 5.});

} // namespace timestream

#endif // PORT_TIMESTREAM_SENSITIVITY_H
