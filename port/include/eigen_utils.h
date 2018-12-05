#ifndef EIGEN_UTILS_H
#define EIGEN_UTILS_H

#include <Eigen/Core>

namespace Eigen {

using VectorXI = Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>;
using MatrixXI = Eigen::Matrix<Eigen::Index, Eigen::Dynamic, Eigen::Dynamic>;

}

namespace utils {

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

template <typename _Scalar>
class Interval: public Eigen::Array<_Scalar, 2, 1>
{
public:
    using Scalar = _Scalar;
    using Base = Eigen::Array<_Scalar, 2, 1>;
    inline static const Scalar Inf = std::numeric_limits<Scalar>::infinity();

    Interval(): Base() {
        this->operator[](0) = -Inf;
        this->operator[](1) = Inf;
    }
    Interval(std::initializer_list<Scalar> interval): Interval() {
        if (interval.size() == 0) return;
        if (interval.size() == 2)
        {
           auto it = interval.begin();
           this->operator[](0) = *it;
           ++it;
           this->operator[](1) = *it;
           return;
        }
        throw std::invalid_argument("empty initalize_list or {left, right} is required");
    }

    template<typename OtherDerived>
    Interval(const Eigen::ArrayBase<OtherDerived>& other): Base(other) {}

    template<typename OtherDerived>
    Interval& operator=(const Eigen::ArrayBase<OtherDerived>& other)
    {
        this->Base::operator=(other);
        return *this;
    }
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Scalar& left() {
        return this->x();
    }
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Scalar& right() {
        return this->y();
    }
};

} // namespace utils

#endif /* !EIGEN_UTILS_H */
