#pragma once

#pragma clang diagnostic ignored "-Wreserved-id-macro"
#define _USE_MATH_DEFINES
#include <cmath>

#include <Eigen/Core>
#ifndef EIGEN_NO_DEBUG
#else
#undef EIGEN_NO_DEBUG
#endif // !EIGEN_NO_DEBUG
#define EIGEN_NO_MALLOC

#include <logging.h>

namespace Eigen {

using VectorXI = Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>;
using MatrixXI = Eigen::Matrix<Eigen::Index, Eigen::Dynamic, Eigen::Dynamic>;

} // namespace Eigen

namespace eiu {

namespace internal {
template <typename Derived>
struct has_storage
    : std::is_base_of<Eigen::PlainObjectBase<std::decay_t<Derived>>,
                      std::decay_t<Derived>> {};

} // namespace internal

template <typename Derived>
using const_ref = typename Eigen::internal::ref_selector<Derived>::type;

template <typename Derived>
using ref = typename Eigen::internal::ref_selector<Derived>::non_const_type;


template <typename UnaryOp> struct MoveEnabledUnaryOp {
    MoveEnabledUnaryOp(const UnaryOp &func = UnaryOp()) : m_func(func) {}

    template <typename T, typename... Args>
    decltype(auto) operator()(T &&in, Args &&... args) {
        auto logger = logging::createLogger("move_enabled_unary_op", this);
        if constexpr (std::is_lvalue_reference<T>::value ||
                      !internal::has_storage<T>::value) {
            // lvalue ref, either expression or non-expression
            // return expression that builds on top of input
            SPDLOG_LOGGER_TRACE(logger, "passed in lvalue reference");
            return m_func(std::forward<T>(in),
                          std::forward<decltype(args)>(args)...);
            // NOLINTNEXTLINE(readability-else-after-return)
        } else {
            // rvalue ref
            // in this case we need to call the function and update inplace
            // first and move to return
            SPDLOG_LOGGER_TRACE(logger, "passed in rvalue reference");
            in = m_func(std::forward<T>(in),
                        std::forward<decltype(args)>(args)...);
            return std::forward<T>(in);
        }
    }

  protected:
    const UnaryOp &m_func;
};

template <typename T = Eigen::ArrayXd> T *ei_nullptr() {
    return static_cast<T *>(nullptr);
}

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
    }
    // even sized vector -> average the two middle values
    auto max_it = std::max_element(vv.begin(), vv.begin() + n);
    return (*max_it + vv[n]) / 2.0;
}

template <typename _Scalar>
class Interval : public Eigen::Array<_Scalar, 2, 1> {
  public:
    using Scalar = _Scalar;
    using Base = Eigen::Array<_Scalar, 2, 1>;
    inline static const Scalar inf = std::numeric_limits<Scalar>::infinity();

    Interval() : Base() {
        this->operator[](0) = -inf;
        this->operator[](1) = inf;
    }
    Interval(std::initializer_list<Scalar> interval) : Interval() {
        if (interval.size() == 0)
            return;
        if (interval.size() == 2) {
            auto it = interval.begin();
            this->operator[](0) = *it;
            ++it;
            this->operator[](1) = *it;
            return;
        }
        throw std::invalid_argument(
            "empty initalize_list or {left, right} is required");
    }

    template <typename OtherDerived>
    Interval(const Eigen::ArrayBase<OtherDerived> &other) : Base(other) {}

    template <typename OtherDerived>
    Interval &operator=(const Eigen::ArrayBase<OtherDerived> &other) {
        this->Base::operator=(other);
        return *this;
    }
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Scalar &left() { return this->x(); }
    EIGEN_DEVICE_FUNC
    EIGEN_STRONG_INLINE Scalar &right() { return this->y(); }
};

} // namespace eiu