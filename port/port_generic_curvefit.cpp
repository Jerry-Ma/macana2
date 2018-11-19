#include <iostream>

#include <unsupported/Eigen/NonLinearOptimization>
#include "port_generic_curvefit.h"

namespace generic {

Gaussian1D::Gaussian1D(double amplitude, double mean, double stddev)
    : Model<3, 1>({amplitude, mean, stddev}) {}

Gaussian1D::ValueType Gaussian1D::eval(const Gaussian1D::InputType& p, const Gaussian1D::InputDataType& x) const
{
    return p[0] * (-0.5 * (x.array() - p[1]).square() / p[2] / p[2]).exp();
}

Gaussian1D::ValueType Gaussian1D::operator() (const Gaussian1D::InputType& p, const Gaussian1D::InputDataType& x) const
{
    return eval(p, x);
}

Gaussian1D::ValueType Gaussian1D::operator() (const Gaussian1D::InputDataType& x) const
{
    return operator() (this->params, x);
}


Gaussian2D::Gaussian2D(double amplitude, double xmean, double ymean, double xstddev, double ystddev, double theta)
    //             0          1      2      3        4        5
    : Model<6, 2>({amplitude, xmean, ymean, xstddev, ystddev, theta}) {}

Gaussian2D::ValueType Gaussian2D::eval(const Gaussian2D::InputType& p, const Gaussian2D::InputDataType& xy) const
{
    double cost2 = cos(p[5]) * cos(p[5]);
    double sint2 = sin(p[5]) * sin(p[5]);
    double sin2t = sin(2. * p[5]);
    double xstd2 = p[3] * p[3];
    double ystd2 = p[4] * p[4];
    double a = - 0.5 * ((cost2 / xstd2) + (sint2 / ystd2));
    double b = - 0.5 * ((sin2t / xstd2) - (sin2t / ystd2));
    double c = - 0.5 * ((sint2 / xstd2) + (cost2 / ystd2));
    std::cout << "g2d eval" << xy.cols() << ", " << xy.rows() << std::endl;
    return p[0] * (
                    (xy.col(0).array() - p[1]).square() * a +
                    (xy.col(0).array() - p[1]) * (xy.col(1).array() - p[2]) * b +
                    (xy.col(1).array() - p[2]).square() * c
                ).exp();
}

Gaussian2D::DataType Gaussian2D::operator() (
        const Gaussian2D::InputType& p,
        const Gaussian2D::InputDataBasisType& x,
        const Gaussian2D::InputDataBasisType& y) const
{
    return eval(p, meshgrid(x, y));
}

Gaussian2D::DataType Gaussian2D::operator() (
        const Gaussian2D::InputDataBasisType& x,
        const Gaussian2D::InputDataBasisType& y) const
{
    return operator() (this->params, x, y);
}

}  // namespace generic