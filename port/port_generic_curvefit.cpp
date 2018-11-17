#include <iostream>

#include <unsupported/Eigen/NonLinearOptimization>
#include "port_generic_curvefit.h"

namespace generic {

using IType = Model::InputType;
using VType = Model::ValueType;

Model::Model(): _Base() {}

Model::Model(int inputs): _Base(inputs, Dynamic), m_params(inputs) {}

IType& Model::params() {return m_params;}

Gaussian1D::Gaussian1D(double amplitude, double mean, double stddev)
    :Model(3)
{
    m_params << amplitude, mean, stddev;
}

VType Gaussian1D::operator() (const IType& p, const VType& x) const
{
    return p[0] * (-0.5 * (x.array() - p[1]).square() / p[2] / p[2]).exp();
}

VType Gaussian1D::operator() (const VType& x) const
{
    return operator() (m_params, x);
}

}  // namespace generic