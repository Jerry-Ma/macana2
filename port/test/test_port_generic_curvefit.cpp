#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <vector>
#include <memory>
#include <random>

#include <Eigen/Dense>

#include "port_generic_curvefit.h"

namespace {

using namespace generic;
using Eigen::VectorXd;

class CurveFitTest : public ::testing::Test
{
protected:
    CurveFitTest(): e2(sd()) {}
    ~CurveFitTest() override {}
    void SetUp() override {
    }
    void TearDown() override {}

    std::random_device sd;  // seed generator
    std::mt19937 e2;  // rand generator
    std::uniform_real_distribution<double> urand{-0.5, 0.5};

    void sample0_Gaussian1D(VectorXd& xdata, VectorXd& ydata, VectorXd& sigma)
    {
        xdata.resize(5);
        ydata.resize(5);
        sigma.resize(5);
        xdata << -2., -1., 0., 1., 2.;
        ydata << 0., 0.2, 0.5, 0.2, 0.;
        sigma.array() = 0.01;
    }

    Gaussian1D generate_Gaussian1D(const VectorXd& params, const VectorXd& xdata, VectorXd& ydata, VectorXd& sigma)
    {
        Gaussian1D g{params[0], params[1], params[2]};
        ydata.resize(xdata.size());
        sigma.resize(xdata.size());
        VectorXd noise(xdata.size());
        for (int i = 0; i < noise.size(); ++i)
            noise << urand(e2);
        noise *= params[0] * 0.002;
        ydata = g(xdata) + noise;
        sigma.array() = params[0] * 0.002;
        return g;
    }

};

using namespace testing;

/*
MATCHER_P(NearWithPrecision, precision, "") {
    return abs(get<0>(arg) - get<1>(arg)) < precision * 0.5 * abs(get<0>(arg) + get<1>(arg));
}
*/

TEST_F(CurveFitTest, curvefit_eigen3_simple) {

    using namespace generic;

    VectorXd xdata, ydata, sigma;
    // populate with sample data
    sample0_Gaussian1D(xdata, ydata, sigma);

    // create model with initial guess
    Gaussian1D g(1., 0, 1.);

    // run the fit
    Gaussian1D::ValueType _p = curvefit_eigen3(g, g.params(), xdata, ydata, sigma);

    // test
    std::vector<double> p(_p.data(), _p.data() + _p.size());
    std::vector<double> expected_p = {0.501708129, -6.33013581e-11, 0.73154623};
    EXPECT_THAT(p, Pointwise(DoubleNear(1e-5), expected_p));
    // EXPECT_THAT(pp, Pointwise(NearWithPrecision(1e-5), expected_pp));
}

TEST_F(CurveFitTest, curvefit_eigen3_roundtrip) {


    VectorXd xdata, ydata, sigma;
    VectorXd init_p(3);
    init_p << 5., 4., 3.;

    xdata.setLinSpaced(100, 0., 10.);
    Gaussian1D g = generate_Gaussian1D(init_p, xdata, ydata, sigma);

    // run the fit
    Gaussian1D::ValueType _p = curvefit_eigen3(g, g.params(), xdata, ydata, sigma);

    // test
    std::vector<double> p(_p.data(), _p.data() + _p.size());
    std::vector<double> expected_p(init_p.data(), init_p.data() + init_p.size());
    EXPECT_THAT(p, Pointwise(DoubleNear(1e-5), expected_p));
    // EXPECT_THAT(pp, Pointwise(NearWithPrecision(1e-5), expected_pp));
}

}