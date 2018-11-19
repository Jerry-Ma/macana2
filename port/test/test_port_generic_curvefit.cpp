#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <vector>
#include <memory>
#include <random>

// #include <opencv2/opencv.hpp>
// #include "Eigen2CV.h"

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
        Gaussian1D g(params);
        auto size = xdata.size();
        ydata.resize(size);
        sigma.resize(size);
        VectorXd noise(size);
        for (int i = 0; i < noise.size(); ++i)
            noise(i) = urand(e2);
        noise *= params[0] * 0.002;
        ydata = g(xdata) + noise;
        sigma.array() = params[0] * 0.002;
        return g;
    }

    Gaussian2D generate_Gaussian2D(const VectorXd& params, const VectorXd& xdata, const VectorXd& ydata, MatrixXd& zdata, MatrixXd& sigma)
    {
        Gaussian2D g(params);
        auto nx = xdata.size();
        auto ny = ydata.size();
        zdata.resize(ny, nx);
        sigma.resize(ny, nx);
        std::cout << "gen g2d " << zdata.cols() << ", " << zdata.rows() << std::endl;
        MatrixXd noise(ny, nx);
        for (int i = 0; i < ny; ++i)
            for (int j = 0; j < nx; ++j)
                noise(i, j) = urand(e2);
        noise *= params[0] * 0.002;
        std::cout << "gen g2d " << noise.cols() << ", " << noise.rows() << std::endl;
        zdata = g(xdata, ydata) + noise;
        std::cout << "gen g2d " << zdata.cols() << ", " << zdata.rows() << std::endl;
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

TEST_F(CurveFitTest, curvefit_eigen3_gaussian1d_simple) {

    using namespace generic;

    VectorXd xdata, ydata, sigma;
    // populate with sample data
    sample0_Gaussian1D(xdata, ydata, sigma);

    // create model with initial guess
    Gaussian1D g(1., 0, 1.);

    // run the fit
    Gaussian1D g_fit = curvefit_eigen3(g, g.params, xdata, ydata, sigma);
    auto _p = g_fit.params;

    // test
    std::vector<double> p(_p.data(), _p.data() + _p.size());
    std::vector<double> expected_p = {0.501708129, -6.33013581e-11, 0.73154623};
    EXPECT_THAT(p, Pointwise(DoubleNear(1e-5), expected_p));
    // EXPECT_THAT(pp, Pointwise(NearWithPrecision(1e-5), expected_pp));
}

TEST_F(CurveFitTest, curvefit_eigen3_gaussian1d_roundtrip) {


    VectorXd xdata, ydata, sigma;
    VectorXd init_p(3);
    init_p << 5., 4., 3.;

    xdata.setLinSpaced(100, 0., 10.);
    Gaussian1D g = generate_Gaussian1D(init_p, xdata, ydata, sigma);

    // run the fit
    Gaussian1D g_fit = curvefit_eigen3(g, g.params, xdata, ydata, sigma);
    auto _p = g_fit.params;

    // test
    std::vector<double> p(_p.data(), _p.data() + _p.size());
    std::vector<double> expected_p(init_p.data(), init_p.data() + init_p.size());
    EXPECT_THAT(p, Pointwise(DoubleNear(1e-5), expected_p));
    // EXPECT_THAT(pp, Pointwise(NearWithPrecision(1e-5), expected_pp));
}

TEST_F(CurveFitTest, curvefit_eigen3_gaussian2d_roundtrip) {


    VectorXd xdata, ydata;
    MatrixXd zdata, sigma;
    VectorXd init_p(6);
    init_p << 5., 4., 3., 2., 1., 0.;

    xdata.setLinSpaced(80, 0., 8.);
    ydata.setLinSpaced(60, 0., 6.);
    Gaussian2D g = generate_Gaussian2D(init_p, xdata, ydata, zdata, sigma);
    std::cout << xdata.rows() << xdata.cols() << std::endl;
    std::cout << ydata.rows() << ydata.cols() << std::endl;
    std::cout << zdata.rows() << zdata.cols() << std::endl;
    std::cout << sigma.rows() << sigma.cols() << std::endl;

    /*
    cv::Mat _z = octane::eigen2cv(zdata);

    std::cout << _z.rows << _z.cols << std::endl;


    cv::namedWindow("Display Image", cv::WINDOW_AUTOSIZE );
    cv::imshow("Display Image", _z);
    cv::waitKey(0);
    */

    // run the fit
    Gaussian2D g_fit = curvefit_eigen3(
                g, g.params,
                g.meshgrid(xdata, ydata), zdata, sigma);
    auto _p = g_fit.params;

    // test
    std::vector<double> p(_p.data(), _p.data() + _p.size());
    std::vector<double> expected_p(init_p.data(), init_p.data() + init_p.size());
    EXPECT_THAT(p, Pointwise(DoubleNear(1e-5), expected_p));
    // EXPECT_THAT(pp, Pointwise(NearWithPrecision(1e-5), expected_pp));
}
}