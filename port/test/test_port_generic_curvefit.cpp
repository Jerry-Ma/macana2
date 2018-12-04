#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <vector>
#include <memory>
#include <random>

#include <opencv2/opencv.hpp>
#include <Eigen/Core>
#include <opencv2/core/eigen.hpp>

#define CVUI_IMPLEMENTATION
#include <opencv/cvui.h>
#include <opencv/enhancedwindow.h>

#include <port_generic_curvefit.h>

namespace {

using namespace generic;
using Eigen::VectorXd;
using Eigen::MatrixXd;

class CurveFitTest : public ::testing::Test
{
protected:
    CurveFitTest(): e2(sd()) {
        logger = logging::createLogger("test.curvefit", this);
    }
    ~CurveFitTest() override {}
    void SetUp() override {}
    void TearDown() override {}

    std::random_device sd;  // seed generator
    std::mt19937 e2;  // rand generator
    // std::uniform_real_distribution<double> urand{-1., 1.};
    const double PI = std::atan(1.0) * 4;
    std::shared_ptr<spdlog::logger> logger = nullptr;

    template <typename Derived>
    void add_GaussianNoise(double stddev, DenseBase<Derived>& data, DenseBase<Derived>& sigma)
    {
        auto ndist = std::normal_distribution<double>(0, stddev);
        auto ny = data.rows();
        auto nx = data.cols();
        Derived noise;
        noise.resize(ny, nx);
        for (int i = 0; i < ny; ++i)
            for (int j = 0; j < nx; ++j)
                noise(i, j) = ndist(e2);
        data += noise;

        sigma.derived().resize(ny, nx);
        sigma.derived().array() = stddev;
        logger->debug("add gaussian noise of stddev={}", stddev);
        logger->debug("noise{}", logging::pprint(&noise));
        logger->debug("sigma{}", logging::pprint(&sigma));
        logger->debug("data{}", logging::pprint(&data));
    }

    void sample0_Gaussian1D(VectorXd& xdata, VectorXd& ydata, VectorXd& sigma)
    {
        xdata.resize(5);
        ydata.resize(5);
        sigma.resize(5);
        xdata << -2., -1., 0., 1., 2.;
        ydata << 0., 0.2, 0.5, 0.2, 0.;
        sigma.array() = 0.01;
    }

    template <typename Model>
    typename std::enable_if<Model::DimensionsAtCompileTime == 1, Model>::type
    generate_ModelData(const VectorXd& params, const VectorXd& xdata, VectorXd& ydata)
    {
        Model g(params);
        ydata.resize(xdata.size());
        ydata = g(xdata);
        logger->debug("generate {}", g);
        logger->debug("params{}", logging::pprint(&params));
        logger->debug("xdata{}", logging::pprint(&xdata));
        logger->debug("ydata{}", logging::pprint(&ydata));
        return g;
    }
    
    template <typename Model>
    typename std::enable_if<Model::DimensionsAtCompileTime == 2, Model>::type
    generate_ModelData(const VectorXd& params, const VectorXd& xdata, const VectorXd& ydata, MatrixXd& zdata)
    {
        Model g(params);
        auto nx = xdata.size();
        auto ny = ydata.size();
        zdata.resize(ny, nx);
        Map<VectorXd>(zdata.data(), zdata.size()) = g(xdata, ydata);

        logger->debug("generate {}", g);
        logger->debug("params{}", logging::pprint(&params));
        logger->debug("xdata{}", logging::pprint(&xdata));
        logger->debug("ydata{}", logging::pprint(&ydata));
        logger->debug("zdata{}", logging::pprint(&zdata));
        return g;
    }

    void plot2d(const VectorXd& xdata, const VectorXd& ydata, const MatrixXd& zdata)
    {
        cv::Mat _data, data;
        cv::eigen2cv(zdata, _data);
        // stats
        double min, max, mean, stddev;
        cv::Scalar _mean,  _stddev;
        cv::minMaxIdx(_data, &min, &max);
        cv::meanStdDev(_data, _mean, _stddev);
        mean = _mean.val[0];
        stddev = _stddev.val[0];
        cv::normalize(_data, data, 0, 255, cv::NORM_MINMAX, CV_8U);
        cv::cvtColor(data, data, cv::COLOR_GRAY2BGR);
        logger->debug("plot 2d data({}, {})", data.rows, data.cols);

        int minWidth = 800;
        int minHeight = 600;
        int padding = 20;
        int width = std::max(minWidth, data.cols + padding);
        int height = std::max(minHeight, data.rows + padding);
        cv::Mat frame = cv::Mat(height, width, data.type());

        // cut levels
        double low = 0.5;
        double high = 0.995;
        double center = 0.5 * (low + high);

        // settings window using the EnhancedWindow class.
        EnhancedWindow settings(padding, padding, 270, 180, "Limits");
        // render a rectangle on the screen as data viewer
        cv::Rect view(width - padding - data.cols, padding, data.cols, data.rows);

        // init cvui and tell it to create a OpenCV window, i.e. cv::namedWindow(WINDOW_NAME).
        auto wname = fmt::format("{}: {}", logger->name(), "plot2d");
        cvui::init(wname, 20);

        while (true) {
            // reset
            frame = cv::Scalar(49, 52, 49);
            // Render the settings window and its content, if it is not minimized.
            settings.begin(frame);
            if (!settings.isMinimized()) {
                cvui::beginRow();
                    cvui::text("lo");
                    cvui::trackbar(settings.width() - 20, &low, 0., 1.);
                cvui::endRow();
                cvui::beginRow();
                    cvui::text("hi");
                    cvui::trackbar(settings.width() - 20, &high, 0., 1.);
                cvui::endRow();
                cvui::beginRow();
                    cvui::text("ct");
                    cvui::trackbar(settings.width() - 20, &center, 0., 1.);
                cvui::endRow();
            }
            settings.end();

            // update all cvui internal stuff, e.g. handle mouse clicks, and show
            // everything on the screen.
            cvui::rect(frame, view.x, view.y, view.width, view.height, 0xff0000);
            cvui::image(frame, view.x, view.y, data);

            // handle mouseover info
            int status = cvui::iarea(view.x, view.y, view.width, view.height);

            // cvui::iarea() will return the current mouse status:
            //  CLICK: mouse just clicked the interaction are
            //	DOWN: mouse button was pressed on the interaction area, but not released yet.
            //	OVER: mouse cursor is over the interaction area
            //	OUT: mouse cursor is outside the interaction area
            switch (status) {
                case cvui::CLICK:   break;
                case cvui::DOWN:    break;
                case cvui::OVER:    {
                    int x = cvui::mouse().x - view.x;
                    int y = cvui::mouse().y - view.y;
                    cvui::printf(
                            frame, padding, height - padding,
                            fmt::format("mouse pos: ({}, {}) x={} y={} z={} v={}",
                                         x, y, xdata[x], ydata[y], zdata.coeff(y, x), _data.at<double>(y, x)).c_str());
                                    break;
                            }
                case cvui::OUT:     cvui::printf(
                            frame, padding, height - padding,
                            fmt::format("data({}, {}) min={} max={} mean={} med={} stddev={}",
                                        data.rows, data.cols, min, max, mean, "N/A", stddev).c_str());
                                    break;
            }
            cvui::update();
            cvui::imshow(wname, frame);
            // check if ESC was pressed
            if (cv::waitKey(20) == 27) {
                break;
            }
        }
    }
};

using namespace testing;

MATCHER_P(NearWithPrecision, precision, "") {
    return abs(get<0>(arg) - get<1>(arg)) < precision * 0.5 * abs(get<0>(arg) + get<1>(arg));
}

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
    // EXPECT_THAT(p, Pointwise(DoubleNear(1e-5), expected_p));
    EXPECT_THAT(p, Pointwise(NearWithPrecision(1e-5), expected_p));
}

TEST_F(CurveFitTest, curvefit_eigen3_gaussian1d_roundtrip) {


    VectorXd xdata, ydata, sigma;
    VectorXd init_p(3);
    init_p << 5., 4., 3.;

    double nsr = 0.001;
    xdata.setLinSpaced(100, 0., 10.);
    auto g = generate_ModelData<Gaussian1D>(init_p, xdata, ydata);
    add_GaussianNoise(init_p[0] * nsr, ydata, sigma);

    // run the fit
    auto g_fit = curvefit_eigen3(g, g.params.reverse(), xdata, ydata, sigma);
    auto _p = g_fit.params;

    // test
    std::vector<double> p(_p.data(), _p.data() + _p.size());
    std::vector<double> expected_p(init_p.data(), init_p.data() + init_p.size());
    // EXPECT_THAT(p, Pointwise(DoubleNear(1e-5), expected_p));
    EXPECT_THAT(p, Pointwise(NearWithPrecision(nsr), expected_p));
}

TEST_F(CurveFitTest, curvefit_eigen3_gaussian2d_roundtrip) {


    VectorXd xdata, ydata;
    MatrixXd zdata, sigma;
    VectorXd init_p(6);
    init_p << 5., 4., 3., 2., 1., PI / 4;

    xdata.setLinSpaced(80, 0., 8.);
    ydata.setLinSpaced(60, 0., 6.);
    auto g = generate_ModelData<Gaussian2D>(init_p, xdata, ydata, zdata);

    double nsr = 0.001;
    add_GaussianNoise(init_p[0] * nsr, zdata, sigma);
    // plot2d(xdata, ydata, zdata);

    // run the fit
    auto g_fit = curvefit_eigen3(
                g, g.params.reverse(),
                g.meshgrid(xdata, ydata), zdata, sigma);
    auto _p = g_fit.params;

    // test
    std::vector<double> p(_p.data(), _p.data() + _p.size());
    std::vector<double> expected_p(init_p.data(), init_p.data() + init_p.size());
    // EXPECT_THAT(p, Pointwise(DoubleNear(1e-5), expected_p));
    EXPECT_THAT(p, Pointwise(NearWithPrecision(nsr), expected_p));
}

TEST_F(CurveFitTest, curvefit_eigen3_symmetricgaussian2d_roundtrip) {


    VectorXd xdata, ydata;
    MatrixXd zdata, sigma;
    VectorXd init_p(4);
    init_p << 5., 4., 3., 2.;

    xdata.setLinSpaced(80, 0., 8.);
    ydata.setLinSpaced(60, 0., 6.);
    auto g = generate_ModelData<SymmetricGaussian2D>(init_p, xdata, ydata, zdata);
    double nsr = 0.001;
    add_GaussianNoise(init_p[0] * nsr, zdata, sigma);
    // plot2d(xdata, ydata, zdata);

    // run the fit
    auto g_fit = curvefit_eigen3(
                g, g.params.reverse(),
                g.meshgrid(xdata, ydata), zdata, sigma);
    auto _p = g_fit.params;

    // test
    std::vector<double> p(_p.data(), _p.data() + _p.size());
    std::vector<double> expected_p(init_p.data(), init_p.data() + init_p.size());
    // EXPECT_THAT(p, Pointwise(DoubleNear(1e-5), expected_p));
    EXPECT_THAT(p, Pointwise(NearWithPrecision(nsr), expected_p));
}

TEST_F(CurveFitTest, curvefit_ceres_gaussian2d_roundtrip) {


    VectorXd xdata, ydata;
    MatrixXd zdata, sigma;
    VectorXd init_p(6);
    init_p << 5., 4., 3., 2., 1., PI / 4;

    xdata.setLinSpaced(80, 0., 8.);
    ydata.setLinSpaced(60, 0., 6.);
    auto g = generate_ModelData<Gaussian2D>(init_p, xdata, ydata, zdata);
    double nsr = 0.001;
    add_GaussianNoise(init_p[0] * nsr, zdata, sigma);
    // plot2d(xdata, ydata, zdata);

    // run the fit
    auto _p = g.params;
    auto g_fit(g);
    auto xy = g.meshgrid(xdata, ydata);
    // for (int i = 0; i < 7000; ++i)
    curvefit_ceres(
                    g, _p,
                    xy, zdata, sigma);
    _p = g_fit.params;

    // test
    std::vector<double> p(_p.data(), _p.data() + _p.size());
    std::vector<double> expected_p(init_p.data(), init_p.data() + init_p.size());
    // EXPECT_THAT(p, Pointwise(DoubleNear(1e-5), expected_p));
    EXPECT_THAT(p, Pointwise(NearWithPrecision(nsr), expected_p));
}

}  // namespace