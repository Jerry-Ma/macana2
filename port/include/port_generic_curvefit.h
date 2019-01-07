#ifndef PORT_GENERIC_CURVEFIT_H
#define PORT_GENERIC_CURVEFIT_H

#include <ceres/ceres.h>
#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>

#include "logging.h"

namespace generic {

using namespace Eigen;

// A general functor, assuming X inputs and Y outputs
template <typename _Scalar, int NX=Dynamic, int NY=Dynamic>
struct DenseFunctor
{
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    using Scalar = _Scalar;
    using InputType = Matrix<Scalar,InputsAtCompileTime, 1>;
    using ValueType = Matrix<Scalar,ValuesAtCompileTime, 1>;
    constexpr static std::string_view name = "functor";
    template<typename OStream, typename _Functor>
    friend OStream &operator<<(OStream &os, const _Functor& f) {
        return os << f.name << "<" << static_cast<int>(_Functor::InputsAtCompileTime) << ", " << static_cast<int>(_Functor::ValuesAtCompileTime) << ">(" << f.inputs() << ", " << f.values() << ")";
    }

    DenseFunctor(int inputs, int values) : m_inputs(inputs), m_values(values) {}
    DenseFunctor(): DenseFunctor(InputsAtCompileTime, ValuesAtCompileTime) {}

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }

    // operator()  // should be defined in subclasses

protected:
    int m_inputs = InputsAtCompileTime;
    int m_values = ValuesAtCompileTime;

};


// Model is a functor working on double (do we need model of other types?)
// NP -- number of input parameters
// ND -- number of demensions of input data
template <int NP=Dynamic, int ND=Dynamic>
struct Model: DenseFunctor<double, NP, Dynamic>
{
    enum {
        DimensionsAtCompileTime = ND
    };
    using _Base = DenseFunctor<double, NP, Dynamic>;
    using DataType = Matrix<double, Dynamic, Dynamic>;
    using InputDataType = Matrix<double, Dynamic, ND>;
    using InputDataBasisType = Matrix<double, Dynamic, 1>;
    constexpr static std::string_view name = "model";

    // constructors
    // via known size of params
    Model(int inputs): _Base(inputs, Dynamic), params(inputs) {}
    // via copy of params
    Model(const typename _Base::InputType& params): Model(static_cast<int>(params.size())) {this->params=params;}
    // via initializer list of params
    Model(std::initializer_list<double> params): Model(static_cast<int>(params.size())) {
        int i = 0;
        for (auto& p: params) {
            this->params(i) = p;
            ++i;
        }
    }
    Model(): Model(Model::InputsAtCompileTime) {}

    typename Model::InputType params;

    // eval()
    // should be defined to take (ndata, ndim) mesh and return vector of (ndata, 1)

    // meshgrid()
    // maybe be defined to covert input of (ndata_dim1, ... ndata_dimn) to (ndata, ndim) mesh

    // operator()
    // cound be defined to perform eval in nd shapes

    // using DataType = ValueType;
    // cound be defined to make it semantically clear for what data type it works with

    template <typename T=InputDataType>
    typename std::enable_if<ND == 2, T>::type meshgrid(const InputDataBasisType& x, const InputDataBasisType& y) const {
        // column major
        // [x0, y0,] [x1, y0] [x2, y0] ... [xn, y0] [x0, y1] ... [xn, yn]
        const long nx = x.size(), ny = y.size();
        InputDataType xy(nx * ny, 2);
        // map xx [ny, nx] to the first column of xy
        Map<MatrixXd> xx(xy.data(), ny, nx);
        // map yy [ny, nx] to the second column of xy
        Map<MatrixXd> yy(xy.data() + xy.rows(), ny, nx);
        // populate xx such that each row is x
        for (Index i = 0; i < ny; ++i) {
            xx.row(i) = x.transpose();
        }
        // populate yy such that each col is y
        for (Index j = 0; j < nx; ++j) {
            yy.col(j) = y;
        }
        return xy;
    }

    // for applying variable substitution
    typename Model::InputType transform(const typename Model::InputType& p) const {
        return p;
    }
    typename Model::InputType inverseTransform(const typename Model::InputType& p) const {
        return p;
    }

    struct Parameter {
        std::string name = "unnammed";
        bool fixed = false;
        bool bounded = false;
        eiu::Interval<double> bounds = {};
    };
    //std::vector<Parameter> param_settings{};
    // should be populated in models
};


struct Gaussian1D: Model<3, 1> // 3 params, 1 dimen
{
    constexpr static std::string_view name = "gaussian1d";

    using Model<3, 1>::Model; // model constructors
    Gaussian1D(double amplitude=1., double mean=0., double stddev=1.);

    ValueType eval(const InputType& p, const InputDataType& x) const;

    ValueType operator() (const InputType& p, const InputDataType& x) const;
    ValueType operator() (const InputDataType& x) const;
    std::vector<Parameter> param_settings{
        {"amplitude"},
        {"mean"},
        {"stddev"},
    };
};

struct Gaussian2D: Model<6, 2>  // 6 params, 2 dimen
{
    constexpr static std::string_view name = "gaussian2d";

    using Model<6, 2>::Model; // model constructors;
    Gaussian2D(double amplitude=1., double xmean=0., double ymean=0., double xstddev=1., double ystddev=1., double theta=0.);

    // operates on meshgrid xy of shape (ny * nx, 2), return a flat vector
    ValueType eval(const InputType& p, const InputDataType& xy) const;

    // operates on x and y coords separately. return a (ny, nx) matrix
    DataType operator() (
            const InputType& p,
            const InputDataBasisType& x,
            const InputDataBasisType& y) const;
    DataType operator() (
            const InputDataBasisType& x,
            const InputDataBasisType& y) const;

    InputType transform(const InputType& p) const;
    InputType inverseTransform(const InputType& p) const;

    std::vector<Parameter> param_settings{
        {"amplitude"},
        {"xmean"},
        {"ymean"},
        {"xstddev"},
        {"ystddev"},
        {"theta", false, true, {0., M_PI / 2.}},
    };
};

struct SymmetricGaussian2D: Model<4, 2>  // 4 params, 2 dimen
{
    constexpr static std::string_view name = "symmetricgaussian2d";
    using Model<4, 2>::Model; // model constructors;
    SymmetricGaussian2D(double amplitude=1., double xmean=0., double ymean=0., double stddev=1.);

    // operates on meshgrid xy of shape (ny * nx, 2), return a flat vector
    ValueType eval(const InputType& p, const InputDataType& xy) const;

    // operates on x and y coords separately. return a (ny, nx) matrix
    DataType operator() (
            const InputType& p,
            const InputDataBasisType& x,
            const InputDataBasisType& y) const;
    DataType operator() (
            const InputDataBasisType& x,
            const InputDataBasisType& y) const;
    std::vector<Parameter> param_settings{
        {"amplitude"},
        {"xmean"},
        {"ymean"},
        {"stddev"},
    };
};

// Fitter is a functor that matches the data types of the Model.
// Fitter relies on the eval() method
template <typename _Model>
struct Fitter: _Model::_Base
{
    using _Base = typename _Model::_Base;
    using Model = _Model;

    Fitter (const Model* model, int values): _Base(model->inputs(), values), m_model(model) {}
    Fitter (const Model* model): Fitter(model, Fitter::InputsAtCompileTime) {}

    const Model* model() const {return m_model;}

    const Map<const typename Model::InputDataType>* xdata = nullptr;  // set via meshgrid
    const Map<const typename Fitter::ValueType>* ydata = nullptr;
    const Map<const typename Fitter::ValueType>* sigma = nullptr;

    //int operator()(const InputType &x, ValueType& fvec) { }
    // should be defined in derived classes

private:
    const Model* m_model = nullptr;
};

// LSQ Fitter provides concreate method for least-square minimization
template <typename Model>
struct LSQFitter: Fitter<Model>
{
    using _Base = Fitter<Model>;
    using JacobianType = Matrix<typename _Base::Scalar, Dynamic, Dynamic>;
    using _Base::_Base;  // the base constructors

    int operator() (const typename LSQFitter::InputType& tp, typename LSQFitter::ValueType& fvec) const
    {
        // tp is transformed for constraint
        fvec = (this->ydata->array() - this->model()->eval(this->model()->inverseTransform(tp), *this->xdata).array()) / this->sigma->array();
        SPDLOG_TRACE("fit with xdata{}", logging::pprint(*this->xdata));
        SPDLOG_TRACE("         ydata{}", logging::pprint(*this->ydata));
        SPDLOG_TRACE("         sigma{}", logging::pprint(*this->sigma));
        SPDLOG_TRACE("residual{}", logging::pprint(fvec));
        return 0;
    }

    //int df(const InputType &x, JacobianType& fjac) { }
    // should be defined in derived classes if fitting using LevMar algorithm
    // TODO: figure out a place to store fit info
    // int info = 0;
    // int result = 0;
};

template <typename Model>
Model curvefit_eigen3(
                    const Model& model,  // y = f(x)
                    const typename Model::InputType& p,         // initial guess of model parameters
                    const typename Model::InputDataType& xdata,     // x data values, independant variable
                    const typename Model::DataType& ydata,     // y data values, measurments
                    const typename Model::DataType& sigma      // uncertainty
                    )
{
    auto logger = logging::createLogger("curvefit_eigen3", &model);
    SPDLOG_LOGGER_DEBUG(logger, "fit model {} on data{}", model, logging::pprint(xdata));

    using Fitter = LSQFitter<Model>;
    Fitter fitter(&model, ydata.size());
    Map<const typename Model::InputDataType> _x(xdata.data(), xdata.rows(), xdata.cols());
    Map<const typename Fitter::ValueType> _y(ydata.data(), ydata.size());
    Map<const typename Fitter::ValueType> _s(sigma.data(), sigma.size());
    fitter.xdata = &_x;
    fitter.ydata = &_y;
    fitter.sigma = &_s;

    using LevMarLSQ = NumericalDiff<Fitter>;
    LevMarLSQ lmlsq(fitter);

    LevenbergMarquardt<LevMarLSQ, typename Model::Scalar> lm(lmlsq);

    VectorXd pp, tp;
    SPDLOG_LOGGER_DEBUG(logger, "initial params{}", logging::pprint(p));

    tp = model.transform(p);  // use the transformed params to minimize
    int info = lm.minimize(tp);
    pp = model.inverseTransform(tp);  // get the minimized

    SPDLOG_LOGGER_DEBUG(logger, "fitted params{}", logging::pprint(pp));
    SPDLOG_LOGGER_DEBUG(logger, "info={}, nfev={}, njev={}", info, lm.nfev, lm.njev);
    SPDLOG_LOGGER_DEBUG(logger, "fvec.squaredNorm={}", lm.fvec.squaredNorm());
    return Model(pp);;
}

// ceres-solver

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::CauchyLoss;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;

// CeresAutoDiff Fitter provides concreate method for least-square minimization using ceres
template <typename Model>
struct CeresAutoDiffFitter: Fitter<Model>
{
    using _Base = Fitter<Model>;
    using _Base::_Base;  // the base constructors

    template <typename T>
    /*
    bool operator()(const T* const params, T* residual) const
    {
        Map<const typename Model::InputType> p(params, this->inputs());
        Map<typename Model::ValueType> r(residual, this->values());
        r = ((this->ydata->array() - this->model()->eval(p, *this->xdata).array()) / this->sigma->array()).eval();
        SPDLOG_LOGGER_TRACE(this->logger, "fit with xdata{}", logging::pprint(this->xdata));
        SPDLOG_LOGGER_TRACE(this->logger, "         ydata{}", logging::pprint(this->ydata));
        SPDLOG_LOGGER_TRACE(this->logger, "         sigma{}", logging::pprint(this->sigma));
        SPDLOG_LOGGER_TRACE(this->logger, "residual{}", logging::pprint(&r));
        return true;
    }
    */
    bool operator()(const T* const p, T* r) const
    {
        auto cost2 = cos(p[5]) * cos(p[5]);
        auto sint2 = sin(p[5]) * sin(p[5]);
        auto sin2t = sin(2. * p[5]);
        auto xstd2 = p[3] * p[3];
        auto ystd2 = p[4] * p[4];
        auto a = - 0.5 * ((cost2 / xstd2) + (sint2 / ystd2));
        auto b = - 0.5 * ((sin2t / xstd2) - (sin2t / ystd2));
        auto c = - 0.5 * ((sint2 / xstd2) + (cost2 / ystd2));
        for (int i=0; i < this->values(); ++i)
            r[i] =  (
                        this->ydata->coeffRef(i) -
                        p[0] * exp(
                            pow(this->xdata->coeffRef(i, 0) - p[1], 2) * a +
                            (this->xdata->coeffRef(i, 0) - p[1]) * (this->xdata->coeffRef(i, 1) - p[2]) * b +
                            pow(this->xdata->coeffRef(i, 1) - p[2], 2) * c
                            )
                    ) / this->sigma->coeffRef(i);
        return true;
    }

    //int df(const InputType &x, JacobianType& fjac) { }
    // should be defined in derived classes if fitting using LevMar algorithm
    // TODO: figure out a place to store fit info
    // int info = 0;
    // int result = 0;
    std::shared_ptr<Problem> createProblem(double* params)
    {
        std::shared_ptr<Problem> problem = std::make_shared<Problem>();
        problem->AddParameterBlock(params, this->model()->params.size());
        for (int i = 0; i < this->model()->params.size(); ++i)
        {
            typename Model::Parameter p = this->model()->param_settings.at(i);
            if (p.fixed) problem->SetParameterBlockConstant(params);
            if (p.bounded)
            {
                problem->SetParameterLowerBound(params, i, p.bounds.left());
                problem->SetParameterUpperBound(params, i, p.bounds.right());
            }
        }
        return problem;
    }
};

template <typename Model>
Model curvefit_ceres(
                    const Model& model,  // y = f(x)
                    const typename Model::InputType& p,         // initial guess of model parameters
                    const typename Model::InputDataType& xdata,     // x data values, independant variable
                    const typename Model::DataType& ydata,     // y data values, measurments
                    const typename Model::DataType& sigma      // uncertainty
                    )
{
    /*
    auto logger = spdlog::get("curvefit_ceres");
    if (!logger) logger = spdlog::stdout_color_mt("curvefit_ceres");
    logger->set_level(spdlog::level::debug);
    logger->info("fit model {} on data{}", model, logging::pprint(&xdata));
    */

    using Fitter = CeresAutoDiffFitter<Model>;
    Fitter* fitter = new Fitter(&model, ydata.size());
    Map<const typename Model::InputDataType> _x(xdata.data(), xdata.rows(), xdata.cols());
    Map<const typename Fitter::ValueType> _y(ydata.data(), ydata.size());
    Map<const typename Fitter::ValueType> _s(sigma.data(), sigma.size());
    fitter->xdata = &_x;
    fitter->ydata = &_y;
    fitter->sigma = &_s;

    CostFunction* cost_function =
        new AutoDiffCostFunction<Fitter, Fitter::ValuesAtCompileTime, Fitter::InputsAtCompileTime>(fitter, fitter->values());

    // do the fit
    /*
    logger->info("initial params{}", logging::pprint(&p));
    */
    VectorXd pp(p);
    auto problem = fitter->createProblem(pp.data());
    problem->AddResidualBlock(cost_function,
                              new CauchyLoss(0.5),
                              pp.data());

    Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    // options.minimizer_progress_to_stdout = true;
    // options.logging_type = ceres::SILENT;
    Solver::Summary summary;
    Solve(options, problem.get(), &summary);

    /*
    logger->info("{}", summary.BriefReport());
    logger->info("fitted params{}", logging::pprint(&pp));
    */
    return Model(pp);
}

}  // namespace generic
#endif /* !PORT_GENERIC_CURVEFIT_H */
