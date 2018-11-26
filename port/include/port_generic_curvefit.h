#ifndef PORT_GENERIC_CURVEFIT_H
#define PORT_GENERIC_CURVEFIT_H

#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
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
    using InputType = Matrix<Scalar,InputsAtCompileTime,1>;
    using ValueType = Matrix<Scalar,ValuesAtCompileTime,1>;

    DenseFunctor(int inputs, int values) : m_inputs(inputs), m_values(values) {
        setupLogger("functor");
    }
    // default
    DenseFunctor(): DenseFunctor(InputsAtCompileTime, ValuesAtCompileTime) {}

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }

    // operator()
    // should be defined in subclasses
    /*
    template<typename OStream, typename Functor>
    friend OStream &operator<<(OStream &os, const Functor& f)
    {
        return os << f.logger->name();
    }
    */
protected:
    int m_inputs = InputsAtCompileTime;
    int m_values = ValuesAtCompileTime;

    void setupLogger(const std::string_view& name)
    {
        logger = std::make_shared<spdlog::logger>(
                    fmt::format("{}@{:x}", name, reinterpret_cast<std::uintptr_t>(this)),
                    logging::console);
        logger->set_level(spdlog::level::debug);
    }
    std::shared_ptr<spdlog::logger> logger;
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
    // via known size of params
    Model(int inputs): _Base(inputs, Dynamic), params(inputs) {
        Model::setupLogger(Model::name);
    }
    // via copy of params
    Model(const typename _Base::InputType& params): Model(static_cast<int>(params.size())) {this->params=params;}
    // via initializer list of params
    Model(std::initializer_list<double> params): Model(static_cast<int>(params.size()))
    {
        int i = 0;
        for (auto& p: params) {
            this->params(i) = p;
            ++i;
        }
    }
    // default
    Model(): Model(Model::InputsAtCompileTime) {}

    /*
    // via varadic template for cleaner syntax
    template <typename... Ts>
    Model(Ts... xs) : _Base(sizeof...(Ts),  Dynamic), m_params(sizeof...(Ts))
    {
        // only accept double as required by the Model class
        static_assert((std::is_same<Ts, double>::value && ...));
        // populate m_params
        for (int i = 0; const auto& x: {xs...}) {m_params[i++] = x;
    }
    */
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
    typename std::enable_if<ND == 2, T>::type meshgrid(const InputDataBasisType& x, const InputDataBasisType& y) const
    {
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

    template<typename OStream, typename _Model>
    friend OStream &operator<<(OStream &os, const _Model& m)
    {
        return os << m.name << "(NP=" << static_cast<int>(_Model::InputsAtCompileTime) << ",ND=" << static_cast<int>(_Model::DimensionsAtCompileTime) << ")";
    }
};


struct Gaussian1D: Model<3, 1> // 3 params, 1 dimen
{
    constexpr static std::string_view name = "gaussian1d";
    using Model<3, 1>::Model; // model constructors
    Gaussian1D(double amplitude=1., double mean=0., double stddev=1.);

    ValueType eval(const InputType& p, const InputDataType& x) const;
    // no meshgrid needed here

    // convinient methods
    ValueType operator() (const InputType& p, const InputDataType& x) const;
    ValueType operator() (const InputDataType& x) const;
};

struct Gaussian2D: Model<6, 2>  // 6 params, 2 dimen
{
    constexpr static std::string_view name = "gaussian2d";
    using Model<6, 2>::Model; // model constructors;
    Gaussian2D(double amplitude=1., double xmean=0., double ymean=0., double xstddev=1., double ystddev=1., double theta=0.);
    ~Gaussian2D(){}

    // operates on meshgrid xy of shape (ny * nx, 2), return a flat vector
    ValueType eval(const InputType& p, const InputDataType& xy) const;

    // convinient methods
    // operates on x and y coords separately. return a (ny, nx) matrix
    DataType operator() (
            const InputType& p,
            const InputDataBasisType& x,
            const InputDataBasisType& y) const;
    DataType operator() (
            const InputDataBasisType& x,
            const InputDataBasisType& y) const;
};

struct SymmetricGaussian2D: Model<4, 2>  // 4 params, 2 dimen
{
    constexpr static std::string_view name = "symmetricgaussian2d";
    using Model<4, 2>::Model; // model constructors;
    SymmetricGaussian2D(double amplitude=1., double xmean=0., double ymean=0., double stddev=1.);
    ~SymmetricGaussian2D(){}

    // operates on meshgrid xy of shape (ny * nx, 2), return a flat vector
    ValueType eval(const InputType& p, const InputDataType& xy) const;

    // convinient methods
    // operates on x and y coords separately. return a (ny, nx) matrix
    DataType operator() (
            const InputType& p,
            const InputDataBasisType& x,
            const InputDataBasisType& y) const;
    DataType operator() (
            const InputDataBasisType& x,
            const InputDataBasisType& y) const;
};

// Fitter is a functor that matches the data types of the Model.
// Fitter relies on the eval() method
template <typename _Model>
struct Fitter: _Model::_Base
{
    using _Base = typename _Model::_Base;
    using Model = _Model;

    Fitter (const Model* model, int values): _Base(model->inputs(), values), m_model(model) {
        Fitter::setupLogger("fitter");
    }
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
    // using QRSolver = ColPivHouseholderQR<JacobianType>;

    using _Base::_Base;  // the base constructors

    int operator() (const typename LSQFitter::InputType& p, typename LSQFitter::ValueType& fvec, [[maybe_unused]] JacobianType* _j=0) const
    {
        fvec = (this->ydata->array() - this->model()->eval(p, *this->xdata).array()) / this->sigma->array();
        this->logger->debug("fit with xdata{}", logging::pprint(this->xdata));
        this->logger->debug("fit with ydata{}", logging::pprint(this->ydata));
        this->logger->debug("fit with sigma{}", logging::pprint(this->sigma));
        this->logger->debug("residual{}", logging::pprint(&fvec));
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
    auto logger = spdlog::get("curvefit");
    if (!logger) logger = spdlog::stdout_color_mt("curvefit");
    logger->set_level(spdlog::level::debug);
    logger->info("fit model {} on data{}", model, logging::pprint(&xdata));

    using Fitter = LSQFitter<Model>;
    Fitter fitter(&model, xdata.size());
    Map<const typename Model::InputDataType> _x(xdata.data(), xdata.rows(), xdata.cols());
    Map<const typename Fitter::ValueType> _y(ydata.data(), ydata.size());
    Map<const typename Fitter::ValueType> _s(sigma.data(), sigma.size());
    fitter.xdata = &_x;
    fitter.ydata = &_y;
    fitter.sigma = &_s;

    using LevMarLSQ = NumericalDiff<Fitter>;
    LevMarLSQ lmlsq(fitter);

    LevenbergMarquardt<LevMarLSQ, typename Model::Scalar> lm(lmlsq);

    VectorXd pp(p);
    logger->info("initial params{}", logging::pprint(&p));

    int info = lm.minimize(pp);
    logger->debug("fitted params{}", logging::pprint(&pp));
    logger->info("info={}, nfev={}, njev={}", info, lm.nfev, lm.njev);
    logger->info("fvec.squaredNorm={}", lm.fvec.squaredNorm());
    return Model(pp);;
}

}  // namespace generic
#endif /* !PORT_GENERIC_CURVEFIT_H */
