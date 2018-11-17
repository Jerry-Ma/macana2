#ifndef PORT_GENERIC_CURVEFIT_H
#define PORT_GENERIC_CURVEFIT_H

#include <iostream>

#include <unsupported/Eigen/NonLinearOptimization>
#include "port_generic_curvefit.h"

namespace generic {

using namespace Eigen;

// since Eigen is tempalate-based, it might be the best to write template code
// we use template-base inheritance just as Eigen does.

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

    DenseFunctor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    DenseFunctor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }

    // operator()
    // should be defined in subclasses
protected:
    int m_inputs;
    int m_values;
};

// Model is a functor working on double; TODO: do we need model of other types?
struct Model: DenseFunctor<double, Dynamic, Dynamic>
{
    using _Base = DenseFunctor<double, Dynamic, Dynamic>;

    Model();
    Model(int inputs);  // set number of parameters

    InputType& params();

    // operator()
    // should be defined in derived classes
protected:
    InputType m_params;
};

struct Gaussian1D: Model
{
    Gaussian1D(double amplitude=1., double mean=0, double stddev=1);

    Model::ValueType operator() (const typename Model::InputType& p, const typename Model::ValueType& x) const;
    Gaussian1D::ValueType operator() (const typename Model::ValueType& x) const;
};

// Fitter is a functor that matches the data types of the Model.
template <typename _Model>
struct Fitter: _Model::_Base
{
    using _Base = typename _Model::_Base;
    using Model = _Model;

    Fitter (const Model* model): _Base(), m_model(model) {};
    Fitter (const Model* model, int values): _Base(model->inputs(), values), m_model(model) {};

    const typename Fitter::ValueType* xdata = nullptr;
    const typename Fitter::ValueType* ydata = nullptr;
    const typename Fitter::ValueType* sigma = nullptr;

    const Model* model() const {return m_model;}

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
    using JacobianType = Matrix<typename _Base::Scalar, _Base::ValuesAtCompileTime, _Base::InputsAtCompileTime>;
    // using QRSolver = ColPivHouseholderQR<JacobianType>;

    using _Base::_Base;

    int operator() (const typename LSQFitter::InputType& p, typename LSQFitter::ValueType& fvec, JacobianType* _j=0) const
    {
        fvec = (this->ydata->array() - (*this->model())(p, *this->xdata).array()) / this->sigma->array();
        return 0;
    }

    //int df(const InputType &x, JacobianType& fjac) { }
    // should be defined in derived classes if fitting using LevMar algorithm
    int info = 0;
    int result = 0;
};

template <typename Model>
typename Model::InputType curvefit_eigen3(
                    const Model& model,  // y = f(x)
                    const typename Model::InputType& p,         // initial guess of model parameters
                    const typename Model::ValueType& xdata,     // x data values, independant variable
                    const typename Model::ValueType& ydata,     // y data values, measurments
                    const typename Model::ValueType& sigma      // uncertainty
                    )
{
    LSQFitter<Model> fitter(&model, xdata.size());
    fitter.xdata = &xdata;
    fitter.ydata = &ydata;
    fitter.sigma = &sigma;

    using LevMarLSQ = NumericalDiff<LSQFitter<Model>>;
    LevMarLSQ lmlsq(fitter);

    VectorXd pp(p);
    LevenbergMarquardt<LevMarLSQ, typename Model::Scalar> lm(lmlsq);

    int info = lm.minimize(pp);
    printf("info, nfev, njev : %d, %ld, %ld\n", info, lm.nfev, lm.njev);
    printf("fvec.squaredNorm() : %.13g\n", lm.fvec.squaredNorm());

    return pp;
}

}  // namespace generic
#endif /* !PORT_GENERIC_CURVEFIT_H */
