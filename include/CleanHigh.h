#ifndef _CLEAN_HIGH_H
#define _CLEAN_HIGH_H


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "Clean.h"

class CleanHigh : public Clean{
    private:
        int order;			//Set 1 for average, 2 linear, 3 cuadratic template calculation
        MatDoub Trp (MatDoub M) ;
        MatDoub MatrixProduct (MatDoub A, MatDoub B) ;
        MatDoub IdMat (size_t size) ;
        MatDoub CorMat (MatDoub d);
        MatDoub CovMat (MatDoub d);
    
    protected:
        MatDoub Coefficients (MatDoub Template, MatDoub tods) ;
        MatDoub IterTemp (MatDoub coef , MatDoub tods ) ;

    public:
        bool clean();
        CleanHigh(Array*  dataArray, Telescope *telescope);
        void fullMedianSustraction ();
        void fullLinearCorrection ();
        MatDoub Average (MatDoub tods, MatDoub C);
        MatDoub Planar (MatDoub tods, MatDoub C, VecDoub x, VecDoub y) ;
        MatDoub Quadratic (MatDoub tods, MatDoub C, VecDoub x, VecDoub y) ;
        MatDoub BuildTemplate (gsl_matrix *tods, gsl_matrix *C, gsl_matrix *S);
        MatDoub subTemplate (MatDoub tods, VecDoub x, VecDoub y, double crit[2],
                             bool withCov, string method);
        ~CleanHigh();



};

#endif
