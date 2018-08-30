#ifndef _SBSM_h_
#define _SBSM_h_

#include <suitesparse/cs.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_vector.h>

#include "nr3.h"

///SBSM (Sparse B-Spline Matrix).
///This class calculates and hold the base matrix for a B-spline of a given order,
///samples and number of breaks.

class SBSM{
	private:
		size_t order;
		size_t nSamples;
		size_t nBreaks;
		size_t nSpline;

		gsl_vector *time;

		gsl_bspline_workspace *bsw;
		cs *baseMatrix;
		cs *baseMatrix_t;
		cs *btbMatrix;

		void createBaseMatrix();
		void destroyBaseMatrix();

	public:
		cs *getBaseMatrix();
		cs *getBaseMatrix_t();
		SBSM (size_t order, size_t nSamples, size_t nBreaks);
		void resize(size_t order, size_t nSamples, size_t nBreaks);
		VecDoub fitData (double *dataVector, size_t nSamples);
		VecDoub fitData (VecDoub dataVector);
		~SBSM();
};




#endif
