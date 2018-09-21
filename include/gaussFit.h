#ifndef _GAUSSFIT_H_
#define _GAUSSFIT_H_

extern "C" {
  int myfunct_gauss(int m, int n, double *p, double *deviates,
		    double **derivs, void *vars);
}

double mpGaussFit(double* params, double* params_err, double* az, double* el, 
		  double* map, double* sigma, int npts, VecInt &fixme, VecDoub &fixVals, double *iguess=NULL);

#endif
