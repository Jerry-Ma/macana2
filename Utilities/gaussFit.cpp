#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include "nr3.h"
#include "gaussFit.h"
#include "mpfit.h"
#include "vector_utilities.h"

//------------------------ o -------------------------------
//------------------------ o -------------------------------
//------------------------ o -------------------------------
//------------------------ o -------------------------------
//------------------------ o -------------------------------

struct vars_struct {
  double *az;
  double *el;
  double *y;
  double *sigma;
};

//user function for gaussian fit
/*
  int m     - number of data points
  int n     - number of parameters
  double *p - array of n parameters
  double *deviates - array of m deviates to be returned by myfunct()
  double **derivs - used for user-computed derivatives
*/
//int myfunct_gauss(int m, int n, double *p, double *deviates,
//                 double **derivs, void *vars)
//{
//  struct vars_struct *v = (struct vars_struct *) vars;
//  double *az = v->az;
//  double *el = v->el;
//  double *sigma = v->sigma;
//  double *y = v->y;
//
//  // p[0]  -   dc offset
//  // p[1]  -   peak of gaussian
//  // p[2]  -   azimuth sigma
//  // p[3]  -   elevation sigma
//  // p[4]  -   azimuth location of peak
//  // p[5]  -   elevation location of peak
//  // p[6]  -   position angle
//  double x0 = p[0];
//  double x1 = p[1];
//  double x2 = p[2];
//  double x3 = p[3];
//  double x4 = p[4];
//  double x5 = p[5];
//  double x6 = p[6];
//
//  double a = pow(cos(x6)/x2,2) + pow(sin(x6)/x3,2);
//  double b = -sin(2.*x6)/(2.*pow(x2,2)) + sin(2.*x6)/(2.*pow(x3,2));
//  double c = pow(sin(x6)/x2,2) + pow(cos(x6)/x3,2);
//
//  for(int i=0;i<m;i++)
//    {
//      /* Model Yi = x0 + x1 * exp(-0.5*(a(az-x4)^2 +
//	 2*b*(az-x4)*(el-x5)^2 +
//	 c*(el-x5)^2)) */
//      /* where a = (cos(x6)/x2)^2 + (sin(x6)/x3)^2,
//	 b = -sin(2*x6)/(2*x2^2) + sin(2*x6)/(2*x3^2) and
//	 c = (sin(x6)/x2)^2 + (cos(x6)/x3)^2 */
//      double expn = a*pow((az[i]-x4),2) + 2.*b*(az[i]-x4)*(el[i]-x5) +
//	c*pow((el[i]-x5),2);
//      double Yi = x0 + x1*exp(-0.5*expn);
//      deviates[i] = (Yi - y[i])/sigma[i];
//    }
//
//  return 0;
//}


int myfunct_gauss(int m, int n, double *p, double *deviates,
                 double **derivs, void *vars)
{
  struct vars_struct *v = (struct vars_struct *) vars;
  double *az = v->az;
  double *el = v->el;
  double *sigma = v->sigma;
  double *y = v->y;
  
  // p[0]  -   dc offset
  // p[1]  -   peak of gaussian
  // p[2]  -   azimuth sigma
  // p[3]  -   elevation sigma
  // p[4]  -   azimuth location of peak
  // p[5]  -   elevation location of peak
  // p[6]  -   position angle
  double x0 = p[0];
  double x1 = p[1];
  double x2 = p[2];
  double x3 = p[3];
  double x4 = p[4];
  double x5 = p[5];
  double x6 = p[6];
  
  double expn =0;
  double xp = 0;
  double yp = 0;
  double Yi = 0;

  
  for(int i=0;i<m;i++)
    {
	  xp = (az[i]-x4)*cos(x6) - (el[i]-x5)*sin(x6);
	  yp = (az[i]-x4)*sin(x6) + (el[i]-x5)*cos(x6);
	  expn = pow(xp/x2,2.0) +pow(yp/x3,2.0);
	  Yi = x0+x1*exp(-0.5*expn);
	  deviates[i]= (y[i]-Yi)/sigma[i];
    }
  
  return 0;
}



//this is the version to go with the mpfit approach
double mpGaussFit(double* params, double* params_err, double* az, double* el, 
		  double* map, double* sigma, int npts, VecInt &fixme, VecDoub &fixVals, double *iguess)
{
  //setup
  int m = npts;
  int n = 7;
  VecDoub y(m);
  for(int i=0;i<m;i++) y[i]=map[i];
  VecDoub yw(m);
  for(int i=0;i<m;i++) yw[i]=map[i]/sigma[i];
  struct vars_struct vars;
  vars.az = az;
  vars.el = el;
  vars.y = &y[0];
  vars.sigma = sigma;

  //find peak value of map for the initial guess
  int peaki=0;
  double pkw=yw[0];
  double pk=0.;
  for(int i=1;i<m;i++){
    if(yw[i] > pkw){
      pkw = yw[i];
      peaki = i;
      pk = y[i];
    }
  }

  double xall[7] = {0., pk, 6.e-5, 6.e-5, az[peaki], el[peaki], 0.};
  //define the user parameters and the initial guess
  if (iguess)
	  memcpy(xall,iguess,sizeof(double)*7);


  //get the parameter errors
  double perror[7];
  mp_result result;
  memset(&result, 0, sizeof(result));
  result.xerror=perror;

  //constrain the FWHM values to be positive
  mp_par pars[7];
  memset(&pars[0], 0, sizeof(pars));
  pars[1].limited[0] = 1;
  pars[1].limits[0] = 0.;
  pars[2].limited[0] = 1;
  pars[2].limits[0] = 0.;
  pars[3].limited[0] = 1;
  pars[3].limits[0] = 0.;

  //constrain the position angle to be between -pi/4 and pi/4
  //pars[6].limited[0] = 1;
  //pars[6].limited[1] = 1;
  //pars[6].limits[0] = -3.14159/4.;
  //pars[6].limits[1] = 3.14159/2.;
  
  //constrain the pa to be 0.
  pars[6].fixed=1;

  //constrain any other parameters as requested by fixme and fixVals
  for(unsigned int i=0;i<fixme.size();i++){
    if(fixme[i]==1){
      pars[i].fixed=1;
      xall[i] = fixVals[i];
    }
  }

  //call mpfit();
  int status;
  status = mpfit(&myfunct_gauss, m, n, xall, pars, 0, &vars, &result);

  if (status <= 0){
	  cerr<<"GausFit(): Bad fit."<<endl;
  }

  //set the outputs
  for(int i=0;i<n;i++){
    params[i] = xall[i];
    params_err[i] = result.xerror[i];
  }

  //return the chisqperdof
  return result.bestnorm/(result.nfunc-result.nfree);
}
