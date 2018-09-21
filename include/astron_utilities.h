#ifndef _ASTRON_UTILITIES_H_
#define _ASTRON_UTILITIES_H_

#define PI 3.1415926535897932384626433
#define TWO_PI (2.0 * PI)
#define PIO2  (PI / 2.0) 
// Degrees per radian
#define DEG_RAD (360.0 / TWO_PI)
// Seconds per UT day
#define SEC_DAY 86400.0
// One UT second in radians
#define ONESEC (TWO_PI / SEC_DAY)
// Arcseconds in 360 degrees
#define ASEC_CIRC 1296000.0
// Arcseconds in 1 radian
#define RAD_ASEC (TWO_PI / ASEC_CIRC)
// J2000 Standard Epoch
#define TJ2000 2451545.0
// One LST second in radians
#define LST_SEC (1.00273790935 * ONESEC)
// FWHM to sigma convetrsion
#define SIGMA_TO_FWHM 2.35482
#define FWHM_TO_SIGMA 1.0/SIGMA_TO_FWHM


#include "TimePlace.h"


double epoch(int imo, int iday, int year);
bool julday(int imo, int iday, int year, double * juldate, double * mjuldate);
bool juldayFromUTDate(double utDate, double* juldate, double *mjuldate);
bool caldateFromJulDate(double julDate, int* year, int* month, int* day);
void preces(double epoch, double juldate, double tdb, double ra1, double dec1,
	    double *ra2, double *dec2, int mode);
bool absToPhys(double* absRa, double* absDec, 
	       double  centerRa, double centerDec,
	       double* physRa, double* physDec,
	       int nSamples);
bool parallacticAngle(double* hourAngle, double* dec, double latitude,
		      int nSamples, double* paraAngle);
bool azElToRaDec2000(TimePlace *timePlace, double* az, double* el, 
		     double *ra, double *dec, int nSamples);
bool raDecToAzEl(TimePlace *timePlace, double* ra, double* dec, 
		 double* az, double* el, int nSamples);
void coord(double ao, double bo, double ap, double bp, 
	   double a1, double b1, double *a2, double *b2);
bool physToAbs(double* pra, double* pdec, double* cra, double* cdec,
	       double* ara, double* adec, int nSamples);

#endif
