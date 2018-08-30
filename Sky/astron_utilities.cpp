#include <netcdfcpp.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
using namespace std;

#include "nr3.h"
#include "astron_utilities.h"
#include "TimePlace.h"
#include "vector_utilities.h"
#include "eph_manager.h"
#include "novas.h"


double epoch(int imo, int iday, int year)
{
  const int mday[12] = {0,31,59,90,120,151,181,212,243,273,304,334};
  int leap, ifrac, leaps;

  if((year < 1900) || (year > 2100)) {
    return(0.0);
  }

  leap = 0;
  if(year == (4 * (year / 4)))leap = 1;

  ifrac = year / 100;
  ifrac = year - 100 * ifrac;
  leaps = (ifrac - 1) / 4;

  double epoch;
  epoch = (double)year + ((double)(ifrac + leaps + mday[imo-1] + iday))/365.0;
  if(imo > 2) epoch += (double)leap/365.;

  return epoch;
}


bool julday(int imo, int iday, int year, double * juldate, double * mjuldate)

/********************************************************************
C
C THIS ROUTINE CALCULATES THE JULIAN DATE AND MODIFIED JULIAN
C  DATE FOR THE ENTERED MONTH,DAY,YEAR. THE FORMULA USED IS
C  BASED ON THE EXPRESSIONS CONTAINED ON PAGE K2 OF
C  "THE ASTRONOMICAL ALMANAC", 1983 EDITION.
C
C ENTER WITH:
C
C     IMO        (1 -> 12) INTEGER : MONTH
C     IDAY       (1 -> 31) INTEGER : DAY OF THE MONTH
C     YEAR      (1900-2100)INTEGER : YEAR
C
C LEAVE WITH:
C
C     JULDATE     DOUBLE PRECISION : JULIAN DATE AT 0 HOURS UT
C     MJULDATE    DOUBLE PRECISION : MODIFIED JULIAN DATE
C                                    (JD - 2 400 000.5)
C
C*******************************************************************/
{
  const int mday[12] = {0,31,59,90,120,151,181,212,243,273,304,334};
  double base;
  int leap, ifrac, leaps;

  if((year < 1900) || (year > 2100)) {
    *juldate = 0.0;
    *mjuldate = 0.0;
    return(false);
  }

  base = 2415019.5;
  if(year >= 2000)base = 2451544.5;
  if(year == 2000)base -= 1.0;

  leap = 0;
  if(year == (4 * (year / 4)))leap = 1;

  ifrac = year / 100;
  ifrac = year - 100 * ifrac;
  leaps = (ifrac - 1) / 4;

  *juldate = base + 365.0 * ( double) ifrac + (double) leaps + (double) mday[imo-1] 
           + (double) iday;
  if(imo > 2)*juldate += leap;

  *mjuldate = *juldate - 2400000.5;
  return(true);
}



//----------------------------- o ---------------------------------------



bool juldayFromUTDate(double utDate, double* juldate, double *mjuldate)
{
  double base;
  int leap, ifrac, leaps;

  int year = floor(utDate);
  if((year < 1900) || (year > 2100)) {
    *juldate = 0.0;
    *mjuldate = 0.0;
    return(false);
  }

  base = 2415019.5;
  if(year >= 2000)base = 2451544.5;
  if(year == 2000)base -= 1.0;

  leap = 0;
  if(year == (4 * (year / 4)))leap = 1;
  ifrac = year / 100;
  ifrac = year - 100 * ifrac;
  leaps = (ifrac - 1) / 4;

  //kamal says don't include hour minute second
  double days_in_year = 365.0;
  int days = (int)((utDate-year)*days_in_year)+1;

  *juldate = base + 365.0 * ( double) ifrac + (double) leaps + (double) days;
  if(days > 59)*juldate += leap;

  *mjuldate = *juldate - 2400000.5;
  return(true);
}


//----------------------------- o ---------------------------------------


//uses novas to return calendar date
bool caldateFromJulDate(double julDate, int* year, int* month, int* day)
{
  short int y, m, d;
  double h;
  cal_date(julDate, &y, &m, &d, &h);
  *year = (int) y;
  *month = (int) m;
  *day = (int) d;

  return 1;
}


//----------------------------- o ---------------------------------------


void preces(double epoch, double juldate, double tdb, double ra1, double dec1,
	    double *ra2, double *dec2, int mode)

/*********************************************************************
C
C  INPUT:  EPOCH     INITIAL JULIAN EPOCH, E.G. 2000.00  (DP)
C          JULDATE   FINAL DATE
C          TDB       BARYCENTRIC DYNAMICAL TIME ON FINAL DAY, IN RADIANS  (DP)
C          RA1,DEC1  INPUT RA AND DEC, IN RADIANS  (DP)
C
C          MODE      -1 = CALCULATE MEAN PLACE AT EPOCH FROM MEAN PLACE
C                         AT JULDATE
C                    +1 = CALCULATE MEAN PLACE AT JULDATE FROM MEAN PLACE
C                         AT EPOCH
C
C  OUTPUT: RA2,DEC2  OUTPUT RA AND DEC, IN RADIANS (DP)
C
C  CORRECTIONS ARE MADE FOR PRECESSION USING J2000 THEORY
C ********************************************************************/

{
  static double tlast = 1.0;
  static double t0last = 1.0;
  static double zeta0, theta, z, ct, st;
  double t0, t1, t2, t, sr0, cr0, sd0, cd0, cdsr1, cdcr1, sd1;

#if defined (_OPENMP)
#pragma omp threadprivate (tlast,t0last, zeta0, theta, z, ct, st)
#endif

  //  T0 =  TIME FROM 2000.00 TO JULIAN EPOCH IN JULIAN CENTURIES.
  t0 = (epoch - 2000.0) / 100.0;

  //  T1 =  TIME FROM J2000.00 TO JULIAN DATE IN JULIAN CENTURIES.
  t1 = juldate - TJ2000;
  t1 = (t1 + tdb / TWO_PI) / 36525.0;

  if(mode == -1) {

    // Reverse definitions of t0 and t1
    t2 = t0;
    t0 = t1;
    t1 = t2;
  }

  // T = TIME BETWEEN EPOCH AND JULDATE IN JULIAN CENTURIES.
  t = t1 - t0;

  // Only need to do this about once per hour
  if((fabs(t - tlast) > 1.0e-6) || (fabs(t0 - t0last) > 1.0e-6)) {

    zeta0 = t * (2306.2181 + (1.39656 - 0.000139 * t0) * t0
          + t * ((0.30188 - 0.000344 * t0) + t * 0.017998));

    z     = t * (2306.2181 + (1.39656 - 0.000139 * t0) * t0
          + t * ((1.09468 + 0.000066 * t0) + t * 0.018203));

    theta = t * (2004.3109 - (0.85330 + 0.000217 * t0) * t0
          - t * ((0.42665 + 0.000217 * t0) + t * 0.041833));

    zeta0 *= RAD_ASEC;
    z     *= RAD_ASEC;
    theta *= RAD_ASEC;

    ct = cos(theta);
    st = sin(theta);

    tlast = t;
    t0last = t0;
  }

  sr0 = sin(ra1 + zeta0);
  cr0 = cos(ra1 + zeta0);
  sd0 = sin(dec1);
  cd0 = cos(dec1);

  cdsr1 = cd0 * sr0;
  cdcr1 = ct * cd0 * cr0 - st * sd0;
  sd1   = ct * sd0 + st * cd0 * cr0;

  // Precess the coordinates
  *dec2 = asin(sd1);
  *ra2  = atan2(cdsr1, cdcr1) + z;
  if(*ra2 < 0.0)*ra2 += TWO_PI;
}


//----------------------------- o ---------------------------------------


//this is a translation of aztec_abs2phys.pro in the aztec_idl_utilities
//    Takes absolute ra/dec coordinates and a specified tangent point
//    and returns physical coordinates, relative to the tangent point,
//    using tangential (gnomic) projection.
bool absToPhys(double* absRa, double* absDec, 
	       double  centerRa, double centerDec,
	       double* physRa, double* physDec,
	       int nSamples)
{
  //use temporary storage to avoid writing over absRa and absDec
  VecDoub tRa(nSamples);
  VecDoub tDec(nSamples);
  for(int i=0;i<nSamples;i++) tDec[i]=absDec[i];

  //tRa must range from -PI to PI
  for(int i=0;i<nSamples;i++) 
    tRa[i] = (absRa[i] > PI) ? absRa[i]-TWO_PI : absRa[i];

  //same thing for centerRa
  centerRa = (centerRa > PI) ? centerRa-TWO_PI : centerRa;

  //the tangential projection
  VecDoub cosc(nSamples);
  double sCD = sin(centerDec);
  double cCD = cos(centerDec);
  for(int i=0;i<nSamples;i++)
    cosc[i] = sCD*sin(absDec[i]) + cCD*cos(absDec[i])*cos(tRa[i]-centerRa);

  for(int i=0;i<nSamples;i++)
    if(cosc[i]==0.){
      physRa[i] = 0.;
      physDec[i] = 0;
    } else {
      physRa[i] = cos(absDec[i])*sin(tRa[i]-centerRa)/cosc[i];
      physDec[i] = (cCD*sin(absDec[i]) - sCD*cos(absDec[i])*cos(tRa[i]-centerRa))/cosc[i];
    }

   return 1;
}


//----------------------------- o ---------------------------------------


//this is a translation of aztec_phys2abs.pro in the aztec_idl_utilities
//    Takes physical coordinates and a specified tangent point
//    and returns absolute coordinates using tangential (gnomic) projection.
bool physToAbs(double* pra, double* pdec, double* cra, double* cdec,
	       double* ara, double* adec, int nSamples)
{
  for(int i=0;i<nSamples;i++){
    double rho = sqrt(pow(pra[i],2)+pow(pdec[i],2));
    double c = atan(rho);
    if(c == 0.){
      ara[i] = cra[i];
      adec[i] = cdec[i];
    } else {
      double ccwhn0 = cos(c);
      double scwhn0 = sin(c);
      double ccdec = cos(cdec[i]);
      double scdec = sin(cdec[i]);
      double a1;
      double a2;
      a1 = ccwhn0*scdec + pdec[i]*scwhn0*ccdec/rho;
      adec[i] = asin(a1);
      a2 = pra[i]*scwhn0/(rho*ccdec*ccwhn0 - pdec[i]*scdec*scwhn0);
      ara[i] = cra[i] + atan(a2);
    }
  }
    return 1;
}


//----------------------------- o ---------------------------------------

//this is cribbed from Sky.cc but I've removed the nutation and aberation
//correction and replaced it with one from Novas
bool azElToRaDec2000(TimePlace *timePlace,
		     double* az, double* el, 
		     double *ra, double *dec,
		     int nSamples)
{


  //time and place
  double gLat = timePlace->latitude;

  //coordinate transformation
  VecDoub haAct(nSamples);
  VecDoub decAct(nSamples);
  for(int i=0;i<nSamples;i++)
    coord(PI, PIO2 - gLat, 0.0, gLat, az[i], el[i], &haAct[i], &decAct[i]);

  VecDoub raAct(nSamples);
  for(int i=0;i<nSamples;i++){
    raAct[i] = timePlace->lst[i] - haAct[i];
    if(raAct[i] < 0.0) raAct[i] += TWO_PI;
  }
  /*
  //open the ephemeris calculation for NOVAS
  double jd_beg, jd_end;
  short int de_num=0, error;
  if ((error = ephem_open ("Sky/Novas/JPLEPH", &jd_beg,&jd_end,&de_num)) != 0)
    {
      if (error == 1){
	cerr << "JPL ephemeris file not found." << endl;
	exit(1);
      } else {
	cerr << "Error reading JPL ephemeris file header." << endl;
	exit(1);
      }
    }
  else
    {
      cerr << "JPL ephemeris DY" << de_num << " open." << endl;
    }
  */

  const int leap_secs = 33;

  //do the precession and aberration with novas
  for(int i=0;i<nSamples;i++){
    int ms;
    double ra1=raAct[i]*DEG_RAD/360.*24.;
    double dec1=decAct[i]*DEG_RAD;
    double jd_utc, jd_tt;
//#pragma omp critical(novas)
//{
    jd_utc = julian_date(timePlace->year, timePlace->month, 
			 timePlace->day, 
			 timePlace->detUtc[i]+timePlace->timeOffset);
    jd_tt = jd_utc + ((double)leap_secs + 32.184) / 86400.0;

    ms = mean_star(jd_tt, ra1, dec1, 1, &ra[i], &dec[i]);
//}
    if(ms){
      cerr << "azElToRaDec2000(): Problem with mean_star(): ms=" << ms << endl;
      exit(1);
    }
    ra[i] = ra[i]/24.*360./DEG_RAD;
    dec[i] = dec[i]/DEG_RAD;
  } 

  return 1;
}


//----------------------------- o ---------------------------------------

//ra/dec to Az/el using Novas
bool raDecToAzEl(TimePlace *timePlace,
		 double* ra, double* dec, double* az, double* el, int nSamples)
{
//#pragma omp critical (novas)
//{
  //open the ephemeris calculation for NOVAS
  double jd_beg, jd_end;
  short int de_num=0, error;
  string jplpath;
  jplpath.assign("Sky/Novas/JPLEPH");
  char *macanaPath = getenv("AZTEC_MACANA_PATH");
  if (macanaPath)
	  jplpath = string(macanaPath) +"/" + jplpath;


  if ((error = ephem_open (jplpath.c_str(), &jd_beg,&jd_end,&de_num)) != 0)
    {
      if (error == 1){
	cerr << "JPL ephemeris file not found." << endl;
	exit(1);
      } else {
	cerr << "Error reading JPL ephemeris file header." << endl;
	exit(1);
      }
    }
  else
    {
      cerr << "JPL ephemeris DE" << de_num << " open." << endl;
    }


  //set up the on_surface structure
  on_surface geo_loc;
  make_on_surface (timePlace->latitude * DEG_RAD,
		   timePlace->longitude * DEG_RAD,
		   timePlace->elevation,
		   0., 640., &geo_loc);

  //here's the timing
  double jd_utc, jd_tt, jd_ut1, delta_t;
  //these valuse only very accurate for 2008
  const double ut1_utc = -0.387845;
  const int leap_secs = 33;
  
  //storage space for the transformation
  cat_entry star;
  make_cat_entry ("dummy","FK6", 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &star); 

  //now cycle through the ra/dec values
  for(int i=0;i<nSamples;i++){
    //timing
    jd_utc = julian_date(timePlace->year, timePlace->month, 
			 timePlace->day, 
			 timePlace->detUtc[i]+timePlace->timeOffset);
    jd_tt = jd_utc + ((double)leap_secs + 32.184) / 86400.0;
    //the following line has a correction designed to match LMT pointings
    jd_ut1 = jd_utc + ut1_utc / 86400.0 + 0.175/86400.0 - 0.12118220/86400.0;
    delta_t = 32.184 + leap_secs - ut1_utc;

    //set the ra and dec
    star.ra = ra[i]/TWO_PI*24.;
    star.dec = dec[i]*DEG_RAD;

    //topocentric place of star
    short int error;
    double rat, dect;
    if ((error = topo_star (jd_tt,delta_t,&star,&geo_loc, 1, 
			    &rat,&dect)) != 0){
      cerr << "Error %d from topo_star: " << error << endl;
      exit(1);
    }

    //azimuth and zenith distance with NO REFRACTION
    double zd, rar, decr, azd;
    equ2hor (jd_ut1,delta_t,1,0.0,0.0,&geo_loc,rat,dect,0,
	     &zd,&azd,&rar,&decr);

    //finally the azimuth and elvation in radians;
    az[i] = azd / DEG_RAD;
    el[i] = (90.-zd) / DEG_RAD;
  }
//}
  return 1;
}



//----------------------------- o ---------------------------------------

void coord(double ao, double bo, double ap, double bp, 
	   double a1, double b1, double *a2, double *b2)

/*********************************************************************
C    COORDINATE ROTATION USING SPHERICAL TRIANGLES
C
C    AO,BO = INPUT:  ORIGIN OF A2 IN A1,B1 COORDINATE SYSTEM
C    AP,BP = INPUT:  POLE OF B2 IN A1,B1 COORDINATE SYSTEM
C    A1,B1 = INPUT:  COORDINATES IN A1,B1 COORDINATE SYSTEM
C    A2,B2 = OUTPUT: COORDINATES IN A2,B2 COORDINATE SYSTEM
C
C *******************************************************************/

{
  double sbo, cbo, sbp, cbp, sb1, cb1, sb2, cb2, saa, caa, sbb, cbb;
  double sa2, ca2, ta2o2;

  sbo = sin(bo);
  cbo = cos(bo);
  sbp = sin(bp);
  cbp = cos(bp);
  sb1 = sin(b1);
  cb1 = cos(b1);

  sb2 = sbp * sb1 + cbp * cb1 * cos(ap - a1);
  cb2 = sqrt((1.0 - sb2) * (1.0 + sb2));

  *b2 = atan2(sb2, cb2);

  saa = sin(ap - a1) * cb1 / cb2;
  caa = (sb1 - sb2 * sbp) / (cb2 * cbp);

  sbb = sin(ap - ao) * cbo;
  cbb = sbo / cbp;

  sa2 = saa * cbb - caa * sbb;
  ca2 = caa * cbb + saa * sbb;

  if(ca2 > 0.0)ta2o2 = sa2 / (1.0 + ca2);
  else ta2o2 = (1.0 - ca2) / sa2;

  *a2 = 2.0 * atan(ta2o2);
  if(*a2 < 0.0)*a2 += TWO_PI;
}


//----------------------------- o ---------------------------------------



