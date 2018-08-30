#ifndef _CLEAN2DSTRIPE_H_
#define _CLEAN2DSTRIPE_H_

#include <suitesparse/cs.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_spline.h>

#include "Map.h"

class Clean2dStripe{
	private:
		Map *map;

	public:
		Clean2dStripe(Map *map);
		bool removeStripe(int cellSize, double coverage);
		bool fastRemoveStripe(int cellSize, double coverage);
};



class Spline2d{
	//Construnct a two dimentional spline interpolant based on gsl spline funcions
	private:
		double *x;
		double *y;
		double **map;
		size_t nx;
		size_t ny;
		size_t px;
		size_t py;
		int order;
		size_t downsample;
		double dx;
		double dy;

		gsl_bspline_workspace *bswx;

		double **coeff;
		gsl_bspline_workspace *bswy;
		size_t nSplinex;
		size_t nSpliney;
		cs *baseMatrix_x;
		cs *baseMatrix_t_x;
		cs *baseMatrix_y;
		cs *baseMatrix_t_y;

		void createBaseMatrix(gsl_bspline_workspace *bsw, double *time,cs **Matrix, cs **Matrix_t, size_t nSpline, size_t nData);
		void create1DCoeffs(double *data, size_t nData, cs *bMatrix, cs *btbMatrix, size_t nSpline, double *coeffOut);
		void createCoeffs();

	public:

		Spline2d(double *x, size_t nx, double *y, size_t ny, double **map, size_t px, size_t py);
		Spline2d(double *x, size_t nx, double *y, size_t ny, double **map, size_t px, size_t py, size_t downsample);
		~Spline2d();
		void interpolate(double *xi, double *yi, size_t nData, double *dataOut);
		double interpolate_single(double xi, double yi);
};


class Spline2dInterp{
	private:
		double *x;
		double *y;
		double **map;
		size_t npx;
		size_t npy;
		gsl_spline **splineArray;
		gsl_interp_accel **accArray;

	public:
		Spline2dInterp(double *x, size_t nx, double *y, size_t ny, double **map, size_t px, size_t py);
		void interpolate(double *xi, double *yi, size_t nData, double *dataOut);
		double interpolate_single(double xi, double yi);
		~Spline2dInterp();
};

#endif
