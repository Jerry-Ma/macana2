#include <cmath>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit_nlin.h>
#include <suitesparse/cs.h>
#include <gsl/gsl_spline.h>

#include <unistd.h>

#include "Clean2dStripe.h"
#include "Map.h"
#include "vector_utilities.h"

using namespace std;

Clean2dStripe::Clean2dStripe(Map *cmap){
	this->map = cmap;
}


bool Clean2dStripe::fastRemoveStripe(int cellSize, double coverage){
	if (cellSize <= 1){
			cerr<<"Clean2dStripe::removeStripe(): Error. Cell size is less than 1. Skipping Map cleaning"<<endl;
			return false;
		}

		if (coverage <= 0.0){
			coverage = 0.0;
		}

		if (coverage >= 0.7){
			cerr<<"Clean2dStripe::removeStripe(): Error. Cell coverage larger than 0.7. Clipping to 0.7"<<endl;
			coverage = 0.7;
		}
		size_t nx = this->map->getNrows();
		size_t ny = this->map->getNcols();
		size_t ntot = nx * ny;

		double **map2fit = new double *[nx];
		double *weights = new double  [ntot];
		double *xfit = new double [ntot];
		double *x2fit = new double [nx];
		double *yfit = new double [ntot];
		double *fit = new double [ntot];



		for (size_t ix = 0; ix < nx; ix++){
			x2fit[ix] = this->map->getRowCoordsPhys(ix);
			map2fit[ix]= new double [ny];
			for (size_t iy = 0; iy < ny; iy++){
				map2fit[ix][iy] = this->map->image[ix][iy];
				weights[ix*ny+iy] = this->map->weight[ix][iy];
				xfit [ix*ny+iy] = this->map->getRowCoordsPhys(ix);
				yfit [ix*ny+iy] = this->map->getColCoordsPhys(iy);
			}
		}

		double pcut = 1e100;
		size_t pixelCut=0;
		if (coverage > 0.0){
			pcut = coverage*percentile(weights, nx*ny,0.95);
			for (size_t ix=0; ix <nx; ix++)
				for (size_t iy=0; iy<ny; iy++){
					if (weights[ix*ny+iy]<pcut){
						map2fit[ix][iy] = 0.0;
						pixelCut++;
					}
				}
			cerr<<"Cleand2dStripe(): Cut threshold set on "<<pcut<< " number of pixel set to zero: "<<pixelCut<<endl;
		}
		size_t tCell= size_t(cellSize);
		Spline2d *myspline2d  = new Spline2d(x2fit,nx,yfit,ny,map2fit,tCell, tCell,4);
		//Spline2dInterp *myspline2d  = new Spline2dInterp(x2fit,nx,yfit,ny,map2fit,tCell, tCell);
		cout<<"Cleand2dStripe(): 2D-spline interpolator created"<<endl;

		myspline2d->interpolate(xfit,yfit,ntot,fit);

		for (size_t ix = 0; ix < nx; ix++)
			for (size_t iy = 0; iy < ny; iy++)
				this->map->image[ix][iy]-=fit[ix*ny+iy];

		delete myspline2d;

		for (size_t ix=0; ix <nx; ix++)
			delete [] map2fit[ix];
		delete [] map2fit;
		delete [] fit;
		delete [] x2fit;
		delete [] yfit;
		delete [] xfit;
		delete [] weights;
		return true;
}


bool Clean2dStripe::removeStripe(int cellSize, double coverage){

	bool debug = false;

	if (cellSize <= 1){
		cerr<<"Clean2dStripe::removeStripe(): Error. Cell size is less than 1. Skipping Map cleaning"<<endl;
		return false;
	}

	if (coverage <= 0.0){
		coverage = 0.0;
	}

	if (coverage >= 0.7){
		cerr<<"Clean2dStripe::removeStripe(): Error. Cell coverage larger than 0.7. Clipping to 0.7"<<endl;
		coverage = 0.7;
	}
	size_t nx = this->map->getNrows();
	size_t ny = this->map->getNcols();

	double **map2fit = new double * [nx];
	//double mapOut [2+nx*ny];
	double *fit  = new double [2+nx*ny];
	if (!fit){
		cerr<<"Fatal. Not enough memory to compute 2D stripe removal"<<endl;
				exit(-1);
	}
	double *weights  = new double [nx*ny];
	if (!weights){
		cerr<<"Fatal. Not enough memory to compute 2D stripe removal"<<endl;
						exit(-1);
	}
	double *x2fit = new double [nx];
	if (x2fit ==NULL){
		cerr<<"Fatal. Not enough memory to compute 2D stripe removal"<<endl;
		exit(-1);
	}
	double *y2fit = new double [ny];
	if (y2fit == NULL){
		cerr<<"Fatal. Not enough memory to compute 2D stripe removal"<<endl;
		exit(-1);
	}

	//Copy values
	for (size_t ix=0; ix <nx; ix++)
		x2fit[ix] = this->map->getRowCoordsPhys(ix);
	for (size_t iy=0; iy<ny; iy++)
		y2fit[iy] = this->map->getColCoordsPhys(iy);
	for (size_t ix=0; ix <nx; ix++){
		map2fit[ix]= new double [ny];
		for (size_t iy=0; iy<ny; iy++){
			map2fit[ix][iy] = this->map->image[ix][iy];
			weights[ix*ny+iy] = this->map->weight[ix][iy];
		}
	}

	double pcut = 1e100;
	size_t pixelCut=0;
	if (coverage > 0.0){
		pcut = coverage*percentile(weights, nx*ny,0.95);
		for (size_t ix=0; ix <nx; ix++)
				for (size_t iy=0; iy<ny; iy++){
					if (weights[ix*ny+iy]<pcut){
						map2fit[ix][iy] = 0.0;
						pixelCut++;

					}
					//mapOut[2+ix*ny+iy]= map2fit[ix][iy];
				}
	}
//	fit[0]=mapOut[0]= nx;
//	fit[1]=mapOut[1]=ny;
//
//	writeVecOut("Map2fit_fs.txt", mapOut, 2+nx*ny);
	cerr<<"Cleand2dStripe(): Cut threshold set on "<<pcut<< " number of pixel set to zero: "<<pixelCut<<endl;

	size_t tCell= size_t(cellSize);
	//Spline2dInterp *myspline2d  = new Spline2dInterp(x2fit,nx,y2fit,ny,map2fit,tCell, tCell);
	Spline2d *myspline2d  = new Spline2d(x2fit,nx,y2fit,ny,map2fit,tCell, tCell,4);
	cout<<"Cleand2dStripe(): 2D-spline interpolator created"<<endl;
	for (size_t ix=0; ix <nx; ix++){
		for (size_t iy=0; iy<ny; iy++){
			fit[2+ix*ny+iy]= myspline2d->interpolate_single(this->map->getRowCoordsPhys(ix),this->map->getColCoordsPhys(iy));
			this->map->image[ix][iy]-=fit[2+ix*ny+iy];
		}
	}



	if (debug)
		writeVecOut("Fit_fs.txt", fit, 2+nx*ny);
	delete myspline2d;

	for (size_t ix=0; ix <nx; ix++)
		delete [] map2fit[ix];
	delete [] map2fit;
	delete [] fit;
	delete [] x2fit;
	delete [] y2fit;
	delete [] weights;
	return true;
}


Spline2d::Spline2d(double *x, size_t nx, double *y, size_t ny, double **map, size_t px, size_t py){
	this->x = x;
	this->y = y;
	this->map = map;
	this->nx = nx;
	this->ny = ny;

	this->order = 4;
	this->px = px;
	this->py = py;

	this->downsample = 1;
	this->dx = 1;
	this->dy = 1;

	this->nSplinex = this->px +order -2;
	this->nSpliney = this->py + order -2;

	this->bswx = gsl_bspline_alloc(order, this->px);
	gsl_bspline_knots_uniform(this->x[0], this->x[this->nx-1] , this->bswx);
	this->bswy = gsl_bspline_alloc(order, this->py);
	gsl_bspline_knots_uniform(this->y[0], this->y[this->ny-1] , this->bswy);

	this->coeff = NULL;
	this->baseMatrix_x = NULL;
	this->baseMatrix_t_x = NULL;
	this->baseMatrix_y = NULL;
	this->baseMatrix_t_y = NULL;

	this->createBaseMatrix(bswx,this->x, &(this->baseMatrix_x), &(this->baseMatrix_t_x), this->nSplinex, this->nx);
	this->createBaseMatrix(bswy, this->y, &(this->baseMatrix_y), &(this->baseMatrix_t_y), this->nSpliney, this->ny);

	this->createCoeffs();
}

Spline2d::Spline2d(double *x, size_t nx, double *y, size_t ny, double **map, size_t px, size_t py, size_t downsample){
	this->nx = ceil(nx/downsample);
	this->ny = ceil (ny/downsample);
	this->dx = double(nx-1)/double(this->nx-1);
	this->dy = double(ny-1)/double(this->ny-1);

	//this->x = x;
	this->x =  new double [this->nx];
	this->y = new double[this->ny];
	this->map = new double *[this->nx];

	for (size_t i=0; i<this->nx; i++){
		this->x[i]= x[size_t(round(i*dx))];
		this->map[i] = new double [ny];
	}
	for (size_t i=0; i<this->ny; i++)
		this->y[i] = y[size_t(round(i*dy))];

	for (size_t i=0; i<this->nx; i++)
		for (size_t j=0; j<this->ny; j++)
			this->map[i][j]= map[size_t(round(i*dx))][size_t(round(j*dy))];

	//this->map = map;
	this->order = 4;
	this->px = px;
	this->py = py;

	this->nSplinex = this->px +order -2;
	this->nSpliney = this->py + order -2;

	this->bswx = gsl_bspline_alloc(order, this->px);
	gsl_bspline_knots_uniform(this->x[0], this->x[this->nx-1] , this->bswx);
	this->bswy = gsl_bspline_alloc(order, this->py);
	gsl_bspline_knots_uniform(this->y[0], this->y[this->ny-1] , this->bswy);

	this->coeff = NULL;
	this->baseMatrix_x = NULL;
	this->baseMatrix_t_x = NULL;
	this->baseMatrix_y = NULL;
	this->baseMatrix_t_y = NULL;

	this->createBaseMatrix(bswx,this->x, &(this->baseMatrix_x), &(this->baseMatrix_t_x), this->nSplinex, this->nx);
	this->createBaseMatrix(bswy, this->y, &(this->baseMatrix_y), &(this->baseMatrix_t_y), this->nSpliney, this->ny);

	this->createCoeffs();
}

void Spline2d::createBaseMatrix(gsl_bspline_workspace *bsw,double *time,cs **bMatrix, cs **bMatrix_t, size_t nSpline, size_t nData){

	  gsl_vector  *btdata = gsl_vector_alloc(nSpline);
	  long kstart =0;
	  long kend = 0;
	  double c_sample = 0;
	  double tol = 0.0;
	  cs *tmpbMatrix = cs_spalloc(nSpline, nData, 5*nSpline,1,1);
	  cs *tmpbMatrix_t = cs_spalloc(nData, nSpline,5*nSpline,1,1);


	  for (size_t i=0; i<nData; i++){
	    gsl_bspline_eval(time[i], btdata, bsw);
	    kstart = -1;
	    kend = -1;
	    for (size_t k=0; k<nSpline; k++){
	      c_sample = gsl_vector_get(btdata, k);
	      if (kstart == -1 && c_sample > tol){
	    	  kstart = k;
	    	  continue;
	      }
	      if (kend ==-1 && kstart > -1 &&  c_sample <=tol )
	    	  kend = k;
	    }
	    if (kend == -1 && kstart != -1)
	       kend = nSpline;
	    for (long k=kstart; k<kend; k++){
	      c_sample = gsl_vector_get(btdata, k);
	      cs_entry(tmpbMatrix, k, i , c_sample);
	      cs_entry(tmpbMatrix_t, i , k , c_sample);
	    }
	  }
	  (*bMatrix) = cs_compress(tmpbMatrix);
	  (*bMatrix_t) = cs_compress (tmpbMatrix_t);

	  cs_spfree(tmpbMatrix);
	  cs_spfree(tmpbMatrix_t);
	  gsl_vector_free(btdata);
}


void Spline2d::create1DCoeffs(double *data, size_t nData, cs *bMatrix, cs *btbMatrix, size_t nSpline, double *coeffOut){

	for (size_t j=0; j<nSpline; j++)
			coeffOut[j]=0.0;

	cs_gaxpy (bMatrix,data,coeffOut);
	if (!cs_qrsol(3, btbMatrix, coeffOut)){
			cerr<<"Cleand2dStripe():: Unable to find coefficient for row"<<endl;
			exit(-1);
	}

//	for (size_t j=0; j<nSpline;j++)
//		cerr<<coeffOut[j]<<",";
//	cerr<<endl;
}

void Spline2d::createCoeffs(){

	this->coeff = new double* [this->nx];


	cs *btbMatrix = cs_multiply(this->baseMatrix_y, this->baseMatrix_t_y);

	//double *downMap = new double [this->ny];

	for (size_t i=0; i<this->nx; i++){
		this->coeff[i]= new double [this->nSpliney];
//		for (size_t j=0; j<this->ny; j++){
//			downMap[j] = this->map[int(round(i*dx))][int(round(j*dy))];
//		}
		this->create1DCoeffs(this->map[i], this->ny, this->baseMatrix_y, btbMatrix, this->nSpliney, this->coeff[i]);
//		for (size_t j=0; j<this->nSpliney; j++)
//			cerr<<this->coeff[i][j]<<",";
//		cerr<<endl;

	}
//	cin >> a;
	cs_spfree(btbMatrix);
	//delete []downMap;
	cerr<<"Spline2d::createCoeffs(). 1D Coefficients created"<<endl;
}

void Spline2d::interpolate(double *xi, double *yi, size_t nData, double *dataOut){

	double *iSplinex = new double [this->ny];
	double *jValuesx = new double [this->nx];
	double *coeffx = new double [this->nSplinex];
	double *cEval = new double  [this->nx];
	cs *btbMatrix = NULL;

	double xpos=0;
	double ypos=0;
	size_t fpox=0;
	size_t cpox=0;
	size_t fpoy=0;
	size_t cpoy=0;

	btbMatrix = cs_multiply(this->baseMatrix_x, this->baseMatrix_t_x);
	for (size_t i =0; i<nData; i++){
		xpos = (xi[i]-this->x[0])*(this->nx-1)/(this->x[this->nx-1]-this->x[0]);
		fpox = floor(xpos);
		cpox = fpox+1;
		ypos = (yi[i]-this->y[0])*(this->ny-1)/(this->y[this->ny-1]-this->y[0]);
		fpoy = floor(ypos);
		cpoy = fpoy +1;
		for (size_t k=0; k<this->nx; k++){
			for (size_t j=0; j<this->ny; j++)
				iSplinex[j]=0.0;
			cs_gaxpy(this->baseMatrix_t_y, this->coeff[k], iSplinex);
			jValuesx[k]= iSplinex[fpoy];

		}
//		cerr<<"i"<<i<<" x:"<<xi[i]<<" y: "<< y[i]<< " px,py:"<< fpox <<","<<fpoy<<endl;
		this->create1DCoeffs(jValuesx,this->nx,this->baseMatrix_x, btbMatrix, this->nSplinex, coeffx);
//		for (size_t ixx = 0; ixx< this->nSplinex; ixx++)
//			cerr<<coeffx[ixx]<<",";
//		cerr<<endl;
//		sleep(1);

		for (size_t k=0; k<this->nx; k++)
			cEval[k]=0.0;
		cs_gaxpy(this->baseMatrix_t_x, coeffx, cEval);
		dataOut[i]= cEval[fpox];
	}
	cs_spfree(btbMatrix);
}

double Spline2d::interpolate_single(double xi, double yi){
	double result;
	this->interpolate(&xi,&yi,1,&result);
//	for (size_t i=0; i< this->nx; i++){
//		for (size_t j=0; j< this->nSpliney; j++)
//			cerr <<coeff[i][j]<<",";
//		cerr<<endl;
//	}
	return result;
}

Spline2d::~Spline2d(){
	cs_spfree(this->baseMatrix_x);
	cs_spfree(this->baseMatrix_t_x);
	cs_spfree(this->baseMatrix_y);
	cs_spfree(this->baseMatrix_t_y);

	for  (size_t i=0; i<this->nx; i++)
		delete [] this->coeff[i];

	delete [] coeff;

	gsl_bspline_free(this->bswx);
	gsl_bspline_free(this->bswy);
}


Spline2dInterp::Spline2dInterp(double *x, size_t nx, double *y, size_t ny, double **map, size_t px, size_t py){
	this->npx = px +1;
	this->npy = py +1;

	double dx = (double(nx)-1.0)/double(px);
	double dy = (double(ny)-1.0)/double(py);

	this->x = new double[this->npx];
	this->y = new double[this->npy];
	this->map = new double * [this->npx];
	cerr<<"Spline2dInter():  dx = "<<dx<< " X values:  ";
	for (size_t i=0; i<this->npx; i++){
		this->x[i] = x[size_t(round(dx*i))];
		cerr<<this->x[i]<<",";
	}
	cerr<<endl;
	cerr<<"Spline2dInter():  dx = "<<dy<< " X values:  ";
	for (size_t j=0; j<this->npy; j++){
		this->y[j]= y[size_t(round(dy*j))];
		cerr<<this->y[j]<<",";
	}
	cerr<<endl;
	for (size_t i=0; i<this->npx; i++){
		this->map[i] = new double [this->npy];
		for (size_t j=0; j<this->npy; j++){
			this->map[i][j]= map[size_t(round(dx*i))][size_t(round(dy*j))];
		}
	}

	this->splineArray = new gsl_spline* [this->npx];
	this->accArray = new gsl_interp_accel* [this->npx];

	for (size_t i=0; i<npx; i++){
		splineArray[i] = gsl_spline_alloc(gsl_interp_cspline, this->npy);
		accArray[i] = gsl_interp_accel_alloc();
		gsl_spline_init(splineArray[i],this->y,this->map[i], this->npy);
	}
}

void Spline2dInterp::interpolate(double *xi, double *yi, size_t nData, double *dataOut){
	gsl_spline *splinex= NULL;
	gsl_interp_accel *accelx = NULL;
	double xp [this->npx];
	for (size_t i=0; i<nData; i++){
		for (size_t j=0; j<this->npx; j++){
			xp[j]= gsl_spline_eval(this->splineArray[j], yi[i],this->accArray[j]);
		}
		splinex = gsl_spline_alloc(gsl_interp_cspline, this->npx);
		gsl_spline_init(splinex,this->x, xp, this->npx);
		accelx = gsl_interp_accel_alloc();
		dataOut[i]=gsl_spline_eval(splinex, xi[i], accelx);
		gsl_spline_free(splinex);
		gsl_interp_accel_free(accelx);
	}

}

double Spline2dInterp::interpolate_single(double xi, double yi){
	double result;
	this->interpolate(&xi,&yi,1,&result);

	return result;
}

Spline2dInterp::~Spline2dInterp(){
	delete [] x;
	delete [] y;
	for (size_t i=0; i<this->npx; i++){
		delete [] map[i];
		gsl_spline_free(this->splineArray[i]);
		gsl_interp_accel_free(this->accArray[i]);
	}

	delete [] this->splineArray;
	delete [] this->accArray;
}
