#include <netcdfcpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <cstdio>
#include <fftw3.h>
#include <omp.h>
#include <vector>
using namespace std;

#include "nr3.h"
#include "AnalParams.h"
#include "Array.h"
#include "astron_utilities.h"
#include "gaussFit.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include "Map.h"
#include "Telescope.h"
#include "vector_utilities.h"


///Map constructor
/**The Map constructor allocates memory and initializes values
   in the map image to 0.
**/
Map::Map(string name, int nr, int nc, double pixsz, 
	 MatDoub &wt, VecDoub &rcp, VecDoub &ccp)
{
  //initialization
  nrows = nr;
  ncols = nc;
  nPixels = nrows*ncols;

  //the name
  mapName.assign(name);

  //pixels and coordinates
  pixelSize = pixsz;
  rowCoordsPhys.resize(nrows);
  colCoordsPhys.resize(ncols);
  for(int i=0;i<nrows;i++) rowCoordsPhys[i] = rcp[i];
  for(int j=0;j<ncols;j++) colCoordsPhys[j] = ccp[j];
  weight.resize(wt.nrows(),wt.ncols());
  for(int i=0;i<wt.nrows();i++) 
    for(int j=0;j<wt.ncols();j++)  weight[i][j] = wt[i][j];

  //coverage cut
  coverageCut = 0.;
  cutXRange.resize(2);
  cutYRange.resize(2);
  cutXRange[0] = 0;
  cutXRange[1] = nrows-1;
  cutYRange[0] = 0;
  cutYRange[1] = ncols-1;
  weightCut = -99.0;

  //size the image array and fill with zeros
  image.assign(nrows, ncols, 0.);
}


double Map::select(vector<double> input, int index){
 //Partition input based on a selected pivot
    //More on selecting pivot later- it's random
  //you now know the index of the pivot
  //if that is the index you want, return
  //else, if it's larger, recurse on the larger partition
  //if it's smaller, recurse on the smaller partition
  unsigned int pivotIndex = rand() % input.size();
  double pivotValue = input[pivotIndex];
  vector<double> left;
  vector<double> right;
  for(unsigned int x = 0; x < input.size(); x++){
    if(x != pivotIndex){
      if(input[x] > pivotValue){
        right.push_back(input[x]);
      }
      else{
        left.push_back(input[x]);
      }
    }
  }
  if((int) left.size() == index){
    return pivotValue;
  }
  else if((int) left.size() < index){
    return select(right, index - left.size() - 1);
  }
  else{
    return select(left, index);
  }
}


//----------------------------- o ---------------------------------------

///converts between tangential projection coordinates to image index
/** Converts between tangential projection and image index.  This
    is really a two line function but I put in a bunch of error
    checking to avoid dumb mistakes that would lead to segmentation
    faults.
**/
bool Map::raDecPhysToIndex(double ra, double dec, int* irow, int* icol)
{
  if(ra < rowCoordsPhys[0]  || ra > rowCoordsPhys[nrows-1] ){
    cerr << "Map::RaDecPhysToIndex(): ";
    cerr << "Ra for map is out of bounds." << endl;
    cerr << "Ra=" << ra << ", limits are [" << rowCoordsPhys[0] << ",";
    cerr << rowCoordsPhys[nrows-1] << "]." << endl;
    exit(1);
  }

  if(dec < colCoordsPhys[0]  || dec > colCoordsPhys[ncols-1] ){
    cerr << "Map::RaDecPhysToIndex(): ";
    cerr << "Dec for map is out of bounds." << endl;
    cerr << "Dec=" << dec << ", limits are [" << colCoordsPhys[0] << ",";
    cerr << colCoordsPhys[ncols-1] << "]." << endl;
    cerr << "Pixel size: "<< pixelSize <<endl;
    exit(1);
  }

  *irow = floor((ra-rowCoordsPhys[0])/pixelSize);
  *icol = floor((dec-colCoordsPhys[0])/pixelSize);

  //alternate definition
  *irow = ra/pixelSize + (nrows+1.)/2.;
  *icol = dec/pixelSize + (ncols+1.)/2.; 


  if(*irow >= nrows){
    cerr << "Map::RaDecPhysToIndex(): ";
    cerr << "Map row index: " << *irow << " is out of bounds." << endl;
    exit(1);
  }

  if(*icol >= ncols){
    cerr << "Map::RaDecPhysToIndex(): ";
    cerr << "Map column index: " << *icol << " is out of bounds." << endl;
    exit(1);
  }

  return 1;
}


//----------------------------- o ---------------------------------------

///fits image to gaussian located near map center
/** Fits the Map image to a 2-d gaussian using the gsl nonlinear least
    squares fitter as implemented in the Utility gaussFit.  There are
    a few things in here that are hard coded that should be
    generalized.  
    The output best-fit parameters are written into the map's netcdf
    file as attributes to the map that has been fit.
    \todo Generalize hard-coded aspects of this function.
    \todo Add a position angle to the fit so that we can properly
          measure the ellipticity of the fitted gaussian.
**/
double Map::fitToGaussian(VecDoub &pp, VecInt &fixme, VecDoub &fixVals, double *iguess, int deg)
{
  //if map is big, make assumption that gaussian is within 
  //40" of center of map, otherwise skip this
  double rowsize = rowCoordsPhys.size()*pixelSize;
  double colsize = colCoordsPhys.size()*pixelSize;
  double r = (double)(deg)/3600.*TWO_PI/360.;
  int minxi=0;
  int maxxi=nrows;
  int minyi=0;
  int maxyi=ncols;
  double medianRcp = median(rowCoordsPhys);
  double medianCcp = median(colCoordsPhys);

  if(rowsize > 2*r || colsize > 2*r){
    double minx = max(medianRcp-r,rowCoordsPhys[0]);
    double maxx = min(medianRcp+r,rowCoordsPhys[nrows-1]);
    double miny = max(medianCcp-r,colCoordsPhys[0]);
    double maxy = min(medianCcp+r,colCoordsPhys[ncols-1]);
  
    //find corresponding indices for min and max in each dimension
    for(int i=0;i<nrows;i++){
      if(rowCoordsPhys[i] >= minx){
	minxi=i;
	break;
      }
    }
    for(int i=0;i<nrows;i++){
      if(rowCoordsPhys[i] <= maxx){
	maxxi=i;
      }
    }
    for(int i=0;i<ncols;i++){
      if(colCoordsPhys[i] >= miny){
	minyi=i;
	break;
      }
    }
    for(int i=0;i<ncols;i++){
      if(colCoordsPhys[i] <= maxy){
	maxyi=i;
      }
    }
  } else {
    minxi = 0;
    maxxi = rowCoordsPhys.size()-1;
    minyi = 0;
    maxyi = colCoordsPhys.size()-1;
  }

  //make up the input arrays to the fitter
  int nptsrows = maxxi-minxi;
  int nptscols = maxyi-minyi;
  int npts = nptsrows*nptscols;
  VecDoub params(7);
  VecDoub params_err(7);
  VecDoub az(npts);
  VecDoub el(npts);
  VecDoub m(npts);
  VecDoub sigma(npts);
  for(int i=0;i<nptsrows;i++)
    for(int j=0;j<nptscols;j++){
      az[nptscols*i + j] = rowCoordsPhys[i+minxi];
      el[nptscols*i + j] = colCoordsPhys[j+minyi];
      m[nptscols*i + j] = image[i+minxi][j+minyi];
      sigma[nptscols*i + j] = (weight[i+minxi][j+minyi] == 0.) ?
    		  	  	  1. : sqrt(1./weight[i+minxi][j+minyi]);
    }

  //call the fitter for the map
  double chisqpdof;
  chisqpdof = mpGaussFit(&params[0], &params_err[0], &az[0], &el[0], 
			 &m[0], &sigma[0], npts, fixme, fixVals, iguess);


  cerr << mapName << " Gaussian fit parameters: " << endl;
  cerr.precision(4);
  cerr << "   dc offset = " << params[0] << " +/- " << params_err[0]
       << " Jy" << endl;
  cerr << "   amplitude = " << params[1] << " +/- " << params_err[1] 
       << " Jy" << endl;
  cerr << "   FWHM x   = " << params[2]/TWO_PI*360.*3600.*2.3548 
       << " +/- " << params_err[2]/TWO_PI*360.*3600.*2.3548 
       << " arcsec" << endl;
  cerr << "   FWHM y   = " << params[3]/TWO_PI*360.*3600.*2.3548
       << " +/- " << params_err[3]/TWO_PI*360.*3600.*2.3548 
       << " arcsec" << endl;
  cerr << "   offset x  = " << params[4]/TWO_PI*360.*3600.
       << " +/- " << params_err[4]/TWO_PI*360.*3600. 
       << " arcsec" << endl;
  cerr << "   offset y  = " << params[5]/TWO_PI*360.*3600.
       << " +/- " << params_err[5]/TWO_PI*360.*3600. 
       << " arcsec" << endl;
  cerr << "   pos angle =  " << params[6] << " +/- " << params_err[6]
       << " radians" << endl;

  //if mapFile is not the null string then write the results of the fit
  if(mapFile.length() > 0){
#pragma omp critical (dataio)
{
    NcFile ncfid = NcFile(mapFile.c_str(), NcFile::Write);
    NcVar* sigVar= ncfid.get_var(mapName.c_str());
    sigVar->add_att("dc_offset", params[0]);
    sigVar->add_att("dc_offset_err", params_err[0]);
    sigVar->add_att("dc_offset_units", "Jy");
    sigVar->add_att("amplitude", params[1]);
    sigVar->add_att("amplitude_err", params_err[1]);
    sigVar->add_att("amplitude_units", "Jy");
    sigVar->add_att("FWHM_x", params[2]/TWO_PI*360.*3600.*2.3548);
    sigVar->add_att("FWHM_x_err", params_err[2]/TWO_PI*360.*3600.*2.3548);
    sigVar->add_att("FWHM_x_units", "arcseconds");
    sigVar->add_att("FWHM_y", params[3]/TWO_PI*360.*3600.*2.3548);
    sigVar->add_att("FWHM_y_err", params_err[3]/TWO_PI*360.*3600.*2.3548);
    sigVar->add_att("FWHM_y_units", "arcseconds");
    sigVar->add_att("offset_x", params[4]/TWO_PI*360.*3600.);
    sigVar->add_att("offset_x_err", params_err[4]/TWO_PI*360.*3600.);
    sigVar->add_att("offset_x_units", "arcseconds");
    sigVar->add_att("offset_y", params[5]/TWO_PI*360.*3600.);
    sigVar->add_att("offset_y_err", params_err[5]/TWO_PI*360.*3600.);
    sigVar->add_att("offset_y_units", "arcseconds");
    sigVar->add_att("pos_Angle", params[6]);
    sigVar->add_att("pos_Angle_err", params_err[6]);
    sigVar->add_att("pos_Angle_units", "radians");
    ncfid.close();
}
  }

  //copy the parameters to the pp vector before leaving
  pp.resize(14);
  for(int i=0;i<7;i++) pp[i] = params[i];
  for(int i=7;i<14;i++) pp[i] = params_err[i-7];

  return chisqpdof;
}

//----------------------------- o ---------------------------------------

double Map::fitToGaussian()
{
  VecDoub pp(7);
  return fitToGaussian(pp);
}


//----------------------------- o ---------------------------------------

double Map::fitToGaussian(VecDoub &pp)
{
  VecInt fixme(7,0);
  VecDoub fixVals(7,0.);
  return fitToGaussian(pp,fixme,fixVals);
}


//----------------------------- o ---------------------------------------

double Map::fitToGaussian(VecDoub &pp, int deg)
{
  VecInt fixme(7,0);
  VecDoub fixVals(7,0.);
  return fitToGaussian(pp,fixme,fixVals,NULL,deg);
}

//----------------------------- o ---------------------------------------


///finds the coverage cut value in weight
/** Finds the coverage cut value in weight to be used in making a
    boolean map of good coverage and/or getting the low and high
    row and column indices of the good coverage region. This routine
    uses the same algorithm used in the idl utilities to get this
    value. No modification of the image is done.
**/
bool Map::findWeightThresh()
{
  //number of elements in map
  vector<double> og;
  for(int x = 0; x <nrows; x++){
    for(int y = 0; y<ncols; y++){
      if(weight[x][y] > 0.){
	og.push_back(weight[x][y]);
      }
    }
  }
  
  //COVCOV IS THE ORIGINAL UNSORTED VECTOR
  //COVPRM IS THE NEW SORTED VECTOR
  
  //find the point where 25% of nonzero elements have greater weight
  double covlim;
  int covlimi;
  covlimi = 0.75* og.size();

	//COVLIMI IS THE INDEX OF THE VALUE THAT HAS 25% OF NONZERO VALUES ABOVE IT
  covlim = select(og, covlimi);
	//COVLIM IS THE VALUE AT THE INDEX COVLIMI

  double mval;
  double mvali;
  mvali = floor((covlimi+og.size())/2.);
	//MVALI IS NOW THE INDEX THAT IS HALFWAY BEWEEN THE NUMBER OF POINTS AND COVLIMI, THE 75% INDEX
  mval = select(og, mvali);
  //MVAL IS THE VALUE AT INDEX MVALI
	
	//define the weight cut value
  weightCut = coverageCut*mval;
	//SET WEIGHT CUT TO BE COVERAGECUT THAT WAS PREVIOUSLY DEFINED * MVAL
  //done with covcov and covprm
  return 1;
}

//----------------------------- o ---------------------------------------


///finds the column and row indices correspoinding to the coverageCut value
/** Finds the row and column indices (low and high for each) that
    correspond to the previously set coverage cut value. No modification
    of the image is done.
**/
bool Map::setCoverageCutRanges()
{
  //check that findWeightThresh() has already been called
  if(weightCut < 0.0){
    findWeightThresh();
  }

  //find lower row bound
  bool flag = false;
  for(int i=0;i<nrows;i++){
    for(int j=0;j<ncols;j++){
    	if(weight[i][j] >= weightCut){
    		cutXRange[0] = i;
    		flag = true;
    		break;
    	}
    }
	if(flag == true){
		break;
	}
  }

  //find upper row bound
  flag = false;
  for(int i=nrows-1;i>-1;i--){
    for(int j=0;j<ncols;j++){
    	if(weight[i][j] >= weightCut){
    		cutXRange[1] = i;
    		flag = true;
    		break;
    	}
    }
	if(flag == true){
		break;
	}
  }

  //find lower column bound
  flag = false;
  for(int i=0;i<ncols;i++){
    for(int j=cutXRange[0];j<cutXRange[1]+1;j++){
      if(weight[j][i] >= weightCut){
    	  cutYRange[0] = i;
    	  flag = true;
    	  break;
      }
    }
    if(flag == true){
    	break;
    }
  }

  //find upper column bound
  flag = false;
  for(int i=ncols-1;i>-1;i--){
    for(int j=cutXRange[0];j<cutXRange[1]+1;j++){
    	if(weight[j][i] >= weightCut){
    		cutYRange[1] = i;
    		flag = true;
    		break;
    	}
    }
    if(flag == true){
    	break;
    }
  }

  return 1;
}

//----------------------------- o ---------------------------------------


///make boolean map of pixels that are greater than or equal to weightCut
/** Makes a boolean map (MatBool) indicating where the good coverage
    region is with 0 and 1 corresponding to outside good coverage and
    part of good coverage respectively. Done by checking each pixel
    weight if it is greater than or equal to weightCut.
**/
bool Map::makeCovBoolMap()
{
  //check that findWeightThresh() has already been called
  if(weightCut < 0.0){
    findWeightThresh();
  }

  //find where pixels have weight >= weightCut
  coverageBool.resize(nrows, ncols);
  for(int i=0;i<nrows;i++){
    for(int j=0;j<ncols;j++){
      if(weight[i][j] >= weightCut){
	coverageBool[i][j] = 1;
      }
      else{
	coverageBool[i][j] = 0;
      }
    }  
  }

  return 1;
}


//----------------------------- o ---------------------------------------


///calculates the 1d map psd
/** Calculate the 1d map psd using the same method and normalization
    as used in aztec_idl_utilities.
**/
bool Map::calcMapPsd(double covCut)
{
  //make sure we've got up to date coverage cut indices
  coverageCut = covCut;
  findWeightThresh();
  setCoverageCutRanges();

  //make sure our coverage cut map has an even number
  //of rows and columns
  int nx = cutXRange[1]-cutXRange[0]+1;
  int ny = cutYRange[1]-cutYRange[0]+1;
  int cxr0 = cutXRange[0];
  int cyr0 = cutYRange[0];
  int cxr1 = cutXRange[1];
  int cyr1 = cutYRange[1];
  if(nx % 2 == 1){
    cxr1 = cutXRange[1]-1;
    nx--;
  }
  if(ny % 2 == 1){
    cyr1 = cutYRange[1]-1;
    ny--;
  }

  //we will do the fft using fftw
  //intially written with the inefficient way to make bookkeeping simpler
  //here is the memory allocation and the plan setup
  fftw_complex *in;
  fftw_complex *out;
  fftw_plan *p;

  #pragma omp critical (noiseFFT)
  {
	  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny);
	  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny);
	  p = new fftw_plan;
	  *p = fftw_plan_dft_2d(nx, ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
         // cout << nx << " " << ny << " " << in << " " << out << " " << FFTW_FORWARD << " " << FFTW_ESTIMATE << endl;
  }

  //the matrix to get fft'd is cast into vector form in *in;
  int ii,jj,stride,index;
  for(int i=cxr0;i<=cxr1;i++)
    for(int j=cyr0;j<=cyr1;j++){
      ii = i-cxr0;
      jj = j-cyr0;
      stride=cyr1-cyr0+1;
      index = stride*ii+jj;
      in[index][0] = image[i][j];
      in[index][1] = 0.;
    }

  //apply a hanning window
  MatDoub h = hanning(nx, ny);
  for(int i=0;i<nx;i++) for(int j=0;j<ny;j++) in[ny*i+j][0] *= h[i][j];


  //calculate frequencies
  double diffx = rowCoordsPhys[1]-rowCoordsPhys[0];
  double diffy = colCoordsPhys[1]-colCoordsPhys[0];
  double xsize = diffx*nx;
  double ysize = diffy*ny;
  double diffqx = 1./xsize;
  double diffqy = 1./ysize;

  //do the fft and cleanup the plan
  fftw_execute(*p);

  #pragma omp critical (noiseFFT)
  {
	  fftw_destroy_plan(*p);
	  delete p;
  }
  //matching the idl code
  for(int i=0;i<nx*ny;i++){
    out[i][0] *= xsize*ysize/nx/ny;
    out[i][1] *= xsize*ysize/nx/ny;
  }


  //here is the magnitude
  //reuse h to be memory-kind, this is pmfq in the idl code
  for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++)
      h[i][j] = diffqx*diffqy*(pow(out[ny*i+j][0],2)+pow(out[ny*i+j][1],2));


  VecDoub w(nx*ny);
  for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++)
      w[ny*i+j] = h[i][j];


  //free up resources
  #pragma omp critical (noiseFFT)
  {
    	fftw_free(in);
    	fftw_free(out);
  }


  //vectors of frequencies
  VecDoub qx(nx);
  VecDoub qy(ny);
  int shift = nx/2-1;
  for(int i=0;i<nx;i++){
    index = i-shift;
    if(index < 0) index += nx;
    qx[index] = diffqx*(i-(nx/2-1));
  }
  shift = ny/2-1;
  for(int i=0;i<ny;i++){
    index = i-shift;
    if(index < 0) index += ny;
    qy[index] = diffqy*(i-(ny/2-1));
  }

  //shed first row and column of h, qx, qy
  MatDoub pmfq(nx-1,ny-1);
  for(int i=1;i<nx;i++) for(int j=1;j<ny;j++) pmfq[i-1][j-1] = h[i][j];
  for(int i=0;i<nx-1;i++) qx[i] = qx[i+1];
  for(int j=0;j<ny-1;j++) qy[j] = qy[j+1];

  //matrices of frequencies and distances
  MatDoub qmap(nx-1,ny-1);
  MatDoub qsymm(nx-1,ny-1);
  for(int i=1;i<nx;i++)
    for(int j=1;j<ny;j++){
      qmap[i-1][j-1] = sqrt(pow(qx[i],2) + pow(qy[j],2));
      qsymm[i-1][j-1] = qx[i]*qy[j];
    }


  //find max of nx and ny and correspoinding diffq
  int nn;
  double diffq;
  if(nx > ny){
    nn = nx/2+1;
    diffq = diffqx;
  }else{
    nn = ny/2+1;
    diffq = diffqy;
  }

  //generate the final vector of frequencies
  psdFreq.resize(nn);
  for(int i=0;i<nn;i++) psdFreq[i] = diffq*(i+0.5);
  
  //pack up the final vector of psd values
  psd.resize(nn);
  for(int i=0;i<nn;i++){
    int countS=0;
    int countA=0;
    double psdarrS=0.;
    double psdarrA=0.;
    for(int j=0;j<nx-1;j++)
      for(int k=0;k<ny-1;k++){
	if((int) (qmap[j][k] / diffq) == i && qsymm[j][k] >= 0.){
	  countS++;
	  psdarrS += pmfq[j][k];
	}
	if((int) (qmap[j][k] / diffq) == i && qsymm[j][k] < 0.){
	  countA++;
	  psdarrA += pmfq[j][k];
	}
      }
    if(countS != 0) psdarrS /= countS;
    if(countA != 0) psdarrA /= countA;
    psd[i] = min(psdarrS,psdarrA);
  }

  //smooth the psd with a 10-element boxcar filter
  VecDoub tmp(nn);
  smooth_edge_truncate(psd, tmp, 10);
  for(int i=0;i<nn;i++) psd[i]=tmp[i];

  psd2d = pmfq;
  psd2dFreq = qmap;
  //write the results into the netcdf file
  if(mapFile.length() > 0){
#pragma omp critical (dataio)
	  {
		NcFile ncfid = NcFile(mapFile.c_str(), NcFile::Write);

		//create dimension
		string psdDimName = "npsd_";
		psdDimName.append(mapName);
		NcDim* psdDim = ncfid.add_dim(psdDimName.c_str(), nn);

		//define psd variables
		string psdVarName = "psd_";
		string psdFreqName = "psdFreq_";
		psdVarName.append(mapName);
		psdFreqName.append(mapName);
		NcVar *psdVar = ncfid.add_var(psdVarName.c_str(), ncDouble, psdDim);
		NcVar *psdFreqVar = ncfid.add_var(psdFreqName.c_str(), ncDouble, psdDim);

		//put in the values
		psdVar->put(&psd[0], nn);
		psdFreqVar->put(&psdFreq[0], nn);

		//add an attribute to the frequency values to remind us of the units
		psdFreqVar->add_att("units", "1/radians");


		psdDimName = "nxpsd_2d";
		NcDim* psdDim2d_x = ncfid.add_dim(psdDimName.c_str(), nx-1);
		psdDimName = "nypsd_2d";
		NcDim* psdDim2d_y = ncfid.add_dim(psdDimName.c_str(), ny-1);

		psdVarName = "psd_2d";
		psdFreqName = "psdFreq_2d";

		NcVar *psdVar2d = ncfid.add_var(psdVarName.c_str(), ncDouble, psdDim2d_x, psdDim2d_y);
		NcVar *psdFreqVar2d = ncfid.add_var(psdFreqName.c_str(), ncDouble, psdDim2d_x, psdDim2d_y);

		psdVar2d->put(&psd2d[0][0], (nx-1),(ny-1));
		psdFreqVar2d->put(&psd2dFreq[0][0], (nx-1),(ny-1));
		ncfid.close();
	  }
  }    

	
  return 1;
}




//----------------------------- o ---------------------------------------


///generates histogram of image values
/** Generates a histogram of image values for a particular coverage
    cut.  The output histogram is written into the map's netcdf file.
    inputs are:
      - nbins - the number of bins in the output histogram
      - cc - the desired coverage cut (0<cc<1).
 **/
bool Map::calcMapHistogram(int nbins, double cc)
{
  //space for the results
  histBins.resize(nbins);
  histVals.resize(nbins);

  //set the coverage cut ranges
  coverageCut = cc;
  findWeightThresh();
  setCoverageCutRanges();

  //make the coverage cut
  int nx = cutXRange[1]-cutXRange[0]+1;
  int ny = cutYRange[1]-cutYRange[0]+1;
  MatDoub im(nx,ny);
  for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++) 
      im[i][j] = image[cutXRange[0]+i][cutYRange[0]+j];

  //do the histogram
  histogramImage(im, nbins, histBins, histVals);

  //write this into the mapFile if it exists
  if(mapFile.length() > 0){
#pragma omp critical (dataio)
	  {
		NcFile ncfid = NcFile(mapFile.c_str(), NcFile::Write);

		//create dimension
		string histDimName = "nhist_";
		histDimName.append(mapName);
		NcDim* histDim = ncfid.add_dim(histDimName.c_str(), nbins);

		//define hist variables
		string binVarName = "histBins_";
		string valVarName = "histVals_";
		binVarName.append(mapName);
		valVarName.append(mapName);
		NcVar *binVar = ncfid.add_var(binVarName.c_str(), ncDouble, histDim);
		NcVar *valVar = ncfid.add_var(valVarName.c_str(), ncDouble, histDim);

		//put in the values
		binVar->put(&histBins[0], nbins);
		valVar->put(&histVals[0], nbins);
		ncfid.close();
	  }
  }  

  return 1;
}


//----------------------------- o ---------------------------------------


int Map::getNrows()
{
  return nrows;
}


int Map::getNcols()
{
  return ncols;
}


double Map::getPixelSize()
{
  return pixelSize;
}


double Map::getRowCoordsPhys(int i)
{
  return rowCoordsPhys[i];
}


double Map::getColCoordsPhys(int i)
{
  return colCoordsPhys[i];
}


double Map::getCoverageCut()
{
  return coverageCut;
}


bool Map::setCoverageCut(double cc)
{
  coverageCut = cc;
  return 1;
}

int Map::getCutXRangeLow()
{
  return cutXRange[0];
}

int Map::getCutXRangeHigh()
{
  return cutXRange[1];
}

int Map::getCutYRangeLow()
{
  return cutYRange[0];
}

int Map::getCutYRangeHigh()
{
  return cutYRange[1];
}

//----------------------------- o ---------------------------------------

//destructor
Map::~Map()
{

}
