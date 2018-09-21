#include <netcdfcpp.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include <CCfits/CCfits>
using namespace std;
#include "nr3.h"

#define PI 3.1415926535897932384626433
#define TWO_PI (2.0 * PI)
#define PIO2  (PI / 2.0) 
#define DEG_RAD (360.0 / TWO_PI)


bool writeFits(string fitsFile, MatDoub &arr, double xpixsz, double ypixsz,
	       double refC1, double refC2, double refpixC1, double refpixC2) 
{
  //declare axis arrays
  int nrows = arr.nrows();
  int ncols = arr.ncols();
  long naxis = 2;
  long naxes[2] = {nrows, ncols};
  long nelements = naxes[0]*naxes[1];

  // declare auto-pointer to FITS at function scope. Ensures no resources
  // leaked if something fails in dynamic allocation.
  auto_ptr<CCfits::FITS> pFits(0);
  
  //open the fits file and overwrite it if it already exists
  try
    {                
      const std::string fileName(fitsFile);            
      // Create a new FITS object
      pFits.reset(new CCfits::FITS(fileName, DOUBLE_IMG, naxis, naxes));
    }
  catch (CCfits::FITS::CantCreate)
    {
      // ... or not, as the case may be.
      cerr << "Failed to open " << fitsFile << "... aborting" << endl;
      exit(1);
    }

  // the array to be written into the fits file
  // we'll take the transpose and reverse the ra axis too
  std::valarray<double> array(nelements);
  for(int i=0;i<ncols;i++)
    for(int j=0;j<nrows;j++)
      array[nrows*i+j] = arr[nrows-1-j][i];

  // don't forget to adjust the pixel location of the coord ref.
  // The extra +1 here is to account for the 1-based coordinate system of fits files
  refpixC1 = nrows-1-refpixC1+1;
  refpixC2 += 1;

  // I have no idea what this is for
  long fpixel(1);

  //this is what we eventually want
  pFits->pHDU().write(fpixel,nelements,array);

  //add header to fits file
  pFits->pHDU().addKey("SYSTEM","EQUATORIAL(J2000.0)","");
  pFits->pHDU().addKey("EQUINOX",2000.0,"Epoch of reference equinox");
  pFits->pHDU().addKey("CRPIX1", refpixC1, "Reference pixel on axis 1");    
  pFits->pHDU().addKey("CRPIX2", refpixC2, "Reference pixel on axis 2");
  pFits->pHDU().addKey("CD1_1", -xpixsz, "Transformation matrix element");
  pFits->pHDU().addKey("CD2_2", ypixsz, "Transformation matrix element");
  pFits->pHDU().addKey("CD1_2", 0.0, "Transformation matrix element");
  pFits->pHDU().addKey("CD2_1", 0.0, "Transformation matrix element");
  pFits->pHDU().addKey("CRVAL1", refC1, "Value at ref. pixel on axis 1");
  pFits->pHDU().addKey("CRVAL2", refC2, "Value at ref. pixel on axis 2");
  pFits->pHDU().addKey("CTYPE1", "RA---TAN", "Quantity represented by axis 1");
  pFits->pHDU().addKey("CTYPE2", "DEC--TAN", "Quantity represented by axis 2");
  pFits->pHDU().addKey("CRUNIT1", "deg", "");
  pFits->pHDU().addKey("CRUNIT2", "deg", "");
  pFits->pHDU().addKey("CDELT1", -xpixsz, "Transformation matrix element");
  pFits->pHDU().addKey("CDELT2", ypixsz, "Transformation matrix element");

  return 1;
}



int main(int nArgs, char* args[])
{
  //deal with command line inputs
  bool unfilt=0;
  bool postageStamps=0;
  if((nArgs > 3) || (nArgs == 1)){
    cerr << "fitswriter: writes fits files given map nc file input." << endl;
    cerr << "  calling syntax: " << endl;
    cerr << "     ./fitswriter <filename.nc> [-u | --unfiltered]" << endl;
    cerr << "     ./fitswriter <filename.nc> [-p | --postageStamps]" << endl;
    cerr << "  These modes are exclusive." << endl;
    exit(1);
  } else {
    if(nArgs == 3){
      string uft;
      uft.assign(args[2]);
      if((!uft.compare("--unfiltered")) || (!uft.compare("-u"))){
	unfilt=1;
      } else if ((!uft.compare("--postageStamps")) || (!uft.compare("-p"))){
	postageStamps=1;
      } else {
	cerr << "fitswriter: writes fits files given map nc file input." << endl;
	cerr << "  calling syntax: " << endl;
	cerr << "     ./fitswriter <filename.nc> [-u | --unfiltered]" << endl;
	cerr << "     ./fitswriter <filename.nc> [-p | --postageStamps]" << endl;
	cerr << "  These modes are exclusive." << endl;
	exit(1);
      }
    }
  }

  //string for the ncdf filename
  string ncdfFile;
  ncdfFile.assign(args[1]);

  //check that the file exists and fail if it doesn't
  ifstream f(ncdfFile.c_str());
  if (!f.good()) {
    cerr << "file not found." << endl;
    exit(1);
  }   

  //open the netcdf file and fetch the map dims
  NcFile ncfid(ncdfFile.c_str(), NcFile::ReadOnly);
  uint nrows = ncfid.get_dim("nrows")->size();
  uint ncols = ncfid.get_dim("ncols")->size();

  //the coordinate arrays
  VecDoub rcp(nrows);
  VecDoub ccp(ncols);
  NcVar* rcpv = ncfid.get_var("rowCoordsPhys");
  NcVar* ccpv = ncfid.get_var("colCoordsPhys");
  for(uint i=0;i<nrows;i++) rcp[i] = rcpv->as_double(i);
  for(uint j=0;j<ncols;j++) ccp[j] = ccpv->as_double(j);

  //header information
  double xpixsz = abs(rcp[1]-rcp[0]) * DEG_RAD;
  double ypixsz = abs(ccp[1]-ccp[0]) * DEG_RAD;
  NcAtt* mg0Att = ncfid.get_att("MasterGrid[0]");
  NcAtt* mg1Att = ncfid.get_att("MasterGrid[1]");
  double refC1 = mg0Att->as_double(0) * DEG_RAD;
  double refC2 = mg1Att->as_double(0) * DEG_RAD;
  int refpixC2 = ncols/2;
  int refpixC1 = nrows/2;

  //is this a noise file?
  NcError err(NcError::silent_nonfatal);
  NcVar* sVar;
  bool noiseFile=0;
  if((sVar = ncfid.get_var("noise"))) noiseFile=1;

  //signal map to fits
  string sigName;
  if(noiseFile){
    if(unfilt) sigName.assign("noise"); else sigName.assign("filteredNoise");
  } else {
    if(unfilt) sigName.assign("signal"); else sigName.assign("filteredSignal");
  }
  MatDoub arr(nrows,ncols);
  if(!(sVar = ncfid.get_var(sigName.c_str()))){
    cerr << ncdfFile << " does not contain filtered maps.  Try invoking " << endl;
    cerr << "with the -u or --unfiltered option to make fits files of the " << endl;
    cerr << "unfiltered maps." << endl;
    exit(1);
  }
  sVar->get(&arr[0][0], nrows, ncols);
  string sigFitsFile;
  sigFitsFile.assign(ncdfFile);
  if(unfilt){
    sigFitsFile.replace(sigFitsFile.find(".nc"), 3, "_signal_unfilt.fits");
  } else {
    sigFitsFile.replace(sigFitsFile.find(".nc"), 3, "_signal.fits");
  }
  sigFitsFile.insert(0,"!");
  writeFits(sigFitsFile, arr, xpixsz, ypixsz, refC1, refC2, 
	    refpixC1, refpixC2);

  //weight map to fits
  MatDoub wt(nrows,ncols);
  if(unfilt) sigName.assign("weight"); else sigName.assign("filteredWeight");
  sVar = ncfid.get_var(sigName.c_str());
  sVar->get(&wt[0][0], nrows, ncols);
  string wtFitsFile;
  wtFitsFile.assign(ncdfFile);
  if(unfilt){
    wtFitsFile.replace(wtFitsFile.find(".nc"), 3, "_weight_unfilt.fits");
  }else{
    wtFitsFile.replace(wtFitsFile.find(".nc"), 3, "_weight.fits");
  }
  wtFitsFile.insert(0,"!");
  writeFits(wtFitsFile, wt, xpixsz, ypixsz, refC1, refC2, 
	    refpixC1, refpixC2);

  //inttime map to fits
  MatDoub it(nrows,ncols);
  if(unfilt){
    sigName.assign("inttime");
    sVar = ncfid.get_var(sigName.c_str());
    sVar->get(&it[0][0], nrows, ncols);
    string itFitsFile;
    itFitsFile.assign(ncdfFile);
    itFitsFile.replace(itFitsFile.find(".nc"), 3, "_zinttime_unfilt.fits");
    itFitsFile.insert(0,"!");
    writeFits(itFitsFile, it, xpixsz, ypixsz, refC1, refC2, 
	      refpixC1, refpixC2);
  }

  
  //s2n map to fits
  MatDoub s2n(nrows,ncols);
  for(uint i=0;i<nrows;i++) 
    for(uint j=0;j<ncols;j++) 
      s2n[i][j] = sqrt(wt[i][j])*arr[i][j];
  if(!noiseFile){
    string s2nFitsFile;
    s2nFitsFile.assign(ncdfFile);
    if(unfilt){
      s2nFitsFile.replace(s2nFitsFile.find(".nc"), 3, "_sn_unfilt.fits");
    }else{
      s2nFitsFile.replace(s2nFitsFile.find(".nc"), 3, "_sn.fits");
    }
    s2nFitsFile.insert(0,"!");
    writeFits(s2nFitsFile, s2n, xpixsz, ypixsz, refC1, refC2, 
	      refpixC1, refpixC2);
    
    //kernel map to fits
    if(unfilt) sigName.assign("kernel"); else sigName.assign("filteredKernel");
    sVar = ncfid.get_var(sigName.c_str());
    sVar->get(&wt[0][0], nrows, ncols);
    string kerFitsFile;
    kerFitsFile.assign(ncdfFile);
    if(unfilt){
      kerFitsFile.replace(kerFitsFile.find(".nc"), 3, "_psf_unfilt.fits");
    }else{
      kerFitsFile.replace(kerFitsFile.find(".nc"), 3, "_psf.fits");
    }
    kerFitsFile.insert(0,"!");
    writeFits(kerFitsFile, wt, xpixsz, ypixsz, refC1, refC2, 
	      refpixC1, refpixC2);
  }

  if(postageStamps and !noiseFile){
    //check that there are sources
    NcError err(NcError::silent_nonfatal);
    bool hasSources=0;
    if((sVar = ncfid.get_var("pointSource_0"))) hasSources=1;
    if(!hasSources){
      return 0;
    }

    //get the map dimensions
    nrows = ncfid.get_dim("postageStampNRows")->size();
    ncols = ncfid.get_dim("postageStampNCols")->size();    

    //get the number of sources
    int nSources =  ncfid.get_dim("nSources")->size();

    //create a ds9 region file for the main fitsfiles
    ofstream ds9file;
    string regName;
    regName.assign(ncdfFile.c_str());
    regName.replace(regName.find(".nc"),3,"_sources.reg");
    ds9file.open (regName.c_str(), ios::out | ios::trunc); 
    ds9file << "# Region file format: DS9 version 4.0" << endl;
    ds9file << "global color=green font='helvetica 10 normal' ";
    ds9file << "select=1 highlite=1 edit=1 move=1 delete=1 ";
    ds9file << "include=1 fixed=0 source" << endl;
    ds9file << "fk5" << endl;

    //loop through the sources, writing a fits file for each one
    string varname;
    for(int k=0;k<nSources;k++){
      //the coordinate arrays
      rcp.resize(nrows);
      ccp.resize(ncols);

      //fetch row and column coordinates in pysical space
      //variable name is of the form pointSource_k_rcp
      varname.assign("pointSource_");
      stringstream o;
      o << k;
      varname.append(o.str());
      varname.append("_rcp");
      rcpv = ncfid.get_var(varname.c_str());
      varname.replace(varname.find("rcp"),3,"ccp");
      ccpv = ncfid.get_var(varname.c_str());
      for(uint i=0;i<nrows;i++) rcp[i] = rcpv->as_double(i);
      for(uint j=0;j<ncols;j++) ccp[j] = ccpv->as_double(j);

      //header information
      double xpixsz = abs(rcp[1]-rcp[0]) * DEG_RAD;
      double ypixsz = abs(ccp[1]-ccp[0]) * DEG_RAD;

      //mastergrid values are specified as attributes to pointSource_k
      varname.replace(varname.find("_ccp"),4,"");
      sVar = ncfid.get_var(varname.c_str());

      NcAtt* m0;
      NcAtt* m1;
      m0 = sVar->get_att("centerRaAbs");
      m1 = sVar->get_att("centerDecAbs");
      refC1 = m0->as_double(0) * DEG_RAD;
      refC2 = m1->as_double(0) * DEG_RAD;
      refpixC2 = ncols/2;
      refpixC1 = nrows/2;
    
      //map to fits
      arr.resize(nrows,ncols);
      sVar->get(&arr[0][0], nrows, ncols);
      sigFitsFile.assign(ncdfFile);
      varname.insert(0,"_");
      sigFitsFile.replace(sigFitsFile.find(".nc"), 3, varname.append(".fits").c_str());
      sigFitsFile.insert(0,"!");
      writeFits(sigFitsFile, arr, xpixsz, ypixsz, refC1, refC2, 
	      refpixC1, refpixC2);    


      //get the approximate signal to noise of the source in the map
      int rapix = sVar->get_att("centerXPos")->as_int(0);
      int decpix = sVar->get_att("centerYPos")->as_int(0);
      double snpix = s2n[rapix][decpix];

      //grab the last bit of data and put it into the region file.
      double rac = sVar->get_att("raCentroid")->as_double(0)*DEG_RAD;
      double decc = sVar->get_att("decCentroid")->as_double(0)*DEG_RAD;
      double fwhmx = sVar->get_att("FWHM_x")->as_double(0);
      double fwhmy = sVar->get_att("FWHM_y")->as_double(0);
      //double rad = 8.5/2./3600.;
      int oldprec = ds9file.precision();
      //ds9file << "circle(" << setprecision(10) << rac << "," << decc << "," << rad << setprecision(oldprec);
      ds9file << "ellipse(" << setprecision(10) << rac << "," << decc << "," << setprecision(3) << fwhmx/2. << "\"," <<	fwhmy/2. << "\"," << "0" << setprecision(oldprec);
      ds9file << ") ";
      ds9file << "# text = {";
      ds9file << "S"<<k<<": S/N~"<< setprecision(2) << snpix;
      ds9file << "}" << endl;
      ds9file.precision(oldprec);
    }
    ds9file.close();
  }
  return 0;
}

