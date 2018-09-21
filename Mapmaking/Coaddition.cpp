#include <netcdfcpp.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <cstdio>
#include <fftw3.h>
#include <vector>
using namespace std;

#include "nr3.h"
#include "AnalParams.h"
#include "Array.h"
#include "astron_utilities.h"
#include "gaussFit.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include "Coaddition.h"
#include "Telescope.h"
#include "vector_utilities.h"
#include "PointSource.h"

#include "convolution.h"  //trash when Wiener filter is done


///Coaddition constructor
/** There is not much going on in the coaddition constructor, just
    some simple initializations.  All of the real work goes on
    in Coaddition::coaddMaps().
**/
Coaddition::Coaddition(AnalParams* analParams)
{
  //set our ap pointer
  ap = analParams;

  //save typing
  pixelSize = ap->getPixelSize()/3600./360.*TWO_PI;
  double* mpp;
  mpp = ap->getMasterGridJ2000();
  masterGrid.resize(2);
  masterGrid[0] = mpp[0];
  masterGrid[1] = mpp[1];
  sourceSigma = ap->getSourceSigma();
  snglSourceWin = ap->getSnglSourceWin();
  nSources = 0;
  synth = false;

  tKernel = NULL;
  tSignal = NULL;
  tWeight = NULL;


  //initialize nrows and ncols
  nrows=0;
  ncols=0;
}


//----------------------------- o ---------------------------------------

//coaddMaps - coadd many observations together
/** This is the work horse of the coaddition.  Observations are read
    from .nc files and coadded together using the weight maps.  The
    map grid is set by the tangential projection with the masterGrid
    values as the tangent point to the sky.  This may not be the right
    way to do things if the masterGrid value falls far outside of the
    map bounds, this must be reviewed.
    \todo - what happens if masterGrid is outside of map bounds?
**/
bool Coaddition::coaddMaps()
{
  //start by running through map files to determine the coadded
  //maps bounds and to check that all the mastergrids are 
  //identical
  int nFiles = ap->getNFiles();
  VecDoub mg(2);
  double minRowVal=0., minColVal=0.;
  double maxRowVal=0., maxColVal=0.;
  string sourceNameString;
  individualMapsTau.resize(nFiles);
  for(int i=0;i<nFiles;i++){
    NcFile ncfid = NcFile(ap->getMapFileList(i).c_str(), NcFile::ReadOnly);
    NcDim* rowDimVar = ncfid.get_dim("nrows");
    NcDim* colDimVar = ncfid.get_dim("ncols");
    double rowDim = rowDimVar->size();
    double colDim = colDimVar->size();
    NcAtt* mg0Att = ncfid.get_att("MasterGrid[0]");
    NcAtt* mg1Att = ncfid.get_att("MasterGrid[1]");
    NcAtt* sourceName = ncfid.get_att("source");

    individualMapsTau[i] = ncfid.get_att("ArrayAvgTau")->as_double(0);

    //find minimum and maximum row and column values
    NcVar* rcpv = ncfid.get_var("rowCoordsPhys");
    NcVar* ccpv = ncfid.get_var("colCoordsPhys");
    if(i == 0){
      mg[0] = mg0Att->as_double(0);
      mg[1] = mg1Att->as_double(0);
      minRowVal = rcpv->as_double(0);
      maxRowVal = rcpv->as_double(rowDim-1);
      minColVal = ccpv->as_double(0);
      maxColVal = ccpv->as_double(colDim-1);
      char* tmpSourceName = sourceName->as_string(0);
      sourceNameString = tmpSourceName;
      delete [] tmpSourceName;
    }
    else {
      if(mg0Att->as_double(0) != mg[0] || mg1Att->as_double(0) != mg[1]){
    	  cerr << "Mastergrid is not consistent in file ";
    	  cerr << ap->getMapFileList(i) << " ... aborting coadd." << endl;
    	  exit(1);
      }

      if(rcpv->as_double(0) < minRowVal) minRowVal = rcpv->as_double(0);
      if(rcpv->as_double(rowDim-1) > maxRowVal) 
	maxRowVal = rcpv->as_double(rowDim-1);
      if(ccpv->as_double(0) < minColVal) minColVal = ccpv->as_double(0);
      if(ccpv->as_double(colDim-1) > maxColVal) 
	maxColVal = ccpv->as_double(colDim-1);
    }
    delete mg0Att;
    delete mg1Att;
    delete sourceName;
  }

  //set the mastergrid
  this->masterGrid[0] = mg[0];
  this->masterGrid[1] = mg[1];
  this->ap->setSourceName(sourceNameString);

  //explicitly center on mastergrid
  int xminpix = ceil(abs(minRowVal/pixelSize));
  int xmaxpix = ceil(abs(maxRowVal/pixelSize));
  xmaxpix = max(xminpix,xmaxpix);
  nrows = 2.*xmaxpix+4;
  int yminpix = ceil(abs(minColVal/pixelSize));
  int ymaxpix = ceil(abs(maxColVal/pixelSize));
  ymaxpix = max(yminpix,ymaxpix);
  ncols = 2.*ymaxpix+4;
  nPixels = nrows*ncols;

  //physical coordinates: this grid is set up so that physical
  //coordinates are 0 at center of pixel near center of map
  rowCoordsPhys.resize(nrows);
  colCoordsPhys.resize(ncols);
  for(int i=0;i<nrows;i++) rowCoordsPhys[i] = (i-(nrows+1.)/2.)*pixelSize;
  for(int i=0;i<ncols;i++) colCoordsPhys[i] = (i-(ncols+1.)/2.)*pixelSize;

  //sanity check
  if(nrows > 36000 ||
     ncols > 36000 ||
     nPixels > 1.e9){
    cerr << "Map is too big: [" << nrows << ", " << ncols << "]." << endl;
    exit(1);
  }

  //matrices of absolute coordinates
  xCoordsAbs.resize(nrows,ncols);
  yCoordsAbs.resize(nrows,ncols);
  for(int i=0;i<nrows;i++)
    for(int j=0;j<ncols;j++){
      physToAbs(&rowCoordsPhys[i], &colCoordsPhys[j], 
		&masterGrid[0], &masterGrid[1],
		&xCoordsAbs[i][j], &yCoordsAbs[i][j], 1);
    }

  //allocate space for the maps
  MatDoub wtt;
  cerr << "Coaddition(): nrows=" << nrows << ", ncols=" << ncols << endl;
  weight = new Map(string("weight"), nrows, ncols, pixelSize, 
		   wtt, rowCoordsPhys, colCoordsPhys);
  signal = new Map(string("signal"), nrows, ncols, pixelSize, 
		   weight->image, rowCoordsPhys, colCoordsPhys);
  kernel = new Map(string("kernel"), nrows, ncols, pixelSize, 
		   weight->image, rowCoordsPhys, colCoordsPhys);
  inttime = new Map(string("inttime"), nrows, ncols, pixelSize, 
		   weight->image, rowCoordsPhys, colCoordsPhys);


  //if the wiener filter is requested then make space for the filtered maps too
  if(ap->getApplyWienerFilter()){
    filteredWeight = new Map(string("filteredWeight"), nrows, ncols, pixelSize, 
		     wtt, rowCoordsPhys, colCoordsPhys);
    filteredSignal = new Map(string("filteredSignal"), nrows, ncols, pixelSize, 
		     weight->image, rowCoordsPhys, colCoordsPhys);
    filteredKernel = new Map(string("filteredKernel"), nrows, ncols, pixelSize, 
		     weight->image, rowCoordsPhys, colCoordsPhys);
  }

  //now do a simpleminded coaddition
  for(int k=0;k<nFiles;k++){
    NcFile ncfid = NcFile(ap->getMapFileList(k).c_str(), NcFile::ReadOnly);
    double onrows = ncfid.get_dim("nrows")->size();
    double oncols = ncfid.get_dim("ncols")->size();
    NcVar* rcpv = ncfid.get_var("rowCoordsPhys");
    NcVar* ccpv = ncfid.get_var("colCoordsPhys");
    VecDoub rcp(onrows);
    VecDoub ccp(oncols);
    for(int i=0;i<onrows;i++) rcp[i] = rcpv->as_double(i);
    for(int j=0;j<oncols;j++) ccp[j] = ccpv->as_double(j);
    NcVar* osVar = ncfid.get_var("signal");
    NcVar* owVar = ncfid.get_var("weight");
    NcVar* okVar = ncfid.get_var("kernel");
    MatDoub os(onrows,oncols);
    MatDoub ow(onrows,oncols);
    MatDoub ok(onrows,oncols);
    osVar->get(&os[0][0], onrows, oncols);
    owVar->get(&ow[0][0], onrows, oncols);
    okVar->get(&ok[0][0], onrows, oncols);

    //the index deltas
    int deltai = (rcp[0]-rowCoordsPhys[0])/pixelSize;
    int deltaj = (ccp[0]-colCoordsPhys[0])/pixelSize;

    //now loop through the observation maps
    for(int oi=0;oi<onrows;oi++){
      for(int oj=0;oj<oncols;oj++){
	weight->image[oi+deltai][oj+deltaj] += ow[oi][oj];
	signal->image[oi+deltai][oj+deltaj] += ow[oi][oj]*os[oi][oj];
	kernel->image[oi+deltai][oj+deltaj] += ow[oi][oj]*ok[oi][oj];
	inttime->image[oi+deltai][oj+deltaj] += ow[oi][oj]*1./64.;
      }
    }
  }

  //normalization
  for(int i=0;i<nrows;i++)
    for(int j=0;j<ncols;j++){
	signal->image[i][j] = (weight->image[i][j] != 0.) ? 
	  signal->image[i][j]/weight->image[i][j] : 0.;
	kernel->image[i][j] = (weight->image[i][j] != 0.) ? 
	  kernel->image[i][j]/weight->image[i][j] : 0.;
	inttime->image[i][j] = (weight->image[i][j] != 0.) ? 
	  inttime->image[i][j]/weight->image[i][j] : 0.;

	signal->weight[i][j] = weight->image[i][j];
	kernel->weight[i][j] = weight->image[i][j];
	inttime->weight[i][j] = weight->image[i][j];
    }

  signal->calcMapPsd(0.75);

  coaddAvgTau = mean(individualMapsTau);
  //weight->calcMapPsd(ap->getCovCut());

  return 1;
}

bool Coaddition::coaddTimeStreams(){
	//For the time being we assume that map coaddition was already done, so the pixelization grid has already been defined
	size_t nFiles = ap->getNFiles();

	//Use map grid already set in the signal map. If we decide this approach is the right way to go then we need to recalculate map bounds.

    MatDoub wtt(nrows,ncols,0.0);
	cerr << "Coaddition(): Timestream domain coaddition nrows=" << nrows << ", ncols=" << ncols << endl;
	tWeight = new Map(string("weight"), nrows, ncols, pixelSize,
			  wtt, rowCoordsPhys, colCoordsPhys);
	tSignal = new Map(string("signal"), nrows, ncols, pixelSize,
			  weight->image, rowCoordsPhys, colCoordsPhys);
	tKernel = new Map(string("kernel"), nrows, ncols, pixelSize,
			  weight->image, rowCoordsPhys, colCoordsPhys);

	double sx, kx;
	//Now loop into the files
	for(uint i=0;i<nFiles;i++){
	    NcFile ncfid = NcFile(ap->getMapFileList(i).c_str(), NcFile::ReadOnly);
	    NcDim* dimDetectors = ncfid.get_dim("nDetectors");
	    NcDim* dimSamples = ncfid.get_dim("nSamples");
	    NcDim* dimTypes = ncfid.get_dim("types");
	    NcDim* dimScans = ncfid.get_dim("nScans");

	    size_t nDetectors = dimDetectors->size();
	    size_t nSamples = dimSamples->size();
	    size_t bTypes = dimTypes->size();
	    size_t nScans = dimScans->size();

	    //Get data and scans
	    Mat3DDoub fullData (bTypes,nDetectors, nSamples);
	    MatInt scanInfo (2,nScans);

	    NcVar *fData = ncfid.get_var("boloData");
	    NcVar *fScan = ncfid.get_var("scanIndex");

	    if (!fData->get(&fullData[0][0][0],bTypes,nDetectors,nSamples)){
	    	cerr<<"Could not retreive timestream data from file: "<<ap->getMapFileList(i).c_str()<<endl;
	    	exit(-1);
	    }
	    if (!fScan->get(&scanInfo[0][0],2,nScans)){
	    	cerr<<"Could not retreive timestream scan info from file: "<<ap->getMapFileList(i).c_str()<<endl;
	    	exit(-1);
	    }

	    size_t si, ei, scanSamples;
	    double scanWeight, scanMean;

	    int irow, icol;
	    for (size_t iscan =0; iscan < nScans; iscan ++){
	    	si = scanInfo[0][iscan];
	    	ei = scanInfo[1][iscan]+1;
	    	scanSamples = ei-si;

	    	for (size_t ibolo =0; ibolo < nDetectors; ibolo++){
	    		scanMean = mean(&fullData[0][ibolo][si],scanSamples);
	    		scanWeight = 1.0/pow(stddev(&fullData[0][ibolo][si], scanSamples,scanMean),2.0);

	    		//Now find where each sample goes

	    		for (size_t isample=0; isample < scanSamples; isample++){
	    			  if (fullData[1][ibolo][si+isample]){
	    				  tWeight->raDecPhysToIndex(fullData[2][ibolo][si+isample],
	    						  fullData[3][ibolo][si+isample],
	    						   &irow, &icol);

	    				  sx = -1.0*fullData[0][ibolo][si+isample];
	    				  kx = fullData[4][ibolo][si+isample];
	    				  if (sx != sx || kx != kx)
	    					  continue;
	    				  tWeight->image[irow][icol]+= scanWeight;
	    				  tSignal->image[irow][icol]+= sx*scanWeight;
	    				  tKernel->image[irow][icol]+= kx*scanWeight;
	    			  }
	    		}

	    	}
	    }


	}

	//All files done, then normalize

	for (int i=0; i<nrows; i++)
		for (int j=0; j<ncols; j++)
			if (tWeight->image[i][j] != 0.0){
				tSignal->image[i][j]/=tWeight->image[i][j];
				tKernel->image[i][j]/=tWeight->image[i][j];
			}else{
				tSignal->image[i][j]=0.0;
				tKernel->image[i][j]=0.0;
			}


	return 1;
}


//----------------------------- o ---------------------------------------

///writes a set of coadded maps to the ncdf file
/** Writes the coadded images and many attributes to a netcdf file
    specified by AnalParams.
**/
bool Coaddition::writeCoadditionToNcdf()
{
  //create the file
  NcFile ncfid = NcFile(ap->getCoaddOutFile().c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);
  if (!ncfid.is_valid()){
    cerr << "Couldn't open " << ap->getCoaddOutFile() << " for writing.";
    cerr << endl;
    exit(1);
  }

  //if we are here then we had success opening the file for writing
  ncdfFile.assign(ap->getCoaddOutFile());
  
  //create dimensions
  NcDim* rowDim = ncfid.add_dim("nrows", nrows);
  NcDim* colDim = ncfid.add_dim("ncols", ncols);

  //define variables for maps
  NcVar *signalVar = ncfid.add_var("signal", ncDouble, rowDim, colDim);
  NcVar *kernelVar = ncfid.add_var("kernel", ncDouble, rowDim, colDim);
  NcVar *weightVar = ncfid.add_var("weight", ncDouble, rowDim, colDim);
  NcVar *inttimeVar = ncfid.add_var("inttime", ncDouble, rowDim, colDim);
  NcVar *rCPhysVar = ncfid.add_var("rowCoordsPhys", ncDouble, rowDim);
  NcVar *cCPhysVar = ncfid.add_var("colCoordsPhys", ncDouble, colDim);
  NcVar *xCAbsVar = ncfid.add_var("xCoordsAbs", ncDouble, rowDim, colDim);
  NcVar *yCAbsVar = ncfid.add_var("yCoordsAbs", ncDouble, rowDim, colDim);

  if (tWeight){
	  NcVar *tsignalVar = ncfid.add_var("tSignal", ncDouble, rowDim, colDim);
	  NcVar *tkernelVar = ncfid.add_var("tKernel", ncDouble, rowDim, colDim);
	  NcVar *tweightVar = ncfid.add_var("tWeight", ncDouble, rowDim, colDim);

	  tsignalVar->put(&tSignal->image[0][0], nrows, ncols);
	  tkernelVar->put(&tKernel->image[0][0], nrows, ncols);
	  tweightVar->put(&weight->image[0][0], nrows, ncols);
  }



  //and write the maps
  signalVar->put(&signal->image[0][0], nrows, ncols);
  kernelVar->put(&kernel->image[0][0], nrows, ncols);
  weightVar->put(&weight->image[0][0], nrows, ncols);
  inttimeVar->put(&inttime->image[0][0], nrows, ncols);
  rCPhysVar->put(&rowCoordsPhys[0], nrows);
  cCPhysVar->put(&colCoordsPhys[0], ncols);
  xCAbsVar->put(&xCoordsAbs[0][0], nrows, ncols);
  yCAbsVar->put(&yCoordsAbs[0][0], nrows, ncols);


  //get the time and date of this analysis
  time_t rawtime;
  struct tm* timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  string t = asctime(timeinfo);
  t = t.substr(0,t.length()-1);

  //also log all of the analysis parameters used as global atributes
  ncfid.add_att("source",ap->getSourceName().c_str());
  ncfid.add_att("analysisDate", t.c_str());
  ncfid.add_att("despikeSigma", ap->getDespikeSigma());
  ncfid.add_att("lowpassFilterKnee", ap->getLowpassFilterKnee());
  ncfid.add_att("timeOffset", ap->getTimeOffset());
  ncfid.add_att("timeChunk", ap->getTimeChunk());
  ncfid.add_att("neigToCut", ap->getNeigToCut());
  ncfid.add_att("cutStd", ap->getCutStd());
  ncfid.add_att("pixelSize", ap->getPixelSize());
  if (ap->getOrder()>0){
	  ncfid.add_att("bsplineOrder", ap->getOrder());
	  ncfid.add_att("bsplineControlChunk", ap->getControlChunk());
	  ncfid.add_att("bsplineStripe", ap->getCleanStripe());
	  ncfid.add_att("bsplineResampling", ap->getResample());

  }
  ncfid.add_att("approximateWeights", ap->getApproximateWeights());
  ncfid.add_att("MasterGrid[0]",masterGrid[0]);
  ncfid.add_att("MasterGrid[1]",masterGrid[1]);
  ncfid.add_att("noiseMapsPerObs",ap->getNNoiseMapsPerObs());
  ncfid.add_att("azelMap",ap->getAzelMap());
  ncfid.add_att("threadNumber",ap->getNThreads());
  ncfid.add_att("noiseRealizations", ap-> getNRealizations());
  if(ap->getApplyWienerFilter()){
    ncfid.add_att("WienerFilterGaussianTemplate", ap->getGaussianTemplate());
    ncfid.add_att("WienerFilterGaussianTemplateFWHM", ap->getGaussianTemplateFWHM());
    ncfid.add_att("WienerFilterLowpassOnly", ap->getLowpassOnly());
    ncfid.add_att("WienerFilterNormalizeErrors", ap->getNormalizeErrors());
  }
  if(ap->getFindSources()){
    ncfid.add_att("SourceFindingBeamSize",ap->getBeamSize());
    ncfid.add_att("SourceFindingCovCut", ap->getCovCut());
    ncfid.add_att("SourceFindingSnglSourceWindow", ap->getSnglSourceWin());
    ncfid.add_att("SourceFindingNegativeToo", ap->getNegativeToo());
    ncfid.add_att("SourceFindingMapNegative", ap->getMapNegative());
    ncfid.add_att("SourceFindingCentroidSources",ap->getSFCentroidSources());
    ncfid.add_att("SourceFindingFitGaussians", ap->getSFFitGaussians());
  }
  if(ap->getCalcCompleteness()){
    ncfid.add_att("CompletenessNumFluxBins",ap->getNFluxBins());
    ncfid.add_att("CompletenessNumSynthSources",ap->getNSynthSources());
    ncfid.add_att("CompletenessMinFlux", ap->getMinFlux());
    ncfid.add_att("CompletenessMaxFlux", ap->getMaxFlux());
    ncfid.add_att("CompletenessRecovS2N", ap->getRecovS2N());
    ncfid.add_att("CompletenessRecovRadius", ap->getRecovRadius());
  }

  //add the list of raw data files
  string df("dataFile_");
  string d;
  int nFiles = ap->getNFiles();
  for(int i=0;i<nFiles;i++){
    stringstream o;
    o << i;
    d.assign(df);
    d.append(o.str());
    ncfid.add_att(d.c_str(),ap->getMapFileList(i).c_str());
  }

  ncfid.add_att("coaddAvgTau", coaddAvgTau);
  cerr << "Coaddition::writeMapsToNcdf(): Maps written to ";
  cerr << ncdfFile << endl;

  //update the ncdfFile for each Map
  signal->mapFile = ncdfFile;
  weight->mapFile = ncdfFile;
  kernel->mapFile = ncdfFile;


  return 1;
}


//----------------------------- o ---------------------------------------

///writes the set of coadded images to a fits file
/** this still needs to be implemented
    \todo Needs to be implemented.
**/
bool Coaddition::writeCoadditionToFits(string fitsFilename)
{


  return 1;
}


//----------------------------- o ---------------------------------------

///Creates a histogram of the coadded signal map.
/** Create a histogram of the coadded signal map using the signal
    map's built-in histogramming routine.  Inputs are:
      - nbins - the number of bins in the output histogram
      - cc - the desired coverage cut (0<cc<1)
    Alternatively you can use this function without arguments to
    choose the default values of nbins=200 and cc=ap->getCoverageThreshold().
**/
bool Coaddition::histogramSignal(int nbins, double cc)
{
  signal->calcMapHistogram(nbins,cc);
  return 1;
}

bool Coaddition::histogramSignal(double cc)
{
  int nbins = 200.;
  histogramSignal(nbins, cc);
  return 1;
}

bool Coaddition::histogramSignal()
{
  int nbins = 200.;
  double coverage_cut = ap->getCoverageThreshold();
  histogramSignal(nbins, coverage_cut);
  return 1;
}



//----------------------------- o ---------------------------------------

///Creates a histogram of the coadded signal map.
/** Create a histogram of the coadded signal map using the signal
    map's built-in histogramming routine.  Inputs are:
      - nbins - the number of bins in the output histogram
      - cc - the desired coverage cut (0<cc<1)
    Alternatively you can use this function without arguments to
    choose the default values of nbins=200 and cc=ap->getCoverageThreshold().
**/
bool Coaddition::histogramFilteredSignal(int nbins, double cc)
{
  filteredSignal->calcMapHistogram(nbins,cc);
  return 1;
}

bool Coaddition::histogramFilteredSignal(double cc)
{
  int nbins = 200.;
  histogramFilteredSignal(nbins, cc);
  return 1;
}

bool Coaddition::histogramFilteredSignal()
{
  int nbins = 200.;
  double coverage_cut = ap->getCoverageThreshold();
  histogramFilteredSignal(nbins, coverage_cut);
  return 1;
}

//----------------------------- o ---------------------------------------

///writes Wiener filtered maps to netcdf file
/** Writes the filtered coadded images to the coadded netcdf file
    specified by AnalParams.
**/
bool Coaddition::writeFilteredMapsToNcdf()
{
  //open the file
  NcFile ncfid = NcFile(ap->getCoaddOutFile().c_str(), NcFile::Write);
  if (!ncfid.is_valid()){
    cerr << "Couldn't open " << ap->getCoaddOutFile() << " for writing.";
    cerr << endl;
    exit(1);
  }

  //create dimensions
  NcDim* rowDim = ncfid.get_dim("nrows");
  NcDim* colDim = ncfid.get_dim("ncols");

  //define variables for maps
  NcVar *signalVar = ncfid.add_var("filteredSignal", ncDouble, rowDim, colDim);
  NcVar *kernelVar = ncfid.add_var("filteredKernel", ncDouble, rowDim, colDim);
  NcVar *weightVar = ncfid.add_var("filteredWeight", ncDouble, rowDim, colDim);

  //and write the maps
  signalVar->put(&filteredSignal->image[0][0], nrows, ncols);
  kernelVar->put(&filteredKernel->image[0][0], nrows, ncols);
  weightVar->put(&filteredWeight->image[0][0], nrows, ncols);

  cerr << "Coaddition::writeFilteredMapsToNcdf(): Maps written to ";
  cerr << ncdfFile << endl;

  //update the ncdfFile for each Map
  filteredSignal->mapFile = ncdfFile;
  filteredWeight->mapFile = ncdfFile;
  filteredKernel->mapFile = ncdfFile;

  return 1;
}


//----------------------------- o ---------------------------------------

///calls findSources() for completeness simulations
/** Overload of findSources() for completeness simulations.
    Allows us to feed in the signal and S/N maps with synthetic
    source injected to be searched for sources. Also changes
    a couple parameters to comply with analysis parameters
    .xml file for the simulations. Inputs are:
      - signalMap - signal map MatDoub plus injected synthetic source
      - s2nMap - S/N map MatDoub plus injected synthetic source contribution
      - realCoadd - pointer to parent Coaddition class
    \todo - add S/N map creation when we are actually able to pass it here
**/
bool Coaddition::findSources(double* signalMap, double* s2nMap,
			     Coaddition* realCoadd, double maxS2N)
{
  //needs to go when Wiener filter is done
  maxPreS2N = maxS2N;

  //set nrows and ncols
  nrows = realCoadd->getNrows();
  ncols = realCoadd->getNcols();

  //make signal and weight Map objects
  MatDoub wtt;
  filteredWeight = new Map(string("filteredWeight"), nrows, ncols, pixelSize,
			   wtt, realCoadd->rowCoordsPhys,
			   realCoadd->colCoordsPhys);
  filteredSignal = new Map(string("filteredSignal"), nrows, ncols, pixelSize,
			   filteredWeight->image,
			   realCoadd->rowCoordsPhys,
			   realCoadd->colCoordsPhys);
  filteredKernel = new Map(string("filteredKernel"), nrows, ncols, pixelSize,
			   filteredWeight->image,
			   realCoadd->rowCoordsPhys,
			   realCoadd->colCoordsPhys);
  //make signal and S/N maps into ones containing synthetic source
  //set up coverage boolean map
  rowCoordsPhys.resize(nrows);
  colCoordsPhys.resize(ncols);
  xCoordsAbs.resize(nrows, ncols);
  yCoordsAbs.resize(nrows, ncols);
  filteredSignal->coverageBool.resize(nrows, ncols);
  for(int i=0;i<nrows;i++){
    for(int j=0;j<ncols;j++){
       filteredSignal->image[i][j] = signalMap[ncols*i + j];
      //something here for S/N map
      filteredWeight->image[i][j] = realCoadd->filteredWeight->image[i][j];
      filteredSignal->coverageBool[i][j] =
	realCoadd->filteredSignal->coverageBool[i][j];
      xCoordsAbs[i][j] = realCoadd->getXCoordsAbs(i, j);
      yCoordsAbs[i][j] = realCoadd->getYCoordsAbs(i , j);
    }
  }
  for(int i=0;i<nrows;i++){
    for(int j=0;j<ncols;j++){
      filteredSignal->weight[i][j] = realCoadd->filteredWeight->image[i][j];
    }
  }
  for(int i=0;i<nrows;i++){
    rowCoordsPhys[i] = realCoadd->getRowCoordsPhys(i);
  }
  for(int i=0;i<ncols;i++){
    colCoordsPhys[i] = realCoadd->getColCoordsPhys(i);
  }

  //use the recovery S/N threshold
  sourceSigma = ap->getRecovS2N();

  //set synthetic source finding flag
  synth = true;

  //do the source finding
  findSources();

  return 1;
}




//----------------------------- o ---------------------------------------
///finds the potential sources in a coadded map
/** Searches the coadded S/N map for sources with S/N greater
    than or equal to the supplied source sigma. Sources are
    masked to supplied snglSourceWin so sources are not counted
    multiple times. All source locations, fluxes, S/N, etc.
    are stored in an array of PointSource objects.
    \todo - get convolving S/N stuff replaced with Weiner filter results
**/
bool Coaddition::findSources()
{
  //grab appropriate switches
  bool negativeToo = ap->getNegativeToo();
  bool mapNegative = ap->getMapNegative();

  //get good coverage region with Map coverage cut methods
  filteredSignal->setCoverageCut(ap->getCovCut());
  filteredSignal->findWeightThresh();
  filteredSignal->makeCovBoolMap();

  //swap the sign of the signal map if analysis calls for the map negative
  if(mapNegative){
    for(int i=0;i<nrows;i++) 
      for(int j=0;j<ncols;j++){
	filteredSignal->image[i][j] = -filteredSignal->image[i][j];
      }
  }
  VecDoub kpp;
  string resFile = filteredKernel->mapFile;
  filteredKernel->mapFile="";
  filteredKernel->fitToGaussian(kpp);
  filteredKernel->mapFile=resFile;

  //make S/N map and find maximum S/N in good coverage region
  double tempMax = sqrt(filteredWeight->image[0][0])*
                   filteredSignal->image[0][0];
  MatDoub signal2Noise(nrows, ncols);
  signal2Noise[0][0] = tempMax;
  for(int i=1;i<nrows;i++){
    for(int j=1;j<ncols;j++){
      signal2Noise[i][j] = sqrt(filteredWeight->image[i][j])*
	filteredSignal->image[i][j];
      if(filteredSignal->coverageBool[i][j]){
	if(signal2Noise[i][j] > tempMax){
	  tempMax = signal2Noise[i][j];
	}
      }
    }
  }
  if(!synth){
    maxPreS2N = tempMax;
  }

  //find pixels equal or above source sigma
  vector<int> rInd;
  vector<int> cInd;
  if(negativeToo){
    for(int i=0;i<nrows;i++){
      for(int j=0;j<ncols;j++){
    	if(filteredSignal->coverageBool[i][j] == 1){
	  if(abs(signal2Noise[i][j]) >= sourceSigma){
	      rInd.push_back(i);
	      cInd.push_back(j);
	  }
	}
      }
    }
  }
  else{
    for(int i=0;i<nrows;i++){
      for(int j=0;j<ncols;j++){
        if(filteredSignal->coverageBool[i][j] == 1){
	  if(signal2Noise[i][j] >= sourceSigma){
	      rInd.push_back(i);
	      cInd.push_back(j);
	  }
	}
      }
    }
  }

  //did we find any sources?
  if (rInd.size() == 0){
    cerr << "Coaddition::findSources(): ";
    cerr << "No sources found." << endl;
    sources = NULL;
    return 1;
  }

  //make sure source extremum is within good coverage region by
  //searching in index boxes of +/- 1 pixel around hot pixels
  vector<int> rSourceInd;
  vector<int> cSourceInd;
  for(unsigned int i=0;i<rInd.size();i++){
    double extremum;
    if(negativeToo && filteredSignal->image[rInd[i]][cInd[i]] < 0.0){
      extremum = filteredSignal->image[rInd[i]][cInd[i]];
      //find minimum within index box
      for(int j=rInd[i]-1;j<rInd[i]+2;j++){
	for(int k=cInd[i]-1;k<cInd[i]+2;k++){
	  if(filteredSignal->image[j][k] < extremum){
	    extremum = filteredSignal->image[j][k];
          }
        }
      }
    }
    else{
      extremum = filteredSignal->image[rInd[i]][cInd[i]];
      //find maximum within index box
      for(int j=rInd[i]-1;j<rInd[i]+2;j++){
	for(int k=cInd[i]-1;k<cInd[i]+2;k++){
	  if(filteredSignal->image[j][k] > extremum){
	    extremum = filteredSignal->image[j][k];
          }
        }
      }
    }
    //only keep the hot pixel if it is the extremum
    if (filteredSignal->image[rInd[i]][cInd[i]] == extremum){
      rSourceInd.push_back(rInd[i]);
      cSourceInd.push_back(cInd[i]);
    }
  }
  int nRawSources = rSourceInd.size();
  //done with vectors of *all* hot pixels
  rInd.clear();
  cInd.clear();

  //do we still have source candidates left?
  if (nRawSources == 0){
    cerr << "Coaddition::findSources(): ";
    cerr << "No sources found." << endl;
    sources = NULL;
    return 1;
  }

  //find indices of hot pixels close together
  vector<int> rDistInd;
  vector<int> cDistInd;
  for (int i=0;i<nRawSources;i++){
    for (int j=0;j<nRawSources;j++){
      unsigned int rSep = gsl_pow_2(rSourceInd[i] - rSourceInd[j]);
      unsigned int cSep = gsl_pow_2(cSourceInd[i] - cSourceInd[j]);
      double hotDist = sqrt(rSep + cSep);
      if (hotDist <= (snglSourceWin/pixelSize) &&
          hotDist != 0.0){
        rDistInd.push_back(i);
        cDistInd.push_back(j);
      }
    }
  }

  //flag non-maximum hot pixel indices
  if(rDistInd.size() != 0){
    for(unsigned int i=0;i<rDistInd.size();i++){
      if(rSourceInd[rDistInd[i]] == -1 ||
	 cSourceInd[cDistInd[i]] == -1){
          continue;
      }
      double f1 =
	filteredSignal->image[rSourceInd[rDistInd[i]]][cSourceInd[rDistInd[i]]];
      double f2 =
	filteredSignal->image[rSourceInd[cDistInd[i]]][cSourceInd[cDistInd[i]]];
      //determine if same sign and which sign
      if(f1 < 0.0 && f2 < 0.0){
	//negative case
	if(f1 <= f2){
	  rSourceInd[cDistInd[i]] = -1;
	  cSourceInd[cDistInd[i]] = -1;
	}
	else{
	  rSourceInd[rDistInd[i]] = -1;
	  cSourceInd[rDistInd[i]] = -1;
	}
      }
      else{
	//positive case
	if(f1 >= f2){
	  rSourceInd[cDistInd[i]] = -1;
	  cSourceInd[cDistInd[i]] = -1;
	}
	else{
	  rSourceInd[rDistInd[i]] = -1;
	  cSourceInd[rDistInd[i]] = -1;
	}
      }
    }
  }
  rDistInd.clear();
  cDistInd.clear();

  vector<int> rSourceLoc;
  vector<int> cSourceLoc;
  nSources = 0;
  for(int i=0;i<nRawSources;i++){
    if((rSourceInd[i] != -1)
       && (cSourceInd[i] != -1)){
      rSourceLoc.push_back(rSourceInd[i]);
      cSourceLoc.push_back(cSourceInd[i]);
      nSources++;
    }
  }
  //done with flag filled source index vectors
  rSourceInd.clear();
  cSourceInd.clear();

  //do we have any sources?
  if(nSources == 0){
    cerr << "SourceLocate::findSources(): ";
    cerr << "No sources found." << endl;
    sources = NULL;
    return 1;
  }

  //store sources information
  sources = new PointSource[nSources];
  for(int i=0;i<nSources;i++){
    sources[i].sID = i;
    sources[i].nSourcesParentMap = nSources;
    sources[i].centerRaAbs = xCoordsAbs[rSourceLoc[i]][cSourceLoc[i]];
    sources[i].centerDecAbs = yCoordsAbs[rSourceLoc[i]][cSourceLoc[i]];
    sources[i].centerRaPhys = rowCoordsPhys[rSourceLoc[i]];
    sources[i].centerDecPhys = colCoordsPhys[cSourceLoc[i]];
    sources[i].centerXPos = rSourceLoc[i];
    sources[i].centerYPos = cSourceLoc[i];
    sources[i].centerFlux = filteredSignal->image[rSourceLoc[i]][cSourceLoc[i]];
    sources[i].centerNoise = filteredSignal->image[rSourceLoc[i]][cSourceLoc[i]]/
                             signal2Noise[rSourceLoc[i]][cSourceLoc[i]];
    sources[i].centerS2N = signal2Noise[rSourceLoc[i]][cSourceLoc[i]];

    //run PointSource initializers
    sources[i].initialize(ap, &filteredSignal->image[0][0],
			  &filteredSignal->weight[0][0],
			  &rowCoordsPhys[0],
			  &colCoordsPhys[0],
			  &xCoordsAbs[0][0],
			  &yCoordsAbs[0][0],
			  nrows, ncols,
			  pixelSize, ap->getCoaddOutFile().c_str());

    //make the postage stamps
    sources[i].makePostageStamp();

    //centroid sources
    sources[i].centroidSource();

    //fit to gaussian
    cout<<"Using fwhm "<< kpp[2] <<"x" <<kpp[3]<<endl;
    sources[i].fitGaussianToSource(kpp[2], kpp[3]);

    //add the source to the netcdf file
    sources[i].addSourceToNCDF();
  }

  //fix signal map if we searched its negative
  if(mapNegative){
    for(int i=0;i<nrows;i++){
      for(int j=0;j<ncols;j++){
	filteredSignal->image[i][j] = -filteredSignal->image[i][j];
      }
    }
  }

  return 1;
}



//----------------------------- o ---------------------------------------



bool Coaddition::normalizeErrors(double noiseRms, double cov)
{
  //work with stored weight map
  MatDoub myWeight(nrows,ncols,0.);
  for(int i=0;i<nrows;i++)
    for(int j=0;j<ncols;j++)
      myWeight[i][j] = filteredWeight->image[i][j];

  //find the weight threshold
  double myWeightCut = findWeightThreshold(myWeight,cov);

  //following IDL nomenclature
  double mean_sqerr=0.;
  int counter=0;
    for(int i=0;i<nrows;i++)
      for(int j=0;j<ncols;j++){
	if(myWeight[i][j] >= myWeightCut){
	  mean_sqerr += (1./myWeight[i][j]);
	  counter++;
	}
      }
    mean_sqerr /= counter;

    //the renormalization factor
    double nfac = (1./pow(noiseRms,2.))*mean_sqerr;
    cerr << "The renormalization factor is: " << nfac << endl;

    for(int i=0;i<nrows;i++)
      for(int j=0;j<ncols;j++)
	filteredWeight->image[i][j] *= nfac;

    //make sure that the filtered signal map has the right weight associated with it
    for(int i=0;i<nrows;i++)
      for(int j=0;j<ncols;j++)
	filteredSignal->weight[i][j] = filteredWeight->image[i][j];

  return 1;
}





//----------------------------- o ---------------------------------------

int Coaddition::getNrows()
{
  return nrows;
}

int Coaddition::getNcols()
{
  return ncols;
}

double Coaddition::getPixelSize()
{
  return pixelSize;
}

double Coaddition::getRowCoordsPhys(int i)
{
  return rowCoordsPhys[i];
}

double Coaddition::getColCoordsPhys(int i)
{
  return colCoordsPhys[i];
}

double Coaddition::getXCoordsAbs(int i, int j)
{
  return xCoordsAbs[i][j];
}

double Coaddition::getYCoordsAbs(int i, int j)
{
  return yCoordsAbs[i][j];
}

int Coaddition::getNSources()
{
  return nSources;
}


//----------------------------- o ---------------------------------------

Coaddition::~Coaddition()
{

  delete signal;
  delete kernel;
  delete weight;

  
  /*
  if(!synth){
    delete signal;
    delete kernel;
    delete weight;
  }
  delete filteredSignal;
  delete filteredKernel;
  delete filteredWeight;
  if(sources != NULL){
    delete[] sources;
  }
  */
}
