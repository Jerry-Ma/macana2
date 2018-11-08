#include <netcdfcpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <cstdio>
#include <time.h>
#include <fftw3.h>
#include <CCfits/CCfits>
using namespace std;

#include "nr3.h"
#include "AnalParams.h"
#include "Array.h"
#include "astron_utilities.h"
#include "gaussFit.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include "Observation.h"
#include "Telescope.h"
#include "vector_utilities.h"
#include "intarray2bmp.h"


///Observation constructor
/** Just simple initialization going on here.
 **/
Observation::Observation(AnalParams* analParams)
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

	//initialize nrows and ncols
	nrows=0;
	ncols=0;

	atmTemplate = NULL;
        signal = NULL;
        kernel = NULL;
        weight = NULL;
        inttime = NULL; 
	saveTimestreams = ap->getSaveTimestreams();
}


//----------------------------- o ---------------------------------------
///this method generates the beammap for each detector in the passed in Array
/*These beammaps are stored in the beammapSignal and beammapWeight MatDoub
One dimension is the detector number, and the other is the signal or weight values
stored row after row
*/

bool Observation::generateBeammaps(Array* a, Telescope* tel)
{
  this->tel = tel;
  //set up coordinate systems and map size, and intialize maps
  mapGenerationPrep(a);

  a->updateDetectorIndices();
  int* di = a->getDetectorIndices();

  //grab the weights
  MatDoub tmpwt = calculateWeights(a, tel);

  int nDetectors = a->getNDetectors();
  int nScans = tel->scanIndex.ncols();

  MatDoub wtt(nrows, ncols, 0.);
  string mname;
  weight = new Map(mname.assign("w"), nrows, ncols, pixelSize, wtt, rowCoordsPhys, colCoordsPhys);

  //allocate space for the maps
  beammapSignal.resize(nDetectors, nrows*ncols);
  beammapWeight.resize(nDetectors, nrows*ncols);
  beammapIntTime.resize(nDetectors, nrows*ncols);

  for(int i=0;i<nDetectors;i++){
    for(int j=0;j<nrows*ncols;j++){
      beammapSignal[i][j] = 0;
      beammapWeight[i][j] = 0;
    }
  }

  //for each detector, calculate signal and weight values
  for(int i=0;i<nDetectors;i++){
    MatInt times(nrows, ncols, 0.);
    for(int k=0;k<nScans;k++){
      int si=tel->scanIndex[0][k];
      int ei=tel->scanIndex[1][k]+1;
      for(int j=si;j<ei;j++){
        if(a->detectors[di[i]].hSampleFlags[j]){
          int irow;
          int icol;
          weight->raDecPhysToIndex(a->detectors[di[i]].hRa[j],
            a->detectors[di[i]].hDec[j], &irow, &icol);

          beammapWeight[i][irow*ncols+icol] += tmpwt[i][k];
          //cout << beammapWeight[i][irow*ncols+icol] << endl;

          double hx = a->detectors[di[i]].hValues[j];
          if(hx != hx){
            cerr << "NaN detected on file: " << ap->getMapFile() << endl;
            cerr << "tmpwt: " << tmpwt[i][k] << endl;
            cerr << "det: " << a->detectors[di[i]].hValues[j] << endl;
            cerr << "i: " << i << endl;
            cerr << "j: " << j << endl;
            cerr << "k: " << k << endl;
          }

          beammapSignal[i][irow*ncols+icol] += hx;
          beammapIntTime[i][irow*ncols+icol] += 1./64.;
          times[irow][icol]++;
        }
      }
    }

    for(int x=0;x<nrows;x++) for(int y=0;y<ncols;y++){
      int rep = times[x][y];
      if(rep != 0){
        beammapSignal[i][x*ncols + y] = beammapSignal[i][x*ncols + y]/rep;
        beammapSignal[i][x*ncols + y] = -beammapSignal[i][x*ncols + y];
      }
      else{
        beammapSignal[i][x*ncols + y] = 0.;
      }
    }
  }
  return 1;
}


//----------------------------- o ---------------------------------------

///generates coordinates, signal, weight, kernel and noise maps
/**This is the workhorse of Observation.  The signal, weight, kernel
   and five noise maps are all generated here taking an Array and
   Telescope objects as input.  For the point source maps I've tested
   there are some very minor differences between the results generated
   here and those of the idl utilities.  These need to be investigated.
   \todo Find root of minor differences between these maps and idl utils.
 **/
bool Observation::generateMaps(Array* a, Telescope* tel)
{
  this->tel = tel;

  //set up coordinate systems and map size and initialize maps
  mapGenerationPrep(a);

  //allocate space for the maps
  cerr << "Observation(): nrows=" << nrows << ", ncols=" << ncols << endl;
  MatDoub wtt(nrows,ncols,0.);
  string mname;
  weight = new Map(mname.assign("weight"), nrows, ncols, pixelSize,
	   wtt, rowCoordsPhys, colCoordsPhys);
  signal = new Map(mname.assign("signal"), nrows, ncols, pixelSize,
	   weight->image, rowCoordsPhys, colCoordsPhys);
  kernel = new Map(mname.assign("kernel"), nrows, ncols, pixelSize,
	   weight->image, rowCoordsPhys, colCoordsPhys);
  inttime = new Map(mname.assign("inttime"), nrows, ncols, pixelSize,
	   weight->image, rowCoordsPhys, colCoordsPhys);

  cerr << "making noiseMaps with " << ap->getNNoiseMapsPerObs() << " maps." << endl;

  noiseMaps.resize(ap->getNNoiseMapsPerObs(),nrows*ncols);
  for(int k=0;k<ap->getNNoiseMapsPerObs();k++)
    for(int i=0;i<nrows;i++) for(int j=0;j<ncols;j++){
      noiseMaps[k][ncols*i+j]=0.;
    }

    //populate the maps
    a->updateDetectorIndices();
    int* di=a->getDetectorIndices();
    int nDetectors = a->getNDetectors();
    int nScans = tel->scanIndex.ncols();

    if (a->detectors[di[0]].atmTemplate.size() > 0.0){
      atmTemplate = new Map(mname.assign("atmTemplate"), nrows, ncols, pixelSize,
      weight->image, rowCoordsPhys, colCoordsPhys);
    }

    //calculate or set the weights
    MatDoub tmpwt = calculateWeights(a, tel);

    //put together the maps
    for(int k=0;k<nScans;k++){
      VecDoub sn(ap->getNNoiseMapsPerObs());
    for(int kk=0;kk<ap->getNNoiseMapsPerObs();kk++){
	      sn[kk] = (ap->macanaRandom->uniformDeviate(-1.,1.)<0) ? -1. : 1.;
	    }
	    int si=tel->scanIndex[0][k];
	    int ei=tel->scanIndex[1][k]+1;
	    for(int i=0;i<nDetectors;i++){
	      for(int j=si;j<ei;j++){
		if(a->detectors[di[i]].hSampleFlags[j]){
		  //get the row and column index corresponding to the ra and
		  //dec
		  int irow;
		  int icol;
		  weight->raDecPhysToIndex(a->detectors[di[i]].hRa[j],
					   a->detectors[di[i]].hDec[j],
					   &irow, &icol);
		  //weight map
		  weight->image[irow][icol] += tmpwt[i][k];

		  //inttime map
		  inttime->image[irow][icol] += 1./64.;

		  //check for NaN
		  double hx = tmpwt[i][k]*a->detectors[di[i]].hValues[j];
		  double hk = tmpwt[i][k]*a->detectors[di[i]].hKernel[j];

		  double ha = 0.0;
		  if (atmTemplate)
			  ha = tmpwt[i][k]*a->detectors[di[i]].atmTemplate[j];
		  if(hx != hx || hk != hk){
		    cerr << "NaN detected on file: "<<ap->getMapFile() << endl;
		    cerr << "tmpwt: " << tmpwt[i][k] << endl;
		    cerr << "det: " << a->detectors[di[i]].hValues[j] << endl;
		    cerr << "ker: " << a->detectors[di[i]].hKernel[j] << endl;
		    cerr << "  i=" << i << endl;
		    cerr << "  j=" << j << endl;
		    cerr << "  k=" << k << endl;
		    exit(1);
		  }

		  //signal map
		  signal->image[irow][icol] += hx;

		  //kernel map
		  kernel->image[irow][icol] += hk;

		  //noise maps
		  for(int kk=0;kk<ap->getNNoiseMapsPerObs();kk++){
		    noiseMaps[kk][irow*ncols+icol] += sn[kk]*hx;
		  }

		  if (atmTemplate)
			  atmTemplate->image[irow][icol] += ha;
		}
	      }
	    }
	  }

	  //some maps need weight normalization
	  //also invert sign of signal map
	  double wt;

	  double atmpix=0;
	  for(int i=0;i<nrows;i++) for(int j=0;j<ncols;j++){
	      wt = weight->image[i][j];
	      if(wt != 0.){
	    	  //if (atmTemplate)
	    		  //atmpix = atmTemplate->image[i][j];
		signal->image[i][j] = -(signal->image[i][j]-atmpix)/wt;
		signal->weight[i][j] = wt;
		kernel->image[i][j] = kernel->image[i][j]/wt;
		kernel->weight[i][j] = wt;

		for(int kk=0;kk<ap->getNNoiseMapsPerObs();kk++){
		  noiseMaps[kk][i*ncols+j] = noiseMaps[kk][i*ncols+j]/wt;
		}
	      }else{
		signal->image[i][j] = 0.;
		kernel->image[i][j] = 0.;
		for(int kk=0;kk<ap->getNNoiseMapsPerObs();kk++){
		  noiseMaps[kk][i*ncols+j] = 0.;
		}
	      }
	    }

	  return 1;
}


//----------------------------- o ---------------------------------------

///Prepares for map generation. Sets map size and prepares coordinate vectors

bool Observation::mapGenerationPrep(Array* a)
{
    array = a;
    //start with the map dimensions but copy IDL utilities
    //for aligning the pixels in various maps
    double y_min_act = a->getMinY();
    double y_max_act = a->getMaxY();
    double x_min_act = a->getMinX();
    double x_max_act = a->getMaxX();

    //explicitly center on mastergrid
    //require that maps are even dimensioned with a few extra pixels.
    int xminpix = ceil(abs(x_min_act/pixelSize));
    int xmaxpix = ceil(abs(x_max_act/pixelSize));
    xmaxpix = max(xminpix,xmaxpix);
    nrows = 2.*xmaxpix+4;
    int yminpix = ceil(abs(y_min_act/pixelSize));
    int ymaxpix = ceil(abs(y_max_act/pixelSize));
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
      cerr << "Observation: Map is too big: [" << nrows << ", " << ncols << "]." << endl;
      cerr << ap->getDataFile() << endl;
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

    return 1;
}

//----------------------------- o ---------------------------------------

///generates the weight map, depending on whether or not getApproximateWeights is set. 
MatDoub Observation::calculateWeights(Array* a, Telescope* tel)
{
    int nDetectors = a->getNDetectors();
    int nScans = tel->scanIndex.ncols();
    int* di = a->getDetectorIndices();
    MatDoub tmpwt(nDetectors, nScans);
    if(ap->getApproximateWeights()){
      double sens;
      double samprate;
      for(int i=0;i<nDetectors;i++)
        for(int j=0;j<nScans;j++){
          sens = a->detectors[di[i]].getSensitivity();
          samprate = a->detectors[di[i]].getSamplerate();
          tmpwt[i][j] = pow(sqrt(samprate)*sens*1.e-3,-2);
        }
      //check for detectors with stupidly high sensitivities
      double medweight;
      double stdweight;
      VecDoub tmpdet(nDetectors);
      for(int i=0;i<nDetectors;i++) tmpdet[i] = tmpwt[i][0];
      medweight = median(&tmpdet[0],nDetectors);
      stdweight = stddev(tmpdet);
      for(int i=0;i<nDetectors;i++){
        if (tmpdet[i] > 3.*stdweight + medweight){
          for(int j=0;j<nScans;j++){
            tmpwt[i][j] = medweight;
          }
        }
      }
    } else {
      for(int i=0;i<nDetectors;i++)
        for(int j=0;j<nScans;j++){
          tmpwt[i][j] = a->detectors[di[i]].scanWeight[j];
        }
    }

    return tmpwt;
}


bool Observation::writeBeammapsToFits(string ncdfFilename)
{

    // Create a FITS primary array containing a 2-D image               
    // declare axis arrays.    
    long naxis    =   2;      
    long naxes[2] = { ncols, nrows };   
    vector<long> naxes_(naxes, end(naxes));
    
    // declare auto-pointer to FITS at function scope. Ensures no resources
    // leaked if something fails in dynamic allocation.
    unique_ptr<CCfits::FITS> pFits;
      
    try
    {                
        // overwrite existing file if the file already exists.
            
        const string fileName(ncdfFilename);            
        
        // Create a new FITS object, specifying the data type and axes for the primary
        // image. Simultaneously create the corresponding file.
        
        // this image is unsigned short data, demonstrating the cfitsio extension
        // to the FITS standard.
        
        pFits.reset( new CCfits::FITS(fileName, CCfits::Write) );
    }
    catch (CCfits::FITS::CantCreate)
    {
          // ... or not, as the case may be.
          cerr << "Failed to open " << ncdfFilename << "... aborting" << endl;
          exit(1);

    }
    
    // references for clarity.
    
    long nelements(1); 
    
    
    // Find the total size of the array. 
    // this is a little fancier than necessary ( It's only
    // calculating naxes[0]*naxes[1]) but it demonstrates  use of the 
    // C++ standard library accumulate algorithm.
    
    nelements = accumulate(&naxes[0],&naxes[naxis],1,std::multiplies<long>());

    std::valarray<double> arr(nelements);

  string obeamsignal;
  size_t nDetectors = array->getNDetectors();
  for (size_t n=0; n < nDetectors; n++) {
    cerr << "writing signal for map " << n << endl;
    obeamsignal.assign("beammapSignal");
    stringstream o;
    o << n+1;
    obeamsignal.append(o.str());
    CCfits::ExtHDU* imageExt = pFits->addImage(obeamsignal, DOUBLE_IMG, naxes_);

    for(int i=0;i<ncols;i++)
      for(int j=0;j<nrows;j++)
        arr[nrows*i+j] = beammapSignal[n][nrows*i + j];
    //
    // I have no idea what this is for
    long fpixel(1);
    //this is what we eventually want
    imageExt->write(fpixel,nelements,arr);

  }
    return 0;
}

//----------------------------- o ---------------------------------------
//writes the entire set of beammaps to a netcdf file
bool Observation::writeBeammapsToNcdf(string ncdfFilename)
{
  NcFile ncfid = NcFile(ncdfFilename.c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);
  if (!ncfid.is_valid()){
    cerr << "Couldn't open map netcdf file for writing!\n";
    return 0;
  }

  ncdfFile.assign(ncdfFilename);

  //create dimensions
  NcDim* rowDim = ncfid.add_dim("nrows", nrows);
  NcDim* colDim = ncfid.add_dim("ncols", ncols);

  NcVar *rCPhysVar = ncfid.add_var("rowCoordsPhys", ncDouble, rowDim);
  NcVar *cCPhysVar = ncfid.add_var("colCoordsPhys", ncDouble, colDim);
  NcVar *xCAbsVar = ncfid.add_var("xCoordsAbs", ncDouble, rowDim, colDim);
  NcVar *yCAbsVar = ncfid.add_var("yCoordsAbs", ncDouble, rowDim, colDim);

  size_t nDetectors = array->getNDetectors();
  size_t nSamples = array->getNSamples();
  size_t nvars = 5;

  rCPhysVar->put(&rowCoordsPhys[0], nrows);
  cCPhysVar->put(&colCoordsPhys[0], ncols);
  xCAbsVar->put(&xCoordsAbs[0][0], nrows, ncols);
  yCAbsVar->put(&yCoordsAbs[0][0], nrows, ncols);

  /*
  NcDim* dimDetector = ncfid.add_dim("nDetectors", nDetectors);
  NcVar* testvar;
  testvar = ncfid.add_var("testvar", ncDouble, dimDetector, rowDim, colDim);
*/
  string obeamsignal;
  NcVar* nVar;
  for(size_t i=0;i<nDetectors;i++){
    cerr << "writing signal for map " << i << endl;
    obeamsignal.assign("beammapSignal");
    stringstream o;
    o << i;
    obeamsignal.append(o.str());
    nVar = ncfid.add_var(obeamsignal.c_str(), ncDouble, rowDim, colDim);
    nVar->put(&beammapSignal[i][0], nrows, ncols);
  }

  string obeamweight;
  for(size_t i=0;i<nDetectors;i++){
    cerr << "writing weight for map " << i << endl;
    obeamweight.assign("beammapWeight");
    stringstream o;
    o << i;
    obeamweight.append(o.str());
    nVar = ncfid.add_var(obeamweight.c_str(), ncDouble, rowDim, colDim);
    nVar->put(&beammapWeight[i][0], nrows, ncols);
  }

  time_t rawtime;
  struct tm* timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  string t = asctime(timeinfo);
  t = t.substr(0,t.length()-1);
  
  //put in the obsNum and file UTdate and Project ID
  ncfid.add_att("ObsNum",tel->getObsNum());
  ncfid.add_att("UTDate",tel->getUTDate());
  ncfid.add_att("ProjectID",tel->getProjectID().c_str());

  //put in the user offsets
  ncfid.add_att("azUserOffset",tel->getAzUserOff());
  ncfid.add_att("elUserOffset",tel->getElUserOff());

  //put in the M2 offsets
  ncfid.add_att("M2XReq",tel->getM2XReq());
  ncfid.add_att("M2YReq",tel->getM2YReq());
  ncfid.add_att("M2ZReq",tel->getM2ZReq());

  //put in the first three Zernike Coeff
  ncfid.add_att("M1ZernikeC0",tel->getM1ZernikeC0());
  ncfid.add_att("M1ZernikeC1",tel->getM1ZernikeC1());
  ncfid.add_att("M1ZernikeC2",tel->getM1ZernikeC2());

  //also log all of the analysis parameters used as global atributes
  ncfid.add_att("source",ap->getSourceName().c_str());
  ncfid.add_att("analysisDate", t.c_str());
  ncfid.add_att("dataFile", ap->getDataFile());
  ncfid.add_att("despikeSigma", ap->getDespikeSigma());
  ncfid.add_att("lowpassFilterKnee", ap->getLowpassFilterKnee());
  ncfid.add_att("timeOffset", ap->getTimeOffset());
  ncfid.add_att("timeChunk", ap->getTimeChunk());
  ncfid.add_att("neigToCut", ap->getNeigToCut());
  ncfid.add_att("cutStd", ap->getCutStd());
 if (ap->getOrder()>0){
    ncfid.add_att("bsplineOrder", ap->getOrder());
    ncfid.add_att("bsplineControlChunk", ap->getControlChunk());
    ncfid.add_att("bsplineStripe", ap->getCleanStripe());
    ncfid.add_att("bsplineResampling", ap->getResample());

  }

  ncfid.add_att("pixelSize", ap->getPixelSize());
  ncfid.add_att("approximateWeights", ap->getApproximateWeights());
  ncfid.add_att("MasterGrid[0]",masterGrid[0]);
  ncfid.add_att("MasterGrid[1]",masterGrid[1]);
  ncfid.add_att ("ArrayAvgTau", array->getAvgTau());

  cerr << "Observation::writeBeammapsToNcdf(): Maps written to ";
  cerr << ncdfFilename << endl;

  return 1;
}

//----------------------------- o ---------------------------------------

///write all fitParams and Errors of beammaps to NetCDF

bool Observation::writeFitParamsToNcdf(string mapFile, MatDoub &fitParams)
{
  NcFile ncfid = NcFile(mapFile.c_str(), NcFile::Write);
  for(int i=0;i<fitParams.nrows();i++){
    stringstream sstm;
    sstm << "beammapSignal" << i;
    NcVar* sigVar= ncfid.get_var(sstm.str().c_str());
    sigVar->add_att("bolo_name", array->detectorNames[i].c_str());
    sigVar->add_att("dc_offset", fitParams[i][0]);
    sigVar->add_att("dc_offset_err", fitParams[i][7]);
    sigVar->add_att("dc_offset_units", "Jy");
    sigVar->add_att("amplitude", fitParams[i][1]);
    sigVar->add_att("amplitude_err", fitParams[i][8]);
    sigVar->add_att("amplitude_units", "Jy");
    sigVar->add_att("FWHM_x", fitParams[i][2]/TWO_PI*360.*3600.*2.3548);
    sigVar->add_att("FWHM_x_err", fitParams[i][9]/TWO_PI*360.*3600.*2.3548);
    sigVar->add_att("FWHM_x_units", "arcseconds");
    sigVar->add_att("FWHM_y", fitParams[i][3]/TWO_PI*360.*3600.*2.3548);
    sigVar->add_att("FWHM_y_err", fitParams[i][10]/TWO_PI*360.*3600.*2.3548);
    sigVar->add_att("FWHM_y_units", "arcseconds");
    sigVar->add_att("offset_x", fitParams[i][4]/TWO_PI*360.*3600.);
    sigVar->add_att("offset_x_err", fitParams[i][11]/TWO_PI*360.*3600.);
    sigVar->add_att("offset_x_units", "arcseconds");
    sigVar->add_att("offset_y", fitParams[i][5]/TWO_PI*360.*3600.);
    sigVar->add_att("offset_y_err", fitParams[i][12]/TWO_PI*360.*3600.);
    sigVar->add_att("offset_y_units", "arcseconds");
    sigVar->add_att("pos_Angle", fitParams[i][6]);
    sigVar->add_att("pos_Angle_err", fitParams[i][13]);
    sigVar->add_att("pos_Angle_units", "radians");
  }
  ncfid.close();
  return true;
}


//----------------------------- o ---------------------------------------


///writes the entire set of maps to a netcdf file 
bool Observation::writeObservationToNcdf(string ncdfFilename)
{
  //create the file
  NcFile ncfid = NcFile(ncdfFilename.c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);
  if (!ncfid.is_valid()){
    cerr << "Couldn't open map netcdf file for writing!\n";
    return 0;
  }

  //if we are here then we had success opening the file for writing
  ncdfFile.assign(ncdfFilename);
  
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

  //define variables for clean time-streams
  size_t nDetectors = array->getNDetectors();
  size_t nSamples = array->getNSamples();
  size_t nvars = 5;
  NcDim* dimDetectors = ncfid.add_dim("nDetectors", nDetectors);
  NcDim* dimSamples = ncfid.add_dim ("nSamples", nSamples);
  NcDim* dimType = ncfid.add_dim ("types", nvars);
  NcVar* bArray =ncfid.add_var("boloData", ncDouble, dimType, dimDetectors, dimSamples);
  size_t nScans = tel->scanIndex.ncols();
  NcDim* dimScans =ncfid.add_dim("nScans", nScans);
  NcDim* dimScansLimit = ncfid.add_dim("scanLimit", 2);
  NcVar* scansVar = ncfid.add_var("scanIndex",ncInt,dimScansLimit,dimScans);

  //and write the maps
  signalVar->put(&signal->image[0][0], nrows, ncols);
  kernelVar->put(&kernel->image[0][0], nrows, ncols);
  weightVar->put(&weight->image[0][0], nrows, ncols);
  inttimeVar->put(&inttime->image[0][0], nrows, ncols);
  rCPhysVar->put(&rowCoordsPhys[0], nrows);
  cCPhysVar->put(&colCoordsPhys[0], ncols);
  xCAbsVar->put(&xCoordsAbs[0][0], nrows, ncols);
  yCAbsVar->put(&yCoordsAbs[0][0], nrows, ncols);

  //the noise maps are done individually
  int nNoiseMaps =  ap->getNNoiseMapsPerObs();
  string onoise;
  NcVar* nVar;
  for(int i=0;i<nNoiseMaps;i++){
    onoise.assign("noise");
    stringstream o;
    o << i+1;
    onoise.append(o.str());
    nVar = ncfid.add_var(onoise.c_str(), ncDouble, rowDim, colDim);
    nVar->put(&noiseMaps[i][0], nrows, ncols);
  }

  //add timestream data only if Cottinham method is invoked
  if(ap->getOrder() > 0){
    Mat3DDoub boloData (nvars, nDetectors, nSamples);
    int * di = array->getDetectorIndices();
    char buff [100];
    for (size_t ibolo=0; ibolo< nDetectors; ibolo++){
      for (size_t iSample=0; iSample <nSamples; iSample++){
	boloData[0][ibolo][iSample] = array->detectors[di[ibolo]].hValues[iSample];
	boloData[1][ibolo][iSample] = (double)array->detectors[di[ibolo]].hSampleFlags[iSample];
	boloData[2][ibolo][iSample] = array->detectors[di[ibolo]].hRa[iSample];
	boloData[3][ibolo][iSample] = array->detectors[di[ibolo]].hDec[iSample];
	boloData[4][ibolo][iSample] = array->detectors[di[ibolo]].hKernel[iSample];
      }
      sprintf (buff, "NameBolo%lu", ibolo);
      bArray->add_att(buff,array->detectors[di[ibolo]].getName().c_str());
    }
    bArray->put(&boloData[0][0][0], nvars,nDetectors,nSamples);
    scansVar->put(&tel->scanIndex[0][0], 2, nScans);
  }

  //get the time and date of this analysis
  time_t rawtime;
  struct tm* timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  string t = asctime(timeinfo);
  t = t.substr(0,t.length()-1);

  //put in the obsNum and file UTdate and Project ID
  ncfid.add_att("ObsNum",tel->getObsNum());
  ncfid.add_att("UTDate",tel->getUTDate());
  ncfid.add_att("ProjectID",tel->getProjectID().c_str());

  //put in the user offsets
  ncfid.add_att("azUserOffset",tel->getAzUserOff());
  ncfid.add_att("elUserOffset",tel->getElUserOff());

  //put in the M2 offsets
  ncfid.add_att("M2XReq",tel->getM2XReq());
  ncfid.add_att("M2YReq",tel->getM2YReq());
  ncfid.add_att("M2ZReq",tel->getM2ZReq());

  //put in the first three Zernike Coeff
  ncfid.add_att("M1ZernikeC0",tel->getM1ZernikeC0());
  ncfid.add_att("M1ZernikeC1",tel->getM1ZernikeC1());
  ncfid.add_att("M1ZernikeC2",tel->getM1ZernikeC2());

  //also log all of the analysis parameters used as global atributes
  ncfid.add_att("source",ap->getSourceName().c_str());
  ncfid.add_att("analysisDate", t.c_str());
  ncfid.add_att("dataFile", ap->getDataFile()); 
  ncfid.add_att("despikeSigma", ap->getDespikeSigma());
  ncfid.add_att("lowpassFilterKnee", ap->getLowpassFilterKnee());
  ncfid.add_att("timeOffset", ap->getTimeOffset());
  ncfid.add_att("timeChunk", ap->getTimeChunk());
  ncfid.add_att("neigToCut", ap->getNeigToCut());
  ncfid.add_att("cutStd", ap->getCutStd());
 if (ap->getOrder()>0){
	  ncfid.add_att("bsplineOrder", ap->getOrder());
	  ncfid.add_att("bsplineControlChunk", ap->getControlChunk());
	  ncfid.add_att("bsplineStripe", ap->getCleanStripe());
	  ncfid.add_att("bsplineResampling", ap->getResample());

  }
  ncfid.add_att("pixelSize", ap->getPixelSize());
  ncfid.add_att("approximateWeights", ap->getApproximateWeights());
  ncfid.add_att("MasterGrid[0]",masterGrid[0]);
  ncfid.add_att("MasterGrid[1]",masterGrid[1]);
  ncfid.add_att ("ArrayAvgTau", array->getAvgTau());
 
  cerr << "Observation::writeMapsToNcdf(): Maps written to ";
  cerr << ncdfFilename << endl;

  //set the netcdf filename for all the maps in the obs
  signal->mapFile = ncdfFile;
  weight->mapFile = ncdfFile;
  inttime->mapFile = ncdfFile;
  kernel->mapFile = ncdfFile;

  //free the memory in the large noiseMaps array
  cerr << "Observation::Deleting noiseMaps array." << endl;
  noiseMaps.resize(0,0);


  //For some reason the following code was copied in here a second
  //time.  It's not clear why since the arrays are the same as above.
  //I'll comment it out on 6/1/17 and see if anyone complains. - GW
  /*
	//and write the maps
	signalVar->put(&signal->image[0][0], nrows, ncols);
	kernelVar->put(&kernel->image[0][0], nrows, ncols);
	weightVar->put(&weight->image[0][0], nrows, ncols);
	inttimeVar->put(&inttime->image[0][0], nrows, ncols);
	rCPhysVar->put(&rowCoordsPhys[0], nrows);
	cCPhysVar->put(&colCoordsPhys[0], ncols);
	xCAbsVar->put(&xCoordsAbs[0][0], nrows, ncols);
	yCAbsVar->put(&yCoordsAbs[0][0], nrows, ncols);
  */
	//if (atmTemplate)
		//atmMap->put(&atmTemplate->image[0][0], nrows, ncols);
	if (saveTimestreams == true){
		size_t nDetectors = array->getNDetectors();
		size_t nSamples = array->getNSamples();
		size_t nvars = 5;
		NcDim* dimDetectors = ncfid.add_dim("nDetectors", nDetectors);
		NcDim* dimSamples = ncfid.add_dim ("nSamples", nSamples);
		NcDim* dimType = ncfid.add_dim ("types", nvars);
		NcVar* bArray =ncfid.add_var("boloData", ncDouble, dimType, dimDetectors, dimSamples);

		size_t nScans = tel->scanIndex.ncols();

		NcDim* dimScans =ncfid.add_dim("nScans", nScans);
		NcDim* dimScansLimit = ncfid.add_dim("scanLimit", 2);
		NcVar* scansVar = ncfid.add_var("scanIndex",ncInt,dimScansLimit,dimScans);
		Mat3DDoub boloData (nvars, nDetectors, nSamples);
		int * di = array->getDetectorIndices();
		char buff [100];
		for (size_t ibolo=0; ibolo< nDetectors; ibolo++){
			for (size_t iSample=0; iSample <nSamples; iSample++){
				boloData[0][ibolo][iSample] = array->detectors[di[ibolo]].hValues[iSample];
				boloData[1][ibolo][iSample] = (double)array->detectors[di[ibolo]].hSampleFlags[iSample];
				boloData[2][ibolo][iSample] = array->detectors[di[ibolo]].hRa[iSample];
				boloData[3][ibolo][iSample] = array->detectors[di[ibolo]].hDec[iSample];
				boloData[4][ibolo][iSample] = array->detectors[di[ibolo]].hKernel[iSample];
			}
			sprintf (buff, "NameBolo%lu", ibolo);
			bArray->add_att(buff,array->detectors[di[ibolo]].getName().c_str());
		}
		bArray->put(&boloData[0][0][0], nvars,nDetectors,nSamples);
		scansVar->put(&tel->scanIndex[0][0], 2, nScans);
	}

	//free the memory in the large noiseMaps array
	cerr << "Observation::Deleting noiseMaps array." << endl;
	noiseMaps.resize(0,0);

	return 1;
}


//----------------------------- o ---------------------------------------


// TODO to be implemented
bool Observation::writeObservationToFits(string fitsFilename)
{

    (void) fitsFilename;
	return 1;
}


//----------------------------- o ---------------------------------------

///writes the (signal only) map to an 8-bit bitmap file
///(intended primarily for display in the M&C window)
bool Observation::writeObservationToBmp(string bmpFilename)
{

	//plan is to use the intarray2bmp header to write out the image
	int w = nrows;
	int h = ncols;
	int *image = new int[w*h];
	int min=0;
	int max=256;
	double minsig,maxsig;

	//find the max and min of the signal map
	maxsig = signal->image[0][0];
	minsig = signal->image[0][0];
	for(int i=0;i<nrows;i++)
		for(int j=0;j<ncols;j++){
			if(signal->image[i][j] > maxsig) maxsig = signal->image[i][j];
			if(signal->image[i][j] < minsig) minsig = signal->image[i][j];
		}

	for(int i=0;i<w;i++)
		for(int j=0;j<h;j++){
			image[i+j*w] = ((signal->image[i][j]-minsig)/(maxsig-minsig))*255;
		}

	//write out the image
	intarray2bmp::intarray2bmp(bmpFilename.c_str(),image,h,w,min,max);

	return 1;
}


//----------------------------- o ---------------------------------------

///writes beammaps to bitmap format for quick, visual checks
bool Observation::writeBeammapsToBmp(string bmpFilename, int nDetectors)
{
  for(int i = 0; i < nDetectors; i++){
    int w = nrows;
    int h = ncols;
    int *image = new int[w*h];
    int min = 0;
    int max = 256;
    double minsig = 0;
    double maxsig = 0;

    maxsig = beammapSignal[i][0];
    minsig = beammapSignal[i][0];
    for(int j = 0; j < w*h; j++){
      if(beammapSignal[i][j] > maxsig) maxsig = beammapSignal[i][j];
      if(beammapSignal[i][j] < minsig) minsig = beammapSignal[i][j];
    }

    MatDoub switcher(nrows, ncols);
    for(int j = 0; j < nrows; j++){
      for(int k = 0; k < ncols; k++){
        switcher[j][k] = beammapSignal[i][j*ncols + k];
      }
    }

    for(int j = 0; j < nrows; j++){
      for(int k = 0; k < ncols; k++){
        image[j + k*nrows] = ((switcher[j][ncols - k - 1] - minsig)/(maxsig-minsig))*255;
      }
    }

    stringstream ss;
    ss << i << bmpFilename;

    intarray2bmp::intarray2bmp(ss.str().c_str(), image, h, w, min, max);
  }
  return 1;
}

//----------------------------- o ---------------------------------------


///make a histogram of the signal map
/** Makes a histogram of the signal map using the following inputs:
     - nbins - number of bins in the output histogram
     - cc - the desired coverage cut
 **/
bool Observation::histogramSignal(int nbins, double cc)
{
	signal->calcMapHistogram(nbins, cc);
	return 1;
}

bool Observation::histogramSignal(double cc)
{
	int nbins = 200.;
	histogramSignal(nbins, cc);
	return 1;
}

bool Observation::histogramSignal()
{
  int nbins = 200.;
  double coverage_cut = ap->getCoverageThreshold();
  histogramSignal(nbins, coverage_cut);
  return 1;
}


//----------------------------- o ---------------------------------------


///make the psd of the signal map
/**Makes a psd of the signal map using the following inputs:
    - cc - coverage cut (0<cc<1)
 **/
bool Observation::signalMapPsd(double cc)
{
	signal->calcMapPsd(cc);
	return 1;
}

bool Observation::signalMapPsd()
{
  double coverage_cut = ap->getCoverageThreshold();
  signal->calcMapPsd(coverage_cut);
  return 1;
}


//----------------------------- o ---------------------------------------


//destructor
Observation::~Observation()
{
  if (weight)
    delete weight;
  if (signal)
    delete signal;
  if (kernel)
    delete kernel;
  if (inttime)
    delete inttime;
  if (atmTemplate)
    delete atmTemplate;
}
