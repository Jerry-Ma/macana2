#include <netcdfcpp.h>
using namespace std;

#include "PointSource.h"
#include "astron_utilities.h"
#include "gaussFit.h"

///PointSource constructor
PointSource::PointSource()
{
  //initialize identification
  sID=0;
  nSourcesParentMap=0;

  //initialize location members
  centerRaAbs = 0.0;
  centerDecAbs = 0.0;
  centerRaPhys = 0.0;
  centerDecPhys = 0.0;
  centerXPos = 0;
  centerYPos = 0;
  raCentroid = 0.0;
  decCentroid= 0.0;
  raPhysCentroid = 0.0;
  decPhysCentroid = 0.0;
  xCentroid = 0.0;
  yCentroid = 0.0;

  //initialize source peak pixel property members
  centerFlux = 0.0;
  centerNoise = 0.0;
  centerS2N = 0.0;

  //storage
  parentMapFile = "";
  pixelSize = 0.0;
}

//----------------------------- o ---------------------------------------

///constructor that sets up pointers to parent map
/** Constructor that stores addresses of parent map vectors and matrices
    along with the analysis parameters pointer and some variables.
    Inputs are:
      - mapSignal - MatDoub containing the parent signal map
      - mapWeight - MatDoub containing the parent weight map
      - mapRowCoordsPhys - VecDoub containing physical
                           (delta-source) row coordinates of
                           parent map
      - mapColCoordsPhys - VecDoub containing physical
                           (delta-source) column coordinates of
                           parent map
      - mapXCoordsAbs - MatDoub containing absolute row coordinates
                        of parent map
      - mapYCoordsAbs - MatDoub containing absolute column
                        coordinates of parent map
      - mapNRows - number of rows in parent map
      - mapNCols - number of columns in parent map
      - pixSize - pixel size in radians
      - parentMFile - c string of parent map filename source found in

**/
bool PointSource::initialize(AnalParams* analParams,
                             double* mSignal,
			     double* mWeight,
			     double* mRowCoordsPhys,
			     double* mColCoordsPhys,
			     double* mXCoordsAbs,
			     double* mYCoordsAbs,
			     int mNRows, int mNCols,
			     double pixSize, std::string parentMFile)
{
  //set up the analysis parameters pointer
  ap = analParams;

  //set up MatDoub and VecDoub pointers
  mapSignal = mSignal;
  mapWeight = mWeight;
  mapRowCoordsPhys = mRowCoordsPhys;
  mapColCoordsPhys = mColCoordsPhys;
  mapXCoordsAbs = mXCoordsAbs;
  mapYCoordsAbs = mYCoordsAbs;

  //simply store some other variables
  mapNRows = mNRows;
  mapNCols = mNCols;
  pixelSize = pixSize;
  centroidWin = ap->getBeamSize()/2.0;
  parentMapFile = parentMFile;

  return 1;
}

//----------------------------- o ---------------------------------------

///store postage stamp of source found with SourceLocate::findSources
/** Stores signal and weight postage stamp for point source
    found with SourceLocate::findSources. Pixels outside of parent
    map are de-weighted (weight=0.0), have signal flagged as -9999
    and all coordinates flagged 9999. Forces postage stamp to be
    an odd number of pixels on a side so source pixel is always
    in center.
 **/
bool PointSource::makePostageStamp()
{
  //pull side length
  int sideLength = ap->getPsSideLength();
  //force sideLength to be odd so source pixel is at center
  if(sideLength%2 == 0){
    sideLength += 1;
  }

  //puts equal pixels on either side of source forcing it
  //to be at center
  int rowRadius = sideLength/2;
  int colRadius = sideLength/2;

  //keeps this general so we can easily switch to having
  //length and width specified by input
  int rowLength = sideLength;
  int colLength = sideLength;

  //get upper and lower indices for rows and columns
  int rowLInd = centerXPos - rowRadius;
  int rowUInd = centerXPos + rowRadius;
  int colLInd = centerYPos - colRadius;
  int colUInd = centerYPos + colRadius;

  //prepare arrays to store values
  MatDoub signalPS(rowLength, colLength);
  MatDoub weightPS(rowLength, colLength);
  VecDoub rowCoordsPhys(rowLength);
  VecDoub colCoordsPhys(colLength);
  xCoordsAbs.resize(rowLength, colLength);
  yCoordsAbs.resize(rowLength, colLength);

  MatDoub rawPS(rowLength, colLength);
  int rCount = 0;
  //loop over row and column ranges to populate stamps
  for(int j=rowLInd;j<rowUInd+1;j++){
    int cCount = 0;
    for(int k=colLInd;k<colUInd+1;k++){
      //de-weight pixel and flag signal and coords if outside map
      if(j < 0 ||
         j > (mapNRows-1) ||
         k < 0 ||
	 k > (mapNCols-1)){
	 signalPS[rCount][cCount] = -9999;
	 weightPS[rCount][cCount] = 0.0;
	 rowCoordsPhys[rCount] = 9999;
	 colCoordsPhys[cCount] = 9999;
	 xCoordsAbs[rCount][cCount] = 9999;
	 yCoordsAbs[rCount][cCount] = 9999;
      }
      //otherwise take from coadded map and place in PointSource
      else{
	signalPS[rCount][cCount] = mapSignal[mapNCols*j + k];
	weightPS[rCount][cCount] = mapWeight[mapNCols*j + k];
        rowCoordsPhys[rCount] = mapRowCoordsPhys[j];
        colCoordsPhys[cCount] = mapColCoordsPhys[k];
	xCoordsAbs[rCount][cCount] = mapXCoordsAbs[mapNCols*j + k];
	yCoordsAbs[rCount][cCount] = mapYCoordsAbs[mapNCols*j + k];
      }
      cCount++;
    }
    rCount++;
  }

  //store nRows and nCols of postage stamp
  int nRowsPS = (rowUInd - rowLInd) + 1;
  int nColsPS = (colUInd - colLInd) + 1;

  //create the Map object
  string mname;
  postageStamp = new Map(mname.assign("postageStampSignal"), nRowsPS, nColsPS, pixelSize, 
		   weightPS, rowCoordsPhys, colCoordsPhys);
  for(int i=0;i<nRowsPS;i++)
    for(int j=0;j<nColsPS;j++){
      postageStamp->image[i][j] = signalPS[i][j];
    }

  return 1;
}

//----------------------------- o ---------------------------------------


///centroid source postage stamp
/** Centroid source postage stamp found with findSources weighted
    by flux^2 and stores the absolute, physical and fractional index
    coordinates in PointSource objects.
**/
bool PointSource::centroidSource()
{
  int nRowsPS = postageStamp->getNrows();
  int nColsPS = postageStamp->getNcols();

  //is centroidWin positive?
  if(centroidWin < 0.0){
    cerr << "PointSource::centroidSource(): ";
    cerr << "Negative centroid window supplied, exiting" << endl;
    exit(1);
  }

  //have postage stamps been made?
  if(nRowsPS == 0 && nColsPS == 0){
    cerr << "PointSource::centroidSource(): ";
    cerr << "Postage stamps have not been made yet, ";
    cerr << "making them first before centroiding." << endl;
    makePostageStamp();
  }

  //find pixels within centroidWin of source pixel
  vector <int> rCentInd;
  vector <int> cCentInd;
  for(int j=0;j<nRowsPS;j++){
    for(int k=0;k<nColsPS;k++){
      double r = sqrt(pow(postageStamp->getRowCoordsPhys(j) - centerRaPhys, 2) +
                      pow(postageStamp->getColCoordsPhys(k) - centerDecPhys, 2));
      if(r <= centroidWin){
        rCentInd.push_back(j);
	cCentInd.push_back(k);
      }
    }
  }
  int nPix = rCentInd.size();

  //calculate flux^2 weighted centroid within centroid window
  double decCent;
  double raCent;
  double rowWeightCent = 0.0;
  double colWeightCent = 0.0;
  double fluxWeight = 0.0;
  for(int k=0;k<nPix;k++){
    rowWeightCent += postageStamp->getRowCoordsPhys(rCentInd[k])*
                     pow(postageStamp->image[rCentInd[k]][cCentInd[k]], 2);
    colWeightCent += postageStamp->getColCoordsPhys(cCentInd[k])*
                     pow(postageStamp->image[rCentInd[k]][cCentInd[k]], 2);
    fluxWeight += pow(postageStamp->image[rCentInd[k]][cCentInd[k]], 2);
  }
  decCent = colWeightCent/fluxWeight;
  raCent = rowWeightCent/fluxWeight;

  //store absolute, physical and index centroided coordinates
  double decOff = decCent - centerDecPhys;
  double decPixOff = decOff/postageStamp->getPixelSize();
  decCentroid = centerDecAbs + decOff;

  double raOff = (raCent - centerRaPhys)/cos(decCentroid);
  double raPixOff = raOff/postageStamp->getPixelSize();
  raCentroid = centerRaAbs + raOff;
  decPhysCentroid = decCent;
  raPhysCentroid = raCent;
  yCentroid = centerYPos + decPixOff;
  xCentroid = centerXPos + raPixOff;

  return 1;
}

//----------------------------- o ---------------------------------------

///fits 2-D Gaussian to source postage stamp using the method from the Map class
/**Fits a 2-D Gaussian to source postage stamp using Map method
   and stores the fit parameters into the gaussParams vector
   data member in PointSource. Fit parameters in their order in
   gaussParams are:
     1) dc offset
     2) peak amplitude
     3) Full Width at Half Maximum in x direction (radians)
     4) Full Width at Half Maximum in y direction (radians)
     5) fitted x position (radians from center of stamp)
     6) fitted y position (radians from center of stamp)
     7) position angle (radians counterclockwise from x-axis)
**/

bool PointSource::fitGaussianToSource(double fwhmx, double fwhmy)
{
  //have postage stamps been made?
  if(postageStamp->getNrows() == 0 && postageStamp->getNcols() == 0){
    cerr << "PointSource::fitGaussianToSource(): ";
    cerr << "Postage stamps have not been made yet, ";
    cerr << "making them first before fitting Gaussian." << endl;
    makePostageStamp();
  }

  double tsigma = 8.0;
  if (ap->getObservatory()=="ASTE")
	  tsigma = 30.0;
  else if (ap->getObservatory() =="JCMT")
	  tsigma = 18.0;

  tsigma *= FWHM_TO_SIGMA * RAD_ASEC;


  //fix dc offset to 0.
  VecInt fixme(7,0);
  VecDoub fixVals(7,0.);
  VecDoub iguess(7,0.0);

  fixme[0]=1;
  fwhmx=0;
  fwhmy=0;
//  cout<<"Assuming a pointsource template of "<< fwhmx *SIGMA_TO_FWHM / RAD_ASEC << ","<<fwhmy *SIGMA_TO_FWHM /RAD_ASEC<<endl;
  if (fwhmx > 0.0){
	  fixVals[2]=fwhmx;
	  fixme[2] = 1;
	  iguess[2] = fwhmx;
  } else
	  iguess[2] = tsigma;
  if (fwhmy > 0.0){
  	  fixVals[3]=fwhmy;
  	  fixme[3] = 1;
  	 iguess[3] = fwhmy;
  } else
	 iguess[3] = tsigma;

  size_t xc = this->postageStamp->image.nrows()/2;
  size_t yc = this->postageStamp->image.ncols()/2;
  iguess[1] = this->postageStamp->image[xc][yc]/2.0;
  iguess[4] = postageStamp->getRowCoordsPhys(xc);
  iguess[5] = postageStamp->getColCoordsPhys(yc);

  cout<<"Using initial guess: ";
  for (size_t igg=0; igg<7; igg++)
	  cout<< iguess[igg]<<", ";
  cout<<endl;

  redchisq=postageStamp->fitToGaussian(gaussParams,fixme,fixVals, &iguess[0]);

  return 1;
}


//----------------------------- o ---------------------------------------

///adds PointSource object to the parent coadded map NCDF file
/** Adds the PointSource object to the parent coadded map NCDF
    file. Writes only the signal array defining the stamp as
    a variable, writes the number of rows and columns of stamp
    as dimensions and writes all data members of the object
    as attributes of the signal array variable (excluding the
    vectors and arrays of coordinates). Input is:
**/

bool PointSource::addSourceToNCDF()
{

  //open coadded file to append sources to
  NcFile ncfid = NcFile(parentMapFile.c_str(), NcFile::Write);
  if (!ncfid.is_valid()){
    cerr << "PointSource::addSourceToNCDF(): ";
    cerr << "Couldn't open parent map file for writing,";
    cerr << " exiting." << endl;
    exit(1);
  }

  //make signal variable name
  string varName;
  ostringstream convert;
  convert << sID;
  varName = convert.str();
  varName = "pointSource_" + varName;

  //create or retrieve dimensions
  NcDim* nSourceDim;
  NcDim* rowDim;
  NcDim* colDim;
  int nrows = postageStamp->getNrows();
  int ncols = postageStamp->getNcols();
  if(varName == "pointSource_0"){
    nSourceDim = ncfid.add_dim("nSources", nSourcesParentMap);
    rowDim = ncfid.add_dim("postageStampNRows", nrows);
    colDim = ncfid.add_dim("postageStampNCols", ncols);
  }
  else{
    rowDim = ncfid.get_dim("postageStampNRows");
    colDim = ncfid.get_dim("postageStampNCols");
  }

  //write in the coordinates
  string rcpvname(varName.c_str());
  string ccpvname(varName.c_str());
  rcpvname.append("_rcp");
  ccpvname.append("_ccp");
  NcVar *rCPhysVar = ncfid.add_var(rcpvname.c_str(), ncDouble, rowDim);
  NcVar *cCPhysVar = ncfid.add_var(ccpvname.c_str(), ncDouble, colDim);
  VecDoub rcp(nrows);
  for(int i=0;i<nrows;i++) rcp[i] = postageStamp->getRowCoordsPhys(i);
  rCPhysVar->put(&rcp[0], nrows);
  VecDoub ccp(ncols);
  for(int i=0;i<ncols;i++) ccp[i] = postageStamp->getColCoordsPhys(i);
  cCPhysVar->put(&ccp[0], ncols);

  //define map variable
  NcVar *signalVar = ncfid.add_var(varName.c_str(),
                                   ncDouble, rowDim, colDim);

  //add attributes to signal variable
    //center pixel positions
  signalVar->add_att("centerRaAbs", centerRaAbs);
  signalVar->add_att("centerDecAbs", centerDecAbs);
  signalVar->add_att("centerRaPhys", centerRaPhys);
  signalVar->add_att("centerDecPhys", centerDecPhys);
  signalVar->add_att("centerXPos", centerXPos);
  signalVar->add_att("centerYPos", centerYPos);

  //centroided positions
  signalVar->add_att("raCentroid", raCentroid);
  signalVar->add_att("decCentroid", decCentroid);
  signalVar->add_att("raPhysCentroid", raPhysCentroid);
  signalVar->add_att("decPhysCentroid", decPhysCentroid);
  signalVar->add_att("xCentroid", xCentroid);
  signalVar->add_att("yCentroid", yCentroid);

  //center values
  signalVar->add_att("centerFlux", centerFlux);
  signalVar->add_att("centerNoise", centerNoise);
  signalVar->add_att("centerS2N", centerS2N);

  //Gaussian fit parameters
  VecDoub gaussParams_err(7,0.);
  for(int i=0;i<7;i++) gaussParams_err[i] = gaussParams[i+7];
  signalVar->add_att("dc_offset", gaussParams[0]);
  signalVar->add_att("dc_offset_err", gaussParams_err[0]);
  signalVar->add_att("dc_offset_units", "Jy");
  signalVar->add_att("amplitude", gaussParams[1]);
  signalVar->add_att("amplitude_err", gaussParams_err[1]);
  signalVar->add_att("amplitude_units", "Jy");
  signalVar->add_att("FWHM_x", gaussParams[2]/TWO_PI*360.*3600.*2.3548);
  signalVar->add_att("FWHM_x_err", gaussParams_err[2]/TWO_PI*360.*3600.*2.3548);
  signalVar->add_att("FWHM_x_units", "arcseconds");
  signalVar->add_att("FWHM_y", gaussParams[3]/TWO_PI*360.*3600.*2.3548);
  signalVar->add_att("FWHM_y_err", gaussParams_err[3]/TWO_PI*360.*3600.*2.3548);
  signalVar->add_att("FWHM_y_units", "arcseconds");
  signalVar->add_att("offset_x", gaussParams[4]/TWO_PI*360.*3600.);
  signalVar->add_att("offset_x_err", gaussParams_err[4]/TWO_PI*360.*3600.);
  signalVar->add_att("offset_x_units", "arcseconds");
  signalVar->add_att("offset_y", gaussParams[5]/TWO_PI*360.*3600.);
  signalVar->add_att("offset_y_err", gaussParams_err[5]/TWO_PI*360.*3600.);
  signalVar->add_att("offset_y_units", "arcseconds");
  signalVar->add_att("pos_Angle", gaussParams[6]);
  signalVar->add_att("pos_Angle_err", gaussParams_err[6]);
  signalVar->add_att("pos_Angle_units", "radians");
  signalVar->add_att("redChisq",redchisq);
  
  //storage
  signalVar->add_att("pixelSize", postageStamp->getPixelSize());
  signalVar->add_att("centroidWin", centroidWin);
  signalVar->add_att("parentMapFile", parentMapFile.c_str());

  //write the signal postage stamp
  signalVar->put(&postageStamp->image[0][0], postageStamp->getNrows(), postageStamp->getNcols());

  return 1;
}

//----------------------------- o ---------------------------------------

PointSource::~PointSource()
{
  //destructor
}
