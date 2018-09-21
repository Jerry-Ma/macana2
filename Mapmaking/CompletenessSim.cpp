#include <netcdfcpp.h>
#include <string>
using namespace std;

#include "nr3.h"
#include "AnalParams.h"
#include "astron_utilities.h"
#include "Coaddition.h"
#include "CompletenessSim.h"
#include "BinomialStats.h"
#include "vector_utilities.h"

//CompletenessSim constructor
CompletenessSim::CompletenessSim(AnalParams* analParams,
				   Coaddition* coaddition)
{
  //set our analysis parameters and coaddition pointers
  ap = analParams;
  cmap = coaddition;

  //grab some parameters
  nBins = ap->getNFluxBins();
  nSources = ap->getNSynthSources();
}


//----------------------------- o ---------------------------------------

//calculates great circle distances
/** Method to calculate the great circle angular separation of two
    points using the exact form of the Haversine formula. (source:
    http://mathforum.org/library/drmath/view/51879.html) Returns
    separation in radians. Inputs are:
      - ra1 - right ascension or longitude of point 1 in degrees
      - dec1 - declination or latitude of point 1 in degrees
      - ra2 - right ascension or longitude of point 2 in degrees
      - dec2 - declination or latitude of point 2 in degrees
**/
double CompletenessSim::gCirc(double ra1, double dec1,
			       double ra2, double dec2)
{
  //use the Haversine formula to find the separation
  double deltRa = ra2 - ra1;
  double deltDec = dec2 - dec1;
  double arg = pow(sin(deltDec/2.0), 2) +
               cos(dec1)*cos(dec2)*pow(sin(deltRa/2.0), 2);
  double sep = 2.0*atan2(sqrt(arg), sqrt(1.0 - arg));

  return sep;
}


//----------------------------- o ---------------------------------------

//calculates map completeness by injecting one source at a time and
//trying to recover it
/**Simulations method to calculate completeness of a coadded map by
   injecting copies of the kernel map, scaled to a known peak
   flux, one at a time into the coadded map and attempting to
   recover them using Coaddition::findSources(double* signalMap,
   double* s2nMap, Coaddition* realCoadd).
**/
//***WARNING!
//            The completeness calculations done here use BOOSTED
//            fluxes! In other words, the fluxes here have *NOT*
//            been de-boosted! The code to do the de-boosting has
//            not been written yet.
//***WARNING!
/// todo - deal with kernel copy ratty edges in a more automatic way,
///        right now it's just hard-coded to remove the first and last
///        85 pixels in both X and Y
///      - add making a S/N map with the injected source that gets fed
///        into overloaded findSources()
bool CompletenessSim::calcCompleteness()
{
  //warn the people that this uses BOOSTED FLUX
  cerr << "*****************" << endl;
  cerr << "****!WARNING!****" << endl;
  cerr << "Fluxes used to calculate completeness are *NOT* de-boosted!" << endl;
  cerr << "****!WARNING!****" << endl;
  cerr << "*****************" << endl;

  //grab some parameters
  double minFlux = ap->getMinFlux();
  double maxFlux = ap->getMaxFlux();
  double recovRadius = ap->getRecovRadius();

  int nRows = cmap->getNrows();
  int nCols = cmap->getNcols();

  //make boolean map indicating where pixels are
  //outside recovRadius of the real sources in addition
  //to inside good coverage region
  MatBool realRecovBool(nRows, nCols);
  for(int i=0;i<nRows;i++){
    for(int j=0;j<nCols;j++){
      if(!cmap->filteredSignal->coverageBool[i][j]){
	realRecovBool[i][j] = 0;
      }
      else{
	int flag = 0;
        for(int k=0;k<cmap->getNSources();k++){
          //use real source centroid coords or not?
	  double dist = (ap->getSFCentroidSources()) ?
                        gCirc(cmap->getXCoordsAbs(i, j),
			      cmap->getYCoordsAbs(i, j),
			      cmap->sources[k].raCentroid,
			      cmap->sources[k].decCentroid) :
                        gCirc(cmap->getXCoordsAbs(i, j),
			      cmap->getYCoordsAbs(i, j),
			      cmap->sources[k].centerRaAbs,
			      cmap->sources[k].centerDecAbs);
	  if(dist  >= recovRadius) flag++;
	}
	if(flag == cmap->getNSources()) realRecovBool[i][j] = 1;
	else realRecovBool[i][j] = 0;
      }
    }
  }

  //populate synthetic source input locations
  bool done = false;
  int count = 0;
  synthRowInput.resize(nBins*nSources);
  synthColInput.resize(nBins*nSources);
  synthRaInput.resize(nBins, nSources);
  synthDecInput.resize(nBins, nSources);
  //grab good coverage intervals to reduce sampling area
  cmap->filteredSignal->setCoverageCutRanges();
  int xRange[2] = {cmap->filteredSignal->getCutXRangeLow(),
		   cmap->filteredSignal->getCutXRangeHigh()};
  int yRange[2] = {cmap->filteredSignal->getCutYRangeLow(),
		   cmap->filteredSignal->getCutYRangeHigh()};
  //randomly draw locations until we have all we need
  while(!done){
    int tempRow = floor(ap->macanaRandom->uniformDeviate(xRange[0], xRange[1]) + 0.5);
    int tempCol = floor(ap->macanaRandom->uniformDeviate(yRange[0], yRange[1]) + 0.5);
    if(realRecovBool[tempRow][tempCol]){
      synthRowInput[count] = tempRow;
      synthColInput[count] = tempCol;
      synthRaInput[count/nSources][count%nSources] =
	cmap->getXCoordsAbs(tempRow, tempCol);
      synthDecInput[count/nSources][count%nSources] = 
	cmap->getYCoordsAbs(tempRow, tempCol);
      count++;
    }
    if(count == nBins*nSources) break;
  }

  //prepare kernel copy for injecting by normalizing
  //to its peak
  /**
  double max = cmap->filteredKernel->image[0][0];
  int currentMaxR = 0;
  int currentMaxC = 0;
  MatDoub kernel2Inject(nRows, nCols);
  for(int i=0;i<nRows;i++){
    for(int j=0;j<nCols;j++){
      //trim ratty edges
      if(cmap->filteredSignal->coverageBool[i][j]){
	kernel2Inject[i][j] = cmap->filteredKernel->image[i][j];
      }
      else kernel2Inject[i][j] = 0.0;
      //find and remember where peak (center) is for shifting
      if(kernel2Inject[i][j] > max){
	max = kernel2Inject[i][j];
	currentMaxR = i;
	currentMaxC = j;
      }
    }
  }
  for(int i=0;i<nRows;i++){
    for(int j=0;j<nCols;j++){
      kernel2Inject[i][j] /= max;
    }
  }
  **/
  //find and remember where kernel peak (center) is for shifting
  MatDoub kernel2Inject(nRows, nCols);
  double max = cmap->filteredKernel->image[0][0];
  int currentMaxR = 0;
  int currentMaxC = 0;
  for(int i=0;i<nRows;i++){
    for(int j=0;j<nCols;j++){
      //trim ratty edges
      if(cmap->filteredSignal->coverageBool[i][j]){
	kernel2Inject[i][j] = cmap->filteredKernel->image[i][j];
      }
      else kernel2Inject[i][j] = 0.0;
      if(kernel2Inject[i][j] > max){
	max = kernel2Inject[i][j];
	currentMaxR = i;
	currentMaxC = j;
      }
    }
  }

  //expand our input flux vector
  inputFluxes.resize(nBins);
  for(int i=0;i<nBins;i++){
    inputFluxes[i] = (i/(nBins-1.0))*(maxFlux - minFlux) + minFlux;
  }

  //expand output matrices
  statusKey.resize(nBins, nSources);
  synthFluxes.resize(nBins, nSources);
  synthS2N.resize(nBins, nSources);
  synthNoises.resize(nBins, nSources);
  synthRas.resize(nBins, nSources);
  synthDecs.resize(nBins, nSources);

  //loop over nBins and nSources injecting the kernel and
  //checking if we recover it
  for(int i=0;i<nBins;i++){
    for(int j=0;j<nSources;j++){
      //put the kernel at random input location and remember where
      shift(kernel2Inject,
	    (currentMaxR - synthRowInput[nSources*i + j]),
	    (currentMaxC - synthColInput[nSources*i + j]));
      currentMaxR = synthRowInput[nSources*i + j];
      currentMaxC = synthColInput[nSources*i + j];

      //make map to search as signal map copy plus scaled kernel copy
      MatDoub map2Search(nRows, nCols);
      //need to add making a S/N map corresponding to including injected source
      for(int k=0;k<nRows;k++){
        for(int l=0;l<nCols;l++){
	  map2Search[k][l] = cmap->filteredSignal->image[k][l] +
	    inputFluxes[i]*1.0e-3*kernel2Inject[k][l];
	}
      }

      //instantiate a Coaddition object to search for sources
      Coaddition synthCoadd(ap);
      //need to eventually pass S/N map made above as second argument
      synthCoadd.findSources(&map2Search[0][0], &map2Search[0][0], cmap,
			     cmap->maxPreS2N);
      for(int k=0;k<synthCoadd.getNSources();k++){
	synthCoadd.sources[k].centroidSource();
      }

      //check synthetic source list for found sources w/in
      //recovery radius of injection point
      int nInRecov = 0;
      vector<int> inRecovI;
      for(int k=0;k<synthCoadd.getNSources();k++){
	double dist = gCirc(cmap->getXCoordsAbs(currentMaxR, currentMaxC),
			    cmap->getYCoordsAbs(currentMaxR, currentMaxC),
			    synthCoadd.sources[k].raCentroid,
			    synthCoadd.sources[k].decCentroid);
	if(dist <= recovRadius){
	  nInRecov++;
	  inRecovI.push_back(k);
	}
      }
      //none close enough to injection point, no detection
      if(nInRecov == 0){
	statusKey[i][j] = 0;
	synthFluxes[i][j] = -99.0;
	synthS2N[i][j] = -99.0;
	synthNoises[i][j] = -99.0;
	synthRas[i][j] = -99.0;
	synthDecs[i][j] = -99.0;
      }
      //one close enough, check if w/in recovery radius of real source
      if(nInRecov == 1){
	bool isReal = false;
	for(int k=0;k<cmap->getNSources();k++){
	  //use real source centroid coords or not?
	  double dist = (ap->getSFCentroidSources()) ?
	                gCirc(cmap->sources[k].raCentroid,
			      cmap->sources[k].decCentroid,
			      synthCoadd.sources[inRecovI[0]].raCentroid,
			      synthCoadd.sources[inRecovI[0]].decCentroid) :
                        gCirc(cmap->sources[k].centerRaAbs,
			      cmap->sources[k].centerDecAbs,
			      synthCoadd.sources[inRecovI[0]].raCentroid,
			      synthCoadd.sources[inRecovI[0]].decCentroid);
	  //just one real source close enough and we can't
	  //distinguish real from synthetic, inconclusive
	  if(dist <= recovRadius){
	    statusKey[i][j] = 2;
	    synthFluxes[i][j] = -99.0;
	    synthS2N[i][j] = -99.0;
	    synthNoises[i][j] = -99.0;
	    synthRas[i][j] = -99.0;
	    synthDecs[i][j] = -99.0;
	    isReal = true;
	    break;
	  }
	}
	//no real sources close enough, synthetic source detected
	if(!isReal){
	  statusKey[i][j] = 1;
	  synthFluxes[i][j] = synthCoadd.sources[inRecovI[0]].centerFlux;
	  synthS2N[i][j] = synthCoadd.sources[inRecovI[0]].centerS2N;
	  synthNoises[i][j] = synthCoadd.sources[inRecovI[0]].centerNoise;
	  synthRas[i][j] = synthCoadd.sources[inRecovI[0]].raCentroid;
	  synthDecs[i][j] = synthCoadd.sources[inRecovI[0]].decCentroid;
	}
      }
      //more than one source found w/in recovery radius of injection
      //point, check those if w/in recovery radius of real source
      if(nInRecov > 1){
	bool containsReal = false;
	for(int k=0;k<cmap->getNSources();k++){
	  for(int l=0;l<nInRecov;l++){
	    //use real source centroid coords or not?
	    double dist = (ap->getSFCentroidSources()) ?
	                  gCirc(cmap->sources[k].raCentroid,
			        cmap->sources[k].decCentroid,
			        synthCoadd.sources[inRecovI[l]].raCentroid,
			        synthCoadd.sources[inRecovI[l]].decCentroid) :
                          gCirc(cmap->sources[k].centerRaAbs,
			        cmap->sources[k].centerDecAbs,
			        synthCoadd.sources[inRecovI[l]].raCentroid,
			        synthCoadd.sources[inRecovI[l]].decCentroid);
	    //just one real source close enough and we can't
	    //distinguish real from synthetic, inconclusive
	    if(dist <= recovRadius){
	      statusKey[i][j] = 2;
	      synthFluxes[i][j] = -99.0;
	      synthS2N[i][j] = -99.0;
	      synthNoises[i][j] = -99.0;
	      synthRas[i][j] = -99.0;
	      synthDecs[i][j] = -99.0;
	      containsReal = true;
	      break;
	    }
	  }
	  if(containsReal){
	    break;
	  }
	}
	//no real sources close enough, find closest source and save
	//as detection of synthetic source
	if(!containsReal){
	  double minDist = gCirc(cmap->getXCoordsAbs(currentMaxR, currentMaxC),
			         cmap->getYCoordsAbs(currentMaxR, currentMaxC),
			         synthCoadd.sources[inRecovI[0]].raCentroid,
			         synthCoadd.sources[inRecovI[0]].raCentroid);
	  int minDI = 0;
	  for(int k=1;k<nInRecov;k++){
	    double dist = gCirc(cmap->getXCoordsAbs(currentMaxR, currentMaxC),
			        cmap->getYCoordsAbs(currentMaxR, currentMaxC),
			        synthCoadd.sources[inRecovI[k]].raCentroid,
			        synthCoadd.sources[inRecovI[k]].raCentroid);
	    if(dist < minDist){
	      minDist = dist;
	      minDI = k;
	    }
	  }
	  statusKey[i][j] = 1;
	  synthFluxes[i][j] = synthCoadd.sources[inRecovI[minDI]].centerFlux;
	  synthS2N[i][j] = synthCoadd.sources[inRecovI[minDI]].centerS2N;
	  synthNoises[i][j] = synthCoadd.sources[inRecovI[minDI]].centerNoise;
	  synthRas[i][j] = synthCoadd.sources[inRecovI[minDI]].raCentroid;
	  synthDecs[i][j] = synthCoadd.sources[inRecovI[minDI]].decCentroid;
	}
      }
    }
  }

  //calculate completeness and normal approx of binomial error
  completeness.resize(nBins);
  nConclusive.resize(nBins);
  for(int i=0;i<nBins;i++){
    completeness[i] = 0.0;
    nConclusive[i] = 0.0;
    for(int j=0;j<nSources;j++){
      if(statusKey[i][j] != 2){
        completeness[i] += statusKey[i][j];
	nConclusive[i]++;
      }
    }
  }
  compErr.resize(nBins);
  compErrLow.resize(nBins);
  compErrHigh.resize(nBins);
  //completeness and binomial errors are stored as percentages
  //not fractional form
  for(int i=0;i<nBins;i++){
    double t = completeness[i]/double(nConclusive[i]);
    compErr[i] = sqrt(completeness[i]*(1.0 - t))/nConclusive[i]*100.0;
    completeness[i] = t*100.0;
    BinomialStats errIntervals(nConclusive[i], completeness[i], 0.68);
    errIntervals.calcIntervals();
    compErrLow[i] = errIntervals.getErrLow();
    compErrHigh[i] = errIntervals.getErrHigh();
  }

  //warn the villagers again that this uses BOOSTED FLUX
  cerr << "*****************" << endl;
  cerr << "****!WARNING!****" << endl;
  cerr << "Fluxes used to calculate completeness are *NOT* de-boosted!" << endl;
  cerr << "****!WARNING!****" << endl;
  cerr << "*****************" << endl;

  return 1;
}


//----------------------------- o ---------------------------------------

///adds completeness results to coadded map NCDF file
/**Adds the flux bin values, number of sources injected per bin,
   completeness values and completeness errors resulting from
   simulations.
**/

bool CompletenessSim::addCompletenessToNcdf()
{
  //open coadded file to append sources to
  NcFile ncfid = NcFile(ap->getCoaddOutFile().c_str(), NcFile::Write);
  if (!ncfid.is_valid()){
    cerr << "PointSource::addSourceToNCDF(): ";
    cerr << "Couldn't open parent map file for writing,";
    cerr << " exiting." << endl;
    exit(1);
  }

  //create nBin dimension
  NcDim* nBinDim = ncfid.add_dim("compNFluxBins", nBins);

  //create variables
  NcVar *inputFluxVar = ncfid.add_var("compInFluxes", ncDouble, nBinDim);
  NcVar *compVar = ncfid.add_var("completeness", ncDouble, nBinDim);
  NcVar *compErrVar = ncfid.add_var("compErr", ncDouble, nBinDim); 
  NcVar *compErrLVar = ncfid.add_var("compErrLow", ncDouble, nBinDim);
  NcVar *compErrHVar = ncfid.add_var("compErrHigh", ncDouble, nBinDim);

  //add nSources per bin attribute to inputFluxes variable
  inputFluxVar->add_att("nSourcesPerBin", nSources);

  //write the variables
  inputFluxVar->put(&inputFluxes[0], nBins);
  compVar->put(&completeness[0], nBins);
  compErrVar->put(&compErr[0], nBins);
  compErrLVar->put(&compErrLow[0], nBins);
  compErrHVar->put(&compErrHigh[0], nBins);

  return 1;
}


//----------------------------- o ---------------------------------------

int CompletenessSim::getNBins()
{
  return nBins;
}

int CompletenessSim::getNSources()
{
  return nSources;
}

double CompletenessSim::getInputFlux(int i)
{
  return inputFluxes[i];
}

int CompletenessSim::getStatusKey(int i, int j)
{
  return statusKey[i][j];
}

double CompletenessSim::getCompleteness(int i)
{
  return completeness[i];
}

double CompletenessSim::getCompErr(int i)
{
  return compErr[i];
}

double CompletenessSim::getCompErrL(int i)
{
  return compErrLow[i];
}

double CompletenessSim::getCompErrH(int i)
{
  return compErrHigh[i];
}

int CompletenessSim::getNConclusive(int i)
{
  return nConclusive[i];
}


//----------------------------- o ---------------------------------------

CompletenessSim::~CompletenessSim()
{
  //destructor
}
