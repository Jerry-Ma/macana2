#include <netcdfcpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <cstdio>
#include <fftw3.h>
#include <gsl/gsl_histogram.h>
#include <omp.h>
using namespace std;

#include "nr3.h"
#include "AnalParams.h"
#include "Array.h"
#include "astron_utilities.h"
#include "gaussFit.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include "Coaddition.h"
#include "NoiseRealizations.h"
#include "Telescope.h"
#include "vector_utilities.h"

///NoiseRealizations constructor
/** This constructor primarily manages initiation of NoiseRealizations
    parameters.  Error checking on the paths and filenames should
    probably be added either here or in AnalParams.cpp.
    \todo - add error checking on paths either here or in AnalParams
**/
NoiseRealizations::NoiseRealizations(AnalParams* analParams)
{
  //set our ap pointer
  ap = analParams;

  //the number of noise files to produce
  nNoiseFiles = ap->getNRealizations();

  //the path for the produced noise files
  noisePath = ap->getNoisePath();
  string np = noisePath;
  np.append("noise");

  //a string array for the filenames
  noiseFiles = new string [nNoiseFiles];
  for(int i=0;i<nNoiseFiles;i++){
    noiseFiles[i] = np;
    stringstream o;
    o << i;
    noiseFiles[i].append(o.str());
    noiseFiles[i].append(".nc");
  }

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
}


//----------------------------- o ---------------------------------------

///the workhorse of NoiseRealizations - generates them all
/** This is where we loop through the observation outputs to generate
    the noise realizations.  We pass a pointer to a coaddition so that
    we can copy common items without needing to regenerate them.
    The strategy for generating the noise realizations is:
     - choose 1 of the 5 jacknifed map from each observation
     - coadd them using the corresponding weight maps
     - generate the histogram of each coadded noise map
     - generate the 1-d psd of each coadded noise map
     - write this all to a netcdf file
     - repeat this a bunch of times.
    Note that the idea of 5 jacknifed maps for each observation is 
    hard coded.  This should probably be made into an AnalParam parameter.
    \todo There could be significant savings in the psd generation for
    larger maps by saving the "gsl_plan" for the fft and reusing it 
    for each noise map realization instead of recalculating the plan
    each time.
**/
bool NoiseRealizations::generateNoiseRealizations(Coaddition* cmap)
{
  //most things were calculated to get cmap so we'll steal them
  int nFiles = ap->getNFiles();
  nrows = cmap->getNrows();
  ncols = cmap->getNcols();
  pixelSize = cmap->getPixelSize();

  //physical coordinates
  rowCoordsPhys.resize(nrows);
  colCoordsPhys.resize(ncols);
  for(int i=0;i<nrows;i++) rowCoordsPhys[i] = cmap->getRowCoordsPhys(i);
  for(int i=0;i<ncols;i++) colCoordsPhys[i] = cmap->getColCoordsPhys(i);

  //matrices of absolute coordinates
  xCoordsAbs.resize(nrows,ncols);
  yCoordsAbs.resize(nrows,ncols);
  for(int i=0;i<nrows;i++)
    for(int j=0;j<ncols;j++){
      xCoordsAbs[i][j] = cmap->getXCoordsAbs(i,j);
      yCoordsAbs[i][j] = cmap->getYCoordsAbs(i,j);
    }

  //and the weight map
  //the weights are the same as in cmap but we need them here later
  MatDoub weight;
  weight.assign(nrows,ncols,0.);
  for(int i=0;i<nrows;i++) 
    for(int j=0;j<ncols;j++) weight[i][j] = cmap->weight->image[i][j];
  cerr << "NoiseRealizations(): nrows=" << nrows << ", ncols=" << ncols << endl;

  //initialize the noise storage map
  noise = new Map(string("noise"), nrows, ncols, pixelSize,
		  weight, rowCoordsPhys, colCoordsPhys);

  Map *myNoise=NULL;
  int mynrows= nrows;
  int myncols = ncols;
  double myPixelSize = pixelSize;
  VecDoub myRow = rowCoordsPhys;
  VecDoub myCol = colCoordsPhys;
  //now do the coadditions picking the obs noise map at random
  //write the noise map to a file before moving on
  cerr << "Generating Noise Realizations: " << endl;
  int n;
  GslRandom* ran;
  ran = ap->macanaRandom;

  #pragma omp parallel shared (ran,weight, mynrows, myncols, myPixelSize, myRow, myCol) private (n,myNoise)
  {
  #pragma omp for schedule(dynamic)
  for(int inoise=0;inoise<nNoiseFiles;inoise++){
	myNoise = new Map(string("noise"), mynrows, myncols, myPixelSize,
			  weight, myRow, myCol);
    myNoise->image.assign(mynrows,myncols,0.);
    for(int k=0;k<nFiles;k++){
    	VecDoub rcp(0);
    	VecDoub ccp(0);
    	MatDoub os(0,0);
    	MatDoub ow(0,0);
    	double onrows=0;
    	double oncols=0;
#pragma omp critical (noiseDataIO)
    	{
	  NcFile ncfid = NcFile(ap->getMapFileList(k).c_str(), NcFile::ReadOnly);
	  
	  onrows = ncfid.get_dim("nrows")->size();
	  oncols = ncfid.get_dim("ncols")->size();
	  NcVar* rcpv = ncfid.get_var("rowCoordsPhys");
	  NcVar* ccpv = ncfid.get_var("colCoordsPhys");
	  rcp.resize(onrows);
	  ccp.resize(oncols);
	  for(int i=0;i<onrows;i++) rcp[i] = rcpv->as_double(i);
	  for(int j=0;j<oncols;j++) ccp[j] = ccpv->as_double(j);
	  
	  //randomly choose one of the noise maps
	  int nN = ap->getNNoiseMapsPerObs();
	  n = floor(ran->uniformDeviate(1,nN));
	  if(n == nN) n = nN;
	  string onoise = "noise";
	  stringstream o;
	  o << n;
	  onoise.append(o.str());
	  NcVar* osVar = ncfid.get_var(onoise.c_str());
	  NcVar* owVar = ncfid.get_var("weight");
	  os.resize(onrows,oncols);
	  ow.resize(onrows,oncols);
	  osVar->get(&os[0][0], onrows, oncols);
	  owVar->get(&ow[0][0], onrows, oncols);
    	}
      //the index deltas
      int deltai = (rcp[0]-rowCoordsPhys[0])/pixelSize;
      int deltaj = (ccp[0]-colCoordsPhys[0])/pixelSize;
      
      //now loop through the observation maps
      for(int oi=0;oi<onrows;oi++){
    	  for(int oj=0;oj<oncols;oj++){
    		  myNoise->image[oi+deltai][oj+deltaj] += ow[oi][oj]*os[oi][oj];
    	  }
      }
    }


    //normalization
    for(int i=0;i<mynrows;i++)
      for(int j=0;j<myncols;j++){
    	  myNoise->image[i][j] = (weight[i][j] != 0.) ?
    	  myNoise->image[i][j]/weight[i][j] : 0.;
      }

    //make a histogram of the noise map and add to average
    myNoise->calcMapHistogram(200,ap->getCoverageThreshold());

    //calculate the noise map psd
    myNoise->calcMapPsd(ap->getCoverageThreshold());

    //write the noise map, psd, and histogram out to a file
	#pragma omp critical (noiseDataIO)
    {
    	writeNoiseMapToNcdf(noiseFiles[inoise], myNoise);
    }
    //what a waste of time this is ... sigh
    int tid;
	#if defined(_OPENMP)
    	tid = omp_get_thread_num() + 1;
    #else
    	 tid = 0;
    #endif
    cerr<<"NoiseRealizations("<<tid<<"): Made Noise map: "<<inoise+1<<endl;

    delete myNoise;

  }
  }
  cerr << endl;
  cerr << "Noise maps, psds and histograms written to " << noisePath << endl;
  
  return 1;
}


//----------------------------- o ---------------------------------------

///writes noise maps, histograms, psds and metadata to netcdf file
bool NoiseRealizations::writeNoiseMapToNcdf(string mapFileName, Map *noise)
{
  //create the file
  NcFile ncfid = NcFile(mapFileName.c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);
  if (!ncfid.is_valid()){
    cerr << "Couldn't open map netcdf file for writing!\n";
    return 0;
  }

  //create dimensions
  NcDim* rowDim = ncfid.add_dim("nrows", nrows);
  NcDim* colDim = ncfid.add_dim("ncols", ncols);

  //define variables for maps
  NcVar *noiseVar = ncfid.add_var("noise", ncDouble, rowDim, colDim);
  NcVar *weightVar = ncfid.add_var("weight", ncDouble, rowDim, colDim);
  NcVar *rCPhysVar = ncfid.add_var("rowCoordsPhys", ncDouble, rowDim);
  NcVar *cCPhysVar = ncfid.add_var("colCoordsPhys", ncDouble, colDim);
  NcVar *xCAbsVar = ncfid.add_var("xCoordsAbs", ncDouble, rowDim, colDim);
  NcVar *yCAbsVar = ncfid.add_var("yCoordsAbs", ncDouble, rowDim, colDim);

  //and write the maps
  noiseVar->put(&noise->image[0][0], nrows, ncols);
  weightVar->put(&noise->weight[0][0], nrows, ncols);
  rCPhysVar->put(&rowCoordsPhys[0], nrows);
  cCPhysVar->put(&colCoordsPhys[0], ncols);
  xCAbsVar->put(&xCoordsAbs[0][0], nrows, ncols);
  yCAbsVar->put(&yCoordsAbs[0][0], nrows, ncols);

  //the histogram bins and values
  //create dimension
  NcDim* histDim = ncfid.add_dim("nhist", noise->histBins.size());
  //define hist variables
  NcVar *binVar = ncfid.add_var("histBins", ncDouble, histDim);
  NcVar *valVar = ncfid.add_var("histVals", ncDouble, histDim);
  //put in the values
  binVar->put(&noise->histBins[0], noise->histBins.size());
  valVar->put(&noise->histVals[0], noise->histBins.size());

  //the psd freqs and values
  //create dimension
  NcDim* psdDim = ncfid.add_dim("npsd", noise->psd.size());
  //define psd variables
  NcVar *freqVar = ncfid.add_var("psdFreq", ncDouble, psdDim);
  NcVar *psdVar = ncfid.add_var("psd", ncDouble, psdDim);
  //put in the values
  freqVar->put(&noise->psdFreq[0], noise->psd.size());
  psdVar->put(&noise->psd[0], noise->psd.size());

  string psdDimName = "nxpsd_2d";

  int nx = noise->psd2d.nrows();
  int ny = noise->psd2d.ncols();

  NcDim* psdDim2d_x = ncfid.add_dim(psdDimName.c_str(), nx);
  psdDimName = "nypsd_2d";
  NcDim* psdDim2d_y = ncfid.add_dim(psdDimName.c_str(), ny);

  string psdVarName = "psd_2d";
  string psdFreqName = "psdFreq_2d";

  NcVar *psdVar2d = ncfid.add_var(psdVarName.c_str(), ncDouble, psdDim2d_x, psdDim2d_y);
  NcVar *psdFreqVar2d = ncfid.add_var(psdFreqName.c_str(), ncDouble, psdDim2d_x, psdDim2d_y);

  psdVar2d->put(&noise->psd2d[0][0], nx,ny);
  psdFreqVar2d->put(&noise->psd2dFreq[0][0], nx,ny);


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
  ncfid.add_att("dataFile", ap->getCoaddOutFile().c_str());
  ncfid.add_att("despikeSigma", ap->getDespikeSigma());
  ncfid.add_att("lowpassFilterKnee", ap->getLowpassFilterKnee());
  ncfid.add_att("timeOffset", ap->getTimeOffset());
  ncfid.add_att("timeChunk", ap->getTimeChunk());
  ncfid.add_att("neigToCut", ap->getNeigToCut());
  ncfid.add_att("cutStd", ap->getCutStd());
  ncfid.add_att("pixelSize", ap->getPixelSize());
  ncfid.add_att("approximateWeights", ap->getApproximateWeights());
  ncfid.add_att("MasterGrid[0]",masterGrid[0]);
  ncfid.add_att("MasterGrid[1]",masterGrid[1]);
 
  return 1;
}


//----------------------------- o ---------------------------------------

//makes the average noise histogram from the set of coadded noise maps
/** Uses gsl histogramming to histogram the individual histograms from
    the noise realizations.  This may be too clever and I have not checked
    the output carefully to see if it is correct.  The output is written
    to its own ncdf file.
    \todo Check output to verify that it is indeed the average histogram.
**/
bool NoiseRealizations::makeAverageHistogram(bool filtered)
{
  //are these the filtered map histograms?

  string dimname="nhist";
  string binname="histBins";
  string valname="histVals";
  if(filtered){
    dimname.assign("nhist_filteredNoise");
    binname.assign("histBins_filteredNoise");
    valname.assign("histVals_filteredNoise");
  }


  //loop through the noise histograms and generate an average histogram
  //start by finding the min and max bins
  double binmin=0., binmax=0.;
  VecDoub bins;
  int nhist=0;
  for(int i=0;i<nNoiseFiles;i++){
    NcFile ncfid = NcFile(noiseFiles[i].c_str(), NcFile::ReadOnly);
    nhist = ncfid.get_dim(dimname.c_str())->size();
    NcVar* histbinv = ncfid.get_var(binname.c_str());
    bins.resize(nhist);
    for(int j=0;j<nhist;j++) bins[j] = histbinv->as_double(j);
    if(i == 0){
      binmin = bins[0];
      binmax = bins[0];
    }
    for(int j=0;j<nhist;j++){
      if(bins[j] < binmin) binmin = bins[j];
      if(bins[j] > binmax) binmax = bins[j];
    }
  }

  //allocate memory for the histogram, bins, and values
  histBins.resize(nhist);
  histVals.resize(nhist);
  gsl_histogram *h = gsl_histogram_alloc(nhist);

  //set new bin ranges
  gsl_histogram_set_ranges_uniform(h, binmin, binmax);

  //fill up the histogram by looping back through the files
  for(int i=0;i<nNoiseFiles;i++){
    NcFile ncfid = NcFile(noiseFiles[i].c_str(), NcFile::ReadOnly);
    int nhisti = ncfid.get_dim(dimname.c_str())->size();
    if (nhisti != nhist){
      cerr << "Noise map histogram bins are not consistently sized." << endl;
      exit(1);
    }
    NcVar* histbinv = ncfid.get_var(binname.c_str());
    NcVar* histvalv = ncfid.get_var(valname.c_str());
    for(int j=0;j<nhist;j++){
      gsl_histogram_accumulate(h, histbinv->as_double(j), histvalv->as_double(j));
    }
    cerr << "Histogram for noise map " << i << "\r";
  }
  cerr << endl;
  for(int i=0;i<nhist;i++){
    histBins[i] = h->range[i];
    histVals[i] = (h->bin[i])/nNoiseFiles;
  }

  //free resources
  gsl_histogram_free(h);
  
  //write this to disk
  string fname = ap->getAvgNoiseHistFile();
  NcFile histfid = NcFile(fname.c_str(), NcFile::Replace);
  if (!histfid.is_valid()){
    cerr << "Couldn't open average histogram file for writing!\n";
    return 0;
  }

  //create dimensions
  NcDim* histDim = histfid.add_dim(dimname.c_str(), histVals.size());

  //define variables
  NcVar *valVar = histfid.add_var(valname.c_str(), ncDouble, histDim);
  NcVar *binVar = histfid.add_var(binname.c_str(), ncDouble, histDim);

  //put in the average hist
  valVar->put(&histVals[0], histVals.size());
  binVar->put(&histBins[0], histBins.size());  

  return 1;
}



//----------------------------- o ---------------------------------------

//makes the average noise psd from the set of coadded noise maps
/** Simply reads the individual noise psds from the coadded noise
    map files and averages them.  The output is written to its own
    netcdf file.
    \todo Check output to verify that it is indeed the average psd.
**/
bool NoiseRealizations::makeAveragePsd()
{
  //use the last noise file to initially size the psd vectors
  //psdFreq.assign(noise->psdFreq.size(),0.);
  //psd.assign(noise->psd.size(),0.);
	MatDoub tmpPsd;
	MatDoub tmpPsdFreq;
  //loop through the noise psd and generate an average psd
	for(int i=0;i<nNoiseFiles;i++){
		NcFile ncfid = NcFile(noiseFiles[i].c_str(), NcFile::ReadOnly);
		uint npsd = ncfid.get_dim("npsd")->size();
		uint nxpsd2d = ncfid.get_dim("nxpsd_2d")->size();
		uint nypsd2d = ncfid.get_dim("nypsd_2d")->size();
		if (i==0){
			cerr<< "PSD: "<< nxpsd2d<<","<<nypsd2d;
			psdFreq.assign(npsd,0.);
			psd.assign(npsd,0.);
			psd2d.assign(nxpsd2d,nypsd2d,0.);
			psd2dFreq.assign (nxpsd2d,nypsd2d,0.);
			tmpPsd.assign(nxpsd2d,nypsd2d,0.);
			tmpPsdFreq.assign(nxpsd2d,nypsd2d,0.);
		} else{
			if(npsd != psd.size()){
				cerr << "NoiseRealizations::makeAveragePsd(): ";
				cerr << "variation in psd lengths.  Check noise files." << endl;
				exit(1);
			}
		}
		NcVar* psdv = ncfid.get_var("psd");
		NcVar* psdfv = ncfid.get_var("psdFreq");
		NcVar* psdv2d = ncfid.get_var("psd_2d");
		NcVar* psdf2d = ncfid.get_var("psdFreq_2d");
		for(uint ip=0;ip<npsd;ip++){
			psd[ip] += psdv->as_double(ip);
			psdFreq[ip] = psdfv->as_double(ip);
		}

		psdv2d->get(&tmpPsd[0][0], nxpsd2d,nypsd2d);
		psdf2d->get(&tmpPsdFreq[0][0], nxpsd2d,nypsd2d);

		for (uint ip=0; ip<nxpsd2d; ip++)
			for (uint jp=0; jp<nypsd2d; jp++){
				psd2d[ip][jp]+=  tmpPsd[ip][jp] /nNoiseFiles;
				psd2dFreq[ip][jp]+= tmpPsdFreq[ip][jp]/nNoiseFiles;
			}
	}
	for(uint i=0;i<psd.size();i++) psd[i] /= nNoiseFiles;

	//write this to disk
	string fname = ap->getAvgNoisePsdFile();
	NcFile psdfid = NcFile(fname.c_str(), NcFile::Replace);
	if (!psdfid.is_valid()){
		cerr << "Couldn't open average psd file for writing!\n";
		return 0;
	}

	//create dimensions
	NcDim* psdDim = psdfid.add_dim("npsd", psd.size());

	//define variables
	NcVar *psdVar = psdfid.add_var("psd", ncDouble, psdDim);
	NcVar *psdfVar = psdfid.add_var("psdFreq", ncDouble, psdDim);

	//put in the average psd
	psdfVar->put(&psdFreq[0], psd.size());
	psdVar->put(&psd[0], psd.size());


	uint nxpsd = psd2d.nrows();
	uint nypsd = psd2d.ncols();
	string psdDimName = "nxpsd_2d";
	NcDim* psdDim2d_x = psdfid.add_dim(psdDimName.c_str(), nxpsd);
	psdDimName = "nypsd_2d";
	NcDim* psdDim2d_y = psdfid.add_dim(psdDimName.c_str(), nypsd);

	string psdVarName = "psd_2d";
	string psdFreqName = "psdFreq_2d";

	NcVar *psdVar2d = psdfid.add_var(psdVarName.c_str(), ncDouble, psdDim2d_x, psdDim2d_y);
	NcVar *psdFreqVar2d = psdfid.add_var(psdFreqName.c_str(), ncDouble, psdDim2d_x, psdDim2d_y);

	psdVar2d->put(&psd2d[0][0], nxpsd,nypsd);
	psdFreqVar2d->put(&psd2dFreq[0][0], nxpsd,nypsd);

	return 1;
}


//----------------------------- o ---------------------------------------

bool NoiseRealizations::calculateAverageFilteredRms(double cov)
{
  VecDoub mapRms(nNoiseFiles);

  for(int k=0;k<nNoiseFiles;k++){
    //open the noise file
    NcFile ncfid = NcFile(noiseFiles[k].c_str());
    
    //the input array dims
    int nr = ncfid.get_dim("nrows")->size();
    int nc = ncfid.get_dim("ncols")->size();
    
    //the actual noise matrix
    MatDoub mynoise(nr, nc);
    NcVar* noisev = ncfid.get_var("filteredNoise");
    noisev->get(&mynoise[0][0], nr, nc);

    //the corresponding weight matrix
    MatDoub myweight(nr, nc);
    NcVar* weightv = ncfid.get_var("filteredWeight");
    weightv->get(&myweight[0][0], nr, nc);

    //find the weight threshold corresponding to coverage cut=cov
    double myWeightCut = findWeightThreshold(myweight,cov);

    //find the rms of the pixels with weight greater than myWeightCut
    int counter=0;
    double rms=0.;
    for(int i=0;i<nr;i++)
      for(int j=0;j<nc;j++){
	if(myweight[i][j] > myWeightCut){
	  counter++;
	  rms += pow(mynoise[i][j],2);
	}
      }
    rms /= (counter);
    mapRms[k] = sqrt(rms);

    cerr << "Filtered noise rms for map " << k << " = " << mapRms[k] 
	 << "\r";
  }
  cerr << endl;

  //calculate the average rms
  double rms=0.;
  for(int k=0;k<nNoiseFiles;k++) rms += mapRms[k];
  rms /= nNoiseFiles;
  averageFilteredRms = rms;

  cerr << "Noise Realizations: Average Filtered RMS=" << averageFilteredRms << endl;

  return 1;
}



//----------------------------- o ---------------------------------------

bool NoiseRealizations::normalizeErrors(double cov)
{
  for(int k=0;k<nNoiseFiles;k++){

    //open the noise file
    NcFile ncfid = NcFile(noiseFiles[k].c_str(), NcFile::Write);
    
    //the input array dims
    int nr = ncfid.get_dim("nrows")->size();
    int nc = ncfid.get_dim("ncols")->size();
    
    //the actual noise matrix
    MatDoub myNoise(nr, nc);
    NcVar* noisev = ncfid.get_var("filteredNoise");
    noisev->get(&myNoise[0][0], nr, nc);

    //the corresponding weight matrix
    MatDoub myWeight(nr, nc);
    NcVar* weightv = ncfid.get_var("filteredWeight");
    weightv->get(&myWeight[0][0], nr, nc);

    //find the weight threshold corresponding to coverage cut=cov
    double myWeightCut = findWeightThreshold(myWeight,cov);

    //find the rms of the pixels with weight greater than myWeightCut
    double counter=0;
    double sig_of_map=0.;
    for(int i=0;i<nr;i++)
      for(int j=0;j<nc;j++){
	if(myWeight[i][j] >= myWeightCut){
	  counter++;
	  sig_of_map += pow(myNoise[i][j],2);
	}
      }
    sig_of_map /= (counter-1);
    sig_of_map = sqrt(sig_of_map);

    //find the mean sqerr (mean(1/wt) for good coverage region
    double mean_sqerr=0;
    counter=0.;
    for(int i=0;i<nr;i++)
      for(int j=0;j<nc;j++){
	if(myWeight[i][j] >= myWeightCut){
	  counter++;
	  mean_sqerr += (1./myWeight[i][j]);
	}
      }
    mean_sqerr /= counter;

    //the renormalization factor
    double nfac = (1./pow(sig_of_map,2.))*mean_sqerr;

    for(int i=0;i<nr;i++)
      for(int j=0;j<nc;j++)
	myWeight[i][j] *= nfac;

    //write the updated weight map back into the noise file
    weightv->put(&myWeight[0][0], nr, nc);
  }
  cerr << endl;
  return 1;
}


//----------------------------- o ---------------------------------------

bool NoiseRealizations::calcFilteredNoiseMapsHistogram(double cov)
{
  for(int k=0;k<nNoiseFiles;k++){

    //open the noise file
    NcFile ncfid = NcFile(noiseFiles[k].c_str(), NcFile::Write);
    
    //the input array dims
    int nr = ncfid.get_dim("nrows")->size();
    int nc = ncfid.get_dim("ncols")->size();
    
    //the actual noise matrix
    MatDoub myNoise(nr, nc);
    NcVar* noisev = ncfid.get_var("filteredNoise");
    noisev->get(&myNoise[0][0], nr, nc);

    //the corresponding weight matrix
    MatDoub myWeight(nr, nc);
    NcVar* weightv = ncfid.get_var("filteredWeight");
    weightv->get(&myWeight[0][0], nr, nc);

    //grab the row and column coordinates (physical)
    VecDoub myRow(nr);
    VecDoub myCol(nc);
    NcVar* rcp = ncfid.get_var("rowCoordsPhys");
    rcp->get(&myRow[0],nr);
    NcVar* ccp = ncfid.get_var("colCoordsPhys");
    ccp->get(&myCol[0],nc);

    //pack these into a Map object
    Map *myNoiseMap;
    myNoiseMap = new Map(string("filteredNoise"), nr, nc, pixelSize,
			 myWeight, myRow, myCol);
    for(int i=0;i<nr;i++)
      for(int j=0;j<nc;j++){
	myNoiseMap->image[i][j] = myNoise[i][j];
	myNoiseMap->weight[i][j] = myWeight[i][j];
      }

    //histogram the map and write into mapFile
    myNoiseMap->mapFile = noiseFiles[k].c_str();
    myNoiseMap->calcMapHistogram(200, cov);

   delete myNoiseMap;

  }
  cerr << endl;
  return 1;
}


//----------------------------- o ---------------------------------------

NoiseRealizations::~NoiseRealizations()
{
  delete noise;
  delete [] noiseFiles;
}
