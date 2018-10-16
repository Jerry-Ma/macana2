#include <netcdfcpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <cstdio>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <omp.h>
using namespace std;

#include "nr3.h"
#include "tinyxml2.h"
#include "AnalParams.h"
#include "Array.h"
#include "Detector.h"
#include "vector_utilities.h"
#include "SBSM.h"


///Array constructor
/** The Array constructor initializes all members.
    In addition, it:
    1) pulls in all relevant data from the boloststs.xml file.
    2) identifies the reference pixel
    3) generates the terms for the digital lowpass filter
**/
Array::Array(AnalParams* anal)
{
int samplerate;
#pragma omp critical (dataio)
{
  //initialize filenames and ap
  ap = anal;
  dataFile = ap->getDataFile();
  bolostatsFile = ap->getBolostatsFile();
  string observatory=ap->getObservatory();

  //get ancillary daa from bolostats xml file
  tinyxml2::XMLDocument bolostats;
  bolostats.LoadFile(bolostatsFile);
 
  //todo - add a check that the bolostatsFile exists - otherwise it
  //crashes here with a segmentation fault

  //get nDetectors: ie, count those with goodFlag=1
  const char* tmp;
  tmp = bolostats.FirstChildElement("nBolos")->
    FirstChildElement("value")->GetText();
  nBolos = atol(tmp);

  //initialize list of detector names
  detectorNames = new string[nBolos];

  //count up those with goodFlag=1
  int count=0;
  char dId [50];
  int n;
  bool goodFlag;
  for(int i=0;i<nBolos;i++){
    n = sprintf(dId, "d%d",i);
    //check the vale of goodFlag
    tmp = bolostats.FirstChildElement(dId)->
      FirstChildElement("goodflag")->
      FirstChildElement("value")->GetText();
    goodFlag = (strcmp(tmp,"1") == 0) ? 1 : 0;
    if(goodFlag) count++;
    detectorNames[i] = bolostats.FirstChildElement(dId)->
      FirstChildElement("name")->
      FirstChildElement("value")->GetText();
    detectorNames[i] = detectorNames[i].substr(18,10);
  }
  nDetectors = count;
  cerr << "Array:: Found " << nDetectors;
  cerr << " detectors with goodFlag=1" << endl;

  //go back through and make a list of detector id's with goodFlag=1
  detectorIdsStr = new string[nDetectors];
  detectorIds.resize(nDetectors);
  count=0;
  for(int i=0;i<nBolos;i++){
    n = sprintf(dId, "d%d",i);
    //get the goodFlag
    tmp = bolostats.FirstChildElement(dId)->
      FirstChildElement("goodflag")->
      FirstChildElement("value")->GetText();
    goodFlag = (strcmp(tmp,"1") == 0) ? 1 : 0;
    //if goodFlag=1 then collect detector id
    if(goodFlag){
      stringstream ss;
      ss << dId;
      ss >> detectorIdsStr[count];
      detectorIds[count]=i;
      count++;
    }
  }

  //initialize other attributes
  tau = 0;

  //get ref pix info
  NcFile ncfid(dataFile, NcFile::ReadOnly);
  string refpixfilename;
  bool LMT=0;
  bool ASTE=0;
  char *strRefPixAtt;


  if(observatory.compare("LMT") == 0){
    LMT=1;
    //get the reference pixel from the netcdf file data
    NcVar* refPixVar=ncfid.get_var("Header.AztecBackend.ReferenceChannel");
    strRefPixAtt = refPixVar->as_string(0);
    refpixfilename.assign(strRefPixAtt);
    refpixfilename.erase(refpixfilename.find(" "));
    delete [] strRefPixAtt;
  } else if(observatory.compare("ASTE") == 0){
    ASTE=1;
    //get the reference pixel from the netcdf file header
    NcAtt* refpixAtt = ncfid.get_att("reference_channel");
    char *strRefPixAtt = refpixAtt->as_string(0);
    refpixfilename.assign(strRefPixAtt);
    delete [] strRefPixAtt;
    delete refpixAtt;
  }    

  //find the corresponding id matching the list of good detectors
  //note we stop at first ocurrance due to name degeneracy
  for(int i=0;i<nDetectors;i++){
    n = sprintf(dId, "d%d",i);
    tmp = bolostats.FirstChildElement(dId)->
      FirstChildElement("name")->
      FirstChildElement("value")->GetText();
    string stmp(tmp);
    uint ifind=stmp.find(refpixfilename);
    if(ifind < stmp.size()){
      refPixId=i;
      break;
    }
  }
  refPix = strdup(string(tmp).c_str());
  cerr << "Array:: Reference Pixel: " << refPix;
  cerr << " with id: " << refPixId << endl;


  //create a vector with the indices of the detectors with
  //goodFlag=1.  By definition all of the detectors we have
  //at this point have goodFlag=1
  detectorInd.resize(nDetectors);
  for(int i=0;i<nDetectors;i++) detectorInd[i] = i;


  //get the terms for the lowpass filter for the detectors
  //we need the samplerate from the file but again this is observatory
  //dependent
  string dimname;

  if(LMT) dimname.assign("Data.AztecBackend.time_xlen");
  if(ASTE) dimname.assign("mfpersf");
  NcDim* mfPerSf=ncfid.get_dim(dimname.c_str());
  samplerate = mfPerSf->size();
  ncfid.close();
}
  nFiltTerms=32;   //the number of digital filter terms
  digFiltTerms.resize(2*nFiltTerms+1);
  double nyquist = samplerate/2.;
  double fLow = 0.;
  double fHigh = ap->getLowpassFilterKnee()/nyquist;
  double aGibbs = 50.;
  cerr << "Array:: Creating digital filter params. " << endl;
  digitalFilter(fLow, fHigh, aGibbs, nFiltTerms, &digFiltTerms[0]);
}


//----------------------------- o ---------------------------------------

///Populates the array object with detector values.
/** Array::populate() fetches detector data from the netcdf file.
    Only detectors with goodFlag=1 in the bolostats.xml file are
    included.
**/
bool Array::populate()
{
  detectors = new Detector[nDetectors];
  for(int i=0;i<nDetectors;i++){
//	  #pragma omp critical (dataio)
//	  {
		  detectors[i].initialize(ap, detectorIdsStr[i].c_str());
//	  }
    if(detectors[i].getNSamples() != detectors[0].getNSamples()){
      cerr << "Array::populate(): variation in";
      cerr << "detector nSamples." << endl;
      exit(1);
    }    
  }
  nSamples = detectors[0].getNSamples();
  cerr << "Array::populate(): Populated array with " << nDetectors;
  cerr << " detectors." << endl;
  return 1;
}


//----------------------------- o ---------------------------------------


///Replaces flagged data, which remains flagged as bad, with faked data.
/** Replace flagged data with fake data (defined below).  Note: this fake
    data is NOT used in mapmaking, it remains flagged as bad.
    The faked data is built in a similar way to that of the 
    aztec_idl_utility but uses the array averaged standard deviation
    rather than that detector's standard deviation.
    Overall, the faked data is composed of:
    1) random noise based on calculated variance of unflagged detectors.
    2) The linear baseline of spikey bolometer.
    3) A sky model made from the average of spike free bolometers.
**/
bool Array::fakeFlaggedData(Telescope* telescope)
{
  //create fake data with following signal/noise sources:
  // 1) random noise based on calculated variance
  // 2) The linear baseline of spikey bolometer
  // 3) sky model (average of spike free bolometers

  //our pointer to the vector of good detectors
  int* di=getDetectorIndices();
  
  //do this scan by scan
  int nScans = telescope->scanIndex.ncols();
  for(int k=0;k<nScans;k++){
    int si=telescope->scanIndex[0][k];
    int ei=telescope->scanIndex[1][k];
    bool useAllDet=0;      //use all detectors regardless of spikes
    VecBool hasSpikes(nDetectors);
    for(int i=0;i<nDetectors;i++) hasSpikes[i]=0;    

    //figure out if there are any flag-free detectors
    int nFlagged=0;
    for(int i=0;i<nDetectors;i++){
      for(int j=si;j<=ei;j++){
	if (detectors[di[i]].hSampleFlags[j]==0){
	  nFlagged++;
	  hasSpikes[i]=1;
	  break;
	}
      }
    }
    if(nFlagged > 0.5*nDetectors){
      cerr << "Array::fakeFlaggedData(): scan #" << k;
      cerr << " has fewer than 50% of the detectors flag-free." << endl;
      cerr << "Using entire array (with spikes included) as sky model.";
      cerr << endl;
      useAllDet=1;
    }

    //go through detectors with hSampleFlags=0 and fake data
    for(int i=0;i<nDetectors;i++){
      if(hasSpikes[i]){
	//must assume there are multiple flagged regions
	//flags at beginning of scans are special cases

	//condition sampleflags so that if there is a spike we can make
	//one long flagged or unflagged region.
	//first do spikes from 0 to 1
	for(int j=si+1;j<ei;j++){
	  if(detectors[di[i]].hSampleFlags[j] == 1 &&
	     detectors[di[i]].hSampleFlags[j-1] == 0 &&
	     detectors[di[i]].hSampleFlags[j+1] == 0)
	    detectors[di[i]].hSampleFlags[j] = 0;
	}
	//now spikes from 1 to 0
	for(int j=si+1;j<ei;j++){
	  if(detectors[di[i]].hSampleFlags[j] == 0 &&
	     detectors[di[i]].hSampleFlags[j-1] == 1 &&
	     detectors[di[i]].hSampleFlags[j+1] == 1)
	    detectors[di[i]].hSampleFlags[j] = 1;
	}
	//and the first and last samples 
        detectors[di[i]].hSampleFlags[si] = 
	  detectors[di[i]].hSampleFlags[si+1];
        detectors[di[i]].hSampleFlags[ei] = 
	  detectors[di[i]].hSampleFlags[ei-1];

	//count up the number of flagged regions of data in the scan
	int nFlaggedRegions=0;
	if(detectors[di[i]].hSampleFlags[ei] == 0){
	  nFlaggedRegions++;
	}
	for(int j=si+1;j<=ei;j++)
	  if(detectors[di[i]].hSampleFlags[j]-detectors[di[i]].hSampleFlags[j-1] > 0)
	    nFlaggedRegions++;
	//sanity check: nFlaggedRegions should not equal zero unless above
	//fixes to sampleflags did all that needed to be done.
	if(nFlaggedRegions==0) break;


	//find the start and end index for each flagged region
	VecInt siFlags(nFlaggedRegions,-1);
	VecInt eiFlags(nFlaggedRegions,-1);
	int count=0;
	int j=si;
	while (j<ei){
	  if(detectors[di[i]].hSampleFlags[j] == 0){
	    int jstart=j;
	    int sampcount=0;
	    while(detectors[di[i]].hSampleFlags[j] == 0 && j<=ei){
	      sampcount++;
	      j++;
	    }
	    if(sampcount > 1){
	      	    siFlags[count]=jstart;
		    eiFlags[count] = j-1;
		    count++;
	    } else {
	      j++;
	    }
	  } else {
	    j++;
	  }
	}

	if(count != nFlaggedRegions){
	  cerr << "Scan k=" << k << ", Det i=" << i << endl;
	  cerr << "Array::fakeFlaggedData(): count=" << count;
	  cerr << " but it should be equal to nFlaggedRegions=";
	  cerr << nFlaggedRegions << endl;
	  VecDoub mysf(ei-si+1);
	  for(j=si;j<=ei;j++) mysf[j-si] = (double) detectors[di[i]].hSampleFlags[j];
	  writeVecOut("testing/sf.txt", &mysf[0],mysf.size());
	  exit(1);
	}


	//now loop on the number of flagged regions for the fix
	VecDoub xx(2);
	VecDoub yy(2);
	for(int j=0;j<nFlaggedRegions;j++){

	  //FLAGGED DETECTOR
	  //determine linear baseline for flagged region
	  //but use flat dc level if flagged at endpoints
	  int nFlags = eiFlags[j]-siFlags[j];
	  VecDoub linOffset(nFlags);
	  if(siFlags[j] == si){
	    //first sample in scan is flagged so offset is flat
	    //with the value of the last sample in the flagged region
	    for(int l=0;l<nFlags;l++) 
	      linOffset[l] = detectors[di[i]].hValues[eiFlags[j]+1];
	  } else if(eiFlags[j] == ei){
	    //last sample in scan is flagged so offset is flat
	    //with the value of the first sample in the flagged region
	    for(int l=0;l<nFlags;l++) 
	      linOffset[l] = detectors[di[i]].hValues[siFlags[j]-1];
	  } else {
	    //in this case we linearly interpolate between the before
	    //and after good samples
	    xx[0] = siFlags[j]-1;
	    xx[1] = eiFlags[j]+1;
	    yy[0] = detectors[di[i]].hValues[siFlags[j]-1];
	    yy[1] = detectors[di[i]].hValues[eiFlags[j]+1];
	    VecDoub xLinOffset(nFlags);
	    for(int l=0;l<nFlags;l++) xLinOffset[l]=siFlags[j]+l;
	    interpolateLinear(xx, yy, 2,
			      &xLinOffset[0], &linOffset[0], nFlags);
	  }


	  //ALL NONFLAGGED DETECTORS
	  //do the same thing but for all detectors without spikes
	  //count up spike-free detectors and store their values
	  double detCount=0;
	  for(int ii=0;ii<nDetectors;ii++) if(!hasSpikes[ii]) detCount++;
	  if(useAllDet) detCount=nDetectors;
	  //storage
	  MatDoub det(detCount,nFlags);     //detector values
	  VecDoub res(detCount);            //detector responsivities
	  int c=0;
	  for(int ii=0;ii<nDetectors;ii++){
	    if(!hasSpikes[ii] || useAllDet){
	      for(int l=0;l<nFlags;l++)
		det[c][l] = detectors[ii].hValues[siFlags[j]+l];
	      res[c] = detectors[ii].responsivity;
	      c++;
	    }
	  }
	  //for each of these go through and redo the offset bit
	  MatDoub linOffsetOthers(detCount,nFlags);
	  if(siFlags[j] == si){
	    //first sample in scan is flagged so offset is flat
	    //with the value of the last sample in the flagged region
	    for(int ii=0;ii<detCount;ii++) for(int l=0;l<nFlags;l++) 
	      linOffsetOthers[ii][l] = det[ii][0];
	  } else if(eiFlags[j] == ei){
	    //last sample in scan is flagged so offset is flat
	    //with the value of the first sample in the flagged region
	    for(int ii=0;ii<detCount;ii++) for(int l=0;l<nFlags;l++) 
	      linOffsetOthers[ii][l] = det[ii][nFlags-1];
	  } else {
	    //in this case we linearly interpolate between the before
	    //and after good samples
	    VecDoub xLinOffset(nFlags);
	    VecDoub tmpVec(nFlags);
	    for(int l=0;l<nFlags;l++) xLinOffset[l]=siFlags[j]+l;
	    xx[0] = siFlags[j]-1;
	    xx[1] = eiFlags[j]+1;
	    for(int ii=0;ii<detCount;ii++){
	      yy[0] = det[ii][0];
	      yy[1] = det[ii][nFlags-1];
	      interpolateLinear(xx, yy, 2,
				&xLinOffset[0], &tmpVec[0], nFlags);
	      for(int l=0;l<nFlags;l++) linOffsetOthers[ii][l]=tmpVec[l];
	    }
	  }

	  //subtract off the linear offset from the spike free bolos
	  for(int ii=0;ii<detCount;ii++)
	    for(int l=0;l<nFlags;l++)
	      det[ii][l] -= linOffsetOthers[ii][l];

	  //scale det by responsivities and average to make sky model
	  VecDoub skyModel(nFlags,0.);
	  for(int ii=0;ii<detCount;ii++)
	    for(int l=0;l<nFlags;l++)
	      skyModel[l] += det[ii][l]/res[ii];
	  for(int l=0;l<nFlags;l++) skyModel[l] /= detCount;

	  //find mean standard deviation of sky-model subtracted detectors
	  //this is a different approach than that taken in IDL pipeline
	  //but I think it's a good one considering the PCA to come later.
	  //This is stored in the standard deviation of flag free detectors.
	  VecDoub stdDevFF(detCount,0.);
	  Doub tmpMean;
	  VecDoub tmpVec(nFlags);
	  for(int ii=0;ii<detCount;ii++){
	    tmpMean=0;
	    for(int l=0;l<nFlags;l++) tmpVec[l] = det[ii][l]/res[ii]-skyModel[l];
	    for(int l=0;l<nFlags;l++) tmpMean += tmpVec[l];
	    tmpMean /= nFlags;
	    for(int l=0;l<nFlags;l++) stdDevFF[ii] += pow((tmpVec[l]-tmpMean),2);
	    stdDevFF[ii] = (nFlags ==1.) ? 
	      stdDevFF[ii]/nFlags : stdDevFF[ii]/(nFlags-1.);
	  }

	  Doub meanStdDev=0;
	  for(int ii=0;ii<detCount;ii++) meanStdDev += sqrt(stdDevFF[ii]);
	  meanStdDev /= detCount;

	  //the noiseless fake data is then the sky model plus the
	  //flagged detectors linear offset
	  VecDoub fake(nFlags);
	  for(int l=0;l<nFlags;l++) 
	    fake[l] = skyModel[l]*detectors[di[i]].responsivity + linOffset[l];

	  //add noise to the fake signal
	  meanStdDev *= detectors[di[i]].responsivity;   //put back in volts

	  //replace detector values with fake signal
	  for(int l=0;l<nFlags;l++)
	    detectors[di[i]].hValues[siFlags[j]+l] = fake[l];
	}
      }
    }
  }//index k, looping on scans

  return 1;
}


//----------------------------- o ---------------------------------------



//----------------------------- o ---------------------------------------


/**This method searches the min and max values of all the detectors in 
   the array in order to find the upper and lower bounds on the pointing
   in each coordinate.
**/
bool Array::findMinMaxXY()
{
  //our pointer to the vector of good detectors
  int* di=getDetectorIndices();

  double tmpmxx = detectors[di[0]].getMaxX();
  double tmpmxy = detectors[di[0]].getMaxY();
  double tmpmnx = detectors[di[0]].getMinX();
  double tmpmny = detectors[di[0]].getMinY();

  for(int i=0;i<nDetectors;i++){
    if(detectors[di[i]].getMaxX() > tmpmxx) tmpmxx=detectors[di[i]].getMaxX();
    if(detectors[di[i]].getMinX() < tmpmnx) tmpmnx=detectors[di[i]].getMinX();
    if(detectors[di[i]].getMaxY() > tmpmxy) tmpmxy=detectors[di[i]].getMaxY();
    if(detectors[di[i]].getMinY() < tmpmny) tmpmny=detectors[di[i]].getMinY();
  }

  minX = tmpmnx;
  minY = tmpmny;
  maxX = tmpmxx;
  maxY = tmpmxy;

  return 1;
}


//----------------------------- o ---------------------------------------

///Updates the vector of detector indices with only those with gF=1
/** Updates the vector of detector indices.  We only want to keep
    those with goodFlag=1.  Note that we also need to update
    nDetectors.
**/
bool Array::updateDetectorIndices()
{
  int dcount=0;
  tau = 0;
  for(int i=0;i<nDetectors;i++){
    if(detectors[detectorInd[i]].goodFlag) dcount++;
  }
  
  //throw a message in the event we're throwing out a detector
  if(dcount < nDetectors){
    cerr << "Array::updateDetectorIndices(): ";
    cerr << "Throwing out " << nDetectors-dcount << " detectors ";
    cerr << "with goodFlag != 1." << endl;
  }

  VecInt t(dcount);
  int it=0;
  for(int i=0;i<nDetectors;i++){
    if(detectors[detectorInd[i]].goodFlag){
      t[it]=detectorInd[i];
      it++;
      tau+= detectors[detectorInd[i]].estimatedTau;
    }
  }
  tau/=dcount;
  nDetectors = dcount;
  detectorInd.resize(dcount);
  for(int i=0;i<dcount;i++) detectorInd[i]=t[i];

  updateRefBoloIndex();
  return 1;
}


//----------------------------- o ---------------------------------------


int* Array::getDetectorIndices()
{
  return &detectorInd[0];
}


//----------------------------- o ---------------------------------------


int Array::getNDetectors()
{
  return nDetectors;
}

int Array::getNSamples()
{
  return nSamples;
}

double Array::getMaxX()
{
  return maxX;
}

double Array::getMinX()
{
  return minX;
}

double Array::getMaxY()
{
  return maxY;
}

double Array::getMinY()
{
  return minY;
}

AnalParams* Array::getAp(){
  return ap;
}

bool Array::fakeAtmData(bool addToKernel) {

	double atmFreq= 0.05;
	size_t nb;
	this->updateDetectorIndices();

	int *di = getDetectorIndices();
	nb  = size_t(round(detectors[di[0]].getNSamples()/detectors[di[0]].getSamplerate()*atmFreq));
	SBSM *bspline = new SBSM(4,nSamples, nb);
	for (size_t i=0; i<(size_t)getNDetectors();i++)
		detectors[di[i]].setAtmTemplate(bspline->fitData(detectors[di[i]].hValues));

	delete bspline;

	if (addToKernel){
		for (size_t i=0; i<(size_t)getNDetectors();i++)
			for (size_t j=0; j < (size_t)detectors[di[i]].getNSamples(); j++)
				detectors[di[i]].hKernel[j]+= detectors[di[i]].atmTemplate[j];

	}

	return true;
}

double Array::getAvgTau() {
	return tau;
}

void Array::updateRefBoloIndex() {
	string signalName ="h2b2";
	if (ap->getObservatory() == "LMT")
		signalName = "Data.AztecBackend." + signalName;
	for (size_t i=0; i< detectorInd.size(); i++){
		if (detectors[detectorInd[i]].getName() == signalName){
			refBoloIndex=i;
			break;
		}
	}

}

size_t Array::getRefBoloIndex() {
	return refBoloIndex;
}

//----------------------------- o ---------------------------------------


Array::~Array()
{
  //delete refPix;
  free ((void *) refPix);
  delete [] detectorIdsStr;
  delete [] detectors;
}
