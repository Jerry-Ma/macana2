#include <netcdfcpp.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <gsl/gsl_matrix.h>
#include <omp.h>
#include <ctime>
#include <exception>
#include <sstream>
#include <vector>

using namespace std;

#include "astron_utilities.h"
#include "nr3.h"
#include "Array.h"

#include "Map.h"
#include "Observation.h"
#include "Coaddition.h"
#include "NoiseRealizations.h"
#include "Telescope.h"
#include "TimePlace.h"
#include "Source.h"
#include "AnalParams.h"
#include "vector_utilities.h"
#include "CleanPCA.h"
#include "CleanBspline.h"
#include "CleanSelector.h"
#include "WienerFilter.h"
#include "Clean2dStripe.h"
#include "MapNcFile.h"
#include "SimulatorInserter.h"
#include "Subtractor.h"
#include "nr3.h"

int main(int nArgs, char* args[])
{
  //selects the ap.xml file and stores in apXml
  if(nArgs > 2)
  {
    cerr << "Beammap: " << endl;
    cerr << "  calling syntax ./beammap [analysis.xml] " << endl;
    cerr << " if no xml file is given on input, apbDefault.xml ";
    cerr << " will be used." << endl;;
    return 1;
  }
  string apXml;
  if(nArgs == 2) apXml.assign(args[1]);
  if(nArgs == 1){
    cerr << "Main(): using default analysis xml instructions from ";
    cerr << "apbDefault.xml" << endl;
    apXml.assign("apbDefault.xml");
  }

  //setting up the AnalParams, TimePlace, Array, Source, and Telescope
  AnalParams* ap = new AnalParams(apXml, 1);
  int nFiles = ap->getNFiles();
  
  //begin loop over input files
  for(int fileNum = 0; fileNum < nFiles; fileNum++){  
    cerr << "starting file number " << fileNum << endl;
    ap->setDataFile(fileNum);
    TimePlace* timePlace = new TimePlace(ap);
    Array* array = new Array(ap);
    array->populate();
    Source* source = new Source(ap, timePlace);
    Telescope* telescope = new Telescope(ap, timePlace, source);
    array->updateDetectorIndices();
    int *di = array->getDetectorIndices();

    //setting up physical coordinate system
    telescope->absToPhysEqPointing();
    for(int i = 0; i < array->getNDetectors(); i++){
      array->detectors[di[i]].getPointing(telescope, timePlace, source);
      array->detectors[di[i]].getAzElPointing(telescope);
    }
    array->findMinMaxXY();

    //despiking at detector level, detector by detector
    for(int i = 0; i < array->getNDetectors(); i++){
      array->detectors[di[i]].despike(ap->getDespikeSigma());
    }
    array->updateDetectorIndices();
    di = array->getDetectorIndices();

    //replace flagged data in scans with faked data (useful for pca if needed)
    array->fakeFlaggedData(telescope);

    //make a fake source for psf determination
    for(int i = 0; i <array->getNDetectors(); i++){
      array->detectors[di[i]].makeKernelTimestream(telescope);
    }

    //lowpass the data
    for(int i = 0; i < array->getNDetectors(); i++){
      array->detectors[di[i]].lowpass(&array->digFiltTerms[0], array->nFiltTerms);
    }

    VecBool obsFlags(array->detectors[0].getNSamples());
    for(int i = 0; i < array->detectors[0].getNSamples(); i++){
      obsFlags[i] = 0;
      for(int j = 0; j < array->getNDetectors(); j++){
        if(array->detectors[di[j]].hSampleFlags[i]) obsFlags[i] = 1;
      }
    }
    telescope->checkScanLengths(obsFlags, array->detectors[0].getSamplerate());

    //generate the pointing signals
    for(int i=0;i<array->getNDetectors();i++){
      array->detectors[di[i]].estimateExtinction(array->getAvgTau());
    }

    //create vectors to store future fit parameters for iterative cleaning
    MatDoub fitParams(array->getNDetectors(), 7, 0.);
    MatDoub previousFitParams(array->getNDetectors(), 7, 1.0);

    //create the vector that stores whether to terminate iteration for each detector
    VecDoub needsIteration(array->getNDetectors(), 1.0);

    //store original hValues for iterative cleaning
    int nSamples = array->detectors[di[0]].getNSamples();
    MatDoub originalhVals(array->getNDetectors(), nSamples, 0.);
    for(int i = 0; i < array->getNDetectors(); i++){
      for(int j = 0; j < nSamples; j++){
        originalhVals[i][j] = array->detectors[di[i]].hValues[j];
      } 
    }

    //grab the iteration loop cap and the percent change cutoff
    int iteration = 0;
    int cap = ap->getCleanIterationCap();
    double cutOff = ap->getCleanIterationCutoff();

    //begin the iterative cleaning
    while(iteration < cap){
      cerr << "file number: " << fileNum << endl; 
      cerr << "iteration: " << iteration << endl;
	
      //replace previously cleaned data with original data
      //subtract out previous best guess for signal from original data before cleaning
      // (don't do this on the first iteration)
      if(iteration > 0){
        for(int i = 0; i < array->getNDetectors(); i++){
          for(int j = 0; j < nSamples; j++){
            array->detectors[di[i]].hValues[j] = originalhVals[i][j];
          }
          fitParams[i][1] = -fitParams[i][1];
          array->detectors[di[i]].addGaussian(fitParams, i);
        }
      }

      //cleaning
      Clean* cleaner = CleanSelector::getCleaner(array, telescope);
      cleaner->clean();
      array->updateDetectorIndices();
      di = array->getDetectorIndices();
      delete cleaner;

      if(array->getAvgTau() < 0){
        cerr << "Error opacity is less than 0.0 on file" << endl;
        return 1;
      }

      //add back in previously subtracted out signal
      if(iteration > 0){
        for(int i = 0; i < array->getNDetectors(); i++){
          fitParams[i][1] = -fitParams[i][1];
          array->detectors[di[i]].addGaussian(fitParams, i);
        }
      }

      //following the cleaning, generate the beammaps
      Observation* obs = new Observation(ap);
      obs->generateBeammaps(array, telescope);
      
      //grab each signal and weight map from obs and use fitToGaussian to generate fit parameters
      //on subsequent iterations, only run the fit on detectors that show percent change > cutoff
      int nrows = obs->nrows;
      int ncols = obs->ncols;
      double pixelSize = obs->pixelSize;
      for(int i = 0; i < array->getNDetectors(); i++){
        if(needsIteration[i] == 1){
          VecDoub pp;
          string mname;
          MatDoub sig(nrows, ncols, 0.);
          MatDoub wtt(nrows, ncols, 0.);
          for(int j = 0; j < nrows; j++){
      	    for(int k = 0; k < ncols; k++){
              sig[j][k] = obs->beammapSignal[i][j*ncols+k];
              wtt[j][k] = obs->beammapWeight[i][j*ncols+k];
            }
          }
          stringstream sstm;
          sstm << "beammapSignal" << i;
          Map* signal = new Map(mname.assign(sstm.str()), nrows, ncols, pixelSize, wtt,
                              obs->rowCoordsPhys, obs->colCoordsPhys);
          signal->image = sig;
          signal->weight = wtt;
          signal->mapFile = "";
        
          //the "120" parameter is an optional parameter that sets the size of the area to be fit"
          signal->fitToGaussian(pp, 120);
          for(int j = 0; j < 7; j++){
            previousFitParams[i][j] = fitParams[i][j];
            fitParams[i][j] = pp[j];
          }
          delete signal;
        }
      }
   
      //evaluate the percent change from the previous fit
      int keepGoing = 0;
      for(int i = 0; i < array->getNDetectors(); i++){
        if(needsIteration[i] == 1){
          int detNeedsIter = 0;
          for(int j = 1; j < 6; j++){
            double pctDiff=abs((fitParams[i][j] - previousFitParams[i][j]) /previousFitParams[i][j]);
            if(pctDiff > cutOff){
              cerr << "detector number: " << i << endl;
              cerr << "parameter number: " << j << endl;
              cerr << "previous: " << previousFitParams[i][j] << endl;
              cerr << "current: " << fitParams[i][j] << endl;
              cerr << "percent difference: " << pctDiff * 100 << "%" << endl;
              keepGoing = 1;
              detNeedsIter = 1;
              break;
            }
          }
          
          //keep track of which detectors need another fit due to a significant change in parameters
          needsIteration[i] = detNeedsIter;
        }
      }
      
      //either the cap or the cutoff is met
      //the maps are optionally written to netCDF and the loop terminates
      if(iteration + 1 == cap || keepGoing == 0){
        cerr << "stopping" << endl;
        
        if(ap->getWriteBeammapToNcdf()){
          string mapFile = ap->getOutBeammapNcdf();
          
          cerr << "writing maps to " << mapFile << endl;
          /* obs->writeBeammapsToNcdf(mapFile); */
          obs->writeBeammapsToFits(mapFile);
          
          cerr << "writing fit parameters to " << mapFile << endl;
          obs->writeFitParamsToNcdf(mapFile, fitParams);
        }

        //for the purpose of stopping the terminating loop, but still hitting the delete obs line. 
        iteration = cap;
      }
      else{
        cerr << "iterating" << endl;
      }
      delete obs;
      cerr << "observation successfully deallocated" << endl;
      iteration++;
    }

    //fcf is calculated using the beammapSourceFlux
    for(int i=0;i<fitParams.nrows();i++){
      double amplitude = fitParams[i][1];
      double sourceBrightness = ap->getBeammapSourceFlux();
      double extinction = array->detectors[di[i]].getExtinction();
      double responsivity = array->detectors[di[i]].responsivity;
      double fecGain = array->detectors[di[i]].getFecGain();
      double fcf = sourceBrightness * responsivity * fecGain * 1000 / (amplitude * extinction);
      array->detectors[di[i]].setFcf(fcf);
    }

    //the sensitivity calculation currently suffers from both errors and memory leaks
    for(int i=0;i<array->getNDetectors();i++){
      array->detectors[di[i]].calibrate();
      //array->detectors[di[i]].calculateSensitivity(telescope);
    }

    cerr << "generating bstats file" << endl;
    
    //bstatsFile is the text file to which all of these parameters are written
    ofstream bstatsFile;
    bstatsFile.open(ap->getOutBeammapInfo());
    
    bstatsFile << "Bolo     Az FWHM         El FWHM    	    Az OFF           El OFF          Gain       Sensitivity" << endl;
    bstatsFile << "         (arcsec)        (arcsec)        (arcsec)         (arcsec)        (mJy/nW)   (mJy s^(1/2))" << endl;
    bstatsFile << "**************************************************************************************************" << endl;
    int currDetectorIndex = 0;
    //The offsets are multiplied by -1 for the purpose of matching the IDL utilities
    for(int i = 0; i < array->nBolos; i++){
      bstatsFile << array->detectorNames[i] << "    ";
      if(currDetectorIndex < array->getNDetectors() && 
      array->detectors[di[currDetectorIndex]].getId() == i){
        stringstream sstm;
        bstatsFile << fitParams[currDetectorIndex][2]/TWO_PI*360.*3600.*2.3548 << "          ";
        bstatsFile << fitParams[currDetectorIndex][3]/TWO_PI*360.*3600.*2.3548 << "          ";
        bstatsFile << -1 * fitParams[currDetectorIndex][4]/TWO_PI*360.*3600. << "         ";
        bstatsFile << -1 * fitParams[currDetectorIndex][5]/TWO_PI*360.*3600. << "      ";
        bstatsFile << -1 * array->detectors[di[currDetectorIndex]].getFcf() << "      ";
        //bstatsFile << array->detectors[di[currDetectorIndex]].getSensitivity();
        currDetectorIndex++;
      }
      bstatsFile << endl;

    }
    bstatsFile.close();
    
    //memory deallocation
    delete timePlace;
    delete array;
    delete source;
    delete telescope;
    cerr << "file number " << fileNum << " memory successfully deallocated" << endl; 
    
  }
  //memory deallocation
  delete ap;

  cerr << "finished" << endl;
  return 0;
}

