#include <netcdfcpp.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <gsl/gsl_matrix.h>
#include <omp.h>
#include <ctime>
#include <exception>

using namespace std;

#include "nr3.h"
#include "Array.h"
#include "Detector.h"
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



int main(int nArgs, char* args[])
{
  if(nArgs > 2){
    cerr << "Macanap (Parallel version of macana using OMP): " << endl;
    cerr << "  calling syntax: ./macanap [analysis.xml] " << endl;
    cerr << "  If no xml file is given on input, apDefault.xml ";
    cerr << "  will be used.";
    exit(1);
  }
  time_t start = time(NULL);
  time_t end = 0.0;
  cerr << "Main(): setting analysis parameters." << endl;
  string apXml;
  if(nArgs == 2) apXml.assign(args[1]);
  if(nArgs == 1){
    cerr << "Main(): using default analysis xml instructions from ";
    cerr << "apDefault.xml" << endl;
    apXml.assign("apDefault.xml");
  }
  
  AnalParams* ap = new AnalParams(apXml);
  Array *array=NULL;
  TimePlace *timePlace = NULL;
  Source *source = NULL;
  Telescope *telescope = NULL;
  Observation *obs = NULL;
  AnalParams *tap = NULL;
  string ofile;
  int nFiles = ap->getNFiles();
  int filei=0;
  int tid = -1;
  Clean *cleaner=NULL;
  int *di = NULL;
  int i = 0;
  //Simulate map and object
  MapNcFile *simMap = NULL;
  MapNcFile *subMap = NULL;
  
  if (ap->getSimParams() != NULL){
    simMap = new MapNcFile(ap->getSimParams()->getMapFile());
  }
  if (ap->getDoSubtract()){
    subMap = new MapNcFile(ap->getSubtractFile());
  }
  
  
  
  //Set number of threads for parallel execution
#if defined (_OPENMP)
  omp_set_num_threads(ap->getNThreads());
  //omp_set_nested(1);
#endif
  
  
  cerr << endl;
  cerr << "--------------------- o ------------------------" << endl;
  cerr << endl;
  
  //map individual observations
  if(ap->getMapIndividualObservations())
    {
      
      //Begin the parallel execution.  Non threadsafe operations are
      //walled off with the "omp critical" pragma.
      
#pragma omp parallel shared (ap,nFiles,simMap, subMap, cerr, cout) private(tap,array, timePlace, source, telescope, obs,cleaner, filei, di,i, tid,ofile) default (none)
      {
	
	//cycle through the data files to do the reductions
#pragma omp for schedule(dynamic)
	for(filei=0;filei<nFiles;filei++){
#if defined(_OPENMP)
	  tid = omp_get_thread_num() + 1;
#else
	  tid = 0;
#endif
	  
	  //The following dataio section is walled off since netcdf is not
	  //threadsafe.
	  
	  //set the Analysis Parameters
	  tap = new AnalParams(ap);
	  tap->setDataFile(filei);
	  
	  //create an Array object
	  cerr << "Main("<<tid<<"): Creating an Array object." << endl;
	  array = new Array(tap);
	  
#pragma omp critical (dataio)
	  {
	    cerr << "Main("<<tid<<"): Populating the array with detectors." << endl;
	    array->populate();
	  }
	  //make a TimePlace
	  cerr << "Main("<<tid<<"): Creating the time and place." << endl;
	  
	  timePlace = new TimePlace(tap);
	  
	  //make a Source
	  cerr << "Main("<<tid<<"): Making the source." << endl;
	  source = new Source(tap, timePlace);
	  cerr << "Main("<<tid<<"): Made the source." << endl;
	  double *tmpGrid = tap->getMasterGridJ2000();
	  cout << "Main("<<tid<<"): Source Ra: " <<tmpGrid[0] *180.0/M_PI<< " Dec: "
	       << tmpGrid[1]*180.0/M_PI << endl;
	  
	  
	  //make a telescope
	  cerr << "Main("<<tid<<"): Making a telescope." << endl;
	  telescope=new Telescope(tap, timePlace, source);
	  
	  //set the output file
	  ofile =tap->getMapFile();
	  
	  cerr<< "Main("<<tid<<"):Critical I/O section done. " <<endl;
	  
	  
	  //get a pointer to the good detectors in the array
	  array->updateDetectorIndices();
	  di=array->getDetectorIndices();
	  
	  
	  cerr << "Main("<<tid<<"): Projecting telescope pointing "
	       << "to mastergrid tangent." << endl;
	  telescope->absToPhysEqPointing();
	  cerr << "Main("<<tid<<"): Generating pointing for each detector." << endl;
	  for(i=0;i<array->getNDetectors();i++){
	    array->detectors[di[i]].getPointing(telescope, timePlace, source);
	    array->detectors[di[i]].getAzElPointing(telescope);
	  }
	  
	  cerr << "Main("<<tid<<"): finding map bounds." << endl;
	  array->findMinMaxXY();
	  
	  //do the first round of despiking (steps 1 and 2)
	  cerr << "Main("<<tid<<"): Finding and flagging spikes." << endl;
	  for(i=0;i<array->getNDetectors();i++){
	    array->detectors[di[i]].despike(tap->getDespikeSigma());
	  }
	  array->updateDetectorIndices();
	  di=array->getDetectorIndices();
	  
	  //replace flagged data in scans with faked data (useful for pca if needed)
	  cerr << "Main("<<tid<<"): faking flagged data." << endl;
	  array->fakeFlaggedData(telescope);
	  
	  //make a fake source for psf detemination
	  cerr << "Main("<<tid<<"): making kernel timestreams." << endl;
	  for(i=0;i<array->getNDetectors();i++){
	    array->detectors[di[i]].makeKernelTimestream(telescope);
	  }
	  
	  //lowpass the data
	  cerr << "Main("<<tid<<"): Lowpassing the detector data." << endl;
	  for(i=0;i<array->getNDetectors();i++){
	    array->detectors[di[i]].lowpass(&array->digFiltTerms[0],
					    array->nFiltTerms);
	  }
	  
	  //clean out overflagged scans
	  VecBool obsFlags(array->detectors[0].getNSamples());
	  for (int j=0;j<array->detectors[0].getNSamples();j++){
	    obsFlags[j] = 0;
	    for(i=0;i<array->getNDetectors();i++)
	      if (array->detectors[di[i]].hSampleFlags[j]) obsFlags[j] = 1;
	  }
	  telescope->checkScanLengths(obsFlags,array->detectors[0].getSamplerate());
	  //generate the pointing signals
	  
	  cerr << "Main("<<tid<<"): Estimate average extinction." << endl;
	  for(i=0;i<array->getNDetectors();i++){
	    array->detectors[di[i]].estimateExtinction(array->getAvgTau());
	  }
	  //if we have a simulated map the insert it before cleaning
	  
	  if (!tap->subtractFirst && simMap !=NULL){
	    cout<<"Main("<<tid<<"): Inserting simulated signals into timestreams"<<endl;
	    SimulatorInserter *insert = new SimulatorInserter(simMap, tap->getSimParams());
	    insert->insertIntoArray(array);
	    delete insert;
	    cout<<"Main("<<tid<<"): Simulated signals inserted sucessfully at 0"<<endl;
	  }
	  
	  if (subMap !=NULL){
	    cout<<"Main("<<tid<<"): Subtract input map signals from timestreams"<<endl;
	    Subtractor *subs = new Subtractor(subMap);
	    subs->subtract(array);
	    delete subs;
	    cout<<"Main("<<tid<<"): Input map signals subtracted correctly"<<endl;
	  }
	  
	  if (tap->subtractFirst && simMap !=NULL){
	    cout<<"Main("<<tid<<"): Inserting simulated signals into timestreams"<<endl;
	    SimulatorInserter *insert = new SimulatorInserter(simMap, tap->getSimParams());
	    insert->insertIntoArray(array);
	    delete insert;
	    cout<<"Main("<<tid<<"): Simulated signals inserted sucessfully at 1"<<endl;
	  }
	  
	  
	  cleaner = CleanSelector::getCleaner(array, telescope);
	  cout << "Main("<<tid<<"): Cleaning Process Initiated." << endl;
	  cleaner->clean();
	  cout << "Main("<<tid<<"): Cleaning Done" << endl;
	  array->updateDetectorIndices();
	  di=array->getDetectorIndices();
	  delete cleaner;
	  
	  
	  cerr << "Main("<<tid<<"): Average 225GHz opacity: " << array->getAvgTau() << endl;
	  if (array->getAvgTau()<0.0){
	    cerr<< "Error opacity is less than 0.0 on file:" << tap->getMapFile()<<endl;
	    exit(-1);
	  }
	  
	  //calibrate timestreams
	  cerr << "Main("<<tid<<"): Calibrating detector signals." << endl;
	  for(i=0;i<array->getNDetectors();i++){
	    array->detectors[di[i]].estimateExtinction(array->getAvgTau());
	    array->detectors[di[i]].calibrate();
	  }
	  
	  //calculate detector weights on a scan by scan basis and then
	  //apply the same hack that is done in the IDL utilities by
	  //knocking back absurdly highly weighted detectors.
	  cerr << "Main("<<tid<<"): Generating scan weights." << endl;
	  for(i=0;i<array->getNDetectors();i++){
	    array->detectors[di[i]].calculateScanWeight(telescope);
	  }
	  
	  if (!ap->getApproximateWeights()){
	    int ndet = array->getNDetectors();
	    int nscans = telescope->scanIndex.ncols();
	    double *tmpData = new double [ndet*nscans];
	    double medianWt = 0;
	    for(i=0;i<ndet;i++){
	      for(int j=0;j<nscans;j++){
		tmpData[i*nscans+j] = array->detectors[di[i]].scanWeight[j];
	      }
	    }
	    
	    medianWt = median(tmpData,ndet*nscans);
	    
	    for(i=0;i<ndet;i++){
	      for(int j=0;j<nscans;j++){
		if(array->detectors[di[i]].scanWeight[j] > 2.*medianWt){
		  //cout<<"Bad weight on scan: "<<j<< " bolometer"<<array->detectors[di[i]].getName()<<" filename: "<<tap->getMapFile()<<endl;
		  array->detectors[di[i]].scanWeight[j] = medianWt/2.0;
		}
	      }
	    }
	    delete [] tmpData;
	  }
	  //delete []medianWt;
	  
	  //create the maps
	  cerr << "Main("<<tid<<"): Generating the observation maps." << endl;
	  obs= new Observation(tap);
	  obs->generateMaps(array, telescope);
	  
	  //write obsmaps out to a datafile.  Since netcdf is not threadsafe
	  //this must be walled off.
#pragma omp critical (dataio)
	  {
	    obs->writeObservationToNcdf(ofile);
	  }
	  //fit obs signal map's central region to gaussian
	  if (ap->getAzelMap()!=0){
	    cerr << "Main("<<tid<<"): Fitting obs signal to gaussian." << endl;
	    obs->signal->fitToGaussian();
	  }
	  
	  //calculate the obs signal map psd with coverage cut of 0.9
	  cerr << "Main("<<tid<<"): generating the psd of the obs.signal map."
	       << endl;
	  obs->signalMapPsd(ap->getCoverageThreshold());
	  
	  //histogram the obs signal map with coverage cut=0.9 and default nbins
	  cerr << "Main("<<tid<<"): generating the histogram of "
	       << "the signal obsmap." << endl;
	  obs->histogramSignal(ap->getCoverageThreshold());
	  
	  
	  
	  //this observation is analyzed so clean up all the timestream objects.
	  delete obs;
	  delete telescope;
	  delete source;
	  delete timePlace;
	  delete array;
	  delete tap;
	  cerr <<"Main("<<tid<<"): memory deallocation succeded."<<endl;
	}
	
      }
      
      if (simMap!=NULL)
	delete simMap;
      if (subMap!=NULL)
	delete subMap;
      ap->setDataFile(0);
      
      
      cerr << endl;
      cerr << "--------------------- o ------------------------" << endl;
      cerr << endl;
    }
  
  if(ap->getCoaddObservations()){
    //coadd the maps
    cerr << "Main(): Coadding Maps." << endl;
    Coaddition cmap(ap);
    cmap.coaddMaps();
    cmap.writeCoadditionToNcdf();
    
    //histogram the cmap signal map with coverage cut=0.9
    cerr << "Main(): generating the histogram of the signal cmap." << endl;
    cmap.histogramSignal(ap->getCoverageThreshold());
    
    cerr << endl;
    cerr << "--------------------- o ------------------------" << endl;
    cerr << endl;
    
    //fit signal map to gaussian
    if(ap->getFitCoadditionToGaussian()){
      cerr << "Main(): Fitting cmap to gaussian." << endl;
      cmap.signal->fitToGaussian();
    }
    
    cerr << endl;
    cerr << "--------------------- o ------------------------" << endl;
    cerr << endl;
    
    if(ap->getProduceNoiseMaps()){
      //make the noise maps
      cerr << "Main(): Making the noise realizations." << endl;
      NoiseRealizations noiseMaps(ap);
      noiseMaps.generateNoiseRealizations(&cmap);
      
      cerr << "Main(): Making average noise histogram." << endl;
      noiseMaps.makeAverageHistogram(0);
      
      cerr << "Main(): Making average noise psd." << endl;
      noiseMaps.makeAveragePsd();
      
      cerr << endl;
      cerr << "--------------------- o ------------------------" << endl;
      cerr << endl;
      
      //wiener filter
      if(ap->getApplyWienerFilter()){
    	  cerr << "Main(): apply a wiener filter to the signal map." << endl;
    	  WienerFilter wf(ap, &cmap);
    	  wf.filterCoaddition(&cmap);

    	  cerr << endl;
    	  cerr << "--------------------- o ------------------------" << endl;
    	  cerr << endl;

    	  //wiener filter the noise maps
    	  cerr << "Main(): apply a wiener filter to the noise maps" << endl;
    	  cerr << "Main(): This automatically adds them to the ncdf files." << endl;
    	  wf.filterNoiseMaps(&noiseMaps);

    	  //normalize the errors of the noise maps
    	  cerr << "Main(): normalizing errors of filtered noise maps." << endl;
    	  noiseMaps.normalizeErrors(ap->getCoverageThreshold());

    	  //find the average rms of the filtered noise maps with coverage cut of 0.9
    	  cerr << "Main(): calculating average rms of the filtered noise maps"
    			  << endl;
    	  cerr << "Main(): Using coverage cut of 0.9." << endl;
    	  noiseMaps.calculateAverageFilteredRms(ap->getCoverageThreshold());

    	  //normalize the errors in the filtered signal map with coverage cut of 0.9
    	  cmap.normalizeErrors(noiseMaps.averageFilteredRms, ap->getCoverageThreshold());
    	  cmap.writeFilteredMapsToNcdf();

    	  //calculate the flux histograms of all filtered maps
    	  cerr << "Main(): finding histogram of filtered signal map." << endl;
    	  cmap.histogramFilteredSignal(ap->getCoverageThreshold());
    	  noiseMaps.calcFilteredNoiseMapsHistogram(ap->getCoverageThreshold());

    	  cerr << "Main(): Making average noise histogram (filtered)." << endl;
    	  noiseMaps.makeAverageHistogram(1);

    	  if(ap->getFitCoadditionToGaussian()){
    		  cerr << "Main(): Fitting filtered cmap to gaussian." << endl;
    		  cmap.filteredSignal->fitToGaussian();
    	  }
      } else {
    	  cerr << "No Wiener Filter requested in analysis parameters." << endl;
      }
    }

    //Post Reduction Analysis
    if(ap->getPostReductionAnalysis()){
    	if(ap->getFindSources()){
    		cerr << "Finding sources in coadded map." << endl;
    		cmap.findSources();
    	}
    	if(ap->getCalcCompleteness()){
    		cerr << "Calculating the completeness in the coadded map." << endl;
    	}
    }
  }

  //cleanup
  delete ap;
  
  cerr << "Main(): Finished." << endl;
  end = time(NULL);
  cerr << "Total Reduction time: "<<(double)(end-start)/((double)3600.0)
       <<"hrs."<<endl;
  return 0;
}
