#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <omp.h>

using namespace std;

#include "nr3.h"
#include "CleanPCA.h"
#include "Detector.h"
#include "GslRandom.h"
#include "vector_utilities.h"



CleanPCA::CleanPCA(Array *dataArray,Telescope *telescope) : 
  Clean(dataArray,telescope){

  this->adaptive = false;
  this->damping = 0.0;
  
}



bool CleanPCA::clean(){
  //this cleaning is done scan by scan
  int nScans = telescope->scanIndex.ncols();
  Telescope *tel = telescope;
  //update and get the array of good detectors
  Array *dataArray = this->dataArray; 
  //dataArray->updateDetectorIndices();
  int* di=dataArray->getDetectorIndices();

  double neigToCut = ap->getNeigToCut();
  double cutStd = ap->getCutStd();
  int si=0;
  int ei=0;
  int npts=0;
  int k=0;
  int i=0;
  int j=0;
  int nDetectors = dataArray->getNDetectors();
  
  gsl_matrix* det = NULL;
  gsl_matrix* ker = NULL;
  gsl_matrix* flag = NULL;
  gsl_matrix* denom = NULL;
  gsl_matrix* pcaCorr = NULL;
  
  double mn=0.;
  double mk=0.;

  dataArray->eVectors.resize(nDetectors,nDetectors);
  
   #pragma omp  parallel shared (dataArray, tel, di, nScans, nDetectors, neigToCut, cutStd,cerr)\
			  private (k,i,j, si, ei, npts, det, ker,flag, \
				   mn,mk, denom, pcaCorr) default (none)
  {
#pragma omp for schedule(dynamic)
    for(k=0;k<nScans;k++){
      
      si=tel->scanIndex[0][k];
      ei=tel->scanIndex[1][k]+1;
      npts = ei-si;
      
      
      //allocate a gsl_matrix to store the values
      det = gsl_matrix_alloc(nDetectors,npts);
      ker = gsl_matrix_alloc(nDetectors,npts);
      flag = gsl_matrix_alloc(nDetectors,npts);
      
      VecBool scanFlags (npts);
      double dist;
      //copy the data into them
      for(i=0;i<nDetectors;i++){
	//Get flags
	for(j=si;j<ei;j++){
	  dist = sqrt (pow(dataArray->detectors[di[i]].azElRaPhys[j],2.0)+pow(dataArray->detectors[di[i]].azElDecPhys[j],2.0));
	  if (dist < 0.0/60.0)
	    scanFlags[j-si] = false;
	  else
	    scanFlags[j-si] = true;
	}
	//mn = median(&dataArray->detectors[di[i]].hValues[si], ei-si);
	//mk = median(&dataArray->detectors[di[i]].hKernel[si], ei-si);
	mn = median(&dataArray->detectors[di[i]].hValues[si], &scanFlags[0] , ei-si);
	mk = median(&dataArray->detectors[di[i]].hKernel[si], &scanFlags[0] , ei-si);
	for(j=si;j<ei;j++){
	  gsl_matrix_set(det,i,j-si,dataArray->detectors[di[i]].hValues[j]-mn);
	  gsl_matrix_set(ker,i,j-si,dataArray->detectors[di[i]].hKernel[j]-mk);
	  gsl_matrix_set(flag,i,j-si,dataArray->detectors[di[i]].hSampleFlags[j]);
	}
      }
      
      
      //calculate the denom of the pca_corr
      denom = gsl_matrix_alloc(nDetectors,nDetectors);
      gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.,flag,flag,0.,denom);
      gsl_matrix_add_constant(denom,-1.);
      
      //calculate the the pcaCorr matrix
      pcaCorr = gsl_matrix_alloc(nDetectors,nDetectors);
      gsl_matrix_mul_elements (det,flag);
      gsl_matrix_mul_elements (ker,flag);
      
      gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.,det,det,0.,pcaCorr);
      gsl_matrix_div_elements(pcaCorr,denom);
      gsl_matrix_free(denom);
      gsl_matrix_free(flag);
      
      if (dataArray->detectors[di[0]].atmTemplate.size()>0){
	gsl_matrix * corrMatrix = gsl_matrix_alloc (nDetectors,nDetectors);
	double corTmp=0;
	for (size_t icor = 0; icor < (size_t)nDetectors; icor++){
	  gsl_matrix_set(corrMatrix,icor,icor,1.0);
	  for (size_t jcor = icor+1; jcor <(size_t)nDetectors; jcor++){
	    corTmp =gsl_stats_correlation(&(dataArray->detectors[di[icor]].atmTemplate[si]),1,&(dataArray->detectors[di[jcor]].atmTemplate[si]),1, npts);
	    if (abs(corTmp) > 0.7)
	      corTmp = 0.0;
	    else
						  corTmp = 1.0;
	    gsl_matrix_set(corrMatrix,icor,jcor, corTmp);
	    gsl_matrix_set (corrMatrix,jcor, icor, corTmp);
	  }
	}
	gsl_matrix_mul_elements(pcaCorr,corrMatrix);
	gsl_matrix_free(corrMatrix);
      }
      
      //calculate the eigenvalues and eigenvectors of pcaCorr
      gsl_eigen_symmv_workspace* w=gsl_eigen_symmv_alloc(nDetectors);
      gsl_vector* eVals = gsl_vector_alloc(nDetectors);
      gsl_matrix* eVecs = gsl_matrix_alloc(nDetectors,nDetectors);
      gsl_eigen_symmv(pcaCorr,eVals,eVecs,w);
      gsl_eigen_symmv_free(w);
      gsl_matrix_free(pcaCorr);
      
      //the transpose of the eigenfunction
      gsl_matrix* eFunc = gsl_matrix_alloc(nDetectors,npts);
      gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,eVecs,det,0.,eFunc);
      gsl_matrix* eFunk = gsl_matrix_alloc(nDetectors,npts);
      gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,eVecs,ker,0.,eFunk);
      
      //zero unwanted eFunc and eigenvectors here      
      if(neigToCut > 0){
	if(cutStd > 0){
	  cerr << "Can't have both neigToCut and cutStd non-zero." << endl;
	  exit(1);
	}
	//	double cDamping=0.0;
	//zero the proper values in the eigenVectors
	for(i=0;i<neigToCut;i++) {
	  for(j=0;j<npts;j++){
	    gsl_matrix_set(eFunc,i,j,0.0);
	    gsl_matrix_set(eFunk,i,j,0.0);
	  }
	}
      } else {
	//in this case we'd better have non-zero cutStd
	if(cutStd < 1.){
	  cerr << "cutStd must be greater than 1 ";
	  cerr << "(larger than 2. is recommended)" << endl;
	  exit(1);
	}

	VecDoub eV(nDetectors);
	VecDoub whichGood(nDetectors,1.);
	for(i=0;i<nDetectors;i++){
	  eV[i] = log10(abs(gsl_vector_get(eVals,i)));
	}
	double std = stddev(eV);
	double mev = mean(eV);
	bool keepGoing=1;
	int nKeepLast = nDetectors;
	
	int iterator=0;
	while(keepGoing){
	  //count up number of eigenvalues that pass the cut
	  int count=0;
	  for(i=0;i<nDetectors;i++){
	    if(whichGood[i] > 0){
	      if(abs(eV[i]-mev) > abs(cutStd*std)){
		whichGood[i]=0;
	      } else count++;
	    }
	  }
	  if(count >= nKeepLast){
	    //cerr << "k=" << k << ": Cutting " << nDetectors-count;
	    //cerr << " eigenvectors." << endl;
	    keepGoing=0;
	  } else {
	    //get new mean and stddev for only the good eigenvectors
	    mev=0.;
	    for(i=0;i<nDetectors;i++)
	      if(whichGood[i]>0) mev+=eV[i];
	    mev /= count;
	    std=0.;
	    for(i=0;i<nDetectors;i++)
	      if(whichGood[i]>0) std += (eV[i]-mev)*(eV[i]-mev);
	    std = std/(count-1.);
	    std = sqrt(std);
	    nKeepLast = count;
	  }
	  iterator++;
	}
	
	double cut=mev+cutStd*std;
	cut = pow(10.,cut);
	int cutIndex=0;
	for(i=0;i<nDetectors;i++)
	  if(gsl_vector_get(eVals,i) <= cut){
	    cutIndex=i;
	    break;
	  }
	
	//zero the proper values in the eigenVectors
	for(i=0;i<cutIndex;i++)
	  for(j=0;j<npts;j++){
	    gsl_matrix_set(eFunc,i,j,0.0);
	    gsl_matrix_set(eFunk,i,j,0.0);
	  }
      }
      
      //the output data set is the trans(eigenvects) x eigenfunc
      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,eVecs,eFunc,0.,det);
      gsl_matrix_free(eFunc);
      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,eVecs,eFunk,0.,ker);
      gsl_matrix_free(eFunk);
      
      //save the eigenvectors before deleting them
      
      for(i=0;i<nDetectors;i++)
	for(j=0;j<nDetectors;j++)
	  dataArray->eVectors[i][j] = gsl_matrix_get(eVecs,i,j);
      gsl_matrix_free(eVecs);
      
      //replace the data in the detectors and the kernels
      
      //double debug [nDetectors*npts];
      for(i=0;i<nDetectors;i++){
	for(j=si;j<ei;j++){
	  dataArray->detectors[di[i]].hValues[j] = gsl_matrix_get(det,i,j-si);
	  dataArray->detectors[di[i]].hKernel[j] = gsl_matrix_get(ker,i,j-si);
	  //debug[i*npts+(j-si)] = dataArray->detectors[di[i]].hValues[j];
	}
	
	
      }
      //      if (k==1)
      //    	  writeVecOut("pca_v1.txt",debug, nDetectors*npts);
      
      gsl_vector_free(eVals);
      gsl_matrix_free(det);
      gsl_matrix_free(ker);
    }//iterating on scans, k
  }
  return 1;
  
}

void CleanPCA::setDamping(double damping){
	this->damping = damping;
}

double CleanPCA::getDamping(){
	return this->damping;
}
