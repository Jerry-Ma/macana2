#include "CleanSelector.h"

Clean * CleanSelector::getCleaner(Array *dataArray, Telescope *telescope){
  AnalParams *ap = dataArray->getAp();
  

  if (ap->getCutStd() <= 0 && ap->getNeigToCut() <= 0 && ap->getOrder() <= 0 && ap->getTOrder() <=0){
	printErrorMessage("Cannot determinate desired cleaning method!!!\n\tSet cutStd >0 or neigToCut > 0 for PCA Cleaning.\n\tSet splineOrder >=3 for Cottingham method.");
	exit(-1);
  }

  if (ap->getOrder() > 0 && ap->cutStd>0){
	  printErrorMessage("Options cutStd > 0 and order>0 are  not compatible.");
	  cerr<<ap->cutStd<<endl;
	  exit(-1);
  }

  if (ap->getOrder() > 0){
    if (ap->getOrder()<3){
      printErrorMessage("Incorrect value for spline order. Must be >=3.");
      exit(-1);
    } else{
      cout<<"CleanSelector(): Selected Cottingham Method for Cleaning" <<endl;
      CleanBspline *bspline = new CleanBspline(dataArray,telescope);
      return bspline;
    }
  }
  else if (ap->getTOrder()>0){
	if (ap->getTOrder() >3 ){
		printErrorMessage("Incorrect value for high order atmosphere template cleaning. Must be <=3.");
		exit(-1);
	}else{
		cout<<"CleanSelector(): Selected High Order Atmosphere Template Subtraction for Cleaning" <<endl;
		return new CleanHigh(dataArray,telescope);
	}
  }else{
	  if (ap->getCutStd() > 0 && ap->getNeigToCut() >0){
		  printErrorMessage("Cannot set cutStd > 0 and neigToCut > 0.");
		  exit(-1);
	  }
	  cerr<<"CleanSelector():Selected PCA Cleaning"<<endl;
	  return new CleanPCA(dataArray, telescope);
  }
  return NULL;
}

void CleanSelector::printErrorMessage(const char* message){
	  cerr<<"CleanSelector():Fatal Error:"<<endl;
	  cerr<<"\t-------------------------------------------------------"<<endl;
	  cerr<<"\t"<<message<<endl;
	  cerr<<"\tPlease check analysis parameters xml file"<<endl;
	  cerr<<"\tProgram Aborted"<<endl;
	  cerr<<"\t-------------------------------------------------------"<<endl;
}

Clean2dStripe * CleanSelector::getMapCleaner(Map *cmap, AnalParams *ap){
	if (ap->getOrder() > 0){
		return new Clean2dStripe(cmap);
	}
	return NULL;
}

