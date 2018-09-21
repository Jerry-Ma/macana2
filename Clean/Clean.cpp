using namespace std;

#include "Clean.h"



Clean::Clean(Array *dataArray, Telescope *telescope){
  this->ap = NULL; 
  this->dataArray=dataArray;
  this->telescope= telescope;
  if (this->dataArray != NULL){
    this->ap=dataArray->getAp();
  }
}

