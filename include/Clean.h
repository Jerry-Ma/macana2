#ifndef _CLEAN_H_
#define _CLEAN_H_

#include <netcdfcpp.h>
#include "Array.h"
#include "Telescope.h"
#include "AnalParams.h"

class Clean{
  protected:
    Array *dataArray;
    AnalParams *ap;
    Telescope *telescope;
  public:
    Clean(Array *dataArray,Telescope *telescope);
    virtual bool clean()=0;
    virtual ~Clean(){}
};


#endif