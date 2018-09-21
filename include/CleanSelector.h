#ifndef _CLEAN_SELECTOR_H_
#define _CLEAN_SELECTOR_H_

#include <netcdfcpp.h>

#include "Array.h"
#include "Telescope.h"
#include "Clean.h"
#include "CleanPCA.h"
#include "CleanBspline.h"
#include "CleanHigh.h"
//#include "CleanAPCA.h"

#include "Clean2dStripe.h"
#include "Map.h"

class CleanSelector{
  public:
     static Clean  *getCleaner(Array *dataArray, Telescope *telescope);
     static Clean2dStripe *getMapCleaner(Map *cmap, AnalParams *ap);
  private:
     static void printErrorMessage(const char *message);
};


#endif
