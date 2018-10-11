#ifndef _CLEAN_SELECTOR_H_
#define _CLEAN_SELECTOR_H_

#include <netcdfcpp.h>
#include "Array.h"
#include "Telescope.h"
#include "Map.h"
#include "Clean.h"

class CleanSelector{
  public:
     static Clean  *getCleaner(Array *dataArray, Telescope *telescope);
  private:
     static void printErrorMessage(const char *message);
};


#endif
