/**
 * @file Clean.h
 * @author David Sánchez-Argüelles
 * @brief Head of the hierarchy for the atmosphere removal (a.k.a. Atmosphere Cleaning)
 *
 * This is a virtual class that represent the general structure for the  hierarchy
 * of atmosphere removal techniques. To create a new atmosphere subtraction method you must
 * inherit from this class and modify the CleanSelector factory accordingly
 */

#ifndef _CLEAN_H_
#define _CLEAN_H_

#include <netcdfcpp.h>
#include "Array.h"
#include "Telescope.h"
#include "AnalParams.h"

/**
 *  @brief Head of the atmosphere removal class hierarchy
 */

class Clean{
  protected:
    Array *dataArray;		///Pointer to the bolometer data
    AnalParams *ap;			///Pointer to AP configuration
    Telescope *telescope;   ///Pointer to the telescope data
  public:
    bool nodc;				///If true no median/mean subtraction will be performed to the raw timestream data
    Clean(Array *dataArray,Telescope *telescope);
    virtual bool clean()=0; ///This is the method called by the main thread to carry out the atmosphere subtraction
    virtual ~Clean(){}
};


#endif
