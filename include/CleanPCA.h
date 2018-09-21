#ifndef _CLEAN_PCA_H
#define _CLEAN_PCA_H



#include "Clean.h"
#include "Array.h"
#include "Telescope.h"
#include "AnalParams.h"

class CleanPCA : public Clean {
  private:
    bool adaptive;
    double damping;
  public:
    CleanPCA(Array *dataArray,Telescope *telescope);
    bool clean();
    void setAdaptive(bool adaptive);
    void setDamping(double damping);
    bool getAdaptive();
    double getDamping();
};


#endif
