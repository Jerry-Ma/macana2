#include "port_bandpass.h"

namespace bandpass{

void MacanaBandpass(VecDoub &hValues, const VecDoub &digFiltTerms, VecDoub &storage, const int &th,
                       const size_t &nSamples, const size_t &nTerms, const size_t &nCoef, int &tt,
                       double &tmp){

    //Convolve it!
    for(int t=nTerms;t<th;t++){
        tmp=0.;
        tt = t-nTerms;
        for(int i = 0;i<nCoef;i++){
            tmp += hValues[tt+i]*digFiltTerms[i];
        }
        storage[t] = tmp;
    }

    //Replace old values
    for(int i=nTerms;i<th;i++) {
      hValues[i] = storage[i];
    }

}
} //namespaace
