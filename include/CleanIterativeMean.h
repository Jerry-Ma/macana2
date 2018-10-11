#ifndef _CLEAN_IT_H
#define _CLEAN_IT_H

#include "Clean.h"
#include "nr3.h"

#define SIGN(X) (X<=0.0?-1.0:1.0);

class CleanIterativeMean: public Clean{
private:
	size_t maxIteration;
	size_t maxDelta;

	double linearMeanFit();
	VecDoub getSigns(double *maxCor);
	void signMultiply(VecDoub signs);

public:
	CleanIterativeMean(Array* dataArray, Telescope *telescope);
    ~CleanIterativeMean();
	bool clean();																	////Start cleaning process											////Set the number of component
};


#endif
