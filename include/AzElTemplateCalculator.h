#ifndef _AZELTEMPLATECALCULATOR_H
#define _AZELTEMPLATECALCULATOR_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "nr3.h"



typedef enum  {CONST = 1,LINEAR=3, QUADRATIC = 6, CUBIC = 10} AzElFitType;

class AzElTemplateCalculator{
	private:
		gsl_matrix *data;
		gsl_matrix *flags;
		gsl_vector *azOffsets;
		gsl_vector *elOffsets;
		gsl_matrix *azOffMatrix;
		gsl_matrix *elOffMatrix;
		//gsl_vector *corrCoeffs;
		gsl_matrix *azelTemplate;
		size_t nDetectors;
		size_t nSamples;
		AzElFitType mode;
		bool matrixCoords;
		MatDoub tCoeffs;



	public:
		AzElTemplateCalculator (gsl_matrix *data, gsl_matrix *flags, gsl_vector *azOffsets,gsl_vector *elOffsets);
		AzElTemplateCalculator(gsl_matrix* data,gsl_matrix* flags, gsl_matrix* azOffsets, gsl_matrix* elOffsets);
		AzElFitType calculateMode(size_t numDet);
		void calculateTemplate (gsl_vector *corrCoeffs);
		void calculateTemplate2(gsl_vector *corrCoeffs);
		void decorrelateCoeffs(size_t neig2cut);
		gsl_matrix *getTemplate ();
		void overrideMode(int mode);
		void removeTemplate ();
		void setMode();
		void updateTemplate();
		~AzElTemplateCalculator ();


};

#endif
