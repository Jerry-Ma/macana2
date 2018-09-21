#ifndef _AZELTEMPLATECALCULATOR_H
#define _AZELTEMPLATECALCULATOR_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>



typedef enum  {AVERAGE=1, LINEAR=3, QUADRATIC = 6} AzElFitType;

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

		void setMode();
	public:
		AzElTemplateCalculator (gsl_matrix *data, gsl_matrix *flags, gsl_vector *azOffsets,gsl_vector *elOffsets);
		AzElTemplateCalculator(gsl_matrix* data,gsl_matrix* flags, gsl_matrix* azOffsets, gsl_matrix* elOffsets);
		void calculateTemplate (gsl_vector *corrCoeffs);
		gsl_matrix *getTemplate ();
		void removeTemplate ();
		void overrideMode(AzElFitType mode);
		~AzElTemplateCalculator ();


};

#endif
