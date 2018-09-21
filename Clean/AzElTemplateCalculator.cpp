#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <iostream>
#include "AzElTemplateCalculator.h"

using namespace std;

AzElTemplateCalculator::AzElTemplateCalculator(gsl_matrix* data,
		gsl_matrix* flags, gsl_vector* azOffsets, gsl_vector* elOffsets) {

	this->data  = data;
	this->flags = flags;
	this->azOffsets = azOffsets;
	this->elOffsets = elOffsets;
	//this->corrCoeffs = corrCoeffs;
	nDetectors = data->size1;
	nSamples = data->size2;

	azOffMatrix = NULL;
	elOffMatrix = NULL;
	matrixCoords = false;

	setMode();
	azelTemplate = gsl_matrix_alloc(nDetectors, nSamples);
	gsl_matrix_set_all(azelTemplate,0.0);
}

AzElTemplateCalculator::AzElTemplateCalculator(gsl_matrix* data,
		gsl_matrix* flags, gsl_matrix* azOffsets, gsl_matrix* elOffsets) {

	this->data  = data;
	this->flags = flags;
	this->azOffsets = NULL;
	this->elOffsets = NULL;

	azOffMatrix = azOffsets;
	elOffMatrix = elOffsets;
	//this->corrCoeffs = corrCoeffs;
	nDetectors = data->size1;
	nSamples = data->size2;

	matrixCoords = true;

	setMode();
	azelTemplate = gsl_matrix_alloc(nDetectors, nSamples);
	gsl_matrix_set_all(azelTemplate,0.0);
}

void AzElTemplateCalculator::calculateTemplate(gsl_vector *corrCoeffs) {
	gsl_matrix_set_all(azelTemplate,0.0);

	size_t npars = (size_t) mode;
	size_t ipar = 0;

	gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(nDetectors, npars);
	gsl_matrix *coordMatrix = gsl_matrix_alloc(nDetectors, npars);
	gsl_matrix *cov =gsl_matrix_alloc(npars,npars);
	gsl_vector *tmp = gsl_vector_alloc (nDetectors);;
	gsl_vector *iData = gsl_vector_alloc(nDetectors);
	gsl_vector *coeff = gsl_vector_alloc(npars);
	gsl_vector *fcoeff = gsl_vector_alloc(nDetectors);



	double chisq= 0.0;
	double jaz, jel;
	for (size_t i=0; i<nSamples; i++){
		gsl_vector_set_all(tmp,0.0);
		gsl_vector_memcpy(fcoeff, corrCoeffs);
		for (size_t j=0; j<nDetectors; j++){
			ipar = 0;
			gsl_vector_set(iData, j, gsl_matrix_get(data,j,i));
			if (gsl_matrix_get (flags,j,i)==0)
				gsl_vector_set(fcoeff, j, gsl_vector_get(fcoeff,j)*1e-1);
			if (!matrixCoords){
				jaz = gsl_vector_get(azOffsets,j);
				jel = gsl_vector_get(elOffsets,j);
			}else{
				jaz = gsl_matrix_get(azOffMatrix,j,i);
				jel = gsl_matrix_get (elOffMatrix,j,i);
			}

			gsl_matrix_set(coordMatrix, j,ipar++, 1.0);
			if (mode == LINEAR || mode == QUADRATIC){
				gsl_matrix_set(coordMatrix,j, ipar++, jaz);
				gsl_matrix_set(coordMatrix,j, ipar++, jel);
			}
			if (mode == QUADRATIC){
				gsl_matrix_set(coordMatrix,j, ipar++, jaz*jaz);
				gsl_matrix_set(coordMatrix,j, ipar++, jel*jel);
				gsl_matrix_set(coordMatrix,j, ipar++, jaz*jel);
			}
		}
		gsl_multifit_wlinear(coordMatrix,fcoeff,iData,coeff, cov, &chisq, work);
		gsl_blas_dgemv(CblasNoTrans,1.0,coordMatrix,coeff,0.0,tmp);
		gsl_matrix_set_col(azelTemplate, i,tmp);

	}
	gsl_matrix_free(coordMatrix);
	gsl_matrix_free(cov);
	gsl_vector_free(tmp);
	gsl_vector_free(coeff);
	gsl_vector_free (iData);
	gsl_multifit_linear_free(work);

}

void AzElTemplateCalculator::setMode (){
	if (nDetectors > QUADRATIC)
		mode = QUADRATIC;
	else if (nDetectors >LINEAR)
		mode = LINEAR;
	else if (nDetectors > AVERAGE)
		mode = AVERAGE;
	else{
		cerr<<"Wrong number of detectors. Cannot produce an Az/El template. Imploding"<<endl;
		exit(-1);
	}
}

void AzElTemplateCalculator::overrideMode(AzElFitType mode){
	this->mode=mode;
}

void AzElTemplateCalculator::removeTemplate() {
	if (azelTemplate)
		gsl_matrix_sub(data,azelTemplate);
}

gsl_matrix* AzElTemplateCalculator::getTemplate() {
	return azelTemplate;
}

AzElTemplateCalculator::~AzElTemplateCalculator() {
	if (azelTemplate)
		gsl_matrix_free(azelTemplate);
}
