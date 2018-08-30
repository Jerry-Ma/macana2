#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>


#include <iostream>
#include "AzElTemplateCalculator.h"
#include "vector_utilities.h"
#include "nr3.h"

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
	tCoeffs.resize(0,0);
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

	size_t npars = (size_t) mode-1;
	size_t ipar = 0;

	gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(nDetectors, npars);
	gsl_matrix *coordMatrix = gsl_matrix_alloc(nDetectors, npars);
	gsl_matrix *cov =gsl_matrix_alloc(npars,npars);
	gsl_vector *tmp = gsl_vector_alloc (nDetectors);
	gsl_vector *iData = gsl_vector_alloc(nDetectors);
	gsl_vector *coeff = gsl_vector_alloc(npars);
	gsl_vector *fcoeff = gsl_vector_alloc(nDetectors);

	double downWeight = 1.0e-2 * gsl_vector_min (corrCoeffs);

	tCoeffs.resize(npars,nSamples);

	double chisq= 0.0;
	double jaz, jel;
	for (size_t i=0; i<nSamples; i++){
		gsl_vector_set_all(tmp,0.0);
		gsl_vector_memcpy(fcoeff, corrCoeffs);
		for (size_t j=0; j<nDetectors; j++){
			ipar = 0;
			gsl_vector_set(iData, j, gsl_matrix_get(data,j,i));
			if (gsl_matrix_get (flags,j,i)==0)
				gsl_vector_set(fcoeff, j, downWeight);
			if (!matrixCoords){
				jaz = gsl_vector_get(azOffsets,j);
				jel = gsl_vector_get(elOffsets,j);
			}else{
				jaz = gsl_matrix_get(azOffMatrix,j,i);
				jel = gsl_matrix_get (elOffMatrix,j,i);
			}

			//gsl_matrix_set(coordMatrix, j,ipar++, sqrt(pow(jaz,2)+pow(jel,2)));
			//
			//gsl_matrix_set(coordMatrix, j,ipar++, 1.0);
			if (mode == LINEAR || mode == QUADRATIC || mode == CUBIC){
				gsl_matrix_set(coordMatrix,j, ipar++, jaz);
				gsl_matrix_set(coordMatrix,j, ipar++, jel);
			}
			if (mode == QUADRATIC || mode == CUBIC){
				gsl_matrix_set(coordMatrix,j, ipar++, jaz*jaz);
				gsl_matrix_set(coordMatrix,j, ipar++, jaz*jel);
				gsl_matrix_set(coordMatrix,j, ipar++, jel*jel);
			}
			if (mode == CUBIC){
				gsl_matrix_set(coordMatrix,j, ipar++, jaz*jaz*jaz);
				gsl_matrix_set(coordMatrix,j, ipar++, jaz*jaz*jel);
				gsl_matrix_set(coordMatrix,j, ipar++, jaz*jel*jel);
				gsl_matrix_set(coordMatrix,j, ipar++, jel*jel*jel);
			}

		}
		gsl_multifit_wlinear(coordMatrix,fcoeff,iData,coeff, cov, &chisq, work);
		gsl_blas_dgemv(CblasNoTrans,1.0,coordMatrix,coeff,0.0,tmp);
		gsl_vector_scale(tmp,0.9);
		gsl_matrix_set_col(azelTemplate, i,tmp);
		for (size_t ipar = 0; ipar <npars; ipar++)
			tCoeffs[ipar][i] = gsl_vector_get(coeff,ipar);

	}

//	writeMatOut("azel_pars.txt", tCoeffs);
//	exit(-1);

	gsl_matrix_free(coordMatrix);
	gsl_matrix_free(cov);
	gsl_vector_free(tmp);
	gsl_vector_free(coeff);
	gsl_vector_free (iData);
	gsl_multifit_linear_free(work);

}

void AzElTemplateCalculator::calculateTemplate2(gsl_vector *corrCoeffs) {
	gsl_matrix_set_all(azelTemplate,0.0);

	size_t npars = (size_t) mode;
	size_t ipar = 0;

//	if (mode == LINEAR)
//		npars++;

	gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(nDetectors, npars);
	gsl_matrix *coordMatrix = gsl_matrix_alloc(nDetectors, npars);
	gsl_matrix *cov =gsl_matrix_alloc(npars,npars);
	gsl_vector *tmp = gsl_vector_alloc (nDetectors);
	gsl_vector *iData = gsl_vector_alloc(nDetectors);
	gsl_vector *coeff = gsl_vector_alloc(npars);
	gsl_vector *fcoeff = gsl_vector_alloc(nDetectors);

	double downWeight = 1.0e-6 * gsl_vector_min (corrCoeffs);

	tCoeffs.resize(npars,nSamples);

	double chisq= 0.0;
	double jaz, jel;
	for (size_t i=0; i<nSamples; i++){
		gsl_vector_set_all(tmp,0.0);
		gsl_vector_memcpy(fcoeff, corrCoeffs);

		for (size_t j=0; j<nDetectors; j++){
			ipar = 0;
			size_t dSamples = 0;
			gsl_vector_set(iData, j, gsl_matrix_get(data,j,i));
			if (gsl_matrix_get (flags,j,i)==0){
				gsl_vector_set(fcoeff, j, downWeight);
				dSamples++;
			}
			if (!matrixCoords){
				jaz = gsl_vector_get(azOffsets,j);
				jel = gsl_vector_get(elOffsets,j);
			}else{
				jaz = gsl_matrix_get(azOffMatrix,j,i);
				jel = gsl_matrix_get (elOffMatrix,j,i);
			}

			if (dSamples >= size_t(3.0*nDetectors/4.0)){
				cout<<"Warning to many detectors masked. Reseting weights"<<endl;
				gsl_vector_memcpy(fcoeff, corrCoeffs);
			}

			//gsl_matrix_set(coordMatrix, j,ipar++, sqrt(pow(jaz,2)+pow(jel,2)));
			//
			gsl_matrix_set(coordMatrix, j,ipar++, 1.0);
			if (mode == LINEAR || mode == QUADRATIC || mode == CUBIC){
				gsl_matrix_set(coordMatrix,j, ipar++, jaz);
				gsl_matrix_set(coordMatrix,j, ipar++, jel);
				//
			}
			if (mode == QUADRATIC || mode == CUBIC){
				gsl_matrix_set(coordMatrix,j, ipar++, jaz*jaz);
				gsl_matrix_set(coordMatrix,j, ipar++, jaz*jel);
				gsl_matrix_set(coordMatrix,j, ipar++, jel*jel);
			}
			if (mode == CUBIC){
				gsl_matrix_set(coordMatrix,j, ipar++, jaz*jaz*jaz);
				gsl_matrix_set(coordMatrix,j, ipar++, jaz*jaz*jel);
				gsl_matrix_set(coordMatrix,j, ipar++, jaz*jel*jel);
				gsl_matrix_set(coordMatrix,j, ipar++, jel*jel*jel);
			}

		}
		gsl_multifit_wlinear(coordMatrix,fcoeff,iData,coeff, cov, &chisq, work);
		gsl_blas_dgemv(CblasNoTrans,1.0,coordMatrix,coeff,0.0,tmp);
		gsl_matrix_set_col(azelTemplate, i,tmp);
		for (size_t ipar = 0; ipar <npars; ipar++)
			tCoeffs[ipar][i] = gsl_vector_get(coeff,ipar);

	}

//	writeMatOut("azel_pars.txt", tCoeffs);
//	exit(-1);

	gsl_matrix_free(coordMatrix);
	gsl_matrix_free(cov);
	gsl_vector_free(tmp);
	gsl_vector_free(coeff);
	gsl_vector_free (iData);
	gsl_multifit_linear_free(work);

}

void AzElTemplateCalculator::decorrelateCoeffs(size_t neig2cut){


	for (size_t isample=0; isample<nSamples; isample++)
		tCoeffs[0][isample]*=0.9;



//	size_t nPars = (size_t) mode;
//	gsl_matrix *det = gsl_matrix_alloc (nPars, nSamples);
//	gsl_matrix *pcaCorr = gsl_matrix_alloc(nPars,nPars);
//	gsl_eigen_symmv_workspace* w=gsl_eigen_symmv_alloc(nPars);
//	gsl_vector* eVals = gsl_vector_alloc(nPars);
//	gsl_matrix* eVecs = gsl_matrix_alloc(nPars,nPars);
//	gsl_matrix* eFunc = gsl_matrix_alloc(nPars,nSamples);
//	gsl_matrix *fl = gsl_matrix_alloc (nPars, nSamples);
//	gsl_matrix *denom = gsl_matrix_alloc(nPars,nPars);
//
//	//Copy data
//
//	for (size_t ibolo = 0; ibolo< nPars; ibolo++)
//		for (size_t isample=0; isample<nSamples; isample++){
//			gsl_matrix_set (det, ibolo, isample, tCoeffs[ibolo][isample]);
//			gsl_matrix_set (fl, ibolo, isample, 1.0);
//		}
//
//	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.,fl,fl,0.,denom);
//	gsl_matrix_add_constant(denom,-1.);
//	gsl_matrix_mul_elements (det,fl);
//
//	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.,det,det,0.,pcaCorr);
//	gsl_matrix_div_elements(pcaCorr,denom);
//	gsl_matrix_free(fl);
//	gsl_matrix_free(denom);
//	gsl_eigen_symmv(pcaCorr,eVals,eVecs,w);
//	gsl_eigen_symmv_sort(eVals,eVecs, GSL_EIGEN_SORT_ABS_DESC);
//	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,eVecs,det,0.,eFunc);
//
//	for (size_t ibolo = neig2cut; ibolo< nPars; ibolo++)
//		for (size_t isample=0; isample<nSamples; isample++){
//			gsl_matrix_set (eFunc, ibolo, isample, 0.0);
//		}
//
//	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,eVecs,eFunc,0.,det);
//
//	for (size_t ibolo = 0; ibolo< nPars; ibolo++)
//		for (size_t isample=0; isample<nSamples; isample++)
//			tCoeffs[ibolo][isample]= gsl_matrix_get(det,ibolo,isample);
//
//	gsl_matrix_free (det);
//	gsl_eigen_symmv_free(w);
//	gsl_matrix_free(pcaCorr);
//	gsl_vector_free(eVals);
//	gsl_matrix_free (eVecs);
//	gsl_matrix_free(eFunc);

}


void AzElTemplateCalculator::setMode (){
	this->mode = calculateMode(this->nDetectors);

}

AzElFitType AzElTemplateCalculator::calculateMode(size_t numDet){
		if (numDet > CUBIC)
			return CUBIC;
		else if (numDet > QUADRATIC)
			return QUADRATIC;
		else if (numDet >LINEAR)
			return LINEAR;
		else if (numDet > CONST)
			return CONST;
		else{
			cerr<<"Wrong number of detectors. Cannot produce an Az/El template. Imploding"<<endl;
			exit(-1);
		}
		return LINEAR;
}

void AzElTemplateCalculator::overrideMode(int mode){
	switch (mode){
		case 1:
			this->mode = LINEAR;
			//cout<<"AzELTemplateCaltulator(). Setting model to LINEAR mode."<<endl;
			break;
		case 2:
			this->mode = QUADRATIC;
			//cout<<"AzELTemplateCaltulator(). Setting model to QUADRATIC mode."<<endl;
			break;
		case 3:
			this->mode = CUBIC;
			//cout<<"AzELTemplateCaltulator(). Setting model to CUBIC mode."<<endl;
			break;
		default:
			this->setMode();
			break;
			//cout<<"AzELTemplateCaltulator(). Setting model to AUTO mode."<<endl;
	}
}



void AzElTemplateCalculator::removeTemplate() {
	if (azelTemplate)
		gsl_matrix_sub(data,azelTemplate);
}

gsl_matrix* AzElTemplateCalculator::getTemplate() {
	return azelTemplate;
}

void AzElTemplateCalculator::updateTemplate(){
	double jaz,jel;
	double newtemp;
	size_t ipar;

	//writeMatOut("preDecorrMat.txt", azelTemplate);

	for (size_t i=0; i<nDetectors; i++)
		for (size_t j=0; j<nSamples; j++){
			if (!matrixCoords){
				jaz = gsl_vector_get(azOffsets,i);
				jel = gsl_vector_get(elOffsets,i);
			}else{
				jaz = gsl_matrix_get(azOffMatrix,i,j);
				jel = gsl_matrix_get (elOffMatrix,i,j);
			}
			ipar =0;
			newtemp = tCoeffs[ipar++][j];
			if (mode == LINEAR || mode == QUADRATIC || mode == CUBIC){
				newtemp += tCoeffs[ipar++][j]*jaz;
				newtemp += tCoeffs[ipar++][j]*jel;
			}
			if (mode == QUADRATIC || mode == CUBIC){
				newtemp += tCoeffs[ipar++][j]*jaz*jaz;
				newtemp += tCoeffs[ipar++][j]*jaz*jel;
				newtemp += tCoeffs[ipar++][j]*jel*jel;
			}
			if (mode == CUBIC){
				newtemp += tCoeffs[ipar++][j]*jaz*jaz*jaz;
				newtemp += tCoeffs[ipar++][j]*jaz*jaz*jel;
				newtemp += tCoeffs[ipar++][j]*jaz*jel*jel;
				newtemp += tCoeffs[ipar++][j]*jel*jel*jel;
			}
			gsl_matrix_set(azelTemplate,i,j, newtemp);
		}
	//writeMatOut("postDecorrMat.txt", azelTemplate);
	//exit(-1);
}


AzElTemplateCalculator::~AzElTemplateCalculator() {
	if (azelTemplate)
		gsl_matrix_free(azelTemplate);
}
