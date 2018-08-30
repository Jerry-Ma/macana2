#ifndef _CLEAN_BSPLINE_H
#define _CLEAN_BSPLINE_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit_nlin.h>
#include <suitesparse/cs.h>




#include "Clean.h"


typedef enum  {STRIPE_NONE,STRIPE_PCA, STRIPE_FFT, AZEL_TEMPLATE} CleanBsplineDestriping;

///Perform Cottingham method atmospherical estimation for bolometer timestreams.
///Propouse: Estimate atmospherical emission from bolometer timestreams and subtract it from the original data
///This class uses a B-Spline Basis to fit a B-Spline curve to represent the atmosphere signal.
///All matrix are represented using the sparse matrix library CXSparse. This increases the efficiency of the code
class CleanBspline : public Clean{
  private:
    cs *baseMatrix;							///<Sparse B-Spline Base Matrix
    cs *baseMatrix_t;						///<Sparse B-Spline Base Matrix transpose
    gsl_bspline_workspace *bsw;				///<GSL B-Spline workspace
    gsl_vector *time;						///<Time vector
    bool calibrated;						///<Indicates is time stream is calibrated or not
    int baseMatrixSamples;					///<Number of samples to clean
    int baseMatrixDetectors;				///<Number of detectors to clean

    int **detPerHextant;
    VecInt nDetPerHextant;
    int maxHex;
    double	 bright;
    int tid;
    double resample;
    VecDoub fullMedianValues;
	bool cleanKernel;
    double *fixedData;
    bool *scanStatus;
    size_t totSamples;
    VecDoub atmTemplate;
    size_t currScan;




  protected:  
    ///Use GSL to create the B-Spline Base Matrix, store it in a sparse matrix
    void createBaseMatrix(int nSamples, int nDetectors);
    void createBaseMatrix(int nSamples, int nDetectors, int nCells, VecDoub azOff, VecDoub elOff);
    void createPointingMatrix();
    void destroyBaseMatrix();
    cs* getBaseMatrix();
    //Wrapper function used to call the Array->calibrate method if bolometer data is not calibrated
    void calibrate();
    void removeBadBolos();
    cs *getPMatrix(gsl_vector *ra, gsl_vector *dec, double pixelSize);
    void downSample(gsl_vector **out, double *data, bool *flags, long nSamples, long dowSample);
    bool fixFlags(double *dataVector, bool *flagsVector, int nDetectors, long nSamples);
    void subtractTemplate(double *detector,size_t oSamples, double *aTemplate, size_t nSamples, double increment, char *outName=NULL, bool overwrite = false);
    double *cottingham(double *dataVector, gsl_vector *raVector, gsl_vector *decVector, bool *flags, size_t nDetectors, size_t nSamples, double cleanPixelSize, int refBolo=-1);
    bool cleanScans ();
    void removeLargeScaleResiduals (CleanBsplineDestriping method,double correlate);
    MatDoub linearGainCorrection(double *dataVectorm, size_t nDetectors, size_t nSamples, double *fakeAtm = NULL);
    void cleanStripeFFT(double *data, double *kernel, size_t nDetectors, size_t nSamples, double cutOff, double sampleRate);
    void azelResidualLMT (double *data, double *flags, double *flagdis, double *kvector, size_t nDetectors, size_t nSamples, gsl_vector *azOffsets, gsl_vector *elOffsets, double *sens);
    CleanBsplineDestriping translateStripeMethod();
    bool pca (double *dataVector, double *flags, size_t nSamples, size_t nDetectors, size_t neig2cut);
  public:
    CleanBspline(Array*  dataArray, Telescope *telescope);
    ~CleanBspline();
    bool clean ();

};

#endif
