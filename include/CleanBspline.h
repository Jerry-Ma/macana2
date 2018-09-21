#ifndef _CLEAN_BSPLINE_H
#define _CLEAN_BSPLINE_H
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit_nlin.h>
#include <cs.h>




#include "Clean.h"



typedef enum  {STRIPE_NONE,STRIPE_PCA, STRIPE_FFT, SUBARRAY_SPLINE, STRIPE_SCAN, AZEL_TEMPLATE} CleanBsplineDestriping;

///Perform Cottingham method atmospherical estimation for bolometer timestreams.
///Propouse: Estimate atmospherical emission from bolometer timestreams and subtract it from the original data
///This class uses a B-Spline Basis to fit a B-Spline curve to represent the atmosphere signal.
///All matrix are represented using the sparse matrix library CXSparse. This increases the efficency of the code
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

    double *fixedData;
    bool *scanStatus;
    size_t totSamples;
    VecDoub atmTemplate;
    size_t currScan;



  protected:  
    ///Use GSL to create the B-Spline Base Matrix, store it in a sparse matrix
    void createBaseMatrix(int nSamples, int nDetectors);
    void createPointingMatrix();
    void destroyBaseMatrix();
    void calibrate();
    void removeBadBolos();
    cs *getPMatrix(gsl_vector *ra, gsl_vector *dec, double pixelSize);
    void downSample(gsl_vector **out, double *data, bool *flags, long nSamples, long dowSample);
    bool fixFlags(double *dataVector, bool *flagsVector, int nDetectors, long nSamples);
    void subtractTemplate(double *detector,size_t oSamples, double *aTemplate, size_t nSamples, double increment, char *outName=NULL, bool overwrite = false);
    double *cottingham(double *dataVector, gsl_vector *raVector, gsl_vector *decVector, bool *flags, size_t nDetectors, size_t nSamples, double cleanPixelSize, int refBolo=-1);
    bool cleanScans ();
    double *removeAzElResidual(double *dataVector, gsl_vector *raVector, gsl_vector *decVector, size_t nDetectors, size_t nSamples);
    void removeLargeScaleResiduals (CleanBsplineDestriping method,double correlate);
    MatDoub linearGainCorrection(double *dataVectorm, size_t nDetectors, size_t nSamples, double *fakeAtm = NULL);
    bool removeCorrelations (double *dataVector, double *flagVector, size_t nDetectors, size_t nSamples, double corrFactor, double *azoffset, double *eloffset, gsl_vector *azVector, gsl_vector *elVector);
    bool removeCorrelations2 (double *dataVector, size_t nDetectors, size_t nSamples, double corrFactor, double *azoffset, double *eloffset, double *kVector, double *flagVector);
    void cleanStripeFFT2(double *data, double *kernel, int nDetectors, long nSamples, double cutOff, double cutLow = 0, double sampleRate=64.0);
    void cleanStripeFFT3(double *data, double *kernel, size_t nDetectors, size_t nSamples, double cutOff, double sampleRate);
    bool stripePCA(double *dataVector, double *kernelVector, double *flags, size_t nDetectors, size_t nSamples, size_t neigToCut, bool adaptive = false);
    bool removeScanPattern (double *dataVector, gsl_vector *azVector, gsl_vector *decVector, double *flags, size_t nDetectors, size_t nSamples);
    void azelResidual (double *data, double *flags, size_t nDetectors, size_t nSamples, double *azOffsets, double *elOffsets, double correlate);
    void azelResidual (double *data, double *flags, size_t nDetectors, size_t nSamples, double *azOffsets, double *elOffsets, double *sens);
    void azelResidualLMT (double *data, double *flags, size_t nDetectors, size_t nSamples, gsl_vector *azOffsets, gsl_vector *elOffsets, double *sens);
    void splineResidual (double *data, double *flags, size_t nDetectors, size_t nSamples, gsl_vector  *azVector, gsl_vector  *elVector, double correlate);
    CleanBsplineDestriping translateStripeMethod();
  public:
    CleanBspline(Array*  dataArray, Telescope *telescope);
    ~CleanBspline();
    bool clean ();
    cs* getBaseMatrix();
};

#endif
