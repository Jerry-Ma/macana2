#include <cmath>

#include "CleanHigh.h"
#include "nr3.h"
#include "vector_utilities.h"
#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>

#include <iostream>
#include <iomanip>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>

CleanHigh::CleanHigh(Array* dataArray, Telescope* tel) :
    Clean(dataArray,tel){
        order = dataArray->getAp()->getOrder();
    }

bool CleanHigh::clean(){

    fullMedianSustraction();
    //fullLinearCorrection();
    
	size_t nDetectors = dataArray->getNDetectors();
	size_t nScans = telescope->scanIndex.ncols();
	size_t si, ei;
	size_t nSamples;
	int *di = dataArray->getDetectorIndices(); 	//This variable hold the "alive" bolometers

	MatDoub dataVector(0,0);

    for (size_t iScan = 0; iScan < nScans; iScan++){  //Loop through all the data scans
	//for (size_t iScan = 0; iScan < 1; iScan++){  //Loop through all the data scans
		si=telescope->scanIndex[0][iScan];
		ei=telescope->scanIndex[1][iScan]+1;
		nSamples = ei -si;

		//Copy data from detectors to working variable dataVector
		dataVector.resize(nDetectors,nSamples);   //  <-------------------- size of matrix
        
        for (size_t i=0; i<nDetectors; i++)
			for (size_t j=0; j<nSamples; j++)
				dataVector[i][j]=dataArray->detectors[di[i]].hValues[si+j];

		//Do something with the data, like taking out the median of each detector

		cerr<<"CleanHigh(): Processing Data on Scan "<< iScan<<endl;
//		double meanValue;
//		for (size_t i=0; i<nDetectors; i++){
//			meanValue = mean(&(dataVector[i][0]), nSamples);
//			for (size_t j=0; j<nSamples; j++){
//                dataVector[i][j]-= meanValue;
//            }
//		}

        // MODIFICACION  (
        // Offsets de AzTEC
        VecDoub x(nDetectors), y(nDetectors);
        for(size_t i=0; i<nDetectors; i++){
                x[i] = dataArray->detectors[di[i]].azOffset ;
                y[i] = dataArray->detectors[di[i]].elOffset ;
        }
        // Covariance Matrix
        MatDoub Cov(nDetectors,nDetectors);
        //Cov = CovMat(dataVector);
        Cov = IdMat(nDetectors);
        // Construimos el primer Template
        MatDoub T1(nDetectors,nSamples) ;
        T1 = Quadratic (dataVector, Cov,x,y) ;

        // Limpieza (
        // Substraer el template de los datos
        MatDoub Residuals_1(nDetectors,nSamples);
//        writeMatOut ("D1.txt", dataVector);
        for(size_t i=0 ; i<nDetectors ; i++){
            for(size_t j=0 ; j<nSamples ; j++){
        //        Residuals_1[i][j] = dataVector[i][j] - T1[i][j] ;
            	dataVector[i][j]-=T1[i][j];
            }
        }

//        writeMatOut("T1.txt", T1 );
        Cov = CovMat(dataVector);
        T1 = Quadratic (dataVector, Cov,x,y) ;
        for(size_t i=0 ; i<nDetectors ; i++){
                    for(size_t j=0 ; j<nSamples ; j++){
                //        Residuals_1[i][j] = dataVector[i][j] - T1[i][j] ;
                    	dataVector[i][j]-=T1[i][j];
                    }
                }
//        writeMatOut("T2.txt", T1 );
//        exit (-1);

        // 2do Nivel de limpieza:
//        double criterion[2] = {5 , 0.6} ;
//        // at least [0] bolometers correlated more than [1]
//        bool withCov = 1;
//        string method = "Quadratic";
//        MatDoub T2(nDetectors,nSamples);
//        T2 = subTemplate (Residuals_1, x, y, criterion,withCov,method);
//        MatDoub Residuals_2(nDetectors,nSamples);
//        for(size_t i=0 ; i<nDetectors ; i++){
//            for(size_t j=0 ; j<nSamples ; j++){
//                Residuals_2[i][j] = Residuals_1[i][j] - T2[i][j] ;
//            }
//        }
        // After the atm template subtraction is done move back the data to the Detector object
		for (size_t i=0; i<nDetectors; i++)
					for (size_t j=0; j<nSamples; j++)
						dataArray->detectors[di[i]].hValues[si+j]=dataVector[i][j];
        // Limpieza )
        // MODIFICACION  )
	}
    cout<<"CleanHigh(): done ";
	return true;
}  // end clean


MatDoub CleanHigh::subTemplate (MatDoub tods, VecDoub x, VecDoub y, double crit[2],
                                bool withCov, string method)
    {
    size_t nDetectors = tods.nrows();
    size_t nSamples = tods.ncols();
        
    MatDoub T(nDetectors,nSamples);
    for (size_t i=0 ; i<nDetectors ; i++){
        for (size_t j=0 ; j<nSamples ; j++){
            T[i][j] = 0. ;
        }
    }
    
    MatDoub Cor(nDetectors,nDetectors);
    Cor = CorMat(tods);
    for (size_t nb=0 ; nb<nDetectors ; nb++){
        std::vector<int> ind;
        size_t nb_index = 0;
        // Notice that ind[nb_index] == nb
        for (size_t j=0 ; j<nDetectors ; j++){
            if (Cor[nb][j] > crit[1]) { ind.push_back(j) ;}
            if (j == nb) {nb_index = ind.size()-1 ;}
        }
        size_t subnDetectors = ind.size();
        if ( subnDetectors <= crit[0] ) { continue; }
        else {
            MatDoub dB(subnDetectors , nSamples);
            VecDoub subx(subnDetectors);
            VecDoub suby(subnDetectors);
            for (size_t i=0 ; i<subnDetectors ; i++){
                subx[i] = x[ ind[i] ] ;
                suby[i] = y[ ind[i] ] ;
                for (size_t ns=0 ; ns<nSamples ; ns++){
                    dB[i][ns] = tods[ ind[i] ][ns];
                }
            }
            MatDoub subCor(subnDetectors,subnDetectors);
            subCor = CorMat(dB);
            for (size_t i=0 ; i<subnDetectors ; i++){
                if (subCor[ nb_index ][i] < 0.){
                    for (size_t ns=0 ; ns<nSamples ; ns++){
                        dB[i][ns] *= -1. ;
                    }
                }
            }
            
            // Planar fit uses 3 parameters, for each ns needs 3 data-points
            if(method == "Quadratic" && subnDetectors < 6){method = "Planar"; }
            if(method == "Planar" && subnDetectors < 3){method = "Average"; }
            if(subnDetectors <= 2 ){
                method = "Average";
                withCov = false;
            }
            
            MatDoub subCov(subnDetectors,subnDetectors);
            if (withCov) {subCov = CovMat(dB);}
            else {subCov = IdMat(subnDetectors);}
            
            MatDoub subT(subnDetectors,nSamples);
            if (method == "Average") { subT = CleanHigh::Average(dB, subCov ) ; }
            else if (method == "Planar") { subT = CleanHigh::Planar(dB, subCov, subx, suby ) ;}
            else if (method == "Quadratic"){ subT = CleanHigh::Quadratic(dB, subCov, subx, suby );}
            else {
                cerr<<"\nyou should specify either Average, Planar or Quadratic method"<<endl;
                exit(-1);
            }
            for (size_t ns=0 ; ns<nSamples ; ns++){
                T[nb][ns] = subT[nb_index][ns] ;
            }
        }
    }
    return T;
    
} // subTemplate


void CleanHigh::fullMedianSustraction (){
    size_t nScans = telescope->scanIndex.ncols();
    size_t si=0;
    size_t ei=0;
    size_t nSamples=0;
    
    int *di = dataArray->getDetectorIndices();
    size_t nDetectors = dataArray->getNDetectors();
    
    double scanMedian=0;
    
    for(size_t k=0;k<nScans;k++){
        si=telescope->scanIndex[0][k];
        ei=telescope->scanIndex[1][k]+1;
        nSamples = ei -si;
        
        for (size_t i=0; i<nDetectors; i++){
            scanMedian = median(&dataArray->detectors[di[i]].hValues[si], nSamples);
            for (size_t j=0; j<nSamples; j++)
                dataArray->detectors[di[i]].hValues[si+j]-=scanMedian;
        }
    }
}  // fullMedianSustraction




void CleanHigh::fullLinearCorrection (){
    //Get the full number of scans
    size_t nScans = telescope->scanIndex.ncols();
    size_t totNSamples = 0;
    size_t si=0;
    size_t ei=0;
    size_t nSamples=0;
    
    int *di = dataArray->getDetectorIndices();
    size_t nDetectors = dataArray->getNDetectors();
    
    for(size_t k=0 ; k<nScans ; k++){
        //cerr<<"CleanBspline()::clean. Starting cleaning on scan: "<<k <<" of "<<nScans<<endl;
        si=telescope->scanIndex[0][k];
        ei=telescope->scanIndex[1][k]+1;
        totNSamples +=ei-si;
    }
    
    double *fullAverage = new double [totNSamples];
    double *singleSampleAverage = new double [nDetectors];
    bool *singleSampleFlags = new bool [nDetectors];
    size_t cSample = 0;
    
    for (size_t k=0;k<nScans;k++){
        si=telescope->scanIndex[0][k];
        ei=telescope->scanIndex[1][k]+1;
        nSamples = ei -si;
        for (size_t j=0;j<nSamples; j++){
            for (size_t i=0; i<nDetectors; i++){
                singleSampleAverage[i]=dataArray->detectors[di[i]].hValues[si+j];
                singleSampleFlags[i]=dataArray->detectors[di[i]].hSampleFlags[si+j];
            }
            fullAverage[cSample++]=mean(singleSampleAverage, nDetectors);
        }
    }
    
    delete [] singleSampleAverage;
    delete [] singleSampleFlags;
    double *currentBolo = new double [totNSamples];
    double cv0,cv1,cv2,chisq;
    
    double relOffset;
    double relGain;
    
    for (size_t i=0; i< nDetectors; i++){
        cSample = 0;
        for (size_t k=0; k<nScans; k++){
            si=telescope->scanIndex[0][k];
            ei=telescope->scanIndex[1][k]+1;
            nSamples = ei -si;
            for (size_t j=0; j<nSamples; j++){
                currentBolo[cSample++] = dataArray->detectors[di[i]].hValues[si+j];
            }
        }
        gsl_fit_linear(fullAverage,1,currentBolo,1,nSamples,&relOffset,&relGain,&cv0,&cv1,&cv2,&chisq);
        for (size_t j=0; j< dataArray->detectors[di[i]].hValues.size(); j++)
            dataArray->detectors[di[i]].hValues[j]=dataArray->detectors[di[i]].hValues[j]/relGain-relOffset;
    }
    
    delete [] fullAverage;
    delete [] currentBolo;
}  // end fullLinearCorrection

 
MatDoub CleanHigh::Average (MatDoub tods, MatDoub C){
    /* Toma los TODs y las posiciones y realiza un ajuste planar,
     creando un Template y regresadolo en formato MatDoub
     */
    
    size_t nPars = 1;
    size_t nDetectors = tods.nrows() ;
    size_t nSamples   = tods.ncols() ;
    
    // TODs en formato de gsl
    gsl_matrix * TOD = gsl_matrix_alloc (nDetectors,nSamples);
    gsl_matrix_set_zero (TOD) ;
    for (size_t nb=0 ; nb<nDetectors ; nb++){
        for (size_t ns=0 ; ns<nSamples ; ns++){
            gsl_matrix_set (TOD, nb, ns,   tods[nb][ns]  ) ;
        }
    }
    // reescribir la matriz de Covarianza en formato gsl
    gsl_matrix * Cov = gsl_matrix_alloc (nDetectors,nDetectors);
    gsl_matrix_set_zero (Cov) ;
    for (size_t nb_row=0 ; nb_row<nDetectors ; nb_row++){
        for (size_t nb_col=0 ; nb_col<nDetectors ; nb_col++){
            gsl_matrix_set (Cov, nb_row, nb_col,   C[nb_row][nb_col]  ) ;
        }
    }
    // definir la matriz gsl-S y rellenarla
    gsl_matrix * S = gsl_matrix_alloc (nDetectors, nPars);
    gsl_matrix_set_zero (S) ; // initiate to zero
    for (size_t nb=0 ; nb<nDetectors ; nb++){
        gsl_matrix_set (S, nb, 0, 1.0       ) ;
    }
    
    MatDoub T(nDetectors,nSamples);
    T = BuildTemplate(TOD , Cov , S);

    
    gsl_matrix_free(TOD);gsl_matrix_free(Cov);gsl_matrix_free(S);
    return T ;
} // Average Template


MatDoub CleanHigh::Planar (MatDoub tods, MatDoub C, VecDoub x, VecDoub y){
    /* Toma los TODs y las posiciones y realiza un ajuste planar,
     creando un Template y regresadolo en formato MatDoub
     */
    
    size_t nPars = 3;
    size_t nDetectors = tods.nrows() ;
    size_t nSamples   = tods.ncols() ;
    
    // TODs en formato de gsl
    gsl_matrix * TOD = gsl_matrix_alloc (nDetectors,nSamples);
    gsl_matrix_set_zero (TOD) ;
    for (size_t nb=0 ; nb<nDetectors ; nb++){
        for (size_t ns=0 ; ns<nSamples ; ns++){
            gsl_matrix_set (TOD, nb, ns,   tods[nb][ns]  ) ;
        }
    }
    // reescribir la matriz de Covarianza en formato gsl
    gsl_matrix * Cov = gsl_matrix_alloc (nDetectors,nDetectors);
    gsl_matrix_set_zero (Cov) ;
    for (size_t nb_row=0 ; nb_row<nDetectors ; nb_row++){
        for (size_t nb_col=0 ; nb_col<nDetectors ; nb_col++){
            gsl_matrix_set (Cov, nb_row, nb_col,   C[nb_row][nb_col]  ) ;
        }
    }
    // definir la matriz gsl-S y rellenarla
    gsl_matrix * S = gsl_matrix_alloc (nDetectors, nPars);
    gsl_matrix_set_zero (S) ; // initiate to zero
    for (size_t nb=0 ; nb<nDetectors ; nb++){
        gsl_matrix_set (S, nb, 0, 1.0       ) ;
        gsl_matrix_set (S, nb, 1, x[nb] );
        gsl_matrix_set (S, nb, 2, y[nb] );
    }
    
    MatDoub T(nDetectors,nSamples);
    T = BuildTemplate(TOD , Cov , S);

    
    gsl_matrix_free(TOD);gsl_matrix_free(Cov);gsl_matrix_free(S);
    return T ;
} // Planar Template


MatDoub CleanHigh::Quadratic (MatDoub tods, MatDoub C, VecDoub x, VecDoub y){
    /* Toma los TODs y las posiciones y realiza un ajuste cuadratico,
     creando un Template y regresadolo en formato MatDoub
     */
    
    size_t nPars = 6;
    size_t nDetectors = tods.nrows() ;
    size_t nSamples   = tods.ncols() ;
    
    // TODs en formato de gsl
    gsl_matrix * TOD = gsl_matrix_alloc (nDetectors,nSamples);
    gsl_matrix_set_zero (TOD) ;
    for (size_t nb=0 ; nb<nDetectors ; nb++){
        for (size_t ns=0 ; ns<nSamples ; ns++){
            gsl_matrix_set (TOD, nb, ns,   tods[nb][ns]  ) ;
        }
    }
    // reescribir la matriz de Covarianza en formato gsl
    gsl_matrix * Cov = gsl_matrix_alloc (nDetectors,nDetectors);
    gsl_matrix_set_zero (Cov) ;
    for (size_t nb_row=0 ; nb_row<nDetectors ; nb_row++){
        for (size_t nb_col=0 ; nb_col<nDetectors ; nb_col++){
            gsl_matrix_set (Cov, nb_row, nb_col,   C[nb_row][nb_col]  ) ;
        }
    }
    // definir la matriz gsl-S y rellenarla
    gsl_matrix * S = gsl_matrix_alloc (nDetectors, nPars);
    gsl_matrix_set_zero (S) ; // initiate to zero
    for (size_t nb=0 ; nb<nDetectors ; nb++){
        gsl_matrix_set (S, nb, 0, 1.0       ) ;
        gsl_matrix_set (S, nb, 1, x[nb]     );
        gsl_matrix_set (S, nb, 2, y[nb]     );
        gsl_matrix_set (S, nb, 3, x[nb]*y[nb] );
        gsl_matrix_set (S, nb, 4, x[nb]*x[nb] );
        gsl_matrix_set (S, nb, 5, y[nb]*y[nb] );
    }

    MatDoub T(nDetectors,nSamples);
    T = BuildTemplate(TOD , Cov , S);

    
    gsl_matrix_free(TOD);gsl_matrix_free(Cov);gsl_matrix_free(S);
    return T ;
} // Quadratic Template


MatDoub CleanHigh::BuildTemplate (gsl_matrix *tods, gsl_matrix *C, gsl_matrix *S){
    size_t nDetectors = tods->size1;
    size_t nSamples = tods->size2;
    size_t nPars = S->size2;
    
    // invertir la matriz C
    int entero;
    gsl_matrix * Covinv = gsl_matrix_alloc (nDetectors,nDetectors);
    gsl_permutation * permDetect = gsl_permutation_alloc (nDetectors);
    gsl_linalg_LU_decomp (C, permDetect, &entero);
    gsl_linalg_LU_invert (C, permDetect, Covinv);
    
    // SC = S.T * C.I
    gsl_matrix * SC = gsl_matrix_alloc (nPars, nDetectors);
    gsl_matrix_set_zero (SC) ;
    gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, S, Covinv, 0.0, SC);
    // SCS = S.T * C.I * S
    gsl_matrix * SCS = gsl_matrix_alloc (nPars, nPars);
    gsl_matrix_set_zero (SCS) ;
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, SC, S, 0.0, SCS);
    // SCS.I = (S.T * C.I * S).I
    entero = 0;
    gsl_matrix * SCSinv = gsl_matrix_alloc (nPars,nPars);
    gsl_matrix_set_zero (SCSinv) ;
    gsl_permutation * permPars = gsl_permutation_alloc (nPars);
    gsl_linalg_LU_decomp (SCS, permPars, &entero);
    gsl_linalg_LU_invert (SCS, permPars, SCSinv);
    // SCSS = (S.T * C.I * S).I * S.T
    gsl_matrix * SCSS = gsl_matrix_alloc (nPars, nDetectors);
    gsl_matrix_set_zero (SCSS) ;
    gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, SCSinv, S, 0.0, SCSS);
    // SCSSC = (S.T * C.I * S).I * S.T * C.I
    gsl_matrix * SCSSC = gsl_matrix_alloc (nPars, nDetectors);
    gsl_matrix_set_zero (SCSSC) ;
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, SCSS, Covinv, 0.0, SCSSC);
    // SSCSSC = S * (S.T * C.I * S).I * S.T * C.I
    gsl_matrix * SSCSSC = gsl_matrix_alloc (nDetectors, nDetectors);
    gsl_matrix_set_zero (SSCSSC) ;
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, S, SCSSC, 0.0, SSCSSC);
    // T = SSCSSC * tods
    gsl_matrix * T = gsl_matrix_alloc (nDetectors, nSamples);
    gsl_matrix_set_zero (T) ;
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, SSCSSC, tods, 0.0, T);
    
    // The template in MatDoub format
    MatDoub Template(nDetectors,nSamples) ;
    // Fill the Template in MatDoub
    for (size_t nb=0 ; nb<nDetectors ; nb++){
        for(size_t ns=0 ; ns<nSamples ; ns++){
            Template[nb][ns] = gsl_matrix_get(T,nb,ns) ;
        }
    }


    gsl_matrix_free(Covinv);gsl_matrix_free(SC);gsl_matrix_free(SCS);
    gsl_matrix_free(SCSinv);gsl_matrix_free(SCSS);gsl_matrix_free(SCSSC);
    gsl_matrix_free(SSCSSC);gsl_matrix_free(T);//gsl_matrix_free(tods);

    gsl_permutation_free(permDetect); gsl_permutation_free(permPars);

    return Template;
} // BuildTemplate


MatDoub CleanHigh::IterTemp (MatDoub coef , MatDoub tods ) {
    // Construye un nuevo template Tnew = c * d
    // producto de una iteracion, en terminos de los coeficientes de correlacion y la senal
    size_t nSamples = tods.nrows();
    size_t nDetectors = tods.ncols();
    // TT = tods.T * c / sum(c) [=] (nsamples,ndetector)
    MatDoub Tnew (nSamples,nDetectors) ;
    for (size_t bolcol=0 ; bolcol<nDetectors ; bolcol++){
        
        double csum = 0. ;
        for (size_t bolrow=0 ; bolrow<nDetectors ; bolrow++){
            csum += coef[bolrow][bolcol] ;
        }
        
        for (size_t ns=0 ; ns<nSamples ; ns++){
            Tnew[ns][bolcol] = 0. ; // initialize to zero
            for (size_t k=0 ; k<nDetectors ; k++){
                Tnew[ns][bolcol] += tods[ns][k] * coef[k][bolcol] / csum ;
            }
        }
    }
    return Tnew;
} // IterTemp


MatDoub CleanHigh::Coefficients (MatDoub Template, MatDoub tods){
    // Calcula los coeficientes de correlacion, c = Temp * d.T / Temp^2_n :
    // la correlacion entre el template con la senal de los bolometros.
    size_t nSamples     = tods.nrows();
    size_t nDetectors   = tods.ncols();
    // Primero sacamos la traspuesta de los TODs
    MatDoub todsT(nSamples,nDetectors) ;
    for(size_t R=0 ; R<nSamples ; R++){
        for(size_t C=0 ; C<nDetectors ; C++){
            todsT[R][C] = tods[C][R] ;
        }
    }
    MatDoub coefs(nDetectors,nDetectors);
    // coefs = tods.T * Temp / Temp^2_n
    for (size_t bolcol = 0 ; bolcol<nDetectors ; bolcol++){
        // primero sumamos el factor de normalizacion T^2[nbol][0..nsample]
        double Tsq = 0. ;
        for (size_t ns=0 ; ns<nSamples ; ns++){
            double t = Template[ns][bolcol] ;
            Tsq += t*t ;
        }
        // multiplicamos Temp y TODs dividiendo por el factor de normalizacion
        for (size_t bolrow=0 ; bolrow<nDetectors ; bolrow++){
            coefs[bolrow][bolcol] = 0. ; // initialize coefs to zero
            for (size_t ns=0 ; ns<nSamples ; ns++){
                coefs[bolrow][bolcol] += todsT[bolrow][ns] * Template[ns][bolcol] / Tsq ;
            }
        }
        Tsq = 0. ;
    }
    return coefs ;
    
} // Coefficients

MatDoub CleanHigh::Trp (MatDoub M){
    size_t nCols = M.nrows();
    size_t nRows = M.ncols();
    MatDoub MT(nRows,nCols) ;
    for(size_t R=0 ; R<nRows ; R++){
        for(size_t C=0 ; C<nCols ; C++){
            MT[R][C] = M[C][R] ;
        }
    }
    return MT ;
}

MatDoub CleanHigh::MatrixProduct (MatDoub A, MatDoub B){
    size_t nRows = A.nrows();
    size_t nCols = B.ncols();
    size_t nIntern = A.ncols();
    MatDoub M(nRows,nCols) ;
    if(  nIntern != size_t( B.nrows() )  ){
        cerr<<"Attempted to multiply matrices not aligned"<<endl;
        exit(-1);
    }
    for(size_t R=0 ; R<nRows ; R++){
        for(size_t C=0 ; C<nCols ; C++){
            M[R][C] = 0. ;
            for(size_t k=0 ; k<nIntern ; k++){
                M[R][C] += A[R][k] * B[k][C];
            }
        }
    }
    return M ;
}

MatDoub CleanHigh::IdMat (size_t size){
    MatDoub I(size,size);
    for (size_t i=0 ; i<size ; i++){
        for (size_t j=0 ; j<size ; j++){
            if (i == j) {I[i][j] = 1. ;}
            else {I[i][j] = 0. ;}
        }
    }
    return I;
}

MatDoub CleanHigh::CovMat (MatDoub d){
    size_t nRows = d.nrows();
    size_t nCols = d.ncols();
    MatDoub Cov(nRows,nRows);
    
    MatDoub dT(nCols,nRows);
    dT = Trp(d);
    
    for(size_t i=0 ; i<nRows ; i++){
        for(size_t j=0 ; j<nRows ; j++){
            Cov[i][j] = 0. ;
            for(size_t k=0 ; k<nCols ; k++){
                Cov[i][j] += d[i][k] * dT[k][j] / nCols ;
            }
        }
    }
    return Cov;
}

MatDoub CleanHigh::CorMat (MatDoub d){
    size_t nRows = d.nrows();
    MatDoub Cov(nRows,nRows);
    MatDoub Cor(nRows,nRows);
    
    Cov = CovMat(d);
    for (size_t i=0 ; i<nRows ; i++){
        for (size_t j=0 ; j<nRows ; j++){
            Cor[i][j] = Cov[i][j]/( sqrt(Cov[i][i]*Cov[j][j]) );
        }
    }
    return Cor;
}

CleanHigh::~CleanHigh(){

}

