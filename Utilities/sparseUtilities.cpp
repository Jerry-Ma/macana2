#include <iostream>
#include <omp.h>

#include "sparseUtilities.h"

using namespace std;

CS_INT cs_nonfinite (CS_INT i, CS_INT j, CS_ENTRY aij, void *other)
{
	CS_INT res = finite(aij);
	if (!res){
		cerr<<"Detected a NaN value  on sparse matrix"<<endl;
		exit(-1);
	}
    return (res && aij!=0) ;
}
CS_INT cs_dropnotfinite (cs *A)
{
    return (cs_fkeep (A, &cs_nonfinite, NULL)) ; /* keep all nonzero entries */
}

cs *sparse_invert (cs* matrix, double tol)
{
  cs *inverted =NULL;
  cs *inverted_compress=NULL;
  long size= 0;
  double *unitVector = NULL;

  if (!matrix)
     return NULL;

  if (matrix->n != matrix->m){
    return NULL;
  }

  size = matrix->n;
  inverted = cs_spalloc(size,size,matrix->nzmax, 1,1);

  for (long i = 0; i<size; i++){
    unitVector = new double[size];
    for (int j=0; j<size; j++)
      unitVector[j] = 0.0;
    unitVector[i] = 1.0;
    if (!cs_qrsol(3,matrix,unitVector)){
    	cerr<<"Singular Matrix"<<endl;
    	return NULL;
    }


    for (long j=0; j<size; j++){
      if (!finite(unitVector[j])){
    	  return NULL;
      }
      if (abs(unitVector[j])>tol)
    	  cs_entry(inverted, i, j, unitVector[j]);
    }

    delete [] unitVector;
    unitVector=NULL;
  }
  inverted_compress = cs_compress(inverted); //
  cs_spfree(inverted);

  return inverted_compress;
}


cs *fast_sparse_invert (cs* matrix, double tol)
{
  cs *inverted =NULL;
  cs *inverted_compress=NULL;
  long size= 0;
  double *unitVector = NULL;

  if (!matrix)
     return NULL;

  if (matrix->n != matrix->m){
    return NULL;
  }



  size = matrix->n;
  inverted = cs_spalloc(size,size,matrix->nzmax, 1,1);

  // Now Calculate the QR decomposition of the matrix to save some time:

  css *S = cs_sqr(2,matrix,1);
  csn *N = cs_qr (matrix,S);
  double *x;

//#pragma omp parallel shared (S,N,size,tol, inverted,cerr) private (unitVector,x) default (none)
 //{

//#pragma omp for schedule(dynamic)
  for (long i = 0; i<size; i++){
	unitVector = new double[size];
	x= (double *)cs_malloc(S->m2, sizeof(double));
    for (int j=0; j<size; j++)
      unitVector[j] = 0.0;
    unitVector[i] = 1.0;
    if (!sparseSolve(S,N,x,unitVector,size)){
    	cerr<<"Singular Matrix"<<endl;
    	exit(-1);
    }

    for (long j=0; j<size; j++){
      if (unitVector[j]!=unitVector[j]){
    	  cerr<<"NaN detected on sparse inverse. Imploding"<<endl;
    	  exit(-1);
      }
      if (abs(unitVector[j])>tol)
    	  cs_entry(inverted, i, j, unitVector[j]);
    }
    delete [] unitVector;
    cs_free(x);
  }

 //}
  unitVector=NULL;
  cs_sfree(S);
  cs_nfree(N);

  inverted_compress = cs_compress(inverted); //
  cs_spfree(inverted);

  return inverted_compress;
}

size_t sparseSolve (css *S, csn *N, double *x, double *b, size_t n){
	if (S && N && x){
        cs_ipvec (S->pinv, b, x, n) ;   		/* x(0:m-1) = b(p(0:m-1) */
        for (size_t k = 0 ; k < n ; k++)       /* apply Householder refl. to x */
        {
            cs_happly (N->L, k, N->B [k], x) ;
        }
        cs_usolve (N->U, x) ;           		/* x = R\x */
        cs_ipvec (S->q, x, b, n) ;      		/* b(q(0:n-1)) = x(0:n-1) */
        return 1;
	}
	return 0;
}
