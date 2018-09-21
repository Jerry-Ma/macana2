#ifndef _SPARSE_UTILITIES_H
#define _SPARSE_UTILITIES_H

#include <suitesparse/cs.h>


CS_INT cs_nonfinite (CS_INT i, CS_INT j, CS_ENTRY aij, void *other);
CS_INT cs_dropnotfinite (cs *A);
cs *sparse_invert (cs* matrix, double tol = 0.0);
cs *fast_sparse_invert (cs* matrix, double tol =0.0);
size_t sparseSolve (css *S, csn *N, double *x, double *b, size_t n);

#endif
