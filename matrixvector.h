#include <stdio.h>
#include <stdlib.h>
#include <Accelerate/Accelerate.h>

#if !defined(___matrixvector___)
#define ___matrixvector___

int* AllocPiv(int dim);

double* AllocVec(int dim);

void lapack_LUdecomp(double* A, int* ipiv, int n);

void lapack_linsolve(double* A, double* b, int* ipiv, int n);

double lapack_infnorm(double* A, double* work, int n);

void lapack_inverse(double* A, int* ipiv, double* work, int n);

#endif // ___matrixvector___
