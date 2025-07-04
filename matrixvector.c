#include "matrixvector.h"

int* AllocPiv(int dim)
{
  int* ipiv;

  ipiv = (int*)malloc(dim*sizeof(int));

  /* Check memory */
  if (ipiv == NULL) {
    free(ipiv);
    fprintf(stderr, "ipiv was not allocated!\n");
    exit(EXIT_FAILURE);
  }

  return ipiv; /* Return the first address of the allocated memory */
}

double* AllocVec(int dim)
{
  double* vec;
  vec = (double*)calloc(dim, sizeof(double));

  /* Check memory */
  if (vec == NULL) {
    free(vec);
    fprintf(stderr, "Vector was not allocated!\n");
    exit(EXIT_FAILURE);
  }

  return vec; /* Return the first address of the allocated memory */
}

void lapack_LUdecomp(double* A, int* ipiv, int n) {
  int lda, info;

  lda = n;

  dgetrf_(&n, &n, A, &lda, ipiv, &info);

  if (info != 0) {
    fprintf(stderr, "error in lapack_LUdecomp!\n");
    exit(EXIT_FAILURE);
  }
}

void lapack_linsolve(double* A, double* b, int* ipiv, int n) {
  char trans = 'N';
  int nrhs, lda, ldb, info;

  nrhs = 1;
  lda = n;
  ldb = n;

  dgetrs_(&trans, &n, &nrhs, A, &lda, ipiv, b, &ldb, &info);

  if (info != 0) {
    fprintf(stderr, "error in lapack_linsolve!\n");
    exit(EXIT_FAILURE);
  }
}

double lapack_infnorm(double* A, double* work, int n) {
  char norm = 'I';
  int lda;
  double infnorm;

  lda = n;

  infnorm = dlange_(&norm, &n, &n, A, &lda, work);

  return infnorm;
}

void lapack_inverse(double* A, int* ipiv, double* work, int n) {
  int lda, lwork, info;

  lda = n;
  lwork = n;

  dgetri_(&n, A, &lda, ipiv, work, &lwork, &info);

  if (info != 0) {
    fprintf(stderr, "error in lapack_inverse!\n");
    exit(EXIT_FAILURE);
  }
}
