#include "matrixvector.h"
#include "SE_basis_func.h"

/* kSE1(t, a, b, tau) = k1(t, SE_trans(a, b, tau)) */
double kSE1(double t, double a, double b, double tau)
{
  return t / (1 + exp(-tau));
}

/* kSE2(t, a, b, tau) = k2(t, SE_trans(a, b, tau)) */
double kSE2(double t, double a, double b, double tau)
{
  return t / (1 + exp(-tau));
}

/* gSE(a, b, tau) = g(SE_trans(a, b, tau)) */
double gSE(double a, double b, double tau)
{
  double t = SE_trans(a, b, tau);
  return (2 - t*t*t)*t / 3.0;
}

double u(double t)
{
  return t;
}

double uSEn(double a, double b, double t, double h, int N, double* c_N, int n)
{
  int j;
  double x = SE_trans_inv(a, b, t);
  double ans = 0;

  for (j = N; j > 0; j--) {
    ans += c_N[ j+N] * S( j, h, x);
    ans += c_N[-j+N] * S(-j, h, x);
  }
    ans += c_N[ 0+N] * S( 0, h, x);

    ans += c_N[n]*wa(a, b, t) + c_N[n+1]*wb(a, b, t);

  return ans;
}

double* substitute_xN(double a, double b, double h, int N, int n)
{
  int i;
  double* x_N = AllocVec(n);

  for (i = -N; i <= N; i++) {
    x_N[i+N] = SE_trans(a, b, i*h);
  }

  return x_N;
}

double* substitute_AN(double a, double b, double h, int N, int n, double* x_N)
{
  int i, j;
  double* A_N = AllocVec(n*n); /* Column-major */

  for (i = -N; i <= N; i++) {
    A_N[i+N + (i+N)*n] = 1.0;
  }

  for (j = -N; j<= N; j++) {
    for (i = -N; i <= N; i++) {
      A_N[i+N + (j+N)*n]
        -= h*kSE1(x_N[i+N], a, b, j*h)*SE_trans_div(a, b, j*h)*del1(i, j);
      A_N[i+N + (j+N)*n]
        -= h*kSE2(x_N[i+N], a, b, j*h)*SE_trans_div(a, b, j*h);
    }
  }

  return A_N;
}

double* substitute_fN(double a, double b, double h, int N, int n)
{
  int i;
  double* f_N = AllocVec(n+2);

  for (i = -N; i <= N; i++) {
    f_N[i+N] = gSE(a, b, i*h);
  }

  return f_N;
}

void translate_fN(double h, int N, int n, double* f_N)
{
  int j;
  f_N[n+1] = f_N[n-1];
  f_N[n]   = f_N[0];

  for (j = N; j >= -N; j--) {
    f_N[j+N] = f_N[j+N] - f_N[n]*waSE(j*h) - f_N[n+1]*wbSE(j*h);
  }
}

int main()
{
  double a = 0.0;
  double b = 1.0;
  double d = 3.14;
  double alpha = 1.0;
  double *A_N, *f_N, *x_N;
  int* ipiv;
  int i, n, N;
  double h, err, maxerr, x, norm1, norm2;
  int SAMPLE = 4096;
  double hh = (b - a)/SAMPLE;

  for (N = 5; N <= 200; N += 5) {

    n = 2*N+1;
    h = sqrt(M_PI*d / (alpha * N));

    x_N = substitute_xN(a, b, h, N, n);
    A_N = substitute_AN(a, b, h, N, n, x_N);
    f_N = substitute_fN(a, b, h, N, n);
    ipiv= AllocPiv(n);

    norm1 = lapack_infnorm(A_N, x_N, n);

    lapack_LUdecomp(A_N, ipiv, n);
    lapack_linsolve(A_N, f_N, ipiv, n);

    translate_fN(h, N, n, f_N);
    maxerr = 0;
    for (i = 1; i < SAMPLE; i++) {
      x = a + i*hh;

      err = fabs(u(x) - uSEn(a, b, x, h, N, f_N, n));

      maxerr = fmax(err, maxerr);
    }
    printf("%d\t%e\t", n, maxerr);

    lapack_inverse(A_N, ipiv, x_N, n);
    norm2 = lapack_infnorm(A_N, x_N, n);
    printf("%e\n", norm1*norm2);

    free(ipiv);
    free(f_N);
    free(A_N);
    free(x_N);
  }

  return EXIT_SUCCESS;
}
