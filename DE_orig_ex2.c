#include "matrixvector.h"
#include "DE_basis_func.h"

double beta(double p, double q)
{
  return tgamma(p)*tgamma(q)/tgamma(p+q);
}

/* kDE1(t, a, b, tau) = k1(t, DE_trans(a, b, tau)) */
double kDE1(double t, double a, double b, double tau)
{
//double s = DE_trans(a, b, tau);
//return pow(s, t+0.5);
  return (b - a)*exp(- (t+0.5)*log1p(exp(-M_PI_2*sinh(tau))));
}

/* kDE2(t, a, b, tau) = k2(t, DE_trans(a, b, tau)) */
double kDE2(double t, double a, double b, double tau)
{
//double s = DE_trans(a, b, tau);
//return pow(1-s, t);
  return (b - a)*exp(- t*log1p(exp(M_PI_2*sinh(tau))));
}

/* gDE(a, b, tau) = g(DE_trans(a, b, tau)) */
double gDE(double a, double b, double tau)
{
  double t = DE_trans(a, b, tau);
  return sqrt(t) - pow(t, t+2)/(t+2) - beta(1.5, t+1);
}

double u(double t)
{
  return sqrt(t);
}

double uDEn(double a, double b, double t, double h, int N, double* f_N, int n)
{
  int j;
  double x = DE_trans_inv(a, b, t);
  double ans = 0;

  for (j = N-1; j > 0; j--) {
    ans += f_N[ j+N] * S( j, h, x);
    ans += f_N[-j+N] * S(-j, h, x);
  }
    ans += f_N[ 0+N] * S( 0, h, x);

    ans += f_N[0]*wa(a, b, t) + f_N[n-1]*wb(a, b, t);

  return ans;
}

double* substitute_xN(double a, double b, double h, int N, int n)
{
  int i;
  double* x_N = AllocVec(n);

  for (i = -N; i <= N; i++) {
    x_N[i+N] = DE_trans(a, b, i*h);
  }

  return x_N;
}

double* substitute_AN(double a, double b, double h, int N, int n, double* x_N)
{
  int i, j;
  double* A_N = AllocVec(n*n); /* Column-major */
  double* V_N = AllocVec(n*n); /* Column-major */
  double* p_N = AllocVec(n);
  double* w_N = AllocVec(n);

  for (j = -N; j <= N; j++) {
    for (i = -N; i <= N; i++) {
      V_N[i+N + (j+N)*n]
        = h*kDE1(x_N[i+N], a, b, j*h)*DE_trans_div(a, b, j*h)*del1(i, j)
         +h*kDE2(x_N[i+N], a, b, j*h)*DE_trans_div(a, b, j*h);
    }
  }

    w_N[0]   = 1;
  for (i = -N+1; i <= N-1; i++) {
    w_N[i+N] = waDE(i*h);
  }
    w_N[n-1] = 0;

  /* p_N = V_N w_N */
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1, V_N, n, w_N, 1, 0, p_N, 1);

  /* j = -N */
  j = -N;
    for (i = -N; i <= N; i++) {
      A_N[i+N + (j+N)*n] = w_N[i+N] - p_N[i+N];
    }

    w_N[0]   = 0;
  for (i = -N+1; i <= N-1; i++) {
    w_N[i+N] = wbDE(i*h);
  }
    w_N[n-1] = 1;

  /* p_N = V_N w_N */
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1, V_N, n, w_N, 1, 0, p_N, 1);

  /* j = N */
  j = N;
    for (i = -N; i <= N; i++) {
      A_N[i+N + (j+N)*n] = w_N[i+N] - p_N[i+N];
    }

  /* -N < j < N */
    j = -N;
      V_N[i+N + (j+N)*n]
        = h*kDE1(a, a, b, j*h)*DE_trans_div(a, b, j*h)*0
         +h*kDE2(a, a, b, j*h)*DE_trans_div(a, b, j*h);
    j =  N;
      V_N[i+N + (j+N)*n]
        = h*kDE1(b, a, b, j*h)*DE_trans_div(a, b, j*h)
         +h*kDE2(b, a, b, j*h)*DE_trans_div(a, b, j*h);

    for (i = -N+1; i <= N-1; i++) {
      A_N[i+N + (i+N)*n] = 1.0;
    }

  for (j = -N+1; j <= N-1; j++) {
    for (i = -N; i <= N; i++) {
      A_N[i+N + (j+N)*n] -= V_N[i+N + (j+N)*n];
    }
  }

  free(V_N);
  free(p_N);
  free(w_N);
  return A_N;
}

double* substitute_fN(double a, double b, double h, int N, int n)
{
  int i;
  double* f_N = AllocVec(n);

  for (i = -N; i <= N; i++) {
    f_N[i+N] = gDE(a, b, i*h);
  }

  return f_N;
}

int main()
{
  double a = 0.0;
  double b = 1.0;
  double d = M_PI_4;
  double alpha = 1.0;
  double *A_N, *f_N, *x_N;
  int* ipiv;
  int i, n, N;
  double h, err, maxerr, x, norm1, norm2;
  int SAMPLE = 4096;
  double hh = (b - a)/SAMPLE;

  for (N = 5; N <= 90; N += 5) {

    n = 2*N+1;
    h = log(4*d*N/alpha)/N;

    x_N = substitute_xN(a, b, h, N, n);
    A_N = substitute_AN(a, b, h, N, n, x_N);
    f_N = substitute_fN(a, b, h, N, n);
    ipiv= AllocPiv(n);

    norm1 = lapack_infnorm(A_N, x_N, n);

    lapack_LUdecomp(A_N, ipiv, n);
    lapack_linsolve(A_N, f_N, ipiv, n);

    maxerr = 0;
    for (i = 1; i < SAMPLE; i++) {
      x = a + i*hh;

      err = fabs(u(x) - uDEn(a, b, x, h, N, f_N, n));

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
