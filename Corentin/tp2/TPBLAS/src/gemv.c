#include "mnblas.h"

void mncblas_sgemv (const MNCBLAS_LAYOUT layout,
                 const MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX, const float beta,
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  float* tmp = malloc(M * sizeof(float));
  for (; j < N ; j += incY) {
    tmp[j] = Y[j];
  }
  float res;
  for (; i < M ; i += incX) {
    res = 0.;
    for (j=0; j < N ;j += incY) {
      res += alpha * A[j+N*i] * X[j] + beta * tmp[j];
    }
    Y[i] = res;
  }

  free(tmp);
  return;
}

void mncblas_dgemv (MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  double* tmp = malloc(M * sizeof(double));
  for (; j < N ; j += incY) {
    tmp[j] = Y[j];
  }
  double res;
  for (; i < M ; i += incX) {
    res = 0.;
    for (j=0; j < N ;j += incY) {
      res += alpha * A[j+N*i] * X[j] + beta * tmp[j];
    }
    Y[i] = res;
  }

  free(tmp);
  return;
}

void mncblas_cgemv (MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY)
{
register unsigned int i = 0 ;
register unsigned int j = 0 ;

type_c* tmp = malloc(M * sizeof(type_c));
for (; j < N ; j += incY) {
  tmp[j].re = Y[j].re;
  tmp[j].im = Y[j].im;
}
type_c res;
for (; i < M ; i += incX) {
  res.re = 0.;
  res.im = 0.;
  for (j=0; j < N ;j += incY) {
    res.re += alpha.re * A[j+N*i].re * X[j].re - alpha.re * A[j+N*i].im * X[j].im - alpha.im * A[j+N*i].re * X[j].im - alpha.im * A[j+N*i].im * X[j].re + beta.re * tmp[j].re - beta.im * tmp[j].im;
    res.im += alpha.re * A[j+N*i].re * X[j].im + alpha.re * A[j+N*i].im * X[j].re + alpha.im * A[j+N*i].re * X[j].re - alpha.im * A[j+N*i].im * X[j].im + beta.im * tmp[j].re + beta.re * tmp[j].im;
  }
  Y[i].re = res.re;
  Y[i].im = res.im;
}

free(tmp);
return;
}

void mncblas_zgemv (MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  type_z* tmp = malloc(M * sizeof(type_z));
  for (; j < N ; j += incY) {
    tmp[j].re = Y[j].re;
    tmp[j].im = Y[j].im;
  }
  type_z res;
  for (; i < M ; i += incX) {
    res.re = 0.;
    res.im = 0.;
    for (j=0; j < N ;j += incY) {
      res.re += alpha.re * A[j+N*i].re * X[j].re - alpha.re * A[j+N*i].im * X[j].im - alpha.im * A[j+N*i].re * X[j].im - alpha.im * A[j+N*i].im * X[j].re + beta.re * tmp[j].re - beta.im * tmp[j].im;
      res.im += alpha.re * A[j+N*i].re * X[j].im + alpha.re * A[j+N*i].im * X[j].re + alpha.im * A[j+N*i].re * X[j].re - alpha.im * A[j+N*i].im * X[j].im + beta.im * tmp[j].re + beta.re * tmp[j].im;
    }
    Y[i].re = res.re;
    Y[i].im = res.im;
  }

  free(tmp);
  return;
}
