#include <stdlib.h>
#include "mnblas.h"
#include "complex.h"

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

complex_simple* tmp = malloc(M * sizeof(complex_simple));
for (; j < N ; j += incY) {
  tmp[j].real = Y[j].real;
  tmp[j].imaginary = Y[j].imaginary;
}
complex_simple res;
for (; i < M ; i += incX) {
  res.real = 0.;
  res.imaginary = 0.;
  for (j=0; j < N ;j += incY) {
    res.real += alpha.real * A[j+N*i].real * X[j].real - alpha.real * A[j+N*i].imaginary * X[j].imaginary - alpha.imaginary * A[j+N*i].real * X[j].imaginary - alpha.imaginary * A[j+N*i].imaginary * X[j].real + beta.real * tmp[j].real - beta.imaginary * tmp[j].imaginary;
    res.imaginary += alpha.real * A[j+N*i].real * X[j].imaginary + alpha.real * A[j+N*i].imaginary * X[j].real + alpha.imaginary * A[j+N*i].real * X[j].real - alpha.imaginary * A[j+N*i].imaginary * X[j].imaginary + beta.imaginary * tmp[j].real + beta.real * tmp[j].imaginary;
  }
  Y[i].real = res.real;
  Y[i].imaginary = res.imaginary;
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

  complex_double* tmp = malloc(M * sizeof(complex_double));
  for (; j < N ; j += incY) {
    tmp[j].real = Y[j].real;
    tmp[j].imaginary = Y[j].imaginary;
  }
  complex_double res;
  for (; i < M ; i += incX) {
    res.real = 0.;
    res.imaginary = 0.;
    for (j=0; j < N ;j += incY) {
      res.real += alpha.real * A[j+N*i].real * X[j].real - alpha.real * A[j+N*i].imaginary * X[j].imaginary - alpha.imaginary * A[j+N*i].real * X[j].imaginary - alpha.imaginary * A[j+N*i].imaginary * X[j].real + beta.real * tmp[j].real - beta.imaginary * tmp[j].imaginary;
      res.imaginary += alpha.real * A[j+N*i].real * X[j].imaginary + alpha.real * A[j+N*i].imaginary * X[j].real + alpha.imaginary * A[j+N*i].real * X[j].real - alpha.imaginary * A[j+N*i].imaginary * X[j].imaginary + beta.imaginary * tmp[j].real + beta.real * tmp[j].imaginary;
    }
    Y[i].real = res.real;
    Y[i].imaginary = res.imaginary;
  }

  free(tmp);
  return;
}
