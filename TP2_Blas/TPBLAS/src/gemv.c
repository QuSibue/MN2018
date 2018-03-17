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

struct complex_simple *tmp = malloc(M * sizeof(struct complex_simple));
for (; j < N ; j += incY) {
  tmp[j].real = ((struct complex_simple*)Y)[j].real;
  tmp[j].imaginary = ((struct complex_simple*)Y)[j].imaginary;
}

struct complex_simple res;
for (; i < M ; i += incX) {

    res.real = 0.;
	res.imaginary = 0.;

    for (j=0; j < N ;j += incY) {
      res = addition_cs(res,addition_cs(
								multiplication_cs( 
										multiplication_cs( *((struct complex_simple*)alpha) , ((struct complex_simple*)A)[j+N*i]) 
										,((struct complex_simple*)X)[j] 
								),
								multiplication_cs( *((struct complex_simple*)beta) , tmp[j] )
							) 
			);
    }
    ((struct complex_simple*)Y)[i] = res;
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

struct complex_double *tmp = malloc(M * sizeof(struct complex_double));
  for (; j < N ; j += incY) {
    tmp[j].real = ((struct complex_double*)Y)[j].real;
    tmp[j].imaginary = ((struct complex_double*)Y)[j].imaginary;
  }

struct complex_double res;
for (; i < M ; i += incX) {

    res.real = 0.;
	res.imaginary = 0.;

    for (j=0; j < N ;j += incY) {
      res = addition_cd(res,addition_cd(
								multiplication_cd( 
										multiplication_cd( *((struct complex_double*)alpha) , ((struct complex_double*)A)[j+N*i]) 
										,((struct complex_double*)X)[j] 
								),
								multiplication_cd( *((struct complex_double*)beta) , tmp[j] )
							) 
			);
    }
    ((struct complex_double*)Y)[i] = res;
  }
  free(tmp);
  return;
}
