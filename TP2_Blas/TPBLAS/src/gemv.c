#include <stdlib.h>
#include "mnblas.h"
#include "complex.h"

void mncblas_sgemv (const MNCBLAS_LAYOUT layout,
                 const MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX, const float beta,
                 float *Y, const int incY)
{
  register unsigned int i;
  register unsigned int j;

  float tmp ;
  for (i = 0 ; i < M ; i += incX) {
    tmp = A[i*N] * X[0] ;
    for (j=1; j < N ;j += incY) {
      tmp +=  A[j+N*i] * X[j] ;
    }
    Y[i] =alpha * tmp + beta * Y[i];
  }


  return;
}

void mncblas_dgemv (MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY)
{
	register unsigned int i;
  register unsigned int j;

  double tmp ;
  for (i = 0 ; i < M ; i += incX) {
    tmp = A[i*N] * X[0] ;
    for (j=1; j < N ;j += incY) {
      tmp +=  A[j+N*i] * X[j] ;
    }
    Y[i] =alpha * tmp + beta * Y[i];
  }


  return;
}

void mncblas_cgemv (MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY)
{
	register unsigned int i;
  register unsigned int j;

  struct complex_simple tmp ;
  for (i = 0 ; i < M ; i += incX) {
    tmp = multiplication_cs( ((struct complex_simple*)A)[i*N] , ((struct complex_simple*)X)[0] );
    for (j=1; j < N ;j += incY) {
      tmp =addition_cs(tmp ,  multiplication_cs( ((struct complex_simple*)A)[j+i*N] , ((struct complex_simple*)X)[j] ) );
    }
    ((struct complex_simple*)Y)[i] = addition_cs( multiplication_cs( *((struct complex_simple*)alpha) , tmp )
																									,
																								 multiplication_cs( *((struct complex_simple*)beta) , ((struct complex_simple*)Y)[i] )
																								);
  }

  return;
}

void mncblas_zgemv (MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY)
{
	register unsigned int i;
  register unsigned int j;

  struct complex_double tmp ;
  for (i = 0 ; i < M ; i += incX) {
    tmp = multiplication_cd( ((struct complex_double*)A)[i*N] , ((struct complex_double*)X)[0] );
    for (j=1; j < N ;j += incY) {
      tmp = addition_cd(tmp ,  multiplication_cd( ((struct complex_double*)A)[j+i*N] , ((struct complex_double*)X)[j] ) );
    }
    ((struct complex_double*)Y)[i] = addition_cd( multiplication_cd( *((struct complex_double*)alpha) , tmp )
																									,
																								 multiplication_cd( *((struct complex_double*)beta) , ((struct complex_double*)Y)[i] )
																								);
  }

  return;
}
