#include "mnblas.h"
#include "complex.h"

void mncblas_saxpy(const int N, const float a, const float *X, const int incX,
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      Y [j] = a* X[i] + Y[j] ;
    }

  return ;
}

void mncblas_daxpy(const int N, const double a, const double *X, const int incX,
                 double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      Y [j] = a *X[i] + Y[j] ;
    }

  return ;
}

void mncblas_caxpy(const int N, const void* a, const void *X, const int incX,
                void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      ((struct complex_simple*)Y)[j] = addition_cs( multiplication_cs ( *((struct complex_simple*)a) , ((struct complex_simple*)X)[i] ) , ((struct complex_simple*)Y)[j] );
    }

  return ;
}

void mncblas_zaxpy(const int N, const void* a, const void *X, const int incX,
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      ((struct complex_double*)Y)[j] = addition_cd( multiplication_cd ( *((struct complex_double*)a) , ((struct complex_double*)X)[i] ) , ((struct complex_double*)Y)[j] );
    }

  return ;
}
