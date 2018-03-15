#include "mnblas.h"
#include "complex.h"

void mncblas_scopy(const int N, const float *X, const int incX,
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      Y [j] = X [i] ;
    }

  return ;
}

void mncblas_dcopy(const int N, const double *X, const int incX,
                 double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      Y [j] = X [i] ;
    }

  return ;
}

void mncblas_ccopy(const int N, const void *X, const int incX,
                void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      ((struct complex_simple*)Y)[j].real = ((struct complex_simple*)X)[i].real ;
      ((struct complex_simple*)Y)[j].imaginary = ((struct complex_simple*)X)[i].imaginary ;
    }

  return ;
}

void mncblas_zcopy(const int N, const void *X, const int incX,
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      ((struct complex_double*)Y)[j].real = ((struct complex_double*)X)[i].real ;
      ((struct complex_double*)Y)[j].imaginary = ((struct complex_double*)X)[i].imaginary ;
    }

  return ;
}
