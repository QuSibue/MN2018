#include "mnblas.h"

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

  type_c* vect1 = (type_c*) X;
  type_c* vect2 = (type_c*) Y;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      vect2[j].re = vect1[i].re ;
      vect2[j].im = vect1[i].im ;
    }

  return ;
}

void mncblas_zcopy(const int N, const void *X, const int incX,
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  type_z* vect1 = (type_z*) X;
  type_z* vect2 = (type_z*) Y;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      vect2[j].re = vect1[i].re ;
      vect2[j].im = vect1[i].im ;
    }

  return ;
}
