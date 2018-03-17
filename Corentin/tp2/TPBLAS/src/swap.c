#include "mnblas.h"

void mncblas_sswap(const int N, float *X, const int incX,
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_dswap(const int N, double *X, const int incX,
                 double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register double save ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_cswap(const int N, void *X, const int incX,
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register type_c save ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = ((type_c *)Y)[j] ;
      ((type_c *)Y)[j] = ((type_c *)X)[i] ;
      ((type_c *)X)[i] = save ;
    }

  return ;
}

void mncblas_zswap(const int N, void *X, const int incX,
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register type_z save ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = ((type_z *) Y)[j] ;
      ((type_z *)Y)[j] = ((type_z *) X)[i] ;
      ((type_z *)X)[i] = save ;
    }

  return ;
}
