#include "mnblas.h"
#include <math.h>

float mncblas_sasum(const int N, const float *X,
                const int incX, float *Y, const int incY);
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  float res = 0.0;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      res += abs(X[i]) + abs(Y[j]);
    }

  return res;
}

double mncblas_dasum(const int N, const double *X,
                const int incX, double *Y, const int incY);
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  double res = 0.0;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      res += abs(X[i]) + abs(Y[j]);
    }

  return res;
}

float mncblas_casum(const int N,const void *X,
                const int incX, void *Y, const int incY);
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  float res = 0.0;

  type_c* vect1 = (type_c*) X;
  type_c* vect2 = (type_c*) Y;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      res += abs(vect1[i].re) + abs(vect1[i].im) + abs(vect2[j].re) + abs(vect2[j].im);
    }

  return ;
}

double mncblas_zasum(const int N, const void *X,
                const int incX, void *Y, const int incY);
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  double res = 0.0;

  type_z* vect1 = (type_z*) X;
  type_z* vect2 = (type_z*) Y;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      res += abs(vect1[i].re) + abs(vect1[i].im) + abs(vect2[j].re) + abs(vect2[j].im);
    }

  return ;
}
