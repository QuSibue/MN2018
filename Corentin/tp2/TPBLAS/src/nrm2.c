#include "mnblas.h"
#include <math.h>

float mncblas_snrm2(const int N, const float *X,
                const int incX, float *Y, const int incY);
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  float res = 0.0;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      res += pow(X[i], 2) + pow(Y[j], 2);
    }

  return sqrt(res);
}

double mncblas_dnrm2(const int N, const double *X,
                const int incX, double *Y, const int incY);
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  double res = 0.0;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      res += pow(X[i], 2) + pow(Y[j], 2);
    }

  return sqrt(res);
}

float mncblas_canrm2(const int N,const void *X,
                const int incX, void *Y, const int incY);
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  float res = 0.0;

  type_c* vect1 = (type_c*) X;
  type_c* vect2 = (type_c*) Y;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      res += pow(vect1[i].re, 2) + pow(vect1[i].im, 2) + pow(vect2[j].re, 2) + pow(vect2[j].im, 2);
    }

  return sqrt(res);
}

double mncblas_znrm2(const int N, const void *X,
                const int incX, void *Y, const int incY);
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  double res = 0.0;

  type_z* vect1 = (type_z*) X;
  type_z* vect2 = (type_z*) Y;

  for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      res += pow(vect1[i].re, 2) + pow(vect1[i].im, 2) + pow(vect2[j].re, 2) + pow(vect2[j].im, 2);
    }

  return sqrt(res);
}
