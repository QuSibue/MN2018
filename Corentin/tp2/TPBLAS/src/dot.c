#include "mnblas.h"

float mncblas_sdot(const int N, const float *X, const int incX,
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float dot = 0.0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      dot = dot + X [i] * Y [j] ;
    }

  return dot ;
}

double mncblas_ddot(const int N, const double *X, const int incX,
                 const double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register double dot = 0.0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      dot = dot + X [i] * Y [j] ;
    }

  return dot;
}

void mncblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  type_c* vect1 = (type_c*) X;
  type_c* vect2 = (type_c*) Y;
  type_c* res = (type_c*) dotu;

  res->re = 0.0 ;
  res->im = 0.0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      res->re += X[i].re * Y[j].re ;
      res->re -= X[i].im * Y[j].im ;
      res->im += X[i].im * Y[j].re ;
      res->im += X[i].re * Y[j].im ;
    }

  return dot ;
}

void mncblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  type_c* vect1 = (type_c*) X;
  type_c* vect2 = (type_c*) Y;
  type_c* res = (type_c*) dotc;

  res->re = 0.0 ;
  res->im = 0.0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      res->re  += vect1[i].re * vect2[j].re;
      res->re  += vect1[i].im * vect2[j].im;
      res->im  += vect1[i].re * vect2[j].im;
      res->im  -= vect1[i].im * vect2[j].re;
    }

  return;
}

void mncblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  type_z* vect1 = (type_z*) X;
  type_z* vect2 = (type_z*) Y;
  type_z* res = (type_z*) dotu;

  res->re = 0.0 ;
  res->im = 0.0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      res->re += X[i].re * Y[j].re ;
      res->re -= X[i].im * Y[j].im ;
      res->im += X[i].im * Y[j].re ;
      res->im += X[i].re * Y[j].im ;
    }

  return dot ;
}

void mncblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  type_z* vect1 = (type_z*) X;
  type_z* vect2 = (type_z*) Y;
  type_z* res = (type_z*) dotc;

  res->re = 0.0 ;
  res->im = 0.0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      res->re  += vect1[i].re * vect2[j].re;
      res->re  += vect1[i].im * vect2[j].im;
      res->im  += vect1[i].re * vect2[j].im;
      res->im  -= vect1[i].im * vect2[j].re;
    }

  return;
}
