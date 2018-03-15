#include "mnblas.h"
#include "complex.h"
#include "math.h"

float mncblas_snrm2(const int N,const float *X,const int incX)
{
  register unsigned int i = 0 ;
  float res=0;
  for (; (i < N) ; i += incX)
    {
      res+=X[i]*X[i];
    }
    return sqrt(res);
}

double mncblas_dnrm2(const int N,const double *X,const int incX)
{
  register unsigned int i = 0 ;
  double res=0;
  for (; (i < N) ; i += incX)
    {
      res+=X[i]*X[i];
    }
    return sqrt(res);
}

float mncblas_scnrm2(const int N,const void *X,const int incX)
{
  register unsigned int i = 0 ;
  float res=0;

  for (; (i < N) ; i += incX)
    {
      res +=  ( ((struct complex_simple*)X)[i].real * ((struct complex_simple*)X)[i].real )  +
	            ( ((struct complex_simple*)X)[i].imaginary) * ((struct complex_simple*)X)[i].imaginary;
    }
  return sqrt(res);
}

double  mncblas_dznrm2(const int N,const void *X,const int incX)
{
  register unsigned int i = 0 ;
  double res=0;

  for (; (i < N) ; i += incX)
    {
	res += ( ((struct complex_double*)X)[i].real * ((struct complex_double*)X)[i].real )  +
                     ( ((struct complex_double*)X)[i].imaginary) * ((struct complex_double*)X)[i].imaginary;
    }
  return sqrt(res);
}
