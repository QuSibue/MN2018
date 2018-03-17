#include "mnblas.h"
#include "complex.h"
#include "math.h"

float mncblas_sasum(const int N, const float *X, const int incX)
{
  register unsigned int i = 0 ;
  float resultat=0;
  for (; (i < N) ; i += incX)
    {
      resultat += fabs(X[i]) ;
    }
  return resultat;
}

double mncblas_dasum(const int N, const double *X, const int incX)
{
  register unsigned int i = 0 ;
  double resultat=0;

  for (; (i < N) ; i += incX)
    {
      resultat += fabs(X[i]);
    }

  return resultat;
}

float mncblas_scasum(const int N, const void *X, const int incX)
{
  register unsigned int i = 0 ;
  float resultat = 0;

  for (; (i < N) ; i += incX)
    {
      resultat += fabs(((struct complex_simple*)X)[i].real) + fabs(((struct complex_simple*)X) [i].imaginary );
    }

  return resultat;
}

double mncblas_dzasum(const int N, const void *X, const int incX)
{
  register unsigned int i = 0 ;
  double resultat = 0;

  for (; (i < N); i += incX)
    {
      resultat += fabs( ((struct complex_double*)X)[i].real ) + fabs( ((struct complex_double*)X)[i].imaginary );
    }

  return resultat;

}
