#include <stdio.h>
#include <cblas.h>

#include "mnblas.h"
#include "complex.h"

/*
  Mesure des cycles
*/

#include <x86intrin.h>


#define VECSIZE    1000


//===================================================DEFINITION=============================================================================================//

typedef float vfloat[VECSIZE] ;
typedef double vdouble[VECSIZE];
typedef struct complex_simple vcsimple[VECSIZE];
typedef struct complex_double vcdouble[VECSIZE];


vfloat vec1,vec2;
float resultatf;
vdouble vecd1,vecd2;
double resultatd;
vcsimple veccs1,veccs2;
struct complex_simple resultatcs;
vcdouble veccd1,veccd1;
struct complex_double resultatcd;


//=======================================================================================================================================================//


int main (int argc, char **argv)
{
 unsigned long long start, end ;
 unsigned long long residu ;

 /* Calcul du residu de la mesure */
  start = _rdtsc () ;
  end = _rdtsc () ;
  residu = end - start ;


//====================================================vecteur float===========================================================//
  vector_init (vec1, 1.0) ;
  vector_init (vec2, 3.0) ;

  start = _rdtsc () ;
     cblas_saxpy (VECSIZE,2.0, vec1, 1,vec2,1) ;
  end = _rdtsc () ;

  printf ("cblas_sswap nombre de cycles cblas: %Ld \n", end-start-residu) ;
  vector_print(vec2);

  vector_init (vec2,3.0);

  start = _rdtsc () ;
     mncblas_saxpy (VECSIZE,2.0, vec1, 1,vec2,1) ;
  end = _rdtsc () ;

  printf ("mncblas_sswap: nombre de cycles: %Ld \n", end-start-residu) ;
  vector_print(vec2);

//============================================================================================================================//


//====================================================vecteur double===========================================================//
  vector_init_double(vecd1,1.0);
  vector_init_double(vecd2,3.0);

  start = _rdtsc () ;
     cblas_daxpy (VECSIZE,2.0,vecd1, 1,vecd2,1) ;
  end = _rdtsc () ;

  printf ("cblas_dswapy: nombre de cycles: %Ld \n", end-start-residu) ;
  vector_print_double(vecd2);

  vector_init_double(vecd2,3.0);

  start = _rdtsc () ;
     mncblas_daxpy (VECSIZE,2.0, vecd1, 1,vecd2,1) ;
  end = _rdtsc () ;

  printf ("mncblas_dswapy: nombre de cycles: %Ld \n", end-start-residu) ;
  vector_print_double(vecd2);

//============================================================================================================================//

//====================================================vecteur complex_simple===========================================================//
  struct complex_simple x;
  x.real = 2.0;
  x.imaginary = 1.5;
  vector_init_csimple(veccs1,x);

  x.real = 1.0;
  x.imaginary = 3.0;
  vector_init_csimple(veccs2,x);

  struct complex_simple w;
  w.real = 1.0;
  w.imaginary = -0.5;


  start = _rdtsc () ;
     cblas_caxpy (VECSIZE,w, veccs1, 1,veccs2,1) ;
  end = _rdtsc () ;

  printf ("cblas_cswap: nombre de cycles: %Ld \n", end-start-residu) ;
  vector_print_vcsimple(veccs2);

  vector_init_csimple(veccs2,x);

  start = _rdtsc () ;
     resultatcs = mncblas_caxpy (VECSIZE,w, veccs1, 1,veccs2,1) ;
  end = _rdtsc () ;

  printf ("mncblas_cswap: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f\n",resultatcs);
  vector_print_vcsimple(veccs2);

//=====================================================================================================================================//


//====================================================vecteur complex_double===========================================================//
  struct complex_double y;
  y.real = 2.0;
  y.imaginary = 1.5;
  vector_init_cdouble(veccd1,y);

  y.real = 1.0;
  y.imaginary = 3.0;
  vector_init_cdouble(veccd2,y);

  struct complex_double z;
  z.real = -2.0;
  z.imaginary = 1.5;

  start = _rdtsc () ;
     resultatcd = cblas_dzaxpy (VECSIZE,z, veccd1, 1,veccd2,1) ;
  end = _rdtsc () ;

  printf ("cblas_zswap: nombre de cycles: %Ld \n", end-start-residu) ;
  vector_print_vcdouble(veccd2);

  vector_init_cdouble(veccd2,y);

  start = _rdtsc () ;
     resultatcd = mncblas_dzaxpy (VECSIZE,z, veccd1, 1,veccd2,1) ;
  end = _rdtsc () ;

  printf ("cblas_zswap: nombre de cycles: %Ld \n", end-start-residu) ;
  vector_print_vcdouble(veccd2);
//=====================================================================================================================================//


}
