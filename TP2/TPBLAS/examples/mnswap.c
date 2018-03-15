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


vfloat vec1, vec2 ;
vdouble vecd1,vecd2;
vcsimple veccs1, veccs2;
vcdouble veccd1,veccd2;


//=======================================================================================================================================================//


//===================================================INIT===================================================================================================//

void vector_init (vfloat V, float x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_init_double (vdouble V, double x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V[i] = x ;

  return ;
}


void vector_init_csimple (vcsimple V, struct complex_simple x)
{
  register unsigned int i ;


  for (i = 0; i < VECSIZE; i++){
    V[i].real = x.real ;
    V[i].imaginary = x.imaginary ;
  }


  return ;
}

void vector_init_cdouble (vcdouble V, struct complex_double x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++){
    V[i].real = x.real ;
    V[i].imaginary = x.imaginary ;
  }

  return ;
}

//=========================================================================================================================================================//


//===================================================PRINT==================================================================================================//

void vector_print (vfloat V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;

  return ;
}

void vector_print_double (vdouble V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;

  return ;
}


void vector_print_vcsimple (vcsimple V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f +i%f ", V[i].real,V[i].imaginary) ;
  printf ("\n") ;

  return ;
}

void vector_print_vcdouble (vcdouble V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f +i%f ", V[i].real,V[i].imaginary) ;
  printf ("\n") ;

  return ;
}

//=========================================================================================================================================================//


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

  start = _rdtsc () ;
     cblas_sswap (VECSIZE, vec1, 1, vec2, 1) ;
  end = _rdtsc () ;

  printf ("cblas_sswap nombre de cycles cblas: %Ld \n", end-start-residu) ;



  start = _rdtsc () ;
     mncblas_sswap (VECSIZE, vec1, 1, vec2, 1) ;
  end = _rdtsc () ;


  /*printf ("Vector 2:\n") ;
  vector_print (vec2) ;*/


  printf ("mncblas_sswap: nombre de cycles: %Ld \n", end-start-residu) ;

  start = _rdtsc () ;
     cblas_sswap (VECSIZE, vec1, 1, vec2, 1) ;
  end = _rdtsc () ;

  printf ("cblas_sswap nombre de cycles cblas: %Ld \n", end-start-residu) ;

  /*printf("Vector 2 float :\n");
  vector_print(vec2);*/
//============================================================================================================================//


//====================================================vecteur double===========================================================//
  vector_init_double(vecd1,2.0);
  vector_init_double(vecd2,3.0);
  start = _rdtsc () ;
     mncblas_dswap (VECSIZE, vecd1, 1, vecd2, 1) ;
  end = _rdtsc () ;

  printf ("mncblas_dswapy: nombre de cycles: %Ld \n", end-start-residu) ;

  /*printf("Vector 2 double\n");
  vector_print_double(vecd2);
  vector_print_double(vecd1);*/
//============================================================================================================================//

//====================================================vecteur complex_simple===========================================================//
  struct complex_simple x;
  x.real = 2.0;
  x.imaginary = 3.0;
  struct complex_simple w;
  w.real = 4.0;
  w.imaginary = 5.0;
  vector_init_csimple(veccs1,x);
  vector_init_csimple(veccs2,w);

  start = _rdtsc () ;
     mncblas_cswap (VECSIZE, veccs1, 1, veccs2, 1) ;
  end = _rdtsc () ;

  printf ("mncblas_cswap: nombre de cycles: %Ld \n", end-start-residu) ;

  printf("Vector 2 complex_simple\n");
  vector_print_vcsimple(veccs2);
  vector_print_vcsimple(veccs1);
//=====================================================================================================================================//


//====================================================vecteur complex_double===========================================================//
  struct complex_double y;
  y.real = 4.0;
  y.imaginary = 5.0;
  vector_init_cdouble(veccd1,y);

  start = _rdtsc () ;
     mncblas_zswap (VECSIZE, veccd1, 1, veccd2, 1) ;
  end = _rdtsc () ;

  printf ("mncblas_zswap: nombre de cycles: %Ld \n", end-start-residu) ;

  /*printf("Vector 2 complex_double\n");
  vector_print_vcdouble(veccd2);*/
//=====================================================================================================================================//


}
