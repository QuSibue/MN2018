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
vcdouble veccd1,veccd2;
struct complex_double resultatcd;


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
  vector_init (vec2, 2.0) ;

  start = _rdtsc () ;
     resultatf=blas_sasum (VECSIZE, vec1, 1,vec2,1) ;
  end = _rdtsc () ;

  printf ("cblas_sswap nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat : %f\n",resultatf);


  start = _rdtsc () ;
     resultatf=mncblas_sasum (VECSIZE, vec1, 1,vec2,1) ;
  end = _rdtsc () ;


  /*printf ("Vector 2:\n") ;
  vector_print (vec2) ;*/


  printf ("mncblas_sswap: nombre de cycles: %f \n", end-start-residu) ;
  printf ("resultat : %f \n",resultatf) ;

  /*printf("Vector 2 float :\n");
  vector_print(vec2);*/
//============================================================================================================================//


//====================================================vecteur double===========================================================//
  vector_init_double(vecd1,2.0);
  vector_init_double(vecd2,5.0);
  
  start = _rdtsc () ;
     resultatd = blas_dasum (VECSIZE, vecd1, 1,vecd2,1) ;
  end = _rdtsc () ;

  printf ("cblas_dswapy: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f\n",resultatd);

  start = _rdtsc () ;
     resultatd = mncblas_dasum (VECSIZE, vecd1, 1,vecd2,1) ;
  end = _rdtsc () ;

  printf ("mncblas_dswapy: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f\n",resultatd);

  /*printf("Vector 2 double\n");
  vector_print_double(vecd2);
  vector_print_double(vecd1);*/
//============================================================================================================================//

//====================================================vecteur complex_simple===========================================================//
  struct complex_simple x;
  x.real = 2.0;
  x.imaginary = 3.0;
  vector_init_csimple(veccs1,x);

  start = _rdtsc () ;
     resultatcs = cblas_casum (VECSIZE, veccs1, 1) ;
  end = _rdtsc () ;

  printf ("cblas_cswap: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f %f\n",resultatcs.real,resultat.imaginary);

  start = _rdtsc () ;
     resultatcs = mncblas_casum (VECSIZE, veccs1, 1) ;
  end = _rdtsc () ;

  printf ("mncblas_cswap: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f %f\n",resultatcs.real,resultat.imaginary);


//=====================================================================================================================================//


//====================================================vecteur complex_double===========================================================//
  struct complex_double y;
  y.real = 4.0;
  y.imaginary = 5.0;
  vector_init_cdouble(veccd1,y);

  start = _rdtsc () ;
     resultatcd = cblas_dzasum (VECSIZE, veccd1, 1) ;
  end = _rdtsc () ;

  printf ("cblas_zswap: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat %f %f\n",resultatcd.real,resultatcd.imaginary);

  start = _rdtsc () ;
     resultatcd = mncblas_dzasum (VECSIZE, veccd1, 1) ;
  end = _rdtsc () ;

  printf ("cblas_zswap: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat %f %f\n",resultatcd.real,resultatcd.imaginary);
//=====================================================================================================================================//


}
