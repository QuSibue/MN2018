#include <stdio.h>
#include <cblas.h>

#include "mnblas.h"
#include "complex.h"
#include "fonctions_test.h"
/*
  Mesure des cycles
*/

#include <x86intrin.h>





//===================================================DEFINITION=============================================================================================//

typedef float vfloat[VECSIZE] ;
typedef double vdouble[VECSIZE];
typedef struct complex_simple vcsimple[VECSIZE];
typedef struct complex_double vcdouble[VECSIZE];


vfloat vec1;
float resultatf,resultatcs;
vdouble vecd1;
double resultatd,resultatcd;
vcsimple veccs1;
vcdouble veccd1;



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

  start = _rdtsc () ;
     resultatf=cblas_sasum (VECSIZE, vec1, 1) ;
  end = _rdtsc () ;

  printf ("cblas_sswap nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat : %f\n",resultatf);


  start = _rdtsc () ;
     resultatf=mncblas_sasum (VECSIZE, vec1, 1) ;
  end = _rdtsc () ;


  /*printf ("Vector 2:\n") ;
  vector_print (vec2) ;*/


  printf ("mncblas_sswap: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f \n",resultatf) ;

  /*printf("Vector 2 float :\n");
  vector_print(vec2);*/
//============================================================================================================================//


//====================================================vecteur double===========================================================//
  vector_init_double(vecd1,2.0);

  start = _rdtsc () ;
     resultatd = cblas_dasum (VECSIZE, vecd1, 1) ;
  end = _rdtsc () ;

  printf ("cblas_dswapy: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f\n",resultatd);

  start = _rdtsc () ;
     resultatd = mncblas_dasum (VECSIZE, vecd1, 1) ;
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
     resultatcs = cblas_scasum (VECSIZE, veccs1, 1) ;
  end = _rdtsc () ;

  printf ("cblas_cswap: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f \n",resultatcs);

  start = _rdtsc () ;
     resultatcs = mncblas_scasum (VECSIZE, veccs1, 1) ;
  end = _rdtsc () ;

  printf ("mncblas_cswap: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f \n",resultatcs);


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
  printf ("resultat %f\n",resultatcd);

  start = _rdtsc () ;
     resultatcd = mncblas_dzasum (VECSIZE, veccd1, 1) ;
  end = _rdtsc () ;

  printf ("cblas_zswap: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat %f\n",resultatcd);
//=====================================================================================================================================//


}
