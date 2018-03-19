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


vfloat vec1,blvec1,vec2,blvec2;
float resultatf,resultatcs;
vdouble vecd1,blvecd1,vecd2,blvecd2;
double resultatd,resultatcd;
vcsimple veccs1,veccs2,blveccs1,blveccs2;
vcdouble veccd1,veccd2,blveccd1,blveccd2;
double m_Flops;


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
  vector_init (vec2, 2.0) ;

printf("=========================VECTEUR FLOAT================================\n");

  start = _rdtsc () ;
     resultatf=cblas_sdot (VECSIZE, vec1, 1,vec2,1) ;
  end = _rdtsc () ;

  printf ("cblas_sdot nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat : %f\n",resultatf);
  m_Flops=FLOPS(1,3.4,2*VECSIZE,end-start-residu);
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");


  start = _rdtsc () ;
     resultatf=mncblas_sdot (VECSIZE, vec1, 1,vec2,1) ;
  end = _rdtsc () ;

  printf ("mncblas_sdot: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f \n",resultatf) ;
  m_Flops=FLOPS(1,3.4,2*VECSIZE,end-start-residu);
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  vector_init (vec1, 1.0) ;
  vector_init (vec2, 2.0) ;
  start = _rdtsc () ;
     resultatf=cblas_sdot (VECSIZE, vec1, 1,vec2,1) ;
  end = _rdtsc () ;

  printf ("cblas_sdot nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat : %f\n",resultatf);
  m_Flops=FLOPS(1,3.4,2*VECSIZE,end-start-residu);
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

printf("======================================================================\n\n");
  /*printf("Vector 2 float :\n");
  vector_print(vec2);*/
//============================================================================================================================//


//====================================================vecteur double===========================================================//
  vector_init_double(vecd1,2.0);
  vector_init_double(vecd2,5.0);

printf("=========================VECTEUR DOUBLE================================\n");

  start = _rdtsc () ;
     resultatd = cblas_ddot (VECSIZE, vecd1, 1,vecd2,1) ;
  end = _rdtsc () ;

  printf ("cblas_ddot: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f\n",resultatd);
  m_Flops=FLOPS(1,3.4,2*VECSIZE,end-start-residu);
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  start = _rdtsc () ;
     resultatd = mncblas_ddot (VECSIZE, vecd1, 1,vecd2,1) ;
  end = _rdtsc () ;

  printf ("mncblas_ddot: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f\n",resultatd);
  m_Flops=FLOPS(1,3.4,2*VECSIZE,end-start-residu);
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");


  start = _rdtsc () ;
     resultatd = cblas_ddot (VECSIZE, vecd1, 1,vecd2,1) ;
  end = _rdtsc () ;

  printf ("cblas_ddot: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f\n",resultatd);
  m_Flops=FLOPS(1,3.4,2*VECSIZE,end-start-residu);
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  printf("======================================================================\n\n");

  /*printf("Vector 2 double\n");
  vector_print_double(vecd2);
  vector_print_double(vecd1);*/
//============================================================================================================================//

//====================================================vecteur complex_simple===========================================================//
  struct complex_simple x;
  x.real = 2.0;
  x.imaginary = 3.0;
  vector_init_csimple(veccs1,x);

  struct complex_simple x2;
  x2.real = 1.0;
  x2.imaginary = 5.0;
  vector_init_csimple(veccs2,x2);

  struct complex_simple x3;
  x3.real = 0.0;
  x3.imaginary = 0.0;


printf("=========================VECTEUR COMPLEXE SIMPLE================================\n");

  start = _rdtsc () ;
    cblas_cdotu_sub (VECSIZE, veccs1, 1,veccs2,1,&x3) ;
  end = _rdtsc () ;

  printf ("cblas_cdotu_sub: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f +i%f\n",x3.real,x3.imaginary);
  m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  x3.real = 0.0;
  x3.imaginary = 0.0;

  start = _rdtsc () ;
     mncblas_cdotu_sub (VECSIZE, veccs1, 1,veccs2,1,&x3);
  end = _rdtsc () ;

  printf ("mncblas_cdotu_sub: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f +i%f\n",x3.real,x3.imaginary);
  m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");


  x3.real = 0.0;
  x3.imaginary = 0.0;
  start = _rdtsc () ;
    cblas_cdotu_sub (VECSIZE, veccs1, 1,veccs2,1,&x3) ;
  end = _rdtsc () ;

  printf ("cblas_cdotu_sub: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f +i%f\n",x3.real,x3.imaginary);
  m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  printf("======================================================================================\n\n");


//=====================================================================================================================================//


//====================================================vecteur complex_simple===========================================================//
  struct complex_simple y;
  y.real = 2.0;
  y.imaginary = 3.0;
  vector_init_csimple(veccs1,y);

  struct complex_simple y2;
  y2.real = 1.0;
  y2.imaginary = 4.0;
  vector_init_csimple(veccs2,y2);

  struct complex_simple y3;
  y3.real = 0.0;
  y3.imaginary = 0.0;

printf("=========================VECTEUR COMPLEXE SIMPLE================================\n");
  start = _rdtsc () ;
     cblas_cdotc_sub(VECSIZE, veccs1, 1,veccs2,1,&y3) ;
  end = _rdtsc () ;

  printf ("cblas_cdotc_sub: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f +i%f\n",y3.real,y3.imaginary);
  m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  y3.real = 0.0;
  y3.imaginary = 0.0;

  start = _rdtsc () ;
     mncblas_cdotc_sub (VECSIZE, veccs1, 1,veccs2,1,&y3) ;
  end = _rdtsc () ;

  printf ("mncblas_cdotc_sub: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f +i%f\n",y3.real,y3.imaginary);
  m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");



  y3.real = 0.0;
  y3.imaginary = 0.0;
  start = _rdtsc () ;
     cblas_cdotc_sub(VECSIZE, veccs1, 1,veccs2,1,&y3) ;
  end = _rdtsc () ;

  printf ("cblas_cdotc_sub: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat : %f +i%f\n",y3.real,y3.imaginary);
  m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  printf("======================================================================================\n\n");


//=====================================================================================================================================//


//====================================================vecteur complex_double===========================================================//
  struct complex_double z;
  z.real = 4.0;
  z.imaginary = 5.0;
  vector_init_cdouble(veccd1,z);

  struct complex_double z2;
  z2.real = 4.0;
  z2.imaginary = 5.0;
  vector_init_cdouble(veccd2,z2);

  struct complex_double z3;
  z3.real = 0.0;
  z3.imaginary = 0.0;


printf("=========================VECTEUR COMPLEXE DOUBLE================================\n");

  start = _rdtsc () ;
     cblas_zdotu_sub (VECSIZE, veccs1, 1,veccs2,1,&z3) ;
  end = _rdtsc () ;

  printf ("cblas_zdotu_sub: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat %f +i%f\n",z3.real,z3.imaginary);
  m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  z3.real = 0.0;
  z3.imaginary = 0.0;

  start = _rdtsc () ;
     mncblas_zdotu_sub (VECSIZE, veccs1, 1,veccs2,1,&z3) ;
  end = _rdtsc () ;

  printf ("mncblas_zdotu_sub: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat %f +i%f\n",z3.real,z3.imaginary);
  m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  z3.real = 0.0;
  z3.imaginary = 0.0;

  start = _rdtsc () ;
     cblas_zdotu_sub (VECSIZE, veccs1, 1,veccs2,1,&z3) ;
  end = _rdtsc () ;

  printf ("cblas_zdotu_sub: nombre de cycles: %Ld \n", end-start-residu) ;
  printf ("resultat %f +i%f\n",z3.real,z3.imaginary);
  m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

printf("=============================================================================\n\n");
//=====================================================================================================================================//
struct complex_double z4;
z.real = 3.0;
z.imaginary = 7.0;
vector_init_cdouble(veccd1,z4);

struct complex_double z5;
z2.real = 2.0;
z2.imaginary = 1.0;
vector_init_cdouble(veccd2,z5);

struct complex_double z6;
z3.real = 0.0;
z3.imaginary = 0.0;


printf("=========================VECTEUR COMPLEXE DOUBLE================================\n");

start = _rdtsc () ;
   cblas_zdotc_sub (VECSIZE, veccs1, 1,veccs2,1,&z6) ;
end = _rdtsc () ;

printf ("cblas_zdotc_sub: nombre de cycles: %Ld \n", end-start-residu) ;
printf ("resultat %f +i%f\n",z6.real,z6.imaginary);
m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
printf ("resultat en Gflops : %f\n",m_Flops) ;
printf("\n");


z3.real = 0.0;
z3.imaginary = 0.0;
start = _rdtsc () ;
   mncblas_zdotc_sub (VECSIZE, veccs1, 1,veccs2,1,&z6) ;
end = _rdtsc () ;

printf ("mncblas_zdotc_sub: nombre de cycles: %Ld \n", end-start-residu) ;
printf ("resultat %f +i%f\n",z6.real,z6.imaginary);
m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
printf ("resultat en Gflops : %f\n",m_Flops) ;
printf("\n");

z3.real = 0.0;
z3.imaginary = 0.0;
start = _rdtsc () ;
   cblas_zdotc_sub (VECSIZE, veccs1, 1,veccs2,1,&z6) ;
end = _rdtsc () ;

printf ("cblas_zdotc_sub: nombre de cycles: %Ld \n", end-start-residu) ;
printf ("resultat %f +i%f\n",z6.real,z6.imaginary);
m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
printf ("resultat en Gflops : %f\n",m_Flops) ;
printf("\n");

printf("=============================================================================\n\n");


}
