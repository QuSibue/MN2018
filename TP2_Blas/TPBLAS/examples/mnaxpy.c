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
  vector_init (vec2, 3.0) ;

  vector_init (blvec1, 1.0) ;
  vector_init (blvec2, 3.0) ;

printf("=========================VECTEUR FLOAT================================\n");

  start = _rdtsc () ;
     cblas_saxpy (VECSIZE,2.0, vec1, 1,vec2,1) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,2*VECSIZE,end-start-residu);
  printf ("cblas_saxpy nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  start = _rdtsc () ;
     mncblas_saxpy (VECSIZE,2.0, blvec1, 1,blvec2,1) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,2*VECSIZE,end-start-residu);
  printf ("mncblas_sasum: nombre de cycles: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflop : %f \n",m_Flops) ;
  printf("\n");
  //vector_print(vec2);

  vector_init (vec2, 3.0) ;
  start = _rdtsc () ;
     cblas_saxpy (VECSIZE,2.0, vec1, 1,vec2,1) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,2*VECSIZE,end-start-residu);
  printf ("cblas_saxpy nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  if(comparaisonVecteurFloat(vec2,VECSIZE,blvec2,VECSIZE)){
    printf ("Résultats entre cblas et mnblas identiques\n") ;
  }
  else{
    printf ("Erreurs ! Résultats entre cblas et mnblas différents\n") ;
  }



  printf("======================================================================\n\n");

//============================================================================================================================//


//====================================================vecteur double===========================================================//
  vector_init_double(vecd1,1.0);
  vector_init_double(vecd2,3.0);

  vector_init_double(blvecd1,1.0);
  vector_init_double(blvecd2,3.0);


printf("=========================VECTEUR Double================================\n");

  start = _rdtsc () ;
     cblas_daxpy (VECSIZE,2.0,vecd1, 1,vecd2,1) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,2*VECSIZE,end-start-residu);
  printf ("cblas_daxpy nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");



  start = _rdtsc () ;
     mncblas_daxpy (VECSIZE,2.0, blvecd1, 1,blvecd2,1) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,2*VECSIZE,end-start-residu);
  printf ("mblas_daxpy nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");


  vector_init_double(vecd2,3.0);

  start = _rdtsc () ;
     cblas_daxpy (VECSIZE,2.0,vecd1, 1,vecd2,1) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,2*VECSIZE,end-start-residu);
  printf ("cblas_daxpy nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  if(comparaisonVecteurDouble(vecd2,VECSIZE,blvecd2,VECSIZE)){
    printf ("Résultats entre cblas et mnblas identiques\n") ;
  }
  else{
    printf ("Erreurs ! Résultats entre cblas et mnblas différents\n") ;
  }



  printf("======================================================================\n\n");

//============================================================================================================================//

//====================================================vecteur complex_simple===========================================================//
  struct complex_simple x;
  x.real = 2.0;
  x.imaginary = 1.5;
  vector_init_csimple(veccs1,x);
  vector_init_csimple(blveccs1,x);

  x.real = 1.0;
  x.imaginary = 3.0;
  vector_init_csimple(veccs2,x);
  vector_init_csimple(blveccs2,x);

  struct complex_simple w;
  w.real = 1.0;
  w.imaginary = -0.5;


  start = _rdtsc () ;
     cblas_caxpy (VECSIZE,&w, veccs1, 1,veccs2,1) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
  printf ("cblas_caxpy nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");




  start = _rdtsc () ;
     mncblas_caxpy (VECSIZE,&w,blveccs1,1,blveccs2,1) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
  printf ("mncblas_caxpy nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");



  vector_init_csimple(veccs2,x);
  start = _rdtsc () ;
     cblas_caxpy (VECSIZE,&w, veccs1, 1,veccs2,1) ;
  end = _rdtsc () ;
  //vector_print_vcsimple(veccs2);

  m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
  printf ("cblas_caxpy nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  if(comparaisonVecteurCS(veccs2,VECSIZE,blveccs2,VECSIZE)){
    printf ("Résultats entre cblas et mnblas identiques\n") ;
  }
  else{
    printf ("Erreurs ! Résultats entre cblas et mnblas différents\n") ;
  }



  printf("=============================================================================\n\n");

//=====================================================================================================================================//


//====================================================vecteur complex_double===========================================================//
  struct complex_double y;
  y.real = 2.0;
  y.imaginary = 1.5;
  vector_init_cdouble(veccd1,y);
  vector_init_cdouble(blveccd1,y);

  y.real = 1.0;
  y.imaginary = 3.0;
  vector_init_cdouble(veccd2,y);
  vector_init_cdouble(blveccd2,y);

  struct complex_double z;
  z.real = -2.0;
  z.imaginary = 1.5;

printf("=========================VECTEUR COMPLEXES DOUBLES================================\n");

  start = _rdtsc () ;
     cblas_zaxpy (VECSIZE,&z, veccd1, 1,veccd2,1) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
  printf ("cblas_dzaxpy nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");
  //vector_print_vcdouble(veccd2);

  //vector_init_cdouble(veccd2,y);

  start = _rdtsc () ;
     mncblas_zaxpy (VECSIZE,&z, blveccd1, 1,blveccd2,1) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
  printf ("mncblas_dzaxpy nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");
  //vector_print_vcdouble(veccd2);

  vector_init_cdouble(veccd2,y);
  start = _rdtsc () ;
     cblas_zaxpy (VECSIZE,&z, veccd1, 1,veccd2,1) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,8*VECSIZE,end-start-residu);
  printf ("cblas_dzaxpy nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  if(comparaisonVecteurCD(veccd2,VECSIZE,blveccd2,VECSIZE)){
    printf ("Résultats entre cblas et mnblas identiques\n") ;
  }
  else{
    printf ("Erreurs ! Résultats entre cblas et mnblas différents\n") ;
  }


  printf("=============================================================================\n\n");
//=====================================================================================================================================//


}
