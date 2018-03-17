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
vfloat vec1, vec2,blvec1,blvec2 ;
vdouble vecd1,vecd2,blvecd1,blvecd2;
vcsimple veccs1, veccs2,blveccs1,blveccs2;
vcdouble veccd1,veccd2,blveccd1,blveccd2;
float m_debit;

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
  vector_init (blvec1, 1.0) ;


  printf("=========================VECTEUR FLOAT================================\n");

  start = _rdtsc () ;
     cblas_sswap (VECSIZE, vec1, 1, vec2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("cblas_sswap nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en DEBIT : %f\n",m_debit) ;
  printf("\n");



  start = _rdtsc () ;
     mncblas_sswap (VECSIZE, blvec1, 1, blvec2, 1) ;
  end = _rdtsc () ;


  m_debit=DEBIT(3.4,end-start-residu);
  printf ("mncblas_sswap: nombre de cycles: %Ld \n", end-start-residu) ;
	printf ("resultat en DEBIT : %f \n",m_debit) ;
  printf("\n");

  vector_init (vec1, 1.0) ;

  start = _rdtsc () ;
     cblas_sswap (VECSIZE, vec1, 1, vec2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("cblas_sswap nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en DEBIT : %f\n",m_debit) ;
  printf("\n");

  if(comparaisonVecteurFloat(vec2,VECSIZE,blvec2,VECSIZE)){
    printf ("Résultats entre cblas et mnblas identiques\n") ;
  }
  else{
    printf ("Erreurs ! Résultats entre cblas et mnblas différents\n") ;
  }


  printf("======================================================================\n\n");

  /*printf("Vector 2 float :\n");
  vector_print(vec2);*/
//============================================================================================================================//


//====================================================vecteur double===========================================================//
  vector_init_double(vecd1,2.0);
  vector_init_double(blvecd1,2.0);
  vector_init_double(vecd2,3.0);
  vector_init_double(blvecd2,3.0);

  start = _rdtsc () ;
     cblas_dswap (VECSIZE, vecd1, 1, vecd2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("cblas_dswap nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en debit : %f\n",m_debit) ;
  printf("\n");

  start = _rdtsc () ;
     mncblas_dswap (VECSIZE, blvecd1, 1, blvecd2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("mncblas_dswap nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en debit : %f\n",m_debit) ;
  printf("\n");

  vector_init_double(vecd1,2.0);
  vector_init_double(vecd2,3.0);

  start = _rdtsc () ;
     cblas_dswap (VECSIZE, vecd1, 1, vecd2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("cblas_dswap nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en debit : %f\n",m_debit) ;
  printf("\n");

  if(comparaisonVecteurDouble(vecd2,VECSIZE,blvecd2,VECSIZE)){
    printf ("Résultats entre cblas et mnblas identiques\n") ;
  }
  else{
    printf ("Erreurs ! Résultats entre cblas et mnblas différents\n") ;
  }



  printf("======================================================================\n\n");

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
  vector_init_csimple(blveccs1,x);
  vector_init_csimple(veccs2,w);
  vector_init_csimple(blveccs2,w);

  start = _rdtsc () ;
     cblas_cswap (VECSIZE, veccs1, 1, veccs2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("cblas_cswap nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en debit : %f\n",m_debit) ;
  printf("\n");

  start = _rdtsc () ;
     mncblas_cswap (VECSIZE, blveccs1, 1, blveccs2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("mncblas_cswap nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en debit : %f\n",m_debit) ;
  printf("\n");



  vector_init_csimple(veccs1,x);
  vector_init_csimple(veccs2,w);

  start = _rdtsc () ;
     cblas_cswap (VECSIZE, veccs1, 1, veccs2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("cblas_cswap nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en debit : %f\n",m_debit) ;
  printf("\n");


  if(comparaisonVecteurCS(veccs2,VECSIZE,blveccs2,VECSIZE)){
    printf ("Résultats entre cblas et mnblas identiques\n") ;
  }
  else{
    printf ("Erreurs ! Résultats entre cblas et mnblas différents\n") ;
  }



  printf("=============================================================================\n\n");
  //printf("Vector 2 complex_simple\n");
  //vector_print_vcsimple(veccs2);
  //vector_print_vcsimple(veccs1);
//=====================================================================================================================================//


//====================================================vecteur complex_double===========================================================//
  struct complex_double y;
  y.real = 4.0;
  y.imaginary = 5.0;
  vector_init_cdouble(veccd1,y);
  vector_init_cdouble(blveccd1,y);

  struct complex_double y2;
  y2.real = 3.0;
  y2.imaginary = 9.0;
  vector_init_cdouble(veccd2,y2);
  vector_init_cdouble(blveccd2,y2);


  start = _rdtsc () ;
     cblas_zswap (VECSIZE, veccd1, 1, veccd2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("cblas_zswap nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en debit : %f\n",m_debit) ;
  printf("\n");


  start = _rdtsc () ;
     mncblas_zswap (VECSIZE, blveccd1, 1, blveccd2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("mncblas_zswap nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en debit : %f\n",m_debit) ;
  printf("\n");


  vector_init_cdouble(veccd1,y);
  vector_init_cdouble(veccd2,y2);
  start = _rdtsc () ;
     cblas_zswap (VECSIZE, veccd1, 1, veccd2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("cblas_zswap nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en debit : %f\n",m_debit) ;
  printf("\n");

  if(comparaisonVecteurCD(veccd2,VECSIZE,blveccd2,VECSIZE)){
    printf ("Résultats entre cblas et mnblas identiques\n") ;
  }
  else{
    printf ("Erreurs ! Résultats entre cblas et mnblas différents\n") ;
  }


  printf("=============================================================================\n\n");


  /*printf("Vector 2 complex_double\n");
  vector_print_vcdouble(veccd2);*/
//=====================================================================================================================================//


}
