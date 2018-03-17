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
vfloat vec1,blvec1, vec2,blvec2 ;
vdouble vecd1,blvecd1,vecd2,blvecd2;
vcsimple veccs1,blveccs1,veccs2,blveccs2;
vcdouble veccd1,blveccd1,veccd2,blveccd2;
double m_debit;

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
     cblas_scopy (VECSIZE, vec1, 1, vec2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("cblas_scopy nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en DEBIT : %f\n",m_debit) ;
  printf("\n");

  start = _rdtsc () ;
     mncblas_scopy (VECSIZE, blvec1, 1, blvec2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("mncblas_scopy: nombre de cycles: %Ld \n", end-start-residu) ;
	printf ("resultat en DEBIT : %f \n",m_debit) ;
  printf("\n");


  start = _rdtsc () ;
     cblas_scopy (VECSIZE, vec1, 1, vec2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("cblas_scopy nombre de cycles cblas: %Ld \n", end-start-residu) ;
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

printf("=========================VECTEUR Double================================\n");


  start = _rdtsc () ;
     cblas_dcopy (VECSIZE, vecd1, 1, vecd2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("cblas_dcopy nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en debit : %f\n",m_debit) ;
  printf("\n");

  start = _rdtsc () ;
     mncblas_dcopy (VECSIZE, blvecd1, 1, blvecd2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("mncblas_dcopy nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en debit : %f\n",m_debit) ;
  printf("\n");

  start = _rdtsc () ;
     cblas_dcopy (VECSIZE, vecd1, 1, vecd2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("cblas_dcopy nombre de cycles cblas: %Ld \n", end-start-residu) ;
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
  vector_print_double(vecd2);*/
//============================================================================================================================//

//====================================================vecteur complex_simple===========================================================//
  struct complex_simple x;
  x.real = 2.0;
  x.imaginary = 3.0;
  vector_init_csimple(veccs1,x);
  vector_init_csimple(blveccs1,x);

  start = _rdtsc () ;
     cblas_ccopy (VECSIZE, veccs1, 1, veccs2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("cblas_ccopy nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en debit : %f\n",m_debit) ;
  printf("\n");


  start = _rdtsc () ;
     mncblas_ccopy (VECSIZE, blveccs1, 1, blveccs2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("mncblas_ccopy nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en debit : %f\n",m_debit) ;
  printf("\n");

  start = _rdtsc () ;
     cblas_ccopy (VECSIZE, veccs1, 1, veccs2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("cblas_ccopy nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en debit : %f\n",m_debit) ;
  printf("\n");

  if(comparaisonVecteurCS(veccs2,VECSIZE,blveccs2,VECSIZE)){
    printf ("Résultats entre cblas et mnblas identiques\n") ;
  }
  else{
    printf ("Erreurs ! Résultats entre cblas et mnblas différents\n") ;
  }



  printf("=============================================================================\n\n");

  /*printf("Vector 2 complex_simple\n");
  vector_print_vcsimple(veccs2);*/
//=====================================================================================================================================//


//====================================================vecteur complex_double===========================================================//
  struct complex_double y;
  y.real = 4.0;
  y.imaginary = 5.0;
  vector_init_cdouble(veccd1,y);
  vector_init_cdouble(blveccd1,y);

  start = _rdtsc () ;
     cblas_zcopy (VECSIZE, veccd1, 1, veccd2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("cblas_zcopy nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en debit : %f\n",m_debit) ;
  printf("\n");


  start = _rdtsc () ;
     mncblas_zcopy (VECSIZE, blveccd1, 1, blveccd2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("mncblas_zcopy nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en debit : %f\n",m_debit) ;
  printf("\n");

  start = _rdtsc () ;
     cblas_zcopy (VECSIZE, veccd1, 1, veccd2, 1) ;
  end = _rdtsc () ;

  m_debit=DEBIT(3.4,end-start-residu);
  printf ("cblas_zcopy nombre de cycles cblas: %Ld \n", end-start-residu) ;
  printf ("resultat en debit : %f\n",m_debit) ;
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
