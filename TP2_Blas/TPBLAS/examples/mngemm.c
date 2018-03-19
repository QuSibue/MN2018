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

mfloat vecA,blvecA,vecB,blvecB,vecC,blvecC;
mdouble vecdA,blvecdA,vecdB,blvecdB,vecdC,blvecdC;
mcsimple veccsA,veccsB,blveccsA,blveccsB,blveccsC,veccsC;
struct complex_simple alphacs,betacs;
mcdouble veccd1,veccd2,blveccd1,blveccd2;
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
  vector_Minit (vecA, 1.0) ;
  vector_Minit (vecB, 3.0) ;
  vector_Minit (vecC, 2.0) ;

  vector_Minit (blvecA, 1.0) ;
  vector_Minit (blvecB, 3.0) ;
  vector_Minit (blvecC, 2.0) ;


printf("=========================VECTEUR FLOAT================================\n");

  start = _rdtsc () ;
     cblas_sgemm (MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M,N,K,2.0, vecA,M,vecB,N,1.0,vecC,M) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,2*VECSIZE,end-start-residu);
  printf ("cblas_sgemm nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");
  //vector_print(vecC);

  start = _rdtsc () ;
     mncblas_sgemm (MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M,N,K,2.0, blvecA,M,blvecB,N,1.0,blvecC,M);
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,31*M*N + 61 *M*N,end-start-residu);
  printf ("mncblas_sgemm: nombre de cycles: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflop : %f \n",m_Flops) ;
  printf("\n");
  //vector_print(blvecC);

  vector_Minit (vecC, 2.0) ;
  start = _rdtsc () ;
     cblas_sgemm (MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M,M,K,2.0, vecA,M,vecB,N,1.0,vecC,N) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,31*M*N + 61 *M*N,end-start-residu);
  printf ("cblas_sgemm nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  if(comparaisonVecteurFloat(vecC,M*N,blvecC,M*N)){
    printf ("Résultats entre cblas et mnblas identiques\n") ;
  }
  else{
    printf ("Erreurs ! Résultats entre cblas et mnblas différents\n") ;
  }



  printf("======================================================================\n\n");

  vector_Minit_double (vecdA, 1.0) ;
  vector_Minit_double (vecdB, 3.0) ;
  vector_Minit_double (vecdC, 2.0) ;

  vector_Minit_double (blvecdA, 1.0) ;
  vector_Minit_double (blvecdB, 3.0) ;
  vector_Minit_double (blvecdC, 2.0) ;


printf("=========================VECTEUR DOUBLE================================\n");

  start = _rdtsc () ;
     cblas_dgemm (MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M,N,K,2.0, vecdA,M,vecdB,N,1.0,vecdC,M) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,31*M*N + 61 *M*N,end-start-residu);
  printf ("cblas_sgemm nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  start = _rdtsc () ;
     mncblas_dgemm (MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M,N,K,2.0, blvecdA,M,blvecdB,N,1.0,blvecdC,M);
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,31*M*N + 61 *M*N,end-start-residu);
  printf ("mncblas_sgemm: nombre de cycles: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflop : %f \n",m_Flops) ;
  printf("\n");
  //vector_print(vec2);

  vector_Minit_double (vecdC, 2.0) ;
  start = _rdtsc () ;
     cblas_dgemm (MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M,M,K,2.0, vecdA,M,vecdB,N,1.0,vecdC,N) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,31*M*N + 61 *M*N,end-start-residu);
  printf ("cblas_sgemm nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  if(comparaisonVecteurDouble(vecdC,M*N,blvecdC,M*N)){
    printf ("Résultats entre cblas et mnblas identiques\n") ;
  }
  else{
    printf ("Erreurs ! Résultats entre cblas et mnblas différents\n") ;
  }



  printf("======================================================================\n\n");

  printf("======================================================================\n\n");

  struct complex_simple a;a.real=1.0;a.imaginary=2.0;
  struct complex_simple b;b.real=2.0;b.imaginary=4.0;
  struct complex_simple c;c.real=3.0;b.imaginary=6.0;




  vector_Minit_csimple (veccsA, a) ;
  vector_Minit_csimple (veccsB, b) ;
  vector_Minit_csimple (veccsC, c) ;

  vector_Minit_csimple (blveccsA, a) ;
  vector_Minit_csimple (blveccsB, b) ;
  vector_Minit_csimple (blveccsC, c) ;

  alphacs.real=2.0;
  alphacs.imaginary=1.0;
  betacs.real=1.0;
  betacs.imaginary=3.0;

  printf("=========================VECTEUR COMPLEX_SIMPLE================================\n");

    start = _rdtsc () ;
       cblas_cgemm (MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M,N,K,&alphacs, veccsA,M,veccsB,N,&betacs,veccsC,M) ;
    end = _rdtsc () ;

    m_Flops=FLOPS(1,3.4,31*M*N + 61 *M*N,end-start-residu);
    printf ("cblas_sgemm nombre de cycles cblas: %Ld \n", end-start-residu) ;
  	printf ("resultat en Gflops : %f\n",m_Flops) ;
    printf("\n");

    start = _rdtsc () ;
       mncblas_cgemm (MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M,N,K,&alphacs, blveccsA,M,blveccsB,N,&betacs,blveccsC,M);
    end = _rdtsc () ;

    m_Flops=FLOPS(1,3.4,31*M*N + 61 *M*N,end-start-residu);
    printf ("mncblas_sgemm: nombre de cycles: %Ld \n", end-start-residu) ;
  	printf ("resultat en Gflop : %f \n",m_Flops) ;
    printf("\n");
    //vector_print(vec2);

    vector_Minit_csimple (veccsC, c) ;
    start = _rdtsc () ;
       cblas_cgemm (MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M,M,K,&alphacs, veccsA,M,veccsB,N,&betacs,veccsC,N) ;
    end = _rdtsc () ;

    m_Flops=FLOPS(1,3.4,31*M*N + 61 *M*N,end-start-residu);
    printf ("cblas_sgemm nombre de cycles cblas: %Ld \n", end-start-residu) ;
  	printf ("resultat en Gflops : %f\n",m_Flops) ;
    printf("\n");

    if(comparaisonVecteurCS(veccsC,M*N,blveccsC,M*N)){
      printf ("Résultats entre cblas et mnblas identiques\n") ;
    }
    else{
      printf ("Erreurs ! Résultats entre cblas et mnblas différents\n") ;
    }



    printf("======================================================================\n\n");



}

//============================================================================================================================//
