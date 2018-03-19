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

mfloat mvecA,blmvecA;
vfloat vecA,blvecA;

mdouble vecdA,blvecdA,vecdB,blvecdB,vecdC,blvecdC;
mcsimple veccsA,veccsB,blveccsA,blveccsB,blveccsC,veccsC;
struct complex_simple alphacs,betacs;
mcdouble veccdA,veccdB,blveccdA,blveccdB,blveccdC,veccdC;
struct complex_double alphacd,betacd;
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
  vector_Minit (mvecA, 1.0);
  vector_init ( vecA,2.0);

  vector_Minit (blmvecA, 1.0) ;
  vector_init ( blvecA,2.0);


printf("=========================VECTEUR FLOAT================================\n");

  start = _rdtsc () ;
     cblas_sgemv (MNCblasRowMajor,MNCblasNoTrans,M,N,alpha,A,M,X,incX,beta,Y,incY) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,2*VECSIZE,end-start-residu);
  printf ("cblas_sgemm nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");
  //vector_print(vecC);

  start = _rdtsc () ;
     mncblas_sgemv (MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M,N,K,2.0, blvecA,M,blvecB,N,1.0,blvecC,M);
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,31*M*N + 61 *M*N,end-start-residu);
  printf ("mncblas_sgemm: nombre de cycles: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflop : %f \n",m_Flops) ;
  printf("\n");
  //vector_print(blvecC);

  vector_Minit (vecC, 2.0) ;
  start = _rdtsc () ;
     cblas_sgemv (MNCblasRowMajor,MNCblasNoTrans,MNCblasNoTrans,M,M,K,2.0, vecA,M,vecB,N,1.0,vecC,N) ;
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

}
