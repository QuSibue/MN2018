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
vfloat vecX,blvecX,vecY,blvecY;

mdouble mvecdA,blmvecdA
mcsimple vecdX,vecdY,blvecdX,blvecdY


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
  vector_Minit (mvecdA, 1.0);
  vector_init ( vecdX,2.0);
	vector_init ( vecdY,3.0);

  vector_Minit (blmvecdA, 1.0) ;
  vector_init ( blvecdX,2.0);
	vector_init ( blvecdY,3.0);


printf("=========================VECTEUR FLOAT================================\n");

  start = _rdtsc () ;
     cblas_sgemv (MNCblasRowMajor,MNCblasNoTrans,M,N,1.0,blmvecA,M,blvecX,1,2.0,blvecY,1) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,2*VECSIZE,end-start-residu);
  printf ("cblas_sgemm nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");
  //vector_print(blvecY);

  start = _rdtsc () ;
     mncblas_sgemv(MNCblasRowMajor,MNCblasNoTrans,M,N,1.0,mvecA,M,vecX,1,2.0,vecY,1);
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,31*M*N + 61 *M*N,end-start-residu);
  printf ("mncblas_sgemm: nombre de cycles: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflop : %f \n",m_Flops) ;
  printf("\n");
  vector_print(vecY);

  vector_init (blvecY, 3.0) ;
  start = _rdtsc () ;
     cblas_sgemv (MNCblasRowMajor,MNCblasNoTrans,M,N,1.0,blmvecA,M,blvecX,1,2.0,blvecY,1) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,31*M*N + 61 *M*N,end-start-residu);
  printf ("cblas_sgemm nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");
  //vector_print(blvecY);

  if(comparaisonVecteurFloat(vecY,N,blvecY,N)){
    printf ("Résultats entre cblas et mnblas identiques\n") ;
  }
  else{
    printf ("Erreurs ! Résultats entre cblas et mnblas différents\n") ;
  }


	

//====================================================vecteur double===========================================================//
  vector_Minit_double (mvecA, 1.0);
  vector_init_double ( vecX,2.0);
	vector_init_double ( vecY,3.0);

  vector_Minit_double (blmvecA, 1.0) ;
  vector_init_double ( blvecX,2.0);
	vector_init_double ( blvecY,3.0);


printf("=========================VECTEUR FLOAT================================\n");

  start = _rdtsc () ;
     cblas_sgemv (MNCblasRowMajor,MNCblasNoTrans,M,N,1.0,blmvecdA,M,blvecdX,1,2.0,blvecdY,1) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,2*VECSIZE,end-start-residu);
  printf ("cblas_sgemm nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");

  start = _rdtsc () ;
     mncblas_sgemv(MNCblasRowMajor,MNCblasNoTrans,M,N,1.0,mvecdA,M,vecdX,1,2.0,vecdY,1);
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,31*M*N + 61 *M*N,end-start-residu);
  printf ("mncblas_sgemm: nombre de cycles: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflop : %f \n",m_Flops) ;
  printf("\n");

  vector_init_double (blvecdY, 3.0) ;
  start = _rdtsc () ;
     cblas_sgemv (MNCblasRowMajor,MNCblasNoTrans,M,N,1.0,blmvecdA,M,blvecdX,1,2.0,blvecdY,1) ;
  end = _rdtsc () ;

  m_Flops=FLOPS(1,3.4,31*M*N + 61 *M*N,end-start-residu);
  printf ("cblas_sgemm nombre de cycles cblas: %Ld \n", end-start-residu) ;
	printf ("resultat en Gflops : %f\n",m_Flops) ;
  printf("\n");
  //vector_print(blvecY);

  if(comparaisonVecteurDouble(vecdY,N,blvecdY,N)){
    printf ("Résultats entre cblas et mnblas identiques\n") ;
  }
  else{
    printf ("Erreurs ! Résultats entre cblas et mnblas différents\n") ;
  }




}
