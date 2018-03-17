#include <stdio.h>
#include "fonctions_test.h"



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
double DEBIT(double frequence,long long cycles){
  return ((sizeof(float)*VECSIZE)/(cycles))*frequence;
}

double FLOPS(int coeurs,double frequence,int flop,long long cycles){
	return coeurs*frequence*(flop/cycles);
}

int compare_complex_simple(struct complex_simple *c1,struct complex_simple *c2){
	return ( c1->real == c2->real && c1->imaginary == c2->imaginary );
}

int compare_complex_double(struct complex_double *c1,struct complex_double *c2){
	return ( c1->real == c2->real && c1->imaginary == c2->imaginary );
}


int comparaisonVecteurFloat(float *v1,int taillev1,float *v2,int taillev2){
	if (taillev1 != taillev2){
		return 0;
	}else{
		int i=0;
		while (i<taillev1 && v1[i] == v2[i]){
			i++;
		}
		return i>=taillev1;
	}

}

int comparaisonVecteurDouble(double *v1,int taillev1,double *v2,int taillev2){
	if (taillev1 != taillev2){
		return 0;
	}else{
		int i=0;
		while (i<taillev1 && v1[i] == v2[i]){
			i++;
		}
		return i>=taillev1;
	}

}

int comparaisonVecteurCS(struct complex_simple *v1,int taillev1,struct complex_simple *v2,int taillev2){
	if (taillev1 != taillev2){
		return 0;
	}else{
		int i=0;
		while (i<taillev1 && compare_complex_simple(&v1[i],&v2[i])){
			i++;
		}
		return i>=taillev1;
	}

}

int comparaisonVecteurCD(struct complex_double *v1,int taillev1,struct complex_double *v2,int taillev2){
	if (taillev1 != taillev2){
		return 0;
	}else{
		int i=0;
		while (i<taillev1 && compare_complex_double(&v1[i],&v2[i])){
			i++;
		}
		return i>=taillev1;
	}

}
