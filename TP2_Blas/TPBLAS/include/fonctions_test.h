#ifndef FONCTIONS_TEST_H
#define FONCTIONS_TEST_H
#include "complex.h"


#define VECSIZE  1000
#define M  1000
#define N  1000
#define K  1000



typedef float vfloat[VECSIZE] ;
typedef double vdouble[VECSIZE];
typedef struct complex_simple vcsimple[VECSIZE];
typedef struct complex_double vcdouble[VECSIZE];

typedef float mfloat[M*N] ;
typedef double mdouble[M*N];
typedef struct complex_simple mcsimple[M*N];
typedef struct complex_double mcdouble[M*N];

//=======================================================INIT=========================================================================================//
void vector_init (vfloat V, float x);

void vector_init_double (vdouble V, double x);

void vector_init_csimple (vcsimple V, struct complex_simple x);

void vector_init_cdouble (vcdouble V, struct complex_double x);


void vector_Minit (vfloat V, float x);

void vector_Minit_double (vdouble V, double x);

void vector_Minit_csimple (vcsimple V, struct complex_simple x);

void vector_Minit_cdouble (vcdouble V, struct complex_double x);


//=========================================================================================================================================================//


//===================================================PRINT==================================================================================================//

void vector_print (vfloat V);

void vector_print_double (vdouble V);

void vector_print_vcsimple (vcsimple V);

void vector_print_vcdouble (vcdouble V);


//=========================================================================================================================================================//

//=======================================================MERSURES==========================================================================================//
double DEBIT(double frequence,long long cycles);

double FLOPS(int coeurs,double frequence,int flop,long long cycles);

int compare_complex_simple(struct complex_simple *c1,struct complex_simple *c2);

int compare_complex_double(struct complex_double *c1,struct complex_double *c2);

int comparaisonVecteurFloat(float *v1,int taillev1,float *v2,int taillev2);

int comparaisonVecteurDouble(double *v1,int taillev1,double *v2,int taillev2);

int comparaisonVecteurCS(struct complex_simple *v1,int taillev1,struct complex_simple *v2,int taillev2);

int comparaisonVecteurCD(struct complex_double *v1,int taillev1,struct complex_double *v2,int taillev2);

//LMAO
int pd();

//trololo


//=========================================================================================================================================================//

#endif
