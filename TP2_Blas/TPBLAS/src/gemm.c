#include "mnblas.h"
#include "complex.h"


// Float //


void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc){

	register unsigned int i ;
	register unsigned int j ;
	register unsigned int g ;

	register float temp;

	for (i=0;i<M;i++){
		for (j=0;j<N;j++){
			temp =  A[i*K+0] * B[0*K+j] ;
			for (g=1;g<K;g++){
				temp += A[i*K+g] * B[g*K+j] ;
			}
			C[i*N+j] = alpha * temp + beta * C[i*N+j] ;
		}
	}


}


void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc){

	register unsigned int i ;
	register unsigned int j ;
	register unsigned int g ;

	register double temp;

	for (i=0;i<M;i++){
		for (j=0;j<N;j++){
			temp =  A[i*K+0] * B[0*K+j] ;
			for (g=1;g<K;g++){
				temp += A[i*K+g] * B[g*K+j] ;
			}
			C[i*N+j] = alpha * temp + beta * C[i*N+j] ;
		}
	}


}

void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc){

                   register unsigned int i ;
                 	register unsigned int j ;
                 	register unsigned int g ;

                 	register struct complex_simple  temp;

                 	for (i=0;i<M;i++){
                 		for (j=0;j<N;j++){
                 			temp =  multiplication_cs( ((struct complex_simple*)A)[i*K+0], ((struct complex_simple*)B)[0*K+j] );
                 			for (g=1;g<K;g++){
                 				temp = addition_cs( temp,multiplication_cs ( ((struct complex_simple*)A)[i*K+g] , ((struct complex_simple*)B)[g*K+j] ) );
                 			}
                 			((struct complex_simple*)C)[i*N+j] = addition_cs(
                                              multiplication_cs ( *((struct complex_simple*)alpha) , temp)
                                              ,
                                              multiplication_cs( *((struct complex_simple*)beta) , ((struct complex_simple*)C)[i*N+j] )
                                            );
                 		}
                 	}

}

void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc){

                   register unsigned int i ;
                 	register unsigned int j ;
                 	register unsigned int g ;

                 	register struct complex_double  temp;

                 	for (i=0;i<M;i++){
                 		for (j=0;j<N;j++){
                 			temp =  multiplication_cd( ((struct complex_double*)A)[i*K+0], ((struct complex_double*)B)[0*K+j] );
                 			for (g=1;g<K;g++){
                 				temp = addition_cd( temp,multiplication_cd ( ((struct complex_double*)A)[i*K+g] , ((struct complex_double*)B)[g*K+j]) );
                 			}
                 			((struct complex_double*)C)[i*N+j] = addition_cd(
                                              multiplication_cd ( *((struct complex_double*)alpha) , temp)
                                              ,
                                              multiplication_cd( *((struct complex_double*)beta) , ((struct complex_double*)C)[i*N+j] )
                                            );
                 		}
                 	}


}
