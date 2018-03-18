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

}

void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc){




}
