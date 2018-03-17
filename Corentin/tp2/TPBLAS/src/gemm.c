#include "mnblas.h"

void mncblas_sgemm (MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const float alpha, const float *A,
                 const int lda, const float *B, const int ldb,
                 const float beta, float *C, const int ldc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register unsigned int k = 0 ;

  float* tmp = malloc(M * M * sizeof(float));
  for (; j < M*M ; j += incY) {
    tmp[j] = C[j];
  }
  float res;
  for (; i < M * M ; i += incX) {
    res = 0.;
    for (j=0; j < M ; j += incY) {
      res += alpha * A[j+M*(i/M)] * B[M*j+i%M)] + beta * tmp[i];
    }
    C[i] = res;
  }

  free(tmp);
  return;
}

void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register unsigned int k = 0 ;

  double* tmp = malloc(M * M * sizeof(double));
  for (; j < M*M ; j += incY) {
    tmp[j] = C[j];
  }
  double res;
  for (; i < M * M ; i += incX) {
    res = 0.;
    for (j=0; j < M ; j += incY) {
      res += alpha * A[j+M*(i/M)] * B[M*j+i%M)] + beta * tmp[i];
    }
    C[i] = res;
  }

  free(tmp);
  return;
}

void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc)
{

}

void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
                 MNCBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const void *alpha, const void *A,
                 const int lda, const void *B, const int ldb,
                 const void *beta, void *C, const int ldc)
{

}
