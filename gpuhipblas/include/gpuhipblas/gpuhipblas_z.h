/**
 *
 * @file gpuhipblas_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon GPU CHAMELEON_Complex64_t kernels header
 *
 * @version 1.2.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Loris Lucido
 * @date 2023-01-30
 * @precisions normal z -> c d s
 *
 */
#ifndef _gpuhipblas_z_h_
#define _gpuhipblas_z_h_

/**
 *  Declarations of hip kernels - alphabetical order
 */
int HIP_zgemm(  cham_trans_t transa, cham_trans_t transb, int m, int n, int k, const hipblasDoubleComplex *alpha, const hipblasDoubleComplex *A, int lda, const hipblasDoubleComplex *B, int ldb, const hipblasDoubleComplex *beta, hipblasDoubleComplex *C, int ldc, hipblasHandle_t handle );
int HIP_zhemm(  cham_side_t side, cham_uplo_t uplo, int m, int n, const hipblasDoubleComplex *alpha, const hipblasDoubleComplex *A, int lda, const hipblasDoubleComplex *B, int ldb, const hipblasDoubleComplex *beta, hipblasDoubleComplex *C, int ldc, hipblasHandle_t handle );
int HIP_zher2k( cham_uplo_t uplo, cham_trans_t trans, int n, int k, const hipblasDoubleComplex *alpha, const hipblasDoubleComplex *A, int lda, const hipblasDoubleComplex *B, int ldb, const double *beta, hipblasDoubleComplex *C, int ldc, hipblasHandle_t handle );
int HIP_zherk(  cham_uplo_t uplo, cham_trans_t trans, int n, int k, const double *alpha, const hipblasDoubleComplex *A, int lda, const double *beta, hipblasDoubleComplex *B, int ldb, hipblasHandle_t handle );
int HIP_zsymm(  cham_side_t side, cham_uplo_t uplo, int m, int n, const hipblasDoubleComplex *alpha, const hipblasDoubleComplex *A, int lda, const hipblasDoubleComplex *B, int ldb, const hipblasDoubleComplex *beta, hipblasDoubleComplex *C, int ldc, hipblasHandle_t handle );
int HIP_zsyr2k( cham_uplo_t uplo, cham_trans_t trans, int n, int k, const hipblasDoubleComplex *alpha, const hipblasDoubleComplex *A, int lda, const hipblasDoubleComplex *B, int ldb, const hipblasDoubleComplex *beta, hipblasDoubleComplex *C, int ldc, hipblasHandle_t handle );
int HIP_zsyrk(  cham_uplo_t uplo, cham_trans_t trans, int n, int k, const hipblasDoubleComplex *alpha, const hipblasDoubleComplex *A, int lda, const hipblasDoubleComplex *beta, hipblasDoubleComplex *C, int ldc, hipblasHandle_t handle );
int HIP_ztrmm(  cham_side_t side, cham_uplo_t uplo, cham_trans_t transa, cham_diag_t diag, int m, int n, const hipblasDoubleComplex *alpha, const hipblasDoubleComplex *A, int lda, hipblasDoubleComplex *B, int ldb, hipblasHandle_t handle );
int HIP_ztrsm(  cham_side_t side, cham_uplo_t uplo, cham_trans_t transa, cham_diag_t diag, int m, int n, const hipblasDoubleComplex *alpha, const hipblasDoubleComplex *A, int lda, hipblasDoubleComplex *B, int ldb, hipblasHandle_t handle );
int HIP_zlaswp_get( int m0, int m, int n, int k, const CHAMELEON_Complex64_t *A, int lda, CHAMELEON_Complex64_t       *B, int ldb, const int *perm, hipblasHandle_t handle );
int HIP_zlaswp_set( int m0, int m, int n, int k, const CHAMELEON_Complex64_t *A, int lda, CHAMELEON_Complex64_t *B, int ldb, const int *invp, hipblasHandle_t handle );
int hip_zlaswp_get_zcopy( int m0, int m, int n, int k, const CHAMELEON_Complex64_t *A, int lda, CHAMELEON_Complex64_t       *B, int ldb, const int *perm, hipblasHandle_t handle );
int HIP_zlaswp_set_zcopy( int m0, int m, int n, int k, const CHAMELEON_Complex64_t *A, int lda, CHAMELEON_Complex64_t *B, int ldb, const int *invp, hipblasHandle_t handle );
int HIP_zlacpy(cham_uplo_t uplo, int M, int N,
               const CHAMELEON_Complex64_t *A, int LDA,
               CHAMELEON_Complex64_t *B, int LDB,
               hipblasHandle_t handle );

#endif /* _gpuhipblas_z_h_ */
