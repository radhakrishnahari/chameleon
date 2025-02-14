/**
 *
 * @file hip_hgemm.c
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon hip_hgemm GPU kernel
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 *
 */
#include "gpuhipblas.h"

int
HIP_hgemm( cham_trans_t transa, cham_trans_t transb,
           int m, int n, int k,
           const CHAMELEON_Real16_t *alpha,
           const CHAMELEON_Real16_t *A, int lda,
           const CHAMELEON_Real16_t *B, int ldb,
           const CHAMELEON_Real16_t *beta,
           CHAMELEON_Real16_t *C, int ldc,
           hipblasHandle_t handle )
{
    hipblasStatus_t rc;

    rc = hipblasHgemm( handle,
                       chameleon_hipblas_const(transa), chameleon_hipblas_const(transb),
                       m, n, k,
                       HIPBLAS_VALUE(alpha), A, lda,
                                             B, ldb,
                       HIPBLAS_VALUE(beta),  C, ldc );

    assert( rc == HIPBLAS_STATUS_SUCCESS );
    (void)rc;
    return CHAMELEON_SUCCESS;
}
