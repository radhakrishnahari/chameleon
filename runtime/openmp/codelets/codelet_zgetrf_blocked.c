/**
 *
 * @file openmp/codelet_zgetrf_blocked.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_blocked OpenMP codelets
 *
 * @version 1.3.0
 * @comment Codelets to perform panel factorization with partial pivoting
 *
 * @author Mathieu Faverge
 * @date 2024-03-29
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"

void INSERT_TASK_zgetrf_blocked_diag( const RUNTIME_option_t *options,
                                      int m, int n, int h, int m0, int ib,
                                      CHAM_desc_t *A, int Am, int An,
                                      CHAM_desc_t *U, int Um, int Un,
                                      CHAM_ipiv_t *ipiv )
{
    assert( 0 );
    (void)options;
    (void)m;
    (void)n;
    (void)h;
    (void)m0;
    (void)ib;
    (void)A;
    (void)Am;
    (void)An;
    (void)U;
    (void)Um;
    (void)Un;
    (void)ipiv;
}

void INSERT_TASK_zgetrf_blocked_offdiag( const RUNTIME_option_t *options,
                                         int m, int n, int h, int m0, int ib,
                                         CHAM_desc_t *A, int Am, int An,
                                         CHAM_desc_t *U, int Um, int Un,
                                         CHAM_ipiv_t *ipiv )
{
    assert( 0 );
    (void)options;
    (void)m;
    (void)n;
    (void)h;
    (void)m0;
    (void)ib;
    (void)A;
    (void)Am;
    (void)An;
    (void)U;
    (void)Um;
    (void)Un;
    (void)ipiv;
}

void INSERT_TASK_zgetrf_blocked_trsm( const RUNTIME_option_t *options,
                                      int m, int n, int h, int ib,
                                      CHAM_desc_t *U, int Um, int Un,
                                      CHAM_ipiv_t *ipiv )
{
    assert( 0 );
    (void)options;
    (void)m;
    (void)n;
    (void)h;
    (void)ib;
    (void)U;
    (void)Um;
    (void)Un;
    (void)ipiv;
}
