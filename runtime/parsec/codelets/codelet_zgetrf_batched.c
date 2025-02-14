/**
 *
 * @file starpu/codelet_zgetrf_batched.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpanel batched parsec codelets
 *
 * @version 1.3.0
 * @comment Codelets to perform batched panel factorization with partial pivoting
 *
 * @author Alycia Lisito
 * @date 2024-04-02
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"

void
INSERT_TASK_zgetrf_panel_offdiag_batched( const RUNTIME_option_t *options,
                                          int m, int n, int h, int m0,
                                          void *ws,
                                          CHAM_desc_t *A, int Am, int An,
                                          void **clargs_ptr,
                                          CHAM_ipiv_t *ipiv )
{
    assert( 0 );
    (void)options;
    (void)m;
    (void)n;
    (void)h;
    (void)m0;
    (void)A;
    (void)Am;
    (void)An;
    (void)clargs_ptr;
    (void)ipiv;
}

void
INSERT_TASK_zgetrf_panel_offdiag_batched_flush( const RUNTIME_option_t *options,
                                                CHAM_desc_t *A, int An,
                                                void **clargs_ptr,
                                                CHAM_ipiv_t *ipiv )
{
    assert( 0 );
    (void)options;
    (void)A;
    (void)An;
    (void)clargs_ptr;
    (void)ipiv;
}

void
INSERT_TASK_zgetrf_panel_blocked_batched( const RUNTIME_option_t *options,
                                          int m, int n, int h, int m0,
                                          void *ws,
                                          CHAM_desc_t *A, int Am, int An,
                                          CHAM_desc_t *U, int Um, int Un,
                                          void **clargs_ptr,
                                          CHAM_ipiv_t *ipiv )
{
    assert( 0 );
    (void)options;
    (void)m;
    (void)n;
    (void)h;
    (void)m0;
    (void)ws;
    (void)A;
    (void)Am;
    (void)An;
    (void)U;
    (void)Um;
    (void)Un;
    (void)clargs_ptr;
    (void)ipiv;
}

void
INSERT_TASK_zgetrf_panel_blocked_batched_flush( const RUNTIME_option_t *options,
                                                CHAM_desc_t *A, int An,
                                                CHAM_desc_t *U, int Um, int Un,
                                                void **clargs_ptr,
                                                CHAM_ipiv_t *ipiv )
{
    assert( 0 );
    (void)options;
    (void)A;
    (void)An;
    (void)U;
    (void)Um;
    (void)Un;
    (void)clargs_ptr;
    (void)ipiv;
}
