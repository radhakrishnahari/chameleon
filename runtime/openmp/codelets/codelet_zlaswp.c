/**
 *
 * @file openmp/codelet_zlaswp.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon OpenMP codelets to apply zlaswp on a panel
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matteo Marcos
 * @date 2025-07-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_zlaswp_get( const RUNTIME_option_t *options,
                             cham_side_t side, cham_dir_t dir,
                             int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *U, int Um, int Un )
{
    assert( 0 );
    (void)options;
    (void)side;
    (void)dir;
    (void)m0;
    (void)m;
    (void)n;
    (void)k;
    (void)ipiv;
    (void)ipivk;
    (void)A;
    (void)Am;
    (void)An;
}

void INSERT_TASK_zlaswp_set( const RUNTIME_option_t *options,
                             cham_side_t             side,
                             cham_dir_t              dir,
                             int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
    int         *invp  = NULL; // get invp from ipiv

    assert( A->format & CHAMELEON_TILE_FULLRANK );
    assert( B->format & CHAMELEON_TILE_FULLRANK );

#pragma omp task firstprivate( m0, k, ipiv, A, B ) depend( in:invp ) depend( in:tileA[0] ) depend( inout:tileB[0] )
    {
        TCORE_zlaswp_set( ChamLeft, m0, m, n, k, tileA, tileB, invp );
    }

    (void)options;
    (void)side;
    (void)dir;
}

void INSERT_TASK_zlaswp_ret( const RUNTIME_option_t *options,
                             CHAM_perm_t       *ws, int Wm, int Wn,
                             const CHAM_desc_t *A,  int Am, int An )
{
    assert( 0 );
    (void)options;
    (void)ws;
    (void)Wm;
    (void)Wn;
    (void)A;
    (void)Am;
    (void)An;
}

