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
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_zlaswp_get( const RUNTIME_option_t *options,
                             int m0, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *U, int Um, int Un )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileU = U->get_blktile( U, Um, Un );
    int         *perm  = NULL; // get perm from ipiv

    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileU->format & CHAMELEON_TILE_FULLRANK );

#pragma omp task firstprivate( m0, k, ipiv, tileA, tileU ) depend( in:perm ) depend( in:tileA[0] ) depend( inout:tileU[0] )
    {
        TCORE_zlaswp_get( m0, A->m, A->n, k, tileA, tileU, perm );
    }

    (void)options;
}

void INSERT_TASK_zlaswp_set( const RUNTIME_option_t *options,
                             int m0, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );
    int         *invp  = NULL; // get invp from ipiv

    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );

#pragma omp task firstprivate( m0, k, ipiv, tileA, tileB ) depend( in:invp ) depend( in:tileA[0] ) depend( inout:tileB[0] )
    {
        TCORE_zlaswp_set( m0, A->m, A->n, k, tileA, tileB, invp );
    }

    (void)options;
}
