/**
 *
 * @file pzlaswp.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaswp parallel algorithm
 *
 * @version 1.3.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-03-24
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n)   A,         m, n
#define Wu(m,n)  &(ws->Wu), m, n

/**
 *  Permutation of the panel n at step k
 */
static inline void
chameleon_pzlaswp_panel_permute( struct chameleon_pzgetrf_s *ws,
                                 cham_dir_t                  dir,
                                 CHAM_desc_t                *A,
                                 CHAM_ipiv_t                *ipiv,
                                 int                         k,
                                 int                         n,
                                 RUNTIME_option_t           *options )
{
    int m;
    int tempkm, tempnn;
    int withlacpy;

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempnn = A->get_blkdim( A, n, DIM_n, A->n );

    /* Extract selected rows into U */
    withlacpy = options->withlacpy;
    options->withlacpy = 1;
    INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempnn,
                       A(k, n), Wu(A->myrank, n) );
    options->withlacpy = withlacpy;

    INSERT_TASK_zlaswp_get( options, dir, k*A->mb, tempkm,
                           ipiv, k, A(k, n), Wu(A->myrank, n) );

    for ( m = k + 1; m < A->mt; m++ ) {
        /* Extract selected rows into A(k, n) */
        INSERT_TASK_zlaswp_get( options, dir, m*A->mb, tempkm,
                               ipiv, k, A(m, n), Wu(A->myrank, n) );
        /* Copy rows from A(k,n) into their final position */
        INSERT_TASK_zlaswp_set( options, dir, m*A->mb, tempkm,
                               ipiv, k, A(k, n), A(m, n) );
    }

    INSERT_TASK_zperm_allreduce_row( options, dir, A, Wu(A->myrank, n), ipiv, k, k, n, ws );
}

static inline void
chameleon_pzlaswp_panel( struct chameleon_pzgetrf_s *ws,
                         cham_dir_t                  dir,
                         CHAM_desc_t                *A,
                         CHAM_ipiv_t                *ipiv,
                         int                         k,
                         int                         n,
                         RUNTIME_option_t           *options,
                         RUNTIME_sequence_t         *sequence )
{
    int tempkm, tempnn;

#if defined(CHAMELEON_USE_MPI)
    chameleon_get_proc_involved_in_panelk_2dbc( A, k, n, ws );
    if ( A->myrank == chameleon_getrankof_2d( A, k, k ) ) {
        INSERT_TASK_zperm_allreduce_send_perm( options, dir, ipiv, k, A->myrank, ws->np_involved, ws->proc_involved );
        INSERT_TASK_zperm_allreduce_send_invp_row( options, dir, ipiv, k, A, k, n );
    }
    if ( A->myrank == chameleon_getrankof_2d( A, k, n ) ) {
        INSERT_TASK_zperm_allreduce_send_A( options, A, k, n, A->myrank, ws->np_involved, ws->proc_involved );
    }

    if ( !ws->involved ) {
        return;
    }
#endif

    chameleon_pzlaswp_panel_permute( ws, dir, A, ipiv, k, n, options );

    if ( A->myrank == chameleon_getrankof_2d( A, k, n ) ) {

        tempkm = A->get_blkdim( A, k, DIM_m, A->m );
        tempnn = A->get_blkdim( A, n, DIM_n, A->n );
        INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempnn,
                            Wu(A->myrank, n), A(k, n) );
        RUNTIME_data_flush( sequence, A(k, n) );
    }
}

void
chameleon_pzlaswp( struct chameleon_pzgetrf_s *ws,
                   cham_dir_t                  dir,
                   CHAM_desc_t                *A,
                   CHAM_ipiv_t                *IPIV,
                   RUNTIME_sequence_t         *sequence,
                   RUNTIME_request_t          *request )
{
    CHAM_context_t   *chamctxt;
    RUNTIME_option_t  options;

    int n, k;

    chamctxt = chameleon_context_self();
    if ( sequence->status != CHAMELEON_SUCCESS ) {
        return;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    if ( dir == ChamDirForward ) {
        for ( k = 0; k < IPIV->mt; k++ ) {
            for ( n = 0; n < A->nt; n++ ) {
                options.priority = A->nt-n;

                chameleon_pzlaswp_panel( ws, dir, A, IPIV, k, n, &options, sequence );
            }
            RUNTIME_perm_flushk( sequence, IPIV, k );
        }
    }
    else {
        for ( k = IPIV->mt - 1; k > -1; k-- ) {
            for ( n = 0; n < A->nt; n++ ) {
                options.priority = A->nt-n;
                chameleon_pzlaswp_panel( ws, dir, A, IPIV, k, n, &options, sequence );
            }
            RUNTIME_perm_flushk( sequence, IPIV, k );
        }
    }
    RUNTIME_options_finalize( &options, chamctxt );
}

