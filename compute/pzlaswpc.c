/**
 *
 * @file pzlaswpc.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaswp parallel algorithm for column permutation.
 *
 * @version 1.3.0
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-03-24
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n)   A,       m, n
#define W(m,n)  &(ws->W), m, n

/**
 *  Permutation of the panel n at step k
 */
static inline void
chameleon_pzlaswpc_panel_permute( struct chameleon_pzlaswp_s *ws,
                                  cham_dir_t                  dir,
                                  CHAM_desc_t                *A,
                                  CHAM_ipiv_t                *ipiv,
                                  int                         m,
                                  int                         k,
                                  RUNTIME_option_t           *options )
{
    int                        n;
    int                        tempmm, tempkn;
    int                        withlacpy;

    tempmm = A->get_blkdim( A, m, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );

    /* Extract selected rows into U */
    withlacpy = options->withlacpy;
    options->withlacpy = 1;
    INSERT_TASK_zlacpy( options, ChamUpperLower, tempmm, tempkn,
                        A(m, k), W(m, A->myrank) );
    options->withlacpy = withlacpy;

    INSERT_TASK_zlaswpc_get( options, dir, tempkn, k*A->nb,
                             ipiv, k, A(m, k), W(m, A->myrank) );

    for ( n = k + 1; n < A->nt; n++ ) {
        /* Extract selected rows into A(k, n) */
        INSERT_TASK_zlaswpc_get( options, dir, tempkn, n*A->nb,
                                 ipiv, k, A(m, n), W(m, A->myrank) );
        /* Copy rows from A(k,n) into their final position */
        INSERT_TASK_zlaswpc_set( options, dir, tempkn, n*A->nb,
                                 ipiv, k, A(m, k), A(m, n) );
    }

    if ( ws->allreduce ) {
        INSERT_TASK_zperm_allreduce_col( options, dir, A, W(m, A->myrank), ipiv, k, m, k, ws );
    }
    else {
        INSERT_TASK_zperm_reduce_col( options, dir, A, W(m, A->myrank), ipiv, k, m, k, ws );
    }
}

/**
 *  Permutation of the panel n at step k
 */
static inline void
chameleon_pzlaswpc_panel_permute_batched( struct chameleon_pzlaswp_s *ws,
                                          cham_dir_t                  dir,
                                          CHAM_desc_t                *A,
                                          CHAM_ipiv_t                *ipiv,
                                          int                         m,
                                          int                         k,
                                          RUNTIME_option_t           *options )
{
    int                        n;
    int                        tempmm, tempkn;
    int                        withlacpy;

    void **clargs = malloc( sizeof(char *) );
    *clargs = NULL;

    tempmm = A->get_blkdim( A, m, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );

    /* Extract selected rows into U */
    withlacpy = options->withlacpy;
    options->withlacpy = 1;
    INSERT_TASK_zlacpy( options, ChamUpperLower, tempmm, tempkn,
                        A(m, k), W(m, A->myrank) );
    options->withlacpy = withlacpy;

    INSERT_TASK_zlaswpc_get( options, dir, tempkn, k*A->nb,
                             ipiv, k, A(m, k), W(m, A->myrank) );

    for ( n = k + 1; n < A->nt; n++ ) {
        INSERT_TASK_zlaswpc_batched( options, dir, n*A->nb, tempkn, (void *)ws, ipiv, k,
                                     A(m, n), A(m, k), W(m, A->myrank), clargs );
    }
    INSERT_TASK_zlaswpc_batched_flush( options, dir, ipiv, k, A(m, k), W(m, A->myrank), clargs );

    if ( ws->allreduce ) {
        INSERT_TASK_zperm_allreduce_col( options, dir, A, W(m, A->myrank), ipiv, k, m, k, ws );
    }
    else {
        INSERT_TASK_zperm_reduce_col( options, dir, A, W(m, A->myrank), ipiv, k, m, k, ws );
    }

    free( clargs );
}

static inline void
chameleon_pzlaswpc_panel( struct chameleon_pzlaswp_s *ws,
                          cham_dir_t                  dir,
                          CHAM_desc_t                *A,
                          CHAM_ipiv_t                *ipiv,
                          int                         m,
                          int                         k,
                          RUNTIME_option_t           *options,
                          RUNTIME_sequence_t         *sequence )
{
    CHAM_reduce_t *reduce = &(ws->reduce);
    int            tempmm, tempkn;

#if defined(CHAMELEON_USE_MPI)
    chameleon_get_proc_involved_in_rowpanelk_2dbc( A, m, k, reduce );
    if ( A->myrank == ipiv->get_rankof( ipiv, k, k ) ) {
        INSERT_TASK_zperm_allreduce_send_perm( options, dir, ipiv, k, A->myrank, reduce->np_involved, reduce->proc_involved );
        INSERT_TASK_zperm_allreduce_send_invp_col( options, dir, ipiv, k, A, m, k );
    }
    if ( A->myrank == chameleon_getrankof_2d( A, m, k ) ) {
        INSERT_TASK_zperm_allreduce_send_A( options, A, m, k, A->myrank, reduce->np_involved, reduce->proc_involved );
    }

    if ( !reduce->involved ) {
        return;
    }
#endif

    if ( ws->batch_size_swap == 0 ){
        chameleon_pzlaswpc_panel_permute( ws, dir, A, ipiv, m, k, options );
    }
    else {
        chameleon_pzlaswpc_panel_permute_batched( ws, dir, A, ipiv, m, k, options );
    }

    if ( A->myrank == chameleon_getrankof_2d( A, m, k ) ) {

        tempmm = A->get_blkdim( A, m, DIM_m, A->m );
        tempkn = A->get_blkdim( A, k, DIM_n, A->n );
        INSERT_TASK_zlacpy( options, ChamUpperLower, tempmm, tempkn,
                            W(m, A->myrank), A(m, k) );
        RUNTIME_data_flush( sequence, A(m, k) );
    }
    (void)reduce;
}

void
chameleon_pzlaswpc( struct chameleon_pzlaswp_s *ws,
                    cham_dir_t                  dir,
                    CHAM_desc_t                *A,
                    CHAM_ipiv_t                *IPIV,
                    RUNTIME_sequence_t         *sequence,
                    RUNTIME_request_t          *request )
{
    CHAM_context_t   *chamctxt;
    RUNTIME_option_t  options;

    int m, k;

    chamctxt = chameleon_context_self();
    if ( sequence->status != CHAMELEON_SUCCESS ) {
        return;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    if ( dir == ChamDirForward ) {
        for ( k = 0; k < IPIV->mt; k++ ) {
            for ( m = 0; m < A->mt; m++ ) {
                options.priority = A->mt-m;

                chameleon_pzlaswpc_panel( ws, dir, A, IPIV, m, k, &options, sequence );
            }
            RUNTIME_perm_flushk( sequence, IPIV, k );
        }
    }
    else {
        for ( k = IPIV->mt - 1; k > -1; k-- ) {
            for ( m = 0; m < A->mt; m++ ) {
                options.priority = A->mt-m;
                chameleon_pzlaswpc_panel( ws, dir, A, IPIV, m, k, &options, sequence );
            }
            RUNTIME_perm_flushk( sequence, IPIV, k );
        }
    }
    RUNTIME_options_finalize( &options, chamctxt );
}
