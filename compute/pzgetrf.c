/**
 *
 * @file pzgetrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf parallel algorithm
 *
 * @version 1.3.0
 * @author Omar Zenati
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-10-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include "chameleon/flops.h"

#define A(m,n)   A,                m, n
#define U(m,n)   &(ws->U),         m, n
#define Up(m,n)  &(ws->Up),        m, n
#define Wu(m,n)  &(ws->laswp->Wu), m, n
#define Wl(m,n)  &(ws->Wl),        m, n
#define Ws(m,n)  &(ws->laswp->ws), m, n

/*
 * Function to compute the optimal batch size when generating the graph that
 * covers all the possibilities of pivoting within the factorized panel.
 *
 * @param[in] ws
 *      The workspace data structure associated to the algorithm that holds all
 *      extra information that may be needed for LU factorization
 *
 * @param[in] mt
 *      The number of tiles to factorize within the panel
 *
 * @param[in] mb
 *      The default tile size.
 *
 * @param[in] j
 *      The index of the column to factorize
 *
 * @return The optimal batch size.
 *
 */
static inline int
chameleon_pzgetrf_batch_size( const struct chameleon_pzgetrf_s *ws,
                              int mt, int mb, int j )
{
    if ( ws->batch_adaptive == 0 ) {
        return ( ( j % ws->ib ) != 0 ) ? ws->batch_size_blas2 : ws->batch_size_blas3;
    }

    double flops     = flops_zgetrf_blocked_offdiag( mb, mb, j, ws->ib );
    int    batch_max = chameleon_min( CHAMELEON_BATCH_SIZE, mt );
    int    batch_sze = batch_max;

    /**
     * First solution. (Alycia Lisito)
     *
     * This solution aims at maximizing the batch size.
     */
    if ( j != 0 ) {
        batch_sze = chameleon_min( chameleon_max( ws->flops_min / flops, 1 ), batch_max );
    }

    if ( mt % batch_sze != 0 ) {
        batch_sze = chameleon_min( chameleon_ceil( mt, mt / batch_sze ), batch_max );
    }

    /**
     * Second solution. (Mathieu Faverge)
     *
     * This solution aims at balancing the load rather than increasing the task size.
     */
    /*
     batch_sze = chameleon_ceil( mt,
                                 chameleon_max( chameleon_ceil( mt, batch_max ),
                                                ( mt * flops ) / ws->flops_min ) );
     */
    return batch_sze;
}

/*
 * All the functions below are panel factorization variant.
 * The parameters are:
 *   @param[inout] ws
 *      The data structure associated to the algorithm that holds all extra
 *      information that may be needed for LU factorization
 *
 *   @param[inout] A
 *      The descriptor of the full matrix A (not just the panel)
 *
 *   @param[inout] ipiv
 *      The descriptor of the pivot array associated to A.
 *
 *   @param[in] k
 *      The index of the column to factorize
 *
 *   @param[inout] options
 *      The runtime options data structure to pass through all insert_task calls.
 */
static inline void
chameleon_pzgetrf_panel_facto_percol( struct chameleon_pzgetrf_s *ws,
                                      CHAM_desc_t                *A,
                                      CHAM_ipiv_t                *ipiv,
                                      CHAM_desc_pivot_t          *pivot,
                                      int                         k,
                                      RUNTIME_option_t           *options )
{
    int m, h;
    int tempkm, tempkn, tempmm, minmn;

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );
    minmn  = chameleon_min( tempkm, tempkn );

    /* Update the number of column */
    pivot->n = minmn;

    /*
     * Algorithm per column with pivoting
     */
    for (h=0; h<=minmn; h++){
        INSERT_TASK_zgetrf_percol_diag(
            options,
            tempkm, tempkn, h, k * A->mb,
            A(k, k),
            ipiv, pivot );

        for (m = k+1; m < A->mt; m++) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            INSERT_TASK_zgetrf_percol_offdiag(
                options,
                tempmm, tempkn, h, m * A->mb,
                A(m, k),
                pivot );
        }

        /* Reduce globally (between MPI processes) */
        INSERT_TASK_zipiv_allreduce( options, A, pivot, k, h, tempkn, ws->laswp );
    }

    /* Flush temporary data used for the pivoting */
    INSERT_TASK_ipiv_to_perm( options, k * A->mb, tempkm, minmn, 0, A->m, ipiv, k );
    RUNTIME_pivot_flushk( options->sequence, pivot, A->myrank );
}

/*
 *  Factorization of panel k - dynamic scheduling - batched version / stock
 */
static inline void
chameleon_pzgetrf_panel_facto_percol_batched( struct chameleon_pzgetrf_s *ws,
                                              CHAM_desc_t                *A,
                                              CHAM_ipiv_t                *ipiv,
                                              CHAM_desc_pivot_t          *pivot,
                                              int                         k,
                                              RUNTIME_option_t           *options )
{
    int m, h;
    int tempkm, tempkn, tempmm, minmn;
    void **clargs = malloc( sizeof(char *) );
    memset( clargs, 0, sizeof(char *) );

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );
    minmn  = chameleon_min( tempkm, tempkn );

    /* Update the number of column */
    pivot->n = minmn;

    /*
     * Algorithm per column with pivoting (no recursion)
     */
    /* Iterate on current panel column */
    /* Since index h scales column h-1, we need to iterate up to minmn (included) */
    for ( h = 0; h <= minmn; h++ ) {

        INSERT_TASK_zgetrf_percol_diag( options, tempkm, tempkn, h, k * A->mb, A(k, k), ipiv, pivot );

        for ( m = k+1; m < A->mt; m++ ) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            INSERT_TASK_zgetrf_panel_offdiag_batched( options, tempmm, tempkn, h, m * A->mb,
                                                      (void *)ws, A(m, k), clargs, pivot );
        }
        INSERT_TASK_zgetrf_panel_offdiag_batched_flush( options, A, k, clargs, pivot );

        INSERT_TASK_zipiv_allreduce( options, A, pivot, k, h, tempkn, ws->laswp );
    }

    free( clargs );

    /* Flush temporary data used for the pivoting */
    INSERT_TASK_ipiv_to_perm( options, k * A->mb, tempkm, minmn, 0, A->m, ipiv, k );
    RUNTIME_pivot_flushk( options->sequence, pivot, A->myrank );
}

static inline void
chameleon_pzgetrf_panel_facto_blocked( struct chameleon_pzgetrf_s *ws,
                                       CHAM_desc_t                *A,
                                       CHAM_ipiv_t                *ipiv,
                                       CHAM_desc_pivot_t          *pivot,
                                       int                         k,
                                       RUNTIME_option_t           *options )
{
    const RUNTIME_request_t *request = options->request;
    int m, h, b, nbblock;
    int tempkm, tempkn, tempmm, minmn;

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );
    minmn  = chameleon_min( tempkm, tempkn );

    /* Update the number of column */
    pivot->n = minmn;
    nbblock = chameleon_ceil( minmn, ws->ib );

    /*
     * Algorithm per column with pivoting
     */
    for (b=0; b<nbblock; b++){
        int hmax = b == nbblock-1 ? minmn + 1 - b * ws->ib : ws->ib;

        for (h=0; h<hmax; h++){
            int j = h + b * ws->ib;

            INSERT_TASK_zgetrf_blocked_diag(
                options,
                tempkm, tempkn, j, k * A->mb, ws->ib,
                A(k, k), Up(k, k),
                ipiv, pivot );

            for (m = k+1; m < A->mt; m++) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                INSERT_TASK_zgetrf_blocked_offdiag(
                    options,
                    tempmm, tempkn, j, m * A->mb, ws->ib,
                    A(m, k), Up(k, k),
                    pivot );
            }

            assert( j <= minmn );
            /* Reduce globally (between MPI processes) */
            INSERT_TASK_zipiv_allreduce( options, A, pivot, k, j, tempkn, ws->laswp );

            if ( ( b < (nbblock-1) ) && ( h == hmax-1 ) ) {
                INSERT_TASK_zgetrf_blocked_trsm(
                    options,
                    ws->ib, tempkn, j+1, ws->ib,
                    Up(k, k),
                    pivot );
            }
        }
    }
    chameleon_data_flush( options->sequence, Up(k, k), request->flush );

    /* Flush temporary data used for the pivoting */
    INSERT_TASK_ipiv_to_perm( options, k * A->mb, tempkm, minmn, 0, A->m, ipiv, k );
    RUNTIME_pivot_flushk( options->sequence, pivot, A->myrank );
}

/*
 *  Factorization of panel k - dynamic scheduling - batched version / stock
 */
static inline void
chameleon_pzgetrf_panel_facto_blocked_batched( struct chameleon_pzgetrf_s *ws,
                                               CHAM_desc_t                *A,
                                               CHAM_ipiv_t                *ipiv,
                                               CHAM_desc_pivot_t          *pivot,
                                               int                         k,
                                               RUNTIME_option_t           *options )
{
    int m, h, b, nbblock, hmax, j;
    int tempkm, tempkn, tempmm, minmn;
    void **clargs = malloc( sizeof(char *) );
    memset( clargs, 0, sizeof(char *) );

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );
    minmn  = chameleon_min( tempkm, tempkn );

    /* Update the number of column */
    pivot->n = minmn;
    nbblock = chameleon_ceil( minmn, ws->ib );

    /*
     * Algorithm per column with pivoting (no recursion)
     */
    /* Iterate on current panel column */
    /* Since index h scales column h-1, we need to iterate up to minmn (included) */
    for ( b = 0; b < nbblock; b++ ) {
        hmax = b == nbblock-1 ? minmn + 1 - b * ws->ib : ws->ib;

        for ( h = 0; h < hmax; h++ ) {
            j =  h + b * ws->ib;

            ws->batch_size = chameleon_pzgetrf_batch_size( ws, A->mt - k, A->nb, j );
            for ( m = k; m < A->mt; m++ ) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                INSERT_TASK_zgetrf_panel_blocked_batched( options, tempmm, tempkn, j, m * A->mb,
                                                          (void *)ws, A(m, k), Up(k, k), clargs, ipiv, pivot );
            }
            INSERT_TASK_zgetrf_panel_blocked_batched_flush( options, A, k,
                                                            Up(k, k), clargs, ipiv, pivot );

            assert( j <= minmn );
            /* Reduce globally (between MPI processes) */
            INSERT_TASK_zipiv_allreduce( options, A, pivot, k, j, tempkn, ws->laswp );

            if ( (b < (nbblock-1)) && (h == hmax-1) ) {
                INSERT_TASK_zgetrf_blocked_trsm(
                    options,
                    ws->ib, tempkn, b * ws->ib + hmax, ws->ib,
                    Up(k, k),
                    pivot );
            }
        }
    }

    free( clargs );

    /* Flush temporary data used for the pivoting */
    INSERT_TASK_ipiv_to_perm( options, k * A->mb, tempkm, minmn, 0, A->m, ipiv, k );
    RUNTIME_pivot_flushk( options->sequence, pivot, A->myrank );
}

static inline void
chameleon_pzgetrf_panel_facto( struct chameleon_pzgetrf_s *ws,
                               CHAM_desc_t                *A,
                               CHAM_ipiv_t                *ipiv,
                               CHAM_desc_pivot_t          *pivot,
                               int                         k,
                               RUNTIME_option_t           *options )
{
#if defined(CHAMELEON_USE_MPI)
    CHAM_reduce_t *reduce = &(ws->laswp->reduce);
    chameleon_get_proc_involved_in_panelk_2dbc( A, k, k, reduce );
    if ( !reduce->involved ) {
        return;
    }
#endif

    /* TODO: Should be replaced by a function pointer */
    switch( ws->alg ) {
    case ChamGetrfPPivPerColumn:
        if ( ws->batch_size_blas2 > 0 ) {
            chameleon_pzgetrf_panel_facto_percol_batched( ws, A, ipiv, pivot, k, options );
        }
        else {
            chameleon_pzgetrf_panel_facto_percol( ws, A, ipiv, pivot, k, options );
        }
        break;

    case ChamGetrfPPiv:
    default:
        if ( ws->batch_size_blas2 > 0 ) {
            chameleon_pzgetrf_panel_facto_blocked_batched( ws, A, ipiv, pivot, k, options );
        }
        else {
            chameleon_pzgetrf_panel_facto_blocked( ws, A, ipiv, pivot, k, options );
        }
        break;
    }
}

/**
 *  Permutation of the panel n at step k
 */
static inline void
chameleon_pzgetrf_panel_permute( struct chameleon_pzgetrf_s *ws,
                                 CHAM_desc_t                *A,
                                 CHAM_ipiv_t                *ipiv,
                                 int                         k,
                                 int                         n,
                                 RUNTIME_option_t           *options )
{
    switch( ws->alg ) {
    case ChamGetrfPPiv:
    case ChamGetrfPPivPerColumn:
    {
        int m;
        int tempkm, tempkn, tempnn, tempmm, minmn;
        int withlacpy;

        tempkm = A->get_blkdim( A, k, DIM_m, A->m );
        tempkn = A->get_blkdim( A, k, DIM_n, A->n );
        tempnn = A->get_blkdim( A, n, DIM_n, A->n );
        minmn  = chameleon_min( tempkm, tempkn );

        /* Extract selected rows into U */
        withlacpy = options->withlacpy;
        options->withlacpy = 1;
        INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempnn,
                            A(k, n), Wu(A->myrank, n) );
        options->withlacpy = withlacpy;

        /*
         * perm array is made of size tempkm for the first row especially.
         * Otherwise, the final copy back to the tile may copy only a partial tile
         */
        INSERT_TASK_zlaswp_get( options, ChamLeft, ChamDirForward, k*A->mb, tempkm, tempnn, tempkm,
                                ipiv, k, A(k, n), Wu(A->myrank, n) );

        for(m=k+1; m<A->mt; m++){
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            /* Extract selected rows into A(k, n) */
            INSERT_TASK_zlaswp_get( options, ChamLeft, ChamDirForward, m*A->mb, tempmm, tempnn, minmn,
                                    ipiv, k, A(m, n), Wu(A->myrank, n) );
            /* Copy rows from A(k,n) into their final position */
            INSERT_TASK_zlaswp_set( options, ChamLeft, ChamDirForward, m*A->mb, tempmm, tempnn, minmn,
                                    ipiv, k, A(k, n), A(m, n) );
        }

        INSERT_TASK_zperm_allreduce( options, ChamDirForward, A(k, n), ipiv, k, Wu(A->myrank, n), ws->laswp, A->myrank, n );

        if ( ws->laswp->reduce.np_involved != 1 ) {
            INSERT_TASK_zlaswp_ret( options, Ws(A->myrank, n), Wu(A->myrank, n) );
        }
    }
    break;
    default:
        ;
    }
}

static inline void
chameleon_pzgetrf_panel_permute_batched( struct chameleon_pzgetrf_s *ws,
                                         CHAM_desc_t                *A,
                                         CHAM_ipiv_t                *ipiv,
                                         int                         k,
                                         int                         n,
                                         RUNTIME_option_t           *options )
{
    switch( ws->alg ) {
    case ChamGetrfPPiv:
    case ChamGetrfPPivPerColumn:
    {
        int m;
        int tempkm, tempkn, tempmm, tempnn, minmn;
        int withlacpy;

        void **clargs = malloc( sizeof(char *) );
        *clargs = NULL;

        tempkm = A->get_blkdim( A, k, DIM_m, A->m );
        tempkn = A->get_blkdim( A, k, DIM_n, A->n );
        tempnn = A->get_blkdim( A, n, DIM_n, A->n );
        minmn  = chameleon_min( tempkm, tempkn );

        /* Extract selected rows into U */
        withlacpy = options->withlacpy;
        options->withlacpy = 1;
        INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempnn,
                            A(k, n), Wu(A->myrank, n) );
        options->withlacpy = withlacpy;

        /*
         * perm array is made of size tempkm for the first row especially.
         * Otherwise, the final copy back to the tile may copy only a partial tile
         */
        INSERT_TASK_zlaswp_get( options, ChamLeft, ChamDirForward, k*A->mb, tempkm, tempnn, tempkm,
                                ipiv, k, A(k, n), Wu(A->myrank, n) );

        for(m=k+1; m<A->mt; m++){
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            INSERT_TASK_zlaswp_batched( options, ChamLeft, ChamDirForward, m*A->mb, tempmm, tempnn, minmn, (void *)ws->laswp, ipiv, k,
                                        A(m, n), A(k, n), Wu(A->myrank, n), clargs );
        }
        INSERT_TASK_zlaswp_batched_flush( options, ChamDirForward, ipiv, k, A(k, n), Wu(A->myrank, n), clargs );

        INSERT_TASK_zperm_allreduce( options, ChamDirForward, A(k, n), ipiv, k, Wu(A->myrank, n), ws->laswp, A->myrank, n );

        if ( ws->laswp->reduce.np_involved != 1 ) {
            INSERT_TASK_zlaswp_ret( options, Ws(A->myrank, n), Wu(A->myrank, n) );
        }

        free( clargs );
    }
    break;
    default:
        ;
    }
}

static inline void
chameleon_pzgetrf_panel_permute_forward( struct chameleon_pzgetrf_s *ws,
                                         CHAM_desc_t                *A,
                                         CHAM_ipiv_t                *ipiv,
                                         int                         k,
                                         int                         n,
                                         RUNTIME_option_t           *options,
                                         RUNTIME_sequence_t         *sequence )
{
#if defined(CHAMELEON_USE_MPI)
    CHAM_reduce_t *reduce = &(ws->laswp->reduce);
    chameleon_get_proc_involved_in_panelk_2dbc( A, k, n, reduce );
    if ( A->myrank == chameleon_getrankof_2d( A, k, k ) ) {
        INSERT_TASK_zperm_allreduce_send_perm( options, ChamDirForward, ipiv, k, A->myrank, reduce->np_involved, reduce->proc_involved );
        INSERT_TASK_zperm_allreduce_send_invp_row( options, ChamDirForward, ipiv, k, A, k, n );
    }
    if ( A->myrank == chameleon_getrankof_2d( A, k, n ) ) {
        INSERT_TASK_zperm_allreduce_send_A( options, A, k, n, A->myrank, reduce->np_involved, reduce->proc_involved );
    }

    if ( !reduce->involved ) {
        return;
    }
#endif

    if ( ws->laswp->batch_size_swap > 0 ) {
        chameleon_pzgetrf_panel_permute_batched( ws, A, ipiv, k, n, options );
    }
    else {
        chameleon_pzgetrf_panel_permute( ws, A, ipiv, k, n, options );
    }

    RUNTIME_cpui_flushk( sequence, A->myrank, Ws(A->myrank, n) );
}

static inline void
chameleon_pzgetrf_panel_permute_backward( struct chameleon_pzgetrf_s *ws,
                                          CHAM_desc_t                *A,
                                          CHAM_ipiv_t                *ipiv,
                                          int                         k,
                                          int                         n,
                                          RUNTIME_option_t           *options,
                                          RUNTIME_sequence_t         *sequence )
{
    const RUNTIME_request_t *request = options->request;
    CHAM_reduce_t *reduce = &(ws->laswp->reduce);
    int            tempkm, tempnn;

#if defined(CHAMELEON_USE_MPI)
    chameleon_get_proc_involved_in_panelk_2dbc( A, k, n, reduce );
    if ( A->myrank == chameleon_getrankof_2d( A, k, k ) ) {
        INSERT_TASK_zperm_allreduce_send_perm( options, ChamDirForward, ipiv, k, A->myrank, reduce->np_involved, reduce->proc_involved );
        INSERT_TASK_zperm_allreduce_send_invp_row( options, ChamDirForward, ipiv, k, A, k, n );
    }
    if ( A->myrank == chameleon_getrankof_2d( A, k, n ) ) {
        INSERT_TASK_zperm_allreduce_send_A( options, A, k, n, A->myrank, reduce->np_involved, reduce->proc_involved );
    }

    if ( !reduce->involved ) {
        return;
    }
#endif

    if ( ws->laswp->batch_size_swap > 0 ) {
        chameleon_pzgetrf_panel_permute_batched( ws, A, ipiv, k, n, options );
    }
    else {
        chameleon_pzgetrf_panel_permute( ws, A, ipiv, k, n, options );
    }

    if ( A->myrank == chameleon_getrankof_2d( A, k, n ) ) {

        if ( ws->laswp->reduce.np_involved == 1 ) {
            tempkm = A->get_blkdim( A, k, DIM_m, A->m );
            tempnn = A->get_blkdim( A, n, DIM_n, A->n );
            INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempnn,
                               Wu(A->myrank, n), A(k, n) );
        }
        else {
            INSERT_TASK_zlaswp_ret( options, Ws(A->myrank, n), A(k, n) );
        }
        chameleon_data_flush( sequence, A(k, n), request->flush );
    }
    RUNTIME_cpui_flushk( sequence, A->myrank, Ws(A->myrank, n) );
    (void)reduce;
}

#if defined(CHAMELEON_USE_MPI)
static inline void
chameleon_pzgetrf_panel_update_ws( struct chameleon_pzgetrf_s *ws,
                                   CHAM_desc_t                *A,
                                   int                         k,
                                   RUNTIME_option_t           *options )
{
    CHAM_context_t          *chamctxt = chameleon_context_self();
    const RUNTIME_request_t *request  = options->request;
    int m, n, tempmm, tempkn, tempkm, p, q, involved, np;
    int lookahead = chamctxt->lookahead;
    int P         = chameleon_desc_datadist_get_iparam(A, 0);
    int Q         = chameleon_desc_datadist_get_iparam(A, 1);
    int lq        = (k % lookahead) * Q;
    int myp       = A->myrank / Q;

    tempkn = A->get_blkdim( A, k, DIM_n, A->n );

    if ( k >= ws->ringswitch ) {
        for ( m = k+1; m < A->mt; m++ ) {
            if ( ( m % P ) != myp ) continue;

            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            INSERT_TASK_zlacpy(
                options,
                ChamUpperLower, tempmm, tempkn,
                A( m, k ),
                Wl( m, ( k % Q ) + lq ) );

            for ( q = 1; q < Q; q++ ) {
                INSERT_TASK_zlacpy(
                    options,
                    ChamUpperLower, tempmm, tempkn,
                    Wl( m, ( ( k + q - 1 ) % Q ) + lq ),
                    Wl( m, ( ( k + q )     % Q ) + lq ) );
            }
            chameleon_data_flush( options->sequence, A(m, k), request->flush );
        }
    }
    else {
        for ( m = k+1; m < A->mt; m++ ) {
            if ( ( m % P ) != myp ) continue;

            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            for ( q = 0; q < Q; q++ ) {
                INSERT_TASK_zlacpy(
                    options,
                    ChamUpperLower, tempmm, tempkn,
                    A( m, k ),
                    Wl( m, ( ( k + q )% Q ) + lq ) );
            }
            chameleon_data_flush( options->sequence, A(m, k), request->flush );
        }
    }

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    np = chameleon_desc_datadist_get_iparam(A, 1) * chameleon_desc_datadist_get_iparam(A, 0);

    /* Send Akk for replicated trsm */
    if ( A->myrank == chameleon_getrankof_2d( A, k, k ) ) {
        for ( p = 0; p < np; p++ ) {
            involved = 0;
            for ( n = k+1; n < A->nt; n++ ) {
                if ( chameleon_p_involved_in_panelk_2dbc( A, n, p ) ) {
                    involved = 1;
                    break;
                }
            }
            if ( involved ) {
                INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempkn,
                                    A(k, k), Wu(p, k) );
            }
        }
    }
    else {
        involved = 0;
        for ( n = k+1; n < A->nt; n++ ) {
            if ( chameleon_involved_in_panelk_2dbc( A, n ) ) {
                involved = 1;
                break;
            }
        }
        if ( involved ) {
            INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempkn,
                                A(k, k), Wu(A->myrank, k) );
        }
    }
    chameleon_data_flush( options->sequence, A(k, k), request->flush );
}
#endif

static inline void
chameleon_pzgetrf_panel_update( struct chameleon_pzgetrf_s *ws,
                                CHAM_desc_t                *A,
                                CHAM_ipiv_t                *ipiv,
                                int                         k,
                                int                         n,
                                RUNTIME_option_t           *options,
                                RUNTIME_sequence_t         *sequence )
{
    const CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
    const CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;
    CHAM_context_t             *chamctxt = chameleon_context_self();
    const RUNTIME_request_t    *request  = options->request;
    CHAM_reduce_t              *reduce   = &(ws->laswp->reduce);

    int m, tempkm, tempmm, tempnn;

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempnn = A->get_blkdim( A, n, DIM_n, A->n );

    chameleon_pzgetrf_panel_permute_forward( ws, A, ipiv, k, n, options, sequence );

#if defined(CHAMELEON_USE_MPI)
    int rankAmn;
    int lookahead = chamctxt->lookahead;
    int myq       = A->myrank % chameleon_desc_datadist_get_iparam(A, 1);
    int lq        = (k % lookahead) * chameleon_desc_datadist_get_iparam(A, 1);

    if ( reduce->involved ) {
        INSERT_TASK_ztrsm(
            options,
            ChamLeft, ChamLower, ChamNoTrans, ChamUnit,
            tempkm, tempnn, A->mb,
            zone, Wu(A->myrank, k),
                  Wu(A->myrank, n) );
    }

    for (m = k+1; m < A->mt; m++) {
        tempmm = A->get_blkdim( A, m, DIM_m, A->m );
        rankAmn = A->get_rankof( A, m, n );

        if ( A->myrank == rankAmn ) {
            INSERT_TASK_zgemm(
                options,
                ChamNoTrans, ChamNoTrans,
                tempmm, tempnn, A->mb, A->mb,
                mzone, Wl( m, myq + lq ),
                       Wu( A->myrank, n ),
                zone,  A( m, n ) );
        }
    }

    if ( A->myrank == chameleon_getrankof_2d( A, k, n ) ) {
        INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempnn,
                            Wu(A->myrank, n), A(k, n) );
    }
#else
    INSERT_TASK_ztrsm(
        options,
        ChamLeft, ChamLower, ChamNoTrans, ChamUnit,
        tempkm, tempnn, A->mb,
        zone, A( k, k ),
              Wu( A->myrank, n ) );

    for (m = k+1; m < A->mt; m++) {
        tempmm = A->get_blkdim( A, m, DIM_m, A->m );

        INSERT_TASK_zgemm(
            options,
            ChamNoTrans, ChamNoTrans,
            tempmm, tempnn, A->mb, A->mb,
            mzone, A( m, k ),
                   Wu( A->myrank, n ),
            zone,  A( m, n ) );
    }

    INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempnn,
                        Wu(A->myrank, n), A(k, n) );

#endif

    chameleon_data_flush( options->sequence, Wu(A->myrank, n), request->flush );
    chameleon_data_flush( options->sequence, A(k, n), request->flush );
    (void)reduce;
    (void)chamctxt;
}

/**
 *  Parallel tile LU factorization with no pivoting - dynamic scheduling
 */
void chameleon_pzgetrf( struct chameleon_pzgetrf_s *ws,
                        CHAM_desc_t                *A,
                        CHAM_ipiv_t                *IPIV,
                        RUNTIME_sequence_t         *sequence,
                        RUNTIME_request_t          *request )
{
    CHAM_context_t    *chamctxt;
    RUNTIME_option_t   options;
    CHAM_desc_pivot_t *pivot = &(ws->pivot);

    int k, m, n;
    int min_mnt = chameleon_min( A->mt, A->nt );

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    for (k = 0; k < min_mnt; k++) {
        RUNTIME_iteration_push( chamctxt, k );

        options.priority = A->nt;
        /*
         * Do the panel factorization only if the current proc contributes in the
         * block column k.
         */
        options.forcesub = chameleon_involved_in_panelk_2dbc( A, k );
        if ( chameleon_involved_in_panelk_2dbc( A, k ) ) {
            chameleon_pzgetrf_panel_facto( ws, A, IPIV, pivot, k, &options );
        }
        options.forcesub = 0;

#if defined(CHAMELEON_USE_MPI)
        chameleon_pzgetrf_panel_update_ws( ws, A, k, &options );
#endif

        for (n = k+1; n < A->nt; n++) {
            options.priority = A->nt-n;
            if ( chameleon_involved_in_panelk_2dbc( A, k ) ||
                 chameleon_involved_in_panelk_2dbc( A, n ) )
            {
                chameleon_pzgetrf_panel_update( ws, A, IPIV, k, n, &options, sequence );
            }
        }

        /* Flush panel k */
        for (m = k+1; m < A->mt; m++) {
            chameleon_data_flush( sequence, A(m, k), request->flush );
        }
        chameleon_data_flush( sequence, Wu(A->myrank, k), request->flush );

        RUNTIME_iteration_pop( chamctxt );
    }
    CHAMELEON_Desc_Flush( &(ws->Wl), sequence );
    CHAMELEON_Ipiv_Flush( IPIV, sequence );
    chameleon_pivot_destroy_submit( pivot, sequence );

    /* Backward pivoting */
    for (k = 1; k < min_mnt; k++) {
        for (n = 0; n < k; n++) {
            if ( chameleon_involved_in_panelk_2dbc( A, k ) ||
                 chameleon_involved_in_panelk_2dbc( A, n ) )
            {
                chameleon_pzgetrf_panel_permute_backward( ws, A, IPIV, k, n, &options, sequence );
            }
            chameleon_data_flush( sequence, Wu(A->myrank, n), request->flush );
        }
        RUNTIME_perm_flushk( sequence, IPIV, k );
    }
    CHAMELEON_Desc_Flush( &(ws->laswp->Wu), sequence );

    RUNTIME_options_finalize( &options, chamctxt );
}
