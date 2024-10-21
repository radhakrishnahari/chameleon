/**
 *
 * @file pzgetrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
 * @date 2024-10-18
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n)   A,         m, n
#define U(m,n)   &(ws->U),  m, n
#define Up(m,n)  &(ws->Up), m, n
#define Wu(m,n)  &(ws->Wu), m, n
#define Wl(m,n)  &(ws->Wl), m, n

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
chameleon_pzgetrf_panel_facto_nopiv( struct chameleon_pzgetrf_s *ws,
                                     CHAM_desc_t                *A,
                                     CHAM_ipiv_t                *ipiv,
                                     int                         k,
                                     RUNTIME_option_t           *options )
{
    const CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t) 1.0;
    int m, tempkm, tempkn, tempmm;

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );

    /*
     * Algorithm per block without pivoting
     */
    INSERT_TASK_zgetrf_nopiv(
        options,
        tempkm, tempkn, ws->ib, A->mb,
         A(k, k), 0);

    for (m = k+1; m < A->mt; m++) {
        tempmm = A->get_blkdim( A, m, DIM_m, A->m );
        INSERT_TASK_ztrsm(
            options,
            ChamRight, ChamUpper, ChamNoTrans, ChamNonUnit,
            tempmm, tempkn, A->mb,
            zone, A(k, k),
                  A(m, k) );
    }
}

static inline void
chameleon_pzgetrf_panel_facto_nopiv_percol( struct chameleon_pzgetrf_s *ws,
                                            CHAM_desc_t                *A,
                                            CHAM_ipiv_t                *ipiv,
                                            int                         k,
                                            RUNTIME_option_t           *options )
{
    int m, h;
    int tempkm, tempkn, tempmm, minmn;

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );
    minmn  = chameleon_min( tempkm, tempkn );

    /*
     * Algorithm per column without pivoting
     */
    for(h=0; h<minmn; h++){
        INSERT_TASK_zgetrf_nopiv_percol_diag(
            options, tempkm, tempkn, h,
            A( k, k ), U( k, k ), A->mb * k );

        for (m = k+1; m < A->mt; m++) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            INSERT_TASK_zgetrf_nopiv_percol_trsm(
                options, tempmm, tempkn, h,
                A( m, k ), U( k, k ) );
        }
    }

    RUNTIME_data_flush( options->sequence, U(k, k) );
}

static inline void
chameleon_pzgetrf_panel_facto_percol( struct chameleon_pzgetrf_s *ws,
                                      CHAM_desc_t                *A,
                                      CHAM_ipiv_t                *ipiv,
                                      int                         k,
                                      RUNTIME_option_t           *options )
{
    int m, h;
    int tempkm, tempkn, tempmm, minmn;

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );
    minmn  = chameleon_min( tempkm, tempkn );

    /* Update the number of column */
    ipiv->n = minmn;

    /*
     * Algorithm per column with pivoting
     */
    for (h=0; h<=minmn; h++){
        INSERT_TASK_zgetrf_percol_diag(
            options,
            tempkm, tempkn, h, k * A->mb,
            A(k, k),
            ipiv );

        for (m = k+1; m < A->mt; m++) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            INSERT_TASK_zgetrf_percol_offdiag(
                options,
                tempmm, tempkn, h, m * A->mb,
                A(m, k),
                ipiv );
        }

        /* Reduce globally (between MPI processes) */
        INSERT_TASK_zipiv_allreduce( options, A, ipiv, k, h, tempkn, ws );
    }

    /* Flush temporary data used for the pivoting */
    INSERT_TASK_ipiv_to_perm( options, k * A->mb, tempkm, minmn, ipiv, k );
    RUNTIME_ipiv_flushk( options->sequence, ipiv, A->myrank );
}

/*
 *  Factorization of panel k - dynamic scheduling - batched version / stock
 */
static inline void
chameleon_pzgetrf_panel_facto_percol_batched( struct chameleon_pzgetrf_s *ws,
                                              CHAM_desc_t                *A,
                                              CHAM_ipiv_t                *ipiv,
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
    ipiv->n = minmn;

    /*
     * Algorithm per column with pivoting (no recursion)
     */
    /* Iterate on current panel column */
    /* Since index h scales column h-1, we need to iterate up to minmn (included) */
    for ( h = 0; h <= minmn; h++ ) {

        INSERT_TASK_zgetrf_percol_diag( options, tempkm, tempkn, h, k * A->mb, A(k, k), ipiv );

        for ( m = k+1; m < A->mt; m++ ) {
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            INSERT_TASK_zgetrf_panel_offdiag_batched( options, tempmm, tempkn, h, m * A->mb,
                                                      (void *)ws, A(m, k), clargs, ipiv );
        }
        INSERT_TASK_zgetrf_panel_offdiag_batched_flush( options, A, k, clargs, ipiv );

        INSERT_TASK_zipiv_allreduce( options, A, ipiv, k, h, tempkn, ws );
    }

    free( clargs );

    /* Flush temporary data used for the pivoting */
    INSERT_TASK_ipiv_to_perm( options, k * A->mb, tempkm, minmn, ipiv, k );
    RUNTIME_ipiv_flushk( options->sequence, ipiv, A->myrank );
}

static inline void
chameleon_pzgetrf_panel_facto_blocked( struct chameleon_pzgetrf_s *ws,
                                       CHAM_desc_t                *A,
                                       CHAM_ipiv_t                *ipiv,
                                       int                         k,
                                       RUNTIME_option_t           *options )
{
    int m, h, b, nbblock;
    int tempkm, tempkn, tempmm, minmn;

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempkn = A->get_blkdim( A, k, DIM_n, A->n );
    minmn  = chameleon_min( tempkm, tempkn );

    /* Update the number of column */
    ipiv->n = minmn;
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
                ipiv );

            for (m = k+1; m < A->mt; m++) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                INSERT_TASK_zgetrf_blocked_offdiag(
                    options,
                    tempmm, tempkn, j, m * A->mb, ws->ib,
                    A(m, k), Up(k, k),
                    ipiv );
            }

            assert( j <= minmn );
            /* Reduce globally (between MPI processes) */
            INSERT_TASK_zipiv_allreduce( options, A, ipiv, k, j, tempkn, ws );

            if ( ( b < (nbblock-1) ) && ( h == hmax-1 ) ) {
                INSERT_TASK_zgetrf_blocked_trsm(
                    options,
                    ws->ib, tempkn, j+1, ws->ib,
                    Up(k, k),
                    ipiv );
            }
        }
    }
    RUNTIME_data_flush( options->sequence, Up(k, k) );

    /* Flush temporary data used for the pivoting */
    INSERT_TASK_ipiv_to_perm( options, k * A->mb, tempkm, minmn, ipiv, k );
    RUNTIME_ipiv_flushk( options->sequence, ipiv, A->myrank );
}

/*
 *  Factorization of panel k - dynamic scheduling - batched version / stock
 */
static inline void
chameleon_pzgetrf_panel_facto_blocked_batched( struct chameleon_pzgetrf_s *ws,
                                               CHAM_desc_t                *A,
                                               CHAM_ipiv_t                *ipiv,
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
    ipiv->n = minmn;
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

            for ( m = k; m < A->mt; m++ ) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                INSERT_TASK_zgetrf_panel_blocked_batched( options, tempmm, tempkn, j, m * A->mb,
                                                          (void *)ws, A(m, k), Up(k, k), clargs, ipiv );
            }
            INSERT_TASK_zgetrf_panel_blocked_batched_flush( options, A, k,
                                                            Up(k, k), clargs, ipiv );

            assert( j <= minmn );
            /* Reduce globally (between MPI processes) */
            INSERT_TASK_zipiv_allreduce( options, A, ipiv, k, j, tempkn, ws );

            if ( (b < (nbblock-1)) && (h == hmax-1) ) {
                INSERT_TASK_zgetrf_blocked_trsm(
                    options,
                    ws->ib, tempkn, b * ws->ib + hmax, ws->ib,
                    Up(k, k),
                    ipiv );
            }
        }
    }

    free( clargs );

    /* Flush temporary data used for the pivoting */
    INSERT_TASK_ipiv_to_perm( options, k * A->mb, tempkm, minmn, ipiv, k );
    RUNTIME_ipiv_flushk( options->sequence, ipiv, A->myrank );
}

static inline void
chameleon_pzgetrf_panel_facto( struct chameleon_pzgetrf_s *ws,
                               CHAM_desc_t                *A,
                               CHAM_ipiv_t                *ipiv,
                               int                         k,
                               RUNTIME_option_t           *options )
{
    chameleon_get_proc_involved_in_panelk_2dbc( A, k, k, ws );
    if ( !ws->involved ) {
        return;
    }

    /* TODO: Should be replaced by a function pointer */
    switch( ws->alg ) {
    case ChamGetrfNoPivPerColumn:
        chameleon_pzgetrf_panel_facto_nopiv_percol( ws, A, ipiv, k, options );
        break;

    case ChamGetrfPPivPerColumn:
        if ( ws->batch_size > 0 ) {
            chameleon_pzgetrf_panel_facto_percol_batched( ws, A, ipiv, k, options );
        }
        else {
            chameleon_pzgetrf_panel_facto_percol( ws, A, ipiv, k, options );
        }
        break;

    case ChamGetrfPPiv:
        if ( ws->batch_size > 0 ) {
            chameleon_pzgetrf_panel_facto_blocked_batched( ws, A, ipiv, k, options );
        }
        else {
            chameleon_pzgetrf_panel_facto_blocked( ws, A, ipiv, k, options );
        }
        break;

    case ChamGetrfNoPiv:
    default:
        chameleon_pzgetrf_panel_facto_nopiv( ws, A, ipiv, k, options );
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
        int tempkm, tempkn, tempnn, minmn;
        int withlacpy;

        chameleon_get_proc_involved_in_panelk_2dbc( A, k, n, ws );
        if ( A->myrank == chameleon_getrankof_2d( A, k, k ) ) {
            INSERT_TASK_zperm_allreduce_send_perm( options, ipiv, k, A->myrank, ws->np_involved, ws->proc_involved );
            INSERT_TASK_zperm_allreduce_send_invp( options, ipiv, k, A, k, n );
        }
        if ( A->myrank == chameleon_getrankof_2d( A, k, n ) ) {
            INSERT_TASK_zperm_allreduce_send_A( options, A, k, n, A->myrank, ws->np_involved, ws->proc_involved );
        }

        if ( !ws->involved ) {
            return;
        }

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
        INSERT_TASK_zlaswp_get( options, k*A->mb, tempkm,
                                ipiv, k, A(k, n), Wu(A->myrank, n) );

        for(m=k+1; m<A->mt; m++){
            /* Extract selected rows into A(k, n) */
            INSERT_TASK_zlaswp_get( options, m*A->mb, minmn,
                                    ipiv, k, A(m, n), Wu(A->myrank, n) );
            /* Copy rows from A(k,n) into their final position */
            INSERT_TASK_zlaswp_set( options, m*A->mb, minmn,
                                    ipiv, k, A(k, n), A(m, n) );
        }

        INSERT_TASK_zperm_allreduce( options, A, Wu(A->myrank, n), ipiv, k, k, n, ws );
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
        int tempkm, tempkn, tempnn, minmn;
        int withlacpy;

        chameleon_get_proc_involved_in_panelk_2dbc( A, k, n, ws );
        if ( A->myrank == chameleon_getrankof_2d( A, k, k ) ) {
            INSERT_TASK_zperm_allreduce_send_perm( options, ipiv, k, A->myrank, ws->np_involved, ws->proc_involved );
            INSERT_TASK_zperm_allreduce_send_invp( options, ipiv, k, A, k, n );
        }
        if ( A->myrank == chameleon_getrankof_2d( A, k, n ) ) {
            INSERT_TASK_zperm_allreduce_send_A( options, A, k, n, A->myrank, ws->np_involved, ws->proc_involved );
        }

        if ( !ws->involved ) {
            return;
        }

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
        INSERT_TASK_zlaswp_get( options, k*A->mb, tempkm,
                                ipiv, k, A(k, n), Wu(A->myrank, n) );

        for(m=k+1; m<A->mt; m++){
            INSERT_TASK_zlaswp_batched( options, m*A->mb, minmn, (void *)ws, ipiv, k,
                                        A(m, n), A(k, n), Wu(A->myrank, n), clargs );
        }
        INSERT_TASK_zlaswp_batched_flush( options, ipiv, k, A(k, n), Wu(A->myrank, n), clargs );

        INSERT_TASK_zperm_allreduce( options, A, Wu(A->myrank, n), ipiv, k, k, n, ws );

        free( clargs );
    }
    break;
    default:
        ;
    }
}

static inline void
chameleon_pzgetrf_panel_update_ws( struct chameleon_pzgetrf_s *ws,
                                   CHAM_desc_t                *A,
                                   int                         k,
                                   RUNTIME_option_t           *options )
{
    CHAM_context_t  *chamctxt = chameleon_context_self();
    int m, tempmm, tempkn, q;
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
            RUNTIME_data_flush( options->sequence, A(m, k) );
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
            RUNTIME_data_flush( options->sequence, A(m, k) );
        }
    }
}

static inline void
chameleon_pzgetrf_panel_update( struct chameleon_pzgetrf_s *ws,
                                CHAM_desc_t                *A,
                                CHAM_ipiv_t                *ipiv,
                                int                         k,
                                int                         n,
                                RUNTIME_option_t           *options )
{
    const CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
    const CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;
    CHAM_context_t  *chamctxt = chameleon_context_self();

    int m, tempkm, tempmm, tempnn, rankAmn, p;

    int lookahead = chamctxt->lookahead;
    int myq       = A->myrank % chameleon_desc_datadist_get_iparam(A, 1);
    int lq        = (k % lookahead) * chameleon_desc_datadist_get_iparam(A, 1);

    tempkm = A->get_blkdim( A, k, DIM_m, A->m );
    tempnn = A->get_blkdim( A, n, DIM_n, A->n );

    if ( ws->batch_size > 0 ) {
        chameleon_pzgetrf_panel_permute_batched( ws, A, ipiv, k, n, options );
    }
    else {
        chameleon_pzgetrf_panel_permute( ws, A, ipiv, k, n, options );
    }

    if ( A->myrank == chameleon_getrankof_2d( A, k, k ) ) {
        for ( p = 0; p < ws->np_involved; p++ ) {
            INSERT_TASK_ztrsm(
                options,
                ChamLeft, ChamLower, ChamNoTrans, ChamUnit,
                tempkm, tempnn, A->mb,
                zone, A(k, k),
                      Wu(ws->proc_involved[p], n) );
            RUNTIME_data_flush( options->sequence, Wu(ws->proc_involved[p], n) );
        }
    }
    else if ( ws->involved ) {
        INSERT_TASK_ztrsm(
            options,
            ChamLeft, ChamLower, ChamNoTrans, ChamUnit,
            tempkm, tempnn, A->mb,
            zone, A(k, k),
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

    RUNTIME_data_flush( options->sequence, Wu(A->myrank, n) );
    RUNTIME_data_flush( options->sequence, A(k, k) );
    RUNTIME_data_flush( options->sequence, A(k, n) );
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
    CHAM_context_t  *chamctxt;
    RUNTIME_option_t options;

    int k, m, n, tempkm, tempnn;
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
            chameleon_pzgetrf_panel_facto( ws, A, IPIV, k, &options );
        }
        options.forcesub = 0;

        chameleon_pzgetrf_panel_update_ws( ws, A, k, &options );

        for (n = k+1; n < A->nt; n++) {
            options.priority = A->nt-n;
            if ( chameleon_involved_in_panelk_2dbc( A, k ) ||
                 chameleon_involved_in_panelk_2dbc( A, n ) )
            {
                chameleon_pzgetrf_panel_update( ws, A, IPIV, k, n, &options );
            }
        }

        /* Flush panel k */
        for (m = k; m < A->mt; m++) {
            RUNTIME_data_flush( sequence, A(m, k) );
        }

        RUNTIME_iteration_pop( chamctxt );
    }
    CHAMELEON_Desc_Flush( &(ws->Wl), sequence );

    /* Backward pivoting */
    if ( ws->batch_size > 0 ) {
        for (k = 1; k < min_mnt; k++) {
            for (n = 0; n < k; n++) {
                if ( chameleon_involved_in_panelk_2dbc( A, k ) ||
                    chameleon_involved_in_panelk_2dbc( A, n ) )
                {
                    chameleon_pzgetrf_panel_permute_batched( ws, A, IPIV, k, n, &options );
                    if ( A->myrank == chameleon_getrankof_2d( A, k, n ) ) {
                        tempkm = A->get_blkdim( A, k, DIM_m, A->m );
                        tempnn = A->get_blkdim( A, n, DIM_n, A->n );
                        INSERT_TASK_zlacpy( &options, ChamUpperLower, tempkm, tempnn,
                                            Wu(A->myrank, n), A(k, n) );
                        RUNTIME_data_flush( sequence, A(k, n) );
                    }
                }
                RUNTIME_data_flush( sequence, Wu(A->myrank, n) );
            }
            RUNTIME_perm_flushk( sequence, IPIV, k );
        }
    }
    else {
        for (k = 1; k < min_mnt; k++) {
            for (n = 0; n < k; n++) {
                if ( chameleon_involved_in_panelk_2dbc( A, k ) ||
                    chameleon_involved_in_panelk_2dbc( A, n ) )
                {
                    chameleon_pzgetrf_panel_permute( ws, A, IPIV, k, n, &options );
                    if ( A->myrank == chameleon_getrankof_2d( A, k, n ) ) {
                        tempkm = A->get_blkdim( A, k, DIM_m, A->m );
                        tempnn = A->get_blkdim( A, n, DIM_n, A->n );
                        INSERT_TASK_zlacpy( &options, ChamUpperLower, tempkm, tempnn,
                                            Wu(A->myrank, n), A(k, n) );
                        RUNTIME_data_flush( sequence, A(k, n) );
                    }
                }
                RUNTIME_data_flush( sequence, Wu(A->myrank, n) );
            }
            RUNTIME_perm_flushk( sequence, IPIV, k );
        }
    }
    CHAMELEON_Desc_Flush( &(ws->Wu), sequence );

    /* Initialize IPIV with default values if needed */
    if ( (ws->alg == ChamGetrfNoPivPerColumn) ||
         (ws->alg == ChamGetrfNoPiv ) )
    {
        INSERT_TASK_ipiv_init( &options, IPIV );
    }

    RUNTIME_options_finalize( &options, chamctxt );
}
