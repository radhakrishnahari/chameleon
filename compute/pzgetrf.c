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
 * @date 2024-03-16
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n)  A,        m, n
#define U(m,n)  &(ws->U), m, n
#define Up(m,n)  &(ws->Up), m, n

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

    tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
    tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

    /*
     * Algorithm per block without pivoting
     */
    INSERT_TASK_zgetrf_nopiv(
        options,
        tempkm, tempkn, ws->ib, A->mb,
         A(k, k), 0);

    for (m = k+1; m < A->mt; m++) {
        tempmm = (m == (A->mt - 1)) ? A->m - m * A->mb : A->mb;
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

    tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
    tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
    minmn  = chameleon_min( tempkm, tempkn );

    /*
     * Algorithm per column without pivoting
     */
    for(h=0; h<minmn; h++){
        INSERT_TASK_zgetrf_nopiv_percol_diag(
            options, tempkm, tempkn, h,
            A( k, k ), U( k, k ), A->mb * k );

        for (m = k+1; m < A->mt; m++) {
            tempmm = (m == (A->mt - 1)) ? A->m - m * A->mb : A->mb;
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

    tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
    tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
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
            tempmm = (m == (A->mt - 1)) ? A->m - m * A->mb : A->mb;
            INSERT_TASK_zgetrf_percol_offdiag(
                options,
                tempmm, tempkn, h, m * A->mb,
                A(m, k),
                ipiv );
        }

        /* Reduce globally (between MPI processes) */
        INSERT_TASK_zipiv_allreduce( A, options, ipiv, ws->proc_involved, k, h, tempkn );
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

    tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
    tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
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
            tempmm = (m == (A->mt - 1)) ? A->m - m * A->mb : A->mb;
            INSERT_TASK_zgetrf_panel_offdiag_batched( options, tempmm, tempkn, h, m * A->mb,
                                                      (void *)ws, A(m, k), clargs, ipiv );
        }
        INSERT_TASK_zgetrf_panel_offdiag_batched_flush( options, A, k, clargs, ipiv );

        INSERT_TASK_zipiv_allreduce( A, options, ipiv, ws->proc_involved, k, h, tempkn );
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

    if ( ! ws->involved ) {
        return;
    }

    tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
    tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
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
                tempmm = (m == (A->mt - 1)) ? A->m - m * A->mb : A->mb;
                INSERT_TASK_zgetrf_blocked_offdiag(
                    options,
                    tempmm, tempkn, j, m * A->mb, ws->ib,
                    A(m, k), Up(k, k),
                    ipiv );
            }

            assert( j <= minmn );
            /* Reduce globally (between MPI processes) */
            INSERT_TASK_zipiv_allreduce( A, options, ipiv, ws->proc_involved, k, j, tempkn );

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

    tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
    tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
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
                tempmm = (m == (A->mt - 1)) ? A->m - m * A->mb : A->mb;
                INSERT_TASK_zgetrf_panel_blocked_batched( options, tempmm, tempkn, j, m * A->mb,
                                                          (void *)ws, A(m, k), Up(k, k), clargs, ipiv );
            }
            INSERT_TASK_zgetrf_panel_blocked_batched_flush( options, A, k,
                                                            Up(k, k), clargs, ipiv );

            assert( j <= minmn );
            /* Reduce globally (between MPI processes) */
            INSERT_TASK_zipiv_allreduce( A, options, ipiv, ws->proc_involved, k, j, tempkn );

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
#if defined ( CHAMELEON_USE_MPI )
    int *proc_involved = malloc( sizeof( int ) * chameleon_min( A->p, A->mt - k) );
    int  b;

    /* 2DBC only */
    ws->involved = 0;
    for ( b = k; (b < A->mt) && ((b-k) < A->p); b ++ ) {
        int rank = chameleon_getrankof_2d( A, b, k );
        proc_involved[ b-k ] = rank;
        if ( rank == A->myrank ) {
            ws->involved = 1;
        }
    }
    ws->proc_involved = proc_involved;
    if ( ws->involved == 0 ) {
	free( proc_involved );
        return;
    }
#endif

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
        chameleon_attr_fallthrough;
    default:
        chameleon_pzgetrf_panel_facto_nopiv( ws, A, ipiv, k, options );
    }
#if defined ( CHAMELEON_USE_MPI )
    free( proc_involved );
#endif
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
        chameleon_attr_fallthrough;
    case ChamGetrfPPivPerColumn:
    {
        int m;
        int tempkm, tempkn, tempnn, minmn;

        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
        minmn  = chameleon_min( tempkm, tempkn );

        /* Extract selected rows into U */
        INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempnn,
                            A(k, n), U(k, n) );

        /*
         * perm array is made of size tempkm for the first row especially.
         * Otherwise, the final copy back to the tile may copy only a partial tile
         */
        INSERT_TASK_zlaswp_get( options, k*A->mb, tempkm,
                                ipiv, k, A(k, n), U(k, n) );

        for(m=k+1; m<A->mt; m++){
            /* Extract selected rows into A(k, n) */
            INSERT_TASK_zlaswp_get( options, m*A->mb, minmn,
                                    ipiv, k, A(m, n), U(k, n) );
            /* Copy rows from A(k,n) into their final position */
            INSERT_TASK_zlaswp_set( options, m*A->mb, minmn,
                                    ipiv, k, A(k, n), A(m, n) );
        }

        INSERT_TASK_zlacpy( options, ChamUpperLower, tempkm, tempnn,
                            U(k, n), A(k, n) );

        RUNTIME_data_flush( options->sequence, U(k, n) );
    }
    break;
    default:
        ;
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

    int m, tempkm, tempmm, tempnn;

    tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
    tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

    chameleon_pzgetrf_panel_permute( ws, A, ipiv, k, n, options );

    INSERT_TASK_ztrsm(
        options,
        ChamLeft, ChamLower, ChamNoTrans, ChamUnit,
        tempkm, tempnn, A->mb,
        zone, A(k, k),
              A(k, n) );

    for (m = k+1; m < A->mt; m++) {
        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;

        INSERT_TASK_zgemm(
            options,
            ChamNoTrans, ChamNoTrans,
            tempmm, tempnn, A->mb, A->mb,
            mzone, A(m, k),
                   A(k, n),
            zone,  A(m, n) );
    }

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

    int k, m, n;
    int min_mnt  = chameleon_min( A->mt, A->nt );
    int nb_tasks = 0;

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

        for (n = k+1; n < A->nt; n++) {
            options.priority = A->nt-n;
            chameleon_pzgetrf_panel_update( ws, A, IPIV, k, n, &options );
        }

        /* Flush panel k */
        for (m = k; m < A->mt; m++) {
            RUNTIME_data_flush( sequence, A(m, k) );
        }

        if ( chamctxt->lookahead > 0 ) {
            nb_tasks = RUNTIME_lookahead( chamctxt, k, nb_tasks );
        }
        RUNTIME_iteration_pop( chamctxt );
    }

    /* Backward pivoting */
    for (k = 1; k < min_mnt; k++) {
        for (n = 0; n < k; n++) {
            chameleon_pzgetrf_panel_permute( ws, A, IPIV, k, n, &options );
        }
        RUNTIME_perm_flushk( sequence, IPIV, k );
    }

    /* Initialize IPIV with default values if needed */
    if ( (ws->alg == ChamGetrfNoPivPerColumn) ||
         (ws->alg == ChamGetrfNoPiv ) )
    {
        INSERT_TASK_ipiv_init( &options, IPIV );
    }

    RUNTIME_options_finalize( &options, chamctxt );
}
