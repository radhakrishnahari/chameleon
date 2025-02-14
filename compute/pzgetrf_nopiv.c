/**
 *
 * @file pzgetrf_nopiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_nopiv parallel algorithm
 *
 * @version 1.3.0
 * @author Omar Zenati
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @author Terry Cojean
 * @author Matthieu Kuhn
 * @author Pierre Esterie
 * @date 2025-01-24
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m, n)  A,  m, n
#define WD(m)    WL, m, m
#define WL(m, n) WL, m, n
#define WU(m, n) WU, m, n

/**
 *  Parallel tile LU factorization with no pivoting - dynamic scheduling
 */
void chameleon_pzgetrf_nopiv_generic( CHAM_desc_t        *A,
                                      RUNTIME_sequence_t *sequence,
                                      RUNTIME_request_t  *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int k, m, n, ib;
    int tempkm, tempkn, tempmm, tempnn;

    CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
    CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

    if ( chamctxt->autominmax_enabled && (chamctxt->scheduler == RUNTIME_SCHED_STARPU) ) {
        int lookahead = chamctxt->lookahead;
        int nbtasks_per_step = (A->mt * A->nt) / (chameleon_desc_datadist_get_iparam(A, 0) * chameleon_desc_datadist_get_iparam(A, 1));
        int mintasks = nbtasks_per_step *  lookahead;
        int maxtasks = nbtasks_per_step * (lookahead+1);

        if ( CHAMELEON_Comm_rank() == 0 ) {
            chameleon_warning( "chameleon_pzgetrf_nopiv",
                               "Setting limit for the number of submitted tasks\n" );
        }
        RUNTIME_set_minmax_submitted_tasks( mintasks, maxtasks );
    }

    for (k = 0; k < chameleon_min(A->mt, A->nt); k++) {
        RUNTIME_iteration_push(chamctxt, k);

        tempkm = A->get_blkdim( A, k, DIM_m, A->m );
        tempkn = A->get_blkdim( A, k, DIM_n, A->n );

        options.priority = 2*A->nt - 2*k;
        INSERT_TASK_zgetrf_nopiv(
            &options,
            tempkm, tempkn, ib, A->mb,
            A(k, k), A->mb*k);

        for (m = k+1; m < A->mt; m++) {
            options.priority = 2*A->nt - 2*k - m;
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );
            INSERT_TASK_ztrsm(
                &options,
                ChamRight, ChamUpper, ChamNoTrans, ChamNonUnit,
                tempmm, tempkn, A->mb,
                zone, A(k, k),
                      A(m, k));
        }
        for (n = k+1; n < A->nt; n++) {
            tempnn = A->get_blkdim( A, n, DIM_n, A->n );
            options.priority = 2*A->nt - 2*k - n;
            INSERT_TASK_ztrsm(
                &options,
                ChamLeft, ChamLower, ChamNoTrans, ChamUnit,
                tempkm, tempnn, A->mb,
                zone, A(k, k),
                      A(k, n));

            for (m = k+1; m < A->mt; m++) {
                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                options.priority = 2*A->nt - 2*k  - n - m;
                INSERT_TASK_zgemm(
                    &options,
                    ChamNoTrans, ChamNoTrans,
                    tempmm, tempnn, A->mb, A->mb,
                    mzone, A(m, k),
                           A(k, n),
                    zone,  A(m, n));
            }

            RUNTIME_data_flush( sequence, A(k, n) );
        }

        for (m = k; m < A->mt; m++) {
            RUNTIME_data_flush( sequence, A(m, k) );
        }

        RUNTIME_iteration_pop(chamctxt);
    }

    RUNTIME_options_finalize(&options, chamctxt);
}

void chameleon_pzgetrf_nopiv_ws( CHAM_desc_t        *A,
                                 CHAM_desc_t        *WL,
                                 CHAM_desc_t        *WU,
                                 RUNTIME_sequence_t *sequence,
                                 RUNTIME_request_t  *request )
{
    CHAM_context_t  *chamctxt;
    RUNTIME_option_t options;

    int k, m, n, ib, p, q, lp, lq;
    int tempkm, tempkn, tempmm, tempnn;
    int lookahead, myp, myq;

    CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
    CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;
    lookahead = chamctxt->lookahead;
    myp = A->myrank / chameleon_desc_datadist_get_iparam(A, 1);
    myq = A->myrank % chameleon_desc_datadist_get_iparam(A, 1);

    for (k = 0; k < chameleon_min(A->mt, A->nt); k++) {
        RUNTIME_iteration_push(chamctxt, k);
        lp = (k % lookahead) * chameleon_desc_datadist_get_iparam(A, 0);
        lq = (k % lookahead) * chameleon_desc_datadist_get_iparam(A, 1);

        tempkm = A->get_blkdim( A, k, DIM_m, A->m );
        tempkn = A->get_blkdim( A, k, DIM_n, A->n );

        options.priority = 2*A->nt - 2*k;
        INSERT_TASK_zgetrf_nopiv(
            &options,
            tempkm, tempkn, ib, A->mb,
            A(k, k), A->mb*k);

        /**
         * Broadcast of A(k,k) along rings in both directions
         */
        {
            INSERT_TASK_zlacpy(
                &options,
                ChamUpperLower, tempkm, tempkn,
                A(  k, k ),
                WL( k, (k % chameleon_desc_datadist_get_iparam(A, 1)) + lq ) );
            INSERT_TASK_zlacpy(
                &options,
                ChamUpperLower, tempkm, tempkn,
                A(  k, k ),
                WU( (k % chameleon_desc_datadist_get_iparam(A, 0)) + lp, k ) );

            for ( q=1; q < chameleon_desc_datadist_get_iparam(A, 1); q++ ) {
                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpperLower, tempkm, tempkn,
                    WL( k, ((k+q-1) % chameleon_desc_datadist_get_iparam(A, 1)) + lq ),
                    WL( k, ((k+q)   % chameleon_desc_datadist_get_iparam(A, 1)) + lq ) );
            }

            for ( p=1; p < chameleon_desc_datadist_get_iparam(A, 0); p++ ) {
                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpperLower, tempkm, tempkn,
                    WU( ((k+p-1) % chameleon_desc_datadist_get_iparam(A, 0)) + lp, k ),
                    WU( ((k+p)   % chameleon_desc_datadist_get_iparam(A, 0)) + lp, k ) );
            }
        }
        RUNTIME_data_flush( sequence, A( k, k ) );

        for (m = k+1; m < A->mt; m++) {

            /* Skip the row if you are not involved with */
            if ( m%chameleon_desc_datadist_get_iparam(A, 0) != myp ) {
                continue;
            }

            options.priority = 2*A->nt - 2*k - m;
            tempmm = A->get_blkdim( A, m, DIM_m, A->m );

            assert( A->get_rankof( A, m, k ) == WU->get_rankof( WU, myp + lp, k) );
            INSERT_TASK_ztrsm(
                &options,
                ChamRight, ChamUpper, ChamNoTrans, ChamNonUnit,
                tempmm, tempkn, A->mb,
                zone, WU( myp + lp, k ),
                      A( m, k ) );

            /* Broadcast A(m,k) into temp buffers through a ring */
            {
                assert( A->get_rankof( A, m, k ) == WL->get_rankof( WL,  m, (k % chameleon_desc_datadist_get_iparam(A, 1)) + lq) );
                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpperLower, tempmm, tempkn,
                    A(  m, k ),
                    WL( m, (k % chameleon_desc_datadist_get_iparam(A, 1)) + lq) );

                for ( q=1; q < chameleon_desc_datadist_get_iparam(A, 1); q++ ) {
                    INSERT_TASK_zlacpy(
                        &options,
                        ChamUpperLower, tempmm, tempkn,
                        WL( m, ((k+q-1) % chameleon_desc_datadist_get_iparam(A, 1)) + lq ),
                        WL( m, ((k+q)   % chameleon_desc_datadist_get_iparam(A, 1)) + lq ) );
                }
            }
            RUNTIME_data_flush( sequence, A( m, k ) );
        }

        for (n = k+1; n < A->nt; n++) {

            /* Skip the column if you are not involved with */
            if ( n%chameleon_desc_datadist_get_iparam(A, 1) != myq ) {
                continue;
            }

            tempnn = A->get_blkdim( A, n, DIM_n, A->n );
            options.priority = 2*A->nt - 2*k - n;

            assert( A->get_rankof( A, k, n ) == WL->get_rankof( WL, k, myq+lq) );
            INSERT_TASK_ztrsm(
                &options,
                ChamLeft, ChamLower, ChamNoTrans, ChamUnit,
                tempkm, tempnn, A->mb,
                zone, WL( k, myq + lq ),
                      A( k, n ));

            /* Broadcast A(k,n) into temp buffers through a ring */
            {
                assert( A->get_rankof( A, k, n ) == WU->get_rankof( WU, (k%chameleon_desc_datadist_get_iparam(A, 0)) + lp, n) );
                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpperLower, tempkm, tempnn,
                    A(  k, n ),
                    WU( (k % chameleon_desc_datadist_get_iparam(A, 0)) + lp, n ) );

                for ( p=1; p < chameleon_desc_datadist_get_iparam(A, 0); p++ ) {
                    INSERT_TASK_zlacpy(
                        &options,
                        ChamUpperLower, tempkm, tempnn,
                        WU( ((k+p-1) % chameleon_desc_datadist_get_iparam(A, 0)) + lp, n ),
                        WU( ((k+p)   % chameleon_desc_datadist_get_iparam(A, 0)) + lp, n ) );
                }
            }
            RUNTIME_data_flush( sequence, A( k, n ) );

            for (m = k+1; m < A->mt; m++) {

                /* Skip the row if you are not involved with */
                if ( m%chameleon_desc_datadist_get_iparam(A, 0) != myp ) {
                    continue;
                }

                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                options.priority = 2*A->nt - 2*k  - n - m;

                assert( A->get_rankof( A, m, n ) == WL->get_rankof( WL, m, myq + lq) );
                assert( A->get_rankof( A, m, n ) == WU->get_rankof( WU, myp + lp, n) );

                INSERT_TASK_zgemm(
                    &options,
                    ChamNoTrans, ChamNoTrans,
                    tempmm, tempnn, A->mb, A->mb,
                    mzone, WL( m, myq + lq ),
                           WU( myp + lp, n ),
                    zone,  A( m, n ));
            }
        }
        RUNTIME_iteration_pop( chamctxt );
    }

    CHAMELEON_Desc_Flush( WL, sequence );
    CHAMELEON_Desc_Flush( WU, sequence );

    RUNTIME_options_finalize( &options, chamctxt );
}

void chameleon_pzgetrf_nopiv( struct chameleon_pzgetrf_nopiv_s *ws,
                              CHAM_desc_t                      *A,
                              RUNTIME_sequence_t               *sequence,
                              RUNTIME_request_t                *request )
{
    if ( ws && ws->use_workspace ) {
        chameleon_pzgetrf_nopiv_ws( A, &(ws->WL), &(ws->WU), sequence, request );
    }
    else {
        chameleon_pzgetrf_nopiv_generic( A, sequence, request );
    }
}
