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
    int P = chameleon_desc_datadist_get_iparam( A, 0 );
    int Q = chameleon_desc_datadist_get_iparam( A, 1 );

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
        int nbtasks_per_step = (A->mt * A->nt) / (P * Q);
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

            chameleon_data_flush( sequence, A(k, n), request->flush );
        }

        for (m = k; m < A->mt; m++) {
            chameleon_data_flush( sequence, A(m, k), request->flush );
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

    int            k, m, n, ib, p, q, lp, lq;
    int            tempkm, tempkn, tempmm, tempnn;
    int            lookahead;
    int            rankAmk, rankAkn;
    int            Wu_rank, Wl_rank, Wl_prev_rank, Wu_prev_rank;
    int            P      = chameleon_desc_datadist_get_iparam( A, 0 );
    int            Q      = chameleon_desc_datadist_get_iparam( A, 1 );
    int           *Wl_idx = malloc( sizeof(int) * P * Q );
    int           *Wu_idx = malloc( sizeof(int) * P * Q );
    custom_dist_t *dist   = A->get_rankof_init_arg;
    int            max_q  = ( dist ) ? chameleon_min( dist->dist_n, A->nt ) :
                                       chameleon_min( Q           , A->nt );
    int            max_p  = ( dist ) ? chameleon_min( dist->dist_m, A->mt ) :
                                       chameleon_min( P           , A->mt );

    CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
    CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;
    lookahead = chamctxt->lookahead;

    for ( k = 0; k < chameleon_min(A->mt, A->nt); k++ ) {
        RUNTIME_iteration_push(chamctxt, k);
        lp = (k % lookahead) * P * Q;
        lq = (k % lookahead) * P * Q;

        tempkm = A->get_blkdim( A, k, DIM_m, A->m );
        tempkn = A->get_blkdim( A, k, DIM_n, A->n );
        memset( Wl_idx, 0x00, P * Q * sizeof(int) );
        memset( Wu_idx, 0x00, P * Q * sizeof(int) );

        options.priority = 2*A->nt - 2*k;
        INSERT_TASK_zgetrf_nopiv(
            &options,
            tempkm, tempkn, ib, A->mb,
            A(k, k), A->mb*k);

        /**
         * Broadcast of A(k,k) along rings in both directions
         */
        {
            Wl_rank = A->get_rankof( A, k, k );
            Wu_rank = A->get_rankof( A, k, k );
            INSERT_TASK_zlacpy(
                &options,
                ChamUpperLower, tempkm, tempkn,
                A(  k, k ),
                WL( k, Wl_rank + lq ) );
            INSERT_TASK_zlacpy(
                &options,
                ChamUpperLower, tempkm, tempkn,
                A(  k, k ),
                WU( Wu_rank + lp, k ) );

            Wl_idx[Wl_rank] = 1;
            Wu_idx[Wu_rank] = 1;
            Wl_prev_rank = Wl_rank;
            Wu_prev_rank = Wu_rank;
            for ( q = 0; q < max_q; q++ ) {
                Wl_rank = A->get_rankof( A, k, q );

                /* Skip if already sent to process */
                if ( Wl_idx[Wl_rank] != 0 ) {
                    continue;
                }

                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpperLower, tempkm, tempkn,
                    WL( k, Wl_prev_rank + lq ),
                    WL( k, Wl_rank      + lq ) );
                Wl_idx[Wl_rank] = 1;
                Wl_prev_rank    = Wl_rank;
            }

            for ( p = 0; p < max_p; p++ ) {
                Wu_rank = A->get_rankof( A, p, k );

                /* Skip if already sent to process */
                if ( Wu_idx[Wu_rank] != 0 ) {
                    continue;
                }

                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpperLower, tempkm, tempkn,
                    WU( Wu_prev_rank + lp, k ),
                    WU( Wu_rank + lp, k ) );
                Wu_idx[Wu_rank] = 1;
                Wu_prev_rank    = Wu_rank;
            }
        }
        chameleon_data_flush( sequence, A( k, k ), request->flush );

        for ( m = k + 1; m < A->mt; m++ ) {

            /* Skip the row if you are not involved with */
            if ( !chameleon_involved_in_rowpanelk( A, m ) ) {
                continue;
            }

            options.priority = 2*A->nt - 2*k - m;
            tempmm           = A->get_blkdim( A, m, DIM_m, A->m );
            rankAmk          = A->get_rankof( A, m, k );

            assert( rankAmk == WU->get_rankof( WU, rankAmk + lp, k) );
            INSERT_TASK_ztrsm(
                &options,
                ChamRight, ChamUpper, ChamNoTrans, ChamNonUnit,
                tempmm, tempkn, A->mb,
                zone, WU( rankAmk + lp, k ),
                      A( m, k ) );

            /* Broadcast A(m,k) into temp buffers through a ring */
            {
                memset( Wl_idx, 0x00, P * Q * sizeof(int) );
                Wl_rank = A->get_rankof(A, m, k);

                assert( A->get_rankof( A, m, k ) == WL->get_rankof( WL,  m, Wl_rank + lq) );
                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpperLower, tempmm, tempkn,
                    A(  m, k ),
                    WL( m, Wl_rank + lq) );

                Wl_idx[Wl_rank] = 1;
                Wl_prev_rank    = Wl_rank;

                for ( q = 0; q < max_q; q++ ) {
                    Wl_rank = A->get_rankof( A, m, q );

                    /* Skip if already sent to process */
                    if ( Wl_idx[Wl_rank] != 0 ) {
                        continue;
                    }
                    INSERT_TASK_zlacpy(
                        &options,
                        ChamUpperLower, tempmm, tempkn,
                        WL( m, Wl_prev_rank + lq ),
                        WL( m, Wl_rank + lq ) );
                    Wl_idx[Wl_rank] = 1;
                    Wl_prev_rank    = Wl_rank;
                }
            }
            chameleon_data_flush( sequence, A( m, k ), request->flush );
        }

        for ( n = k + 1; n < A->nt; n++ ) {

            /* Skip the column if you are not involved with */
            if ( !chameleon_involved_in_panelk(A, n) ) {
                continue;
            }

            options.priority = 2*A->nt - 2*k - n;
            tempnn           = A->get_blkdim( A, n, DIM_n, A->n );
            rankAkn          = A->get_rankof( A, k, n );

            assert( rankAkn == WL->get_rankof( WL, k, rankAkn + lq) );
            INSERT_TASK_ztrsm(
                &options,
                ChamLeft, ChamLower, ChamNoTrans, ChamUnit,
                tempkm, tempnn, A->mb,
                zone, WL( k, rankAkn + lq ),
                A( k, n ));

            /* Broadcast A(k,n) into temp buffers through a ring */
            {
                Wu_rank = A->get_rankof( A, k, n );
                assert( rankAkn == WU->get_rankof( WU, Wu_rank + lp, n) );
                memset( Wu_idx, 0x00, P * Q * sizeof(int) );
                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpperLower, tempkm, tempnn,
                    A(  k, n ),
                    WU( Wu_rank + lp, n ) );
                Wu_idx[Wu_rank] = 1;
                Wu_prev_rank    = Wu_rank;

                for ( p=0; p < max_p; p++ ) {
                    Wu_rank = A->get_rankof( A, p, n );

                    /* Skip if already sent to process */
                    if ( Wu_idx[Wu_rank] != 0 ) {
                        continue;
                    }
                    INSERT_TASK_zlacpy(
                        &options,
                        ChamUpperLower, tempkm, tempnn,
                        WU( Wu_prev_rank + lp, n ),
                        WU( Wu_rank      + lp, n ) );
                    Wu_idx[Wu_rank] = 1;
                    Wu_prev_rank    = Wu_rank;
                }
            }
            chameleon_data_flush( sequence, A( k, n ), request->flush );

            for ( m = k+1; m < A->mt; m++ ) {

                /* Skip the row if you are not involved with */
                if ( !chameleon_involved_in_rowpanelk( A, m ) ) {
                    continue;
                }

                tempmm = A->get_blkdim( A, m, DIM_m, A->m );
                options.priority = 2*A->nt - 2*k  - n - m;

                assert( A->get_rankof( A, m, n ) == WL->get_rankof( WL, m, A->myrank + lq) );
                assert( A->get_rankof( A, m, n ) == WU->get_rankof( WU, A->myrank + lp, n) );

                INSERT_TASK_zgemm(
                    &options,
                    ChamNoTrans, ChamNoTrans,
                    tempmm, tempnn, A->mb, A->mb,
                    mzone, WL( m, A->myrank + lq ),
                           WU( A->myrank + lp, n ),
                    zone,  A( m, n ));
            }
        }
        RUNTIME_iteration_pop( chamctxt );
    }
    free( Wl_idx );
    free( Wu_idx );

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
