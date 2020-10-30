/**
 *
 * @file pzunmlq_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmlq_param parallel algorithm
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Raphael Boucherie
 * @date 2020-03-03
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include <stdlib.h>

#define A(m,n) A,  m,  n
#define C(m,n) C,  m,  n
#define T(m,n) T,  m,  n
#define D(m,n) D,  m,  n

/**
 *  Parallel application of Q using tile V - LQ factorization - dynamic scheduling
 */
void chameleon_pzunmlq_param( int genD, const libhqr_tree_t *qrtree,
                              cham_side_t side, cham_trans_t trans,
                              CHAM_desc_t *A, CHAM_desc_t *C,
                              CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t *T;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n, i, p;
    int temppm, temppn, tempmm, tempnn, tempkm,tempkmin;
    int ib, KT, L;
    int node, nbtiles, *tiles;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

    if ( D == NULL ) {
        D    = A;
        genD = 0;
    }

    /*
     * zunmlq  = A->nb * ib
     * ztpmlqt = A->nb * ib
     */
    ws_worker = A->nb * ib;

#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmlq  =     A->nb * ib
     * ztpmlqt = 3 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 3 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    /* Initialisation of tiles */
    tiles = (int*)calloc( qrtree->mt, sizeof(int) );

    KT = chameleon_min( A->mt, A->nt );
    if (side == ChamLeft ) {
        if (trans == ChamNoTrans) {
            /*
             *  ChamLeft / ChamNoTrans
             */
            for (k = 0; k < KT; k++) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkm = k == A->mt - 1 ? A->m - k * A->mb : A->mb;

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    p = qrtree->getm(qrtree, k, i);

                    temppm   = p == C->mt-1 ? C->m - p * C->mb : C->mb;
                    tempkmin = chameleon_min( temppm, tempkm );


                    if ( genD ) {
                        int tempDpn = p == D->nt-1 ? D->n-p*D->nb : D->nb;

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamUpper, tempkmin, tempDpn, A->nb,
                            A(k, p),
                            D(k, p) );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamLower, tempkmin, tempDpn,
                            0., 1.,
                            D(k, p) );
#endif
                    }
                    for (n = 0; n < C->nt; n++) {
                        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
                        INSERT_TASK_zunmlq(
                            &options, side, trans,
                            temppm, tempnn, tempkmin, ib, T->nb,
                            D(k, p),
                            T(k, p),
                            C(p, n));
                    }
                    RUNTIME_data_flush( sequence, D(k, p) );
                    RUNTIME_data_flush( sequence, T(k, p) );
                }

                /* Setting the order of the tiles*/
                nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

                for (i = 0; i < nbtiles; i++) {
                    m = tiles[i];
                    p = qrtree->currpiv(qrtree, k, m);

                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                    if( qrtree->gettype(qrtree, k, m) == LIBHQR_KILLED_BY_TS ) {
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = A->nb;
                        T = TT;
                    }
                    for (n = 0; n < C->nt; n++) {
                        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                        node = C->get_rankof( C, m, n );
                        RUNTIME_data_migrate( sequence, C(p, n), node );
                        RUNTIME_data_migrate( sequence, C(m, n), node );

                        INSERT_TASK_ztpmlqt(
                            &options, side, trans,
                            tempmm, tempnn, tempkm, chameleon_min( L, tempnn ), ib, T->nb,
                            A(k, m),
                            T(k, m),
                            C(p, n),
                            C(m, n));
                    }
                    RUNTIME_data_flush( sequence, A(k, m) );
                    RUNTIME_data_flush( sequence, T(k, m) );
                }

                /* Restore the original location of the tiles */
                for (n = 0; n < C->nt; n++) {
                    RUNTIME_data_migrate( sequence, C(k, n),
                                          C->get_rankof( C, k, n ) );
                }

                RUNTIME_iteration_pop(chamctxt);
            }
        }
        /*
         *  ChamLeft / ChamConjTrans
         */
        else {
            for (k = KT-1; k >= 0; k--) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;

                /* Setting the order of the tiles*/
                nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

                for (i = nbtiles-1; i >= 0; i--) {
                    m = tiles[i];
                    p = qrtree->currpiv(qrtree, k, m);

                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                    if( qrtree->gettype(qrtree, k, m) == LIBHQR_KILLED_BY_TS ) {
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = A->nb;
                        T = TT;
                    }
                    for (n = 0; n < C->nt; n++) {
                        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                        node = C->get_rankof( C, m, n );
                        RUNTIME_data_migrate( sequence, C(p, n), node );
                        RUNTIME_data_migrate( sequence, C(m, n), node );

                        INSERT_TASK_ztpmlqt(
                            &options, side, trans,
                            tempmm, tempnn, tempkm, chameleon_min(L, tempnn), ib, T->nb,
                            A(k, m),
                            T(k, m),
                            C(p, n),
                            C(m, n));
                    }
                    RUNTIME_data_flush( sequence, A(k, m) );
                    RUNTIME_data_flush( sequence, T(k, m) );
                }

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    p = qrtree->getm(qrtree, k, i);

                    temppm   = p == C->mt-1 ? C->m-p*C->mb : C->mb;
                    tempkmin = chameleon_min( temppm, tempkm );


                    if ( genD ) {
                        int tempDpn = p == D->nt-1 ? D->n-p*D->nb : D->nb;

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamUpper, tempkmin, tempDpn, A->nb,
                            A(k, p),
                            D(k, p) );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamLower, tempkmin, tempDpn,
                            0., 1.,
                            D(k, p) );
#endif
                    }

                    for (n = 0; n < C->nt; n++) {
                        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                        RUNTIME_data_migrate( sequence, C(p, n),
                                              C->get_rankof( C, p, n ) );

                        INSERT_TASK_zunmlq(
                            &options, side, trans,
                            temppm, tempnn, tempkmin, ib, T->nb,
                            D(k, p),
                            T(k, p),
                            C(p, n));
                    }

                    RUNTIME_data_flush( sequence, D(k, p) );
                    RUNTIME_data_flush( sequence, T(k, p) );
                }
                RUNTIME_iteration_pop(chamctxt);
            }
        }
    }
    /*
     *  ChamRight / ChamNoTrans
     */
    else {
        if (trans == ChamNoTrans) {
            for (k = KT-1; k >= 0; k--) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;

                /* Setting the order of the tiles*/
                nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

                for (i = nbtiles-1; i >= 0; i--) {
                    n = tiles[i];
                    p = qrtree->currpiv(qrtree, k, n);

                    tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                    if( qrtree->gettype(qrtree, k, n) == LIBHQR_KILLED_BY_TS ) {
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = tempnn;
                        T = TT;
                    }

                    for (m = 0; m < C->mt; m++) {
                        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                        node = C->get_rankof( C, m, n );
                        RUNTIME_data_migrate( sequence, C(m, p), node );
                        RUNTIME_data_migrate( sequence, C(m, n), node );

                        INSERT_TASK_ztpmlqt(
                            &options, side, trans,
                            tempmm, tempnn, tempkm, L, ib, T->nb,
                            A(k, n),
                            T(k, n),
                            C(m, p),
                            C(m, n));
                    }
                    RUNTIME_data_flush( sequence, A(k, n) );
                    RUNTIME_data_flush( sequence, T(k, n) );
                }

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    p = qrtree->getm(qrtree, k, i);

                    temppn   = p == C->nt-1 ? C->n - p * C->nb : C->nb;
                    tempkmin = chameleon_min( temppn, tempkm );

                    if ( genD ) {
                        int tempDpn = p == D->nt-1 ? D->n-p*D->nb : D->nb;

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamUpper, tempkmin, tempDpn, A->nb,
                            A(k, p),
                            D(k, p) );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamLower, tempkmin, tempDpn,
                            0., 1.,
                            D(k, p) );
#endif
                    }

                    for (m = 0; m < C->mt; m++) {
                        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                        RUNTIME_data_migrate( sequence, C(m, p),
                                              C->get_rankof( C, m, p ) );

                        INSERT_TASK_zunmlq(
                            &options, side, trans,
                            tempmm, temppn, tempkmin, ib, T->nb,
                            D(k, p),
                            T(k, p),
                            C(m, p));
                    }
                    RUNTIME_data_flush( sequence, D(k, p) );
                    RUNTIME_data_flush( sequence, T(k, p) );
                }
                RUNTIME_iteration_pop(chamctxt);
            }
        }
        /*
         *  ChamRight / ChamConjTrans
         */
        else {
            for (k = 0; k < KT; k++) {
                RUNTIME_iteration_push(chamctxt, k);

                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;

                T = TS;
                for (i = 0; i < qrtree->getnbgeqrf(qrtree, k); i++) {
                    p = qrtree->getm(qrtree, k, i);

                    temppn   = p == C->nt - 1 ? C->n - p * C->nb : C->nb;
                    tempkmin = chameleon_min( temppn, tempkm );

                    if ( genD ) {
                        int tempDpn = p == D->nt-1 ? D->n-p*D->nb : D->nb;

                        INSERT_TASK_zlacpy(
                            &options,
                            ChamUpper, tempkmin, tempDpn, A->nb,
                            A(k, p),
                            D(k, p) );
#if defined(CHAMELEON_USE_CUDA)
                        INSERT_TASK_zlaset(
                            &options,
                            ChamLower, tempkmin, tempDpn,
                            0., 1.,
                            D(k, p) );
#endif
                    }

                    for (m = 0; m < C->mt; m++) {
                        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                        INSERT_TASK_zunmlq(
                            &options, side, trans,
                            tempmm, temppn, tempkmin, ib, T->nb,
                            D(k, p),
                            T(k, p),
                            C(m, p));
                    }
                    RUNTIME_data_flush( sequence, D(k, p) );
                    RUNTIME_data_flush( sequence, T(k, p) );
                }

                /* Setting the order of tiles */
                nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

                for (i = 0; i < nbtiles; i++) {
                    n = tiles[i];
                    p = qrtree->currpiv(qrtree, k, n);

                    tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                    if( qrtree->gettype(qrtree, k, n) == LIBHQR_KILLED_BY_TS ) {
                        /* TS kernel */
                        L = 0;
                        T = TS;
                    }
                    else {
                        /* TT kernel */
                        L = tempnn;
                        T = TT;
                    }

                    for (m = 0; m < C->mt; m++) {
                        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

                        node = C->get_rankof( C, m, n );
                        RUNTIME_data_migrate( sequence, C(m, p), node );
                        RUNTIME_data_migrate( sequence, C(m, n), node );

                        INSERT_TASK_ztpmlqt(
                            &options, side, trans,
                            tempmm, tempnn, tempkm, L, ib, T->nb,
                            A(k, n),
                            T(k, n),
                            C(m, p),
                            C(m, n));
                    }
                    RUNTIME_data_flush( sequence, A(k, n) );
                    RUNTIME_data_flush( sequence, T(k, n) );
                }

                /* Restore the original location of the tiles */
                for (m = 0; m < C->mt; m++) {
                    RUNTIME_data_migrate( sequence, C( m, k ),
                                          C->get_rankof( C, m, k ) );
                }

                RUNTIME_iteration_pop(chamctxt);
            }
        }
    }

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
