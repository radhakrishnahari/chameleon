/**
 *
 * @file pzungqrrh.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zungqrrh parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2020-03-03
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n)  A, (m), (n)
#define Q(m,n)  Q, (m), (n)
#define T(m,n)  T, (m), (n)
#define T2(m,n) T, (m), (n)+(A->nt)
#define D(m,n)  D, (m), (n)

/**
 *  Parallel construction of Q using tile V (application to identity;
 *  reduction Householder) - dynamic scheduling
 */
void chameleon_pzungqrrh( int genD, int BS,
                          CHAM_desc_t *A, CHAM_desc_t *Q,
                          CHAM_desc_t *T,  CHAM_desc_t *D,
                          RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int K, M, RD, lastRD;
    int tempkn, tempMm, tempnn, tempmm, tempMRDm, tempkmin;
    int ib, node;

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
     * zunmqr  = A->nb * ib
     * ztpmqrt = A->nb * ib
     */
    ws_worker = A->nb * ib;

#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmqr  =     A->nb * ib
     * ztpmqrt = 3 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 3 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    K = chameleon_min(A->mt, A->nt);
    for (k = K-1; k >= 0; k--) {
        RUNTIME_iteration_push(chamctxt, k);

        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        lastRD = 0;
        for (RD = BS; RD < A->mt-k; RD *= 2)
            lastRD = RD;
        for (RD = lastRD; RD >= BS; RD /= 2) {
            for (M = k; M+RD < A->mt; M += 2*RD) {
                tempMRDm = M+RD == A->mt-1 ? A->m-(M+RD)*A->mb : A->mb;
                for (n = k; n < Q->nt; n++) {
                    tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;

                    node = Q->get_rankof( Q, M+RD, n );
                    RUNTIME_data_migrate( sequence, Q(M, n),    node );
                    RUNTIME_data_migrate( sequence, Q(M+RD, n), node );

                    /* TT kernel */
                    INSERT_TASK_ztpmqrt(
                        &options,
                        ChamLeft, ChamNoTrans,
                        tempMRDm, tempnn, tempkn, tempMRDm, ib, T->nb,
                        A (M+RD, k),
                        T2(M+RD, k),
                        Q (M,    n),
                        Q (M+RD, n));
                }

                RUNTIME_data_flush( sequence, A (M+RD, k) );
                RUNTIME_data_flush( sequence, T2(M+RD, k) );
            }
        }
        for (M = k; M < A->mt; M += BS) {
            tempMm   = M == A->mt-1 ? A->m-M*A->mb : A->mb;
            tempkmin = chameleon_min(tempMm, tempkn);
            for (m = chameleon_min(M+BS, A->mt)-1; m > M; m--) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;

                for (n = k; n < Q->nt; n++) {
                    tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;

                    node = Q->get_rankof( Q, m, n );
                    RUNTIME_data_migrate( sequence, Q(M, n), node );
                    RUNTIME_data_migrate( sequence, Q(m, n), node );

                    /* TS kernel */
                    INSERT_TASK_ztpmqrt(
                        &options,
                        ChamLeft, ChamNoTrans,
                        tempmm, tempnn, tempkn, 0, ib, T->nb,
                        A(m, k),
                        T(m, k),
                        Q(M, n),
                        Q(m, n));
                }
                RUNTIME_data_flush( sequence, A(m, k) );
                RUNTIME_data_flush( sequence, T(m, k) );
            }

            if ( genD ) {
                int tempDMm = M == D->mt-1 ? D->m-M*D->mb : D->mb;
                INSERT_TASK_zlacpy(
                    &options,
                    ChamLower, tempDMm, tempkmin, A->nb,
                    A(M, k),
                    D(M, k) );
#if defined(CHAMELEON_USE_CUDA)
                INSERT_TASK_zlaset(
                    &options,
                    ChamUpper, tempDMm, tempkmin,
                    0., 1.,
                    D(M, k) );
#endif
            }
            for (n = k; n < Q->nt; n++) {
                tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;

                /* Restore the original location of the tiles */
                RUNTIME_data_migrate( sequence, Q(M, n),
                                      Q->get_rankof( Q, M, n ) );

                INSERT_TASK_zunmqr(
                    &options,
                    ChamLeft, ChamNoTrans,
                    tempMm, tempnn,
                    tempkmin, ib, T->nb,
                    D(M, k),
                    T(M, k),
                    Q(M, n));
            }
            RUNTIME_data_flush( sequence, D(M, k) );
            RUNTIME_data_flush( sequence, T(M, k) );
        }
        RUNTIME_iteration_pop(chamctxt);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    (void)D;
}
