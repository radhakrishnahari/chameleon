/**
 *
 * @file pzungqr.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zungqr parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2020-03-03
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define Q(m,n) Q,  m,  n
#define T(m,n) T,  m,  n
#define D(k)   D,  k,  k

/**
 *  Parallel construction of Q using tile V (application to identity) - dynamic scheduling
 */
void chameleon_pzungqr( int genD, CHAM_desc_t *A, CHAM_desc_t *Q,
                        CHAM_desc_t *T, CHAM_desc_t *D,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int tempmm, tempnn, tempkmin, tempkm;
    int tempAkm, tempAkn;
    int ib, minMT;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

    if (A->m > A->n) {
        minMT = A->nt;
    } else {
        minMT = A->mt;
    }

    if ( D == NULL ) {
        D    = A;
        genD = 0;
    }

    /*
     * zunmqr  = A->nb * ib
     * ztpmqrt = A->nb * ib
     */
    ws_worker = A->nb * ib;

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmqr  =     A->nb * ib
     * ztpmqrt = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    for (k = minMT-1; k >= 0; k--) {
        RUNTIME_iteration_push(chamctxt, k);

        tempAkm  = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        tempAkn  = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        tempkmin = chameleon_min( tempAkn, tempAkm );
        tempkm   = k == Q->mt-1 ? Q->m-k*Q->mb : Q->mb;
        for (m = Q->mt - 1; m > k; m--) {
            tempmm = m == Q->mt-1 ? Q->m-m*Q->mb : Q->mb;
            for (n = k; n < Q->nt; n++) {
                tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;

                RUNTIME_data_migrate( sequence, Q(k, n),
                                      Q->get_rankof( Q, m, n ) );

                /* TS kernel */
                INSERT_TASK_ztpmqrt(
                    &options,
                    ChamLeft, ChamNoTrans,
                    tempmm, tempnn, tempAkn, 0, ib, T->nb,
                    A(m, k),
                    T(m, k),
                    Q(k, n),
                    Q(m, n));
            }
            RUNTIME_data_flush( sequence, A(m, k) );
            RUNTIME_data_flush( sequence, T(m, k) );
        }

        if ( genD ) {
            int tempDkm = k == D->mt-1 ? D->m-k*D->mb : D->mb;

            INSERT_TASK_zlacpy(
                &options,
                ChamLower, tempDkm, tempkmin, A->nb,
                A(k, k),
                D(k) );
#if defined(CHAMELEON_USE_CUDA)
            INSERT_TASK_zlaset(
                &options,
                ChamUpper, tempDkm, tempkmin,
                0., 1.,
                D(k) );
#endif
        }
        for (n = k; n < Q->nt; n++) {
            tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;

            /* Restore the original location of the tiles */
            RUNTIME_data_migrate( sequence, Q(k, n),
                                  Q->get_rankof( Q, k, n ) );

            INSERT_TASK_zunmqr(
                &options,
                ChamLeft, ChamNoTrans,
                tempkm, tempnn, tempkmin, ib, T->nb,
                D(k),
                T(k, k),
                Q(k, n));
        }
        RUNTIME_data_flush( sequence, D(k) );
        RUNTIME_data_flush( sequence, T(k, k) );

        RUNTIME_iteration_pop(chamctxt);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
