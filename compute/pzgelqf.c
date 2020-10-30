/**
 *
 * @file pzgelqf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgelqf parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2020-03-03
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define T(m,n) T,  m,  n
#define D(k)   D,  k,  k

/**
 *  Parallel tile LQ factorization - dynamic scheduling
 */
void chameleon_pzgelqf( int genD, CHAM_desc_t *A, CHAM_desc_t *T, CHAM_desc_t *D,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int tempkm, tempkn, tempmm, tempnn;
    int ib, minMNT;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

    if (A->m > A->n) {
        minMNT = A->nt;
    } else {
        minMNT = A->mt;
    }

    if ( D == NULL ) {
        D    = A;
        genD = 0;
    }

    /*
     * zgelqt  = A->nb * (ib+1)
     * zunmlq  = A->nb * ib
     * ztplqt  = A->nb * (ib+1)
     * ztpmlqt = A->nb * ib
     */
    ws_worker = A->nb * (ib+1);

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /* Worker space
     *
     * zunmlq  =     A->nb * ib
     * ztpmlqt = 2 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 2 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    for (k = 0; k < minMNT; k++) {
        RUNTIME_iteration_push(chamctxt, k);

        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        INSERT_TASK_zgelqt(
            &options,
            tempkm, tempkn, ib, T->nb,
            A(k, k),
            T(k, k));
        if ( genD ) {
            int tempDkm = k == D->mt-1 ? D->m-k*D->mb : D->mb;
            int tempDkn = k == D->nt-1 ? D->n-k*D->nb : D->nb;
            INSERT_TASK_zlacpy(
                &options,
                ChamUpper, tempDkm, tempDkn, A->nb,
                A(k, k),
                D(k) );
#if defined(CHAMELEON_USE_CUDA)
            INSERT_TASK_zlaset(
                &options,
                ChamLower, tempDkm, tempDkn,
                0., 1.,
                D(k) );
#endif
        }
        for (m = k+1; m < A->mt; m++) {
            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            INSERT_TASK_zunmlq(
                &options,
                ChamRight, ChamConjTrans,
                tempmm, tempkn, tempkn, ib, T->nb,
                D(k),
                T(k, k),
                A(m, k));
        }
        RUNTIME_data_flush( sequence, D(k)    );
        RUNTIME_data_flush( sequence, T(k, k) );

        for (n = k+1; n < A->nt; n++) {
            tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

            RUNTIME_data_migrate( sequence, A(k, k),
                                  A->get_rankof( A, k, n ) );

            /* TS kernel */
            INSERT_TASK_ztplqt(
                &options,
                tempkm, tempnn, 0, ib, T->nb,
                A(k, k),
                A(k, n),
                T(k, n));
            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;

                RUNTIME_data_migrate( sequence, A(m, k),
                                      A->get_rankof( A, m, n ) );

                INSERT_TASK_ztpmlqt(
                    &options,
                    ChamRight, ChamConjTrans,
                    tempmm, tempnn, A->mb, 0, ib, T->nb,
                    A(k, n),
                    T(k, n),
                    A(m, k),
                    A(m, n));
            }
            RUNTIME_data_flush( sequence, A(k, n) );
            RUNTIME_data_flush( sequence, T(k, n) );
        }

        /* Restore the original location of the tiles */
        for (m = k; m < A->mt; m++) {
            RUNTIME_data_migrate( sequence, A(m, k),
                                  A->get_rankof( A, m, k ) );
        }

        RUNTIME_iteration_pop(chamctxt);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    (void)D;
}
