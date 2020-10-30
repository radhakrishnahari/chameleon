/**
 *
 * @file pzgeqrf_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeqrf_param parallel algorithm
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
#include "libhqr.h"

#define A(m,n) A, (m), (n)
#define T(m,n) T, (m), (n)
#define D(m,n) D, (m), (n)

/**
 *  Parallel tile QR factorization (reduction Householder) - dynamic scheduling
 *
 * @param[in] genD
 *         Indicate if the copies of the geqrt tiles must be done to speedup
 *         computations in updates.
 */
void chameleon_pzgeqrf_param( int genD, int K,
                              const libhqr_tree_t *qrtree, CHAM_desc_t *A,
                              CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t *T;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n, i, p;
    int L, nbgeqrt;
    int tempkmin, tempkn, tempnn, tempmm;
    int ib, node, nbtiles, *tiles;

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
     * zgeqrt  = A->nb * (ib+1)
     * zunmqr  = A->nb * ib
     * ztpqrt  = A->nb * (ib+1)
     * ztpmqrt = A->nb * ib
     */
    ws_worker = A->nb * (ib+1);

    /* Allocation of temporary (scratch) working space */
#if defined(CHAMELEON_USE_CUDA)
    /*
     * zunmqr  =     A->nb * ib
     * ztpmqrt = 3 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 3 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    /* Initialisation of temporary tiles array */
    tiles = (int*)calloc(qrtree->mt, sizeof(int));

    for (k = 0; k < K; k++) {
        RUNTIME_iteration_push(chamctxt, k);
        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

        /* The number of geqrt to apply */
        nbgeqrt = qrtree->getnbgeqrf(qrtree, k);

        T = TS;
        for (i = 0; i < nbgeqrt; i++) {
            m = qrtree->getm(qrtree, k, i);
            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            tempkmin = chameleon_min(tempmm, tempkn);

            INSERT_TASK_zgeqrt(
                &options,
                tempmm, tempkn, ib, T->nb,
                A(m, k),
                T(m, k));

            if ( genD ) {
                int tempDmm = m == D->mt-1 ? D->m-m*D->mb : D->mb;
                int tempDkn = k == D->nt-1 ? D->n-k*D->nb : D->nb;

                INSERT_TASK_zlacpy(
                    &options,
                    ChamLower, tempDmm, tempDkn, A->nb,
                    A(m, k),
                    D(m, k) );
#if defined(CHAMELEON_USE_CUDA)
                INSERT_TASK_zlaset(
                    &options,
                    ChamUpper, tempDmm, tempDkn,
                    0., 1.,
                    D(m, k) );
#endif
            }

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
                INSERT_TASK_zunmqr(
                    &options,
                    ChamLeft, ChamConjTrans,
                    tempmm, tempnn, tempkmin, ib, T->nb,
                    D(m, k),
                    T(m, k),
                    A(m, n));
            }
            RUNTIME_data_flush( sequence, D(m, k) );
            RUNTIME_data_flush( sequence, T(m, k) );
        }

        /* Setting the order of the tiles */
        nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

        for (i = 0; i < nbtiles; i++) {
            m = tiles[i];
            p = qrtree->currpiv(qrtree, k, m);

            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;

            if ( qrtree->gettype(qrtree, k, m) == LIBHQR_KILLED_BY_TS ) {
                /* TS kernel */
                T = TS;
                L = 0;
            }
            else {
                /* TT kernel */
                T = TT;
                L = tempmm;
            }

            node = A->get_rankof( A, m, k );
            RUNTIME_data_migrate( sequence, A(p, k), node );
            RUNTIME_data_migrate( sequence, A(m, k), node );

            INSERT_TASK_ztpqrt(
                &options,
                tempmm, tempkn, chameleon_min(L, tempkn), ib, T->nb,
                A(p, k),
                A(m, k),
                T(m, k));

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                node = A->get_rankof( A, m, n );
                RUNTIME_data_migrate( sequence, A(p, n), node );
                RUNTIME_data_migrate( sequence, A(m, n), node );

                INSERT_TASK_ztpmqrt(
                    &options,
                    ChamLeft, ChamConjTrans,
                    tempmm, tempnn, A->nb, L, ib, T->nb,
                    A(m, k),
                    T(m, k),
                    A(p, n),
                    A(m, n));
            }
            RUNTIME_data_flush( sequence, A(m, k) );
            RUNTIME_data_flush( sequence, T(m, k) );
        }

        /* Restore the original location of the tiles */
        for (n = k; n < A->nt; n++) {
            RUNTIME_data_migrate( sequence, A(k, n),
                                  A->get_rankof( A, k, n ) );
        }

        RUNTIME_iteration_pop(chamctxt);
    }

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
