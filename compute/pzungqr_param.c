/**
 *
 * @file pzungqr_param.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zungqr_param parallel algorithm
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

#define A(m,n) A, m, n
#define Q(m,n) Q, m, n
#define T(m,n) T, m, n
#define D(m,n) D, m, n

/**
 *  Parallel construction of Q using tile V (application to identity) - dynamic scheduling
 */
void chameleon_pzungqr_param( int genD, int K,
                              const libhqr_tree_t *qrtree,
                              CHAM_desc_t *A, CHAM_desc_t *Q,
                              CHAM_desc_t *TS, CHAM_desc_t *TT, CHAM_desc_t *D,
                              RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t *T;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n, i, p, L;
    int tempmm, tempnn, tempkmin, tempkn;
    int ib, nbgeqrt, node, nbtiles, *tiles;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

    if (D == NULL) {
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
    /*
     * ztpmqrt = 3 * A->nb * ib
     */
    ws_worker = chameleon_max( ws_worker, ib * A->nb * 3 );
#endif

    ws_worker *= sizeof(CHAMELEON_Complex64_t);
    ws_host   *= sizeof(CHAMELEON_Complex64_t);

    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    /* Initialisation of temporary tiles array */
    tiles = (int*)calloc(qrtree->mt, sizeof(int));

    for (k = K-1; k >=0; k--) {
        RUNTIME_iteration_push(chamctxt, k);

        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

        /* Setting the order of tiles */
        nbtiles = libhqr_walk_stepk( qrtree, k, tiles );

        for (i = nbtiles-1; i >= 0; i--) {
            m = tiles[i];
            p = qrtree->currpiv(qrtree, k, m);

            tempmm = m == Q->mt-1 ? Q->m-m*Q->mb : Q->mb;

            if( qrtree->gettype(qrtree, k, m) == LIBHQR_KILLED_BY_TS ) {
                /* TS kernel */
                T = TS;
                L = 0;
            }
            else {
                /* TT kernel */
                T = TT;
                L = tempmm;
            }

            for (n = k; n < Q->nt; n++) {
                tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;

                node = Q->get_rankof( Q, m, n );
                RUNTIME_data_migrate( sequence, Q(p, n), node );
                RUNTIME_data_migrate( sequence, Q(m, n), node );

                INSERT_TASK_ztpmqrt(
                    &options,
                    ChamLeft, ChamNoTrans,
                    tempmm, tempnn, tempkn, L, ib, T->nb,
                    A(m, k),
                    T(m, k),
                    Q(p, n),
                    Q(m, n));
            }
            RUNTIME_data_flush( sequence, A(m, k) );
            RUNTIME_data_flush( sequence, T(m, k) );
        }

        T = TS;

        /* The number of geqrt to apply */
        nbgeqrt = qrtree->getnbgeqrf(qrtree, k);
        for (i = 0; i < nbgeqrt; i++) {
            m = qrtree->getm(qrtree, k, i);

            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            tempkmin = chameleon_min(tempmm, tempkn);

            if ( genD ) {
                int tempDmm = m == D->mt-1 ? D->m-m*D->mb : D->mb;
                INSERT_TASK_zlacpy(
                    &options,
                    ChamLower, tempDmm, tempkmin, A->nb,
                    A(m, k),
                    D(m, k) );
#if defined(CHAMELEON_USE_CUDA)
                INSERT_TASK_zlaset(
                    &options,
                    ChamUpper, tempDmm, tempkmin,
                    0., 1.,
                    D(m, k) );
#endif
            }

            for (n = k; n < Q->nt; n++) {
                tempnn = n == Q->nt-1 ? Q->n-n*Q->nb : Q->nb;

                /* Restore the original location of the tiles */
                RUNTIME_data_migrate( sequence, Q(m, n),
                                      Q->get_rankof( Q, m, n ) );

                INSERT_TASK_zunmqr(
                    &options,
                    ChamLeft, ChamNoTrans,
                    tempmm, tempnn, tempkmin, ib, T->nb,
                    D(m, k),
                    T(m, k),
                    Q(m, n));
            }
            RUNTIME_data_flush( sequence, D(m, k) );
            RUNTIME_data_flush( sequence, T(m, k) );
        }

        RUNTIME_iteration_pop(chamctxt);
    }

    free(tiles);
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
