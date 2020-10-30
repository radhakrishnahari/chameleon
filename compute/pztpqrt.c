/**
 *
 * @file pztpqrt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2018 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon computational routines
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2020-03-03
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A, (m), (n)
#define B(m,n) B, (m), (n)
#define T(m,n) T, (m), (n)

/**
 *  Parallel tile QR factorization - dynamic scheduling
 */
void chameleon_pztpqrt( int L, CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *T,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    size_t ws_worker = 0;
    size_t ws_host = 0;

    int k, m, n;
    int tempkm, tempkn, tempnn, tempmm, templm;
    int ib;

    /* Dimension of the first column */
    int maxm  = chameleon_max( B->m - L, 1 );
    int maxmt = (maxm % B->mb == 0) ? (maxm / B->mb) : (maxm / B->mb + 1);

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;

    /*
     * ztpqrt  = A->nb * (ib+1)
     * ztpmqrt = A->nb * ib
     */
    ws_worker = A->nb * (ib+1);

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

    for (k = 0; k < A->nt; k++) {
        RUNTIME_iteration_push(chamctxt, k);

        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

        for (m = 0; m < maxmt; m++) {
            tempmm = m == B->mt-1 ? B->m-m*B->mb : B->mb;
            templm = ((L > 0) && (m == maxmt-1)) ? tempmm : 0;
            /* TT kernel */
            INSERT_TASK_ztpqrt(
                &options,
                tempmm, tempkn, templm, ib, T->nb,
                A(k, k),
                B(m, k),
                T(m, k) );

            for (n = k+1; n < B->nt; n++) {
                tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                INSERT_TASK_ztpmqrt(
                    &options,
                    ChamLeft, ChamConjTrans,
                    tempmm, tempnn, tempkm, templm, ib, T->nb,
                    B(m, k),
                    T(m, k),
                    A(k, n),
                    B(m, n) );
            }
        }

        maxmt = chameleon_min( B->mt, maxmt+1 );

        RUNTIME_iteration_pop(chamctxt);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
