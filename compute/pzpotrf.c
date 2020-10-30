/**
 *
 * @file pzpotrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotrf parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2020-03-03
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
/**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 */
void chameleon_pzpotrf(cham_uplo_t uplo, CHAM_desc_t *A,
                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int k, m, n;
    int tempkm, tempmm, tempnn;
    size_t ws_host   = 0;

    CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t) 1.0;
    CHAMELEON_Complex64_t mzone = (CHAMELEON_Complex64_t)-1.0;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    RUNTIME_options_ws_alloc( &options, 0, ws_host );

    /*
     *  ChamLower
     */
    if (uplo == ChamLower) {
        for (k = 0; k < A->mt; k++) {
            RUNTIME_iteration_push(chamctxt, k);

            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;

            options.priority = 2*A->mt - 2*k;
            INSERT_TASK_zpotrf(
                &options,
                ChamLower, tempkm, A->mb,
                A(k, k), A->nb*k);

            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;

                options.priority = 2*A->mt - 2*k - m;
                INSERT_TASK_ztrsm(
                    &options,
                    ChamRight, ChamLower, ChamConjTrans, ChamNonUnit,
                    tempmm, A->mb, A->mb,
                    zone, A(k, k),
                          A(m, k));
            }
            RUNTIME_data_flush( sequence, A(k, k) );

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                options.priority = 2*A->mt - 2*k - n;
                INSERT_TASK_zherk(
                    &options,
                    ChamLower, ChamNoTrans,
                    tempnn, A->nb, A->mb,
                    -1.0, A(n, k),
                     1.0, A(n, n));

                for (m = n+1; m < A->mt; m++) {
                    tempmm = m == A->mt-1 ? A->m - m*A->mb : A->mb;

                    options.priority = 2*A->mt - 2*k - n - m;
                    INSERT_TASK_zgemm(
                        &options,
                        ChamNoTrans, ChamConjTrans,
                        tempmm, tempnn, A->mb, A->mb,
                        mzone, A(m, k),
                               A(n, k),
                        zone,  A(m, n));
                }
                RUNTIME_data_flush( sequence, A(n, k) );
            }
            RUNTIME_iteration_pop(chamctxt);
        }
    }
    /*
     *  ChamUpper
     */
    else {
        for (k = 0; k < A->nt; k++) {
            RUNTIME_iteration_push(chamctxt, k);

            tempkm = k == A->nt-1 ? A->n-k*A->nb : A->nb;

            options.priority = 2*A->nt - 2*k;
            INSERT_TASK_zpotrf(
                &options,
                ChamUpper,
                tempkm, A->mb,
                A(k, k), A->nb*k);

            for (n = k+1; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n - n*A->nb : A->nb;

                options.priority = 2*A->nt - 2*k - n;
                INSERT_TASK_ztrsm(
                    &options,
                    ChamLeft, ChamUpper, ChamConjTrans, ChamNonUnit,
                    A->mb, tempnn, A->mb,
                    zone, A(k, k),
                          A(k, n));
            }
            RUNTIME_data_flush( sequence, A(k, k) );

            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m - m*A->mb : A->mb;

                options.priority = 2*A->nt - 2*k  - m;
                INSERT_TASK_zherk(
                    &options,
                    ChamUpper, ChamConjTrans,
                    tempmm, A->mb, A->mb,
                    -1.0, A(k, m),
                     1.0, A(m, m));

                for (n = m+1; n < A->nt; n++) {
                    tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                    options.priority = 2*A->nt - 2*k - n - m;
                    INSERT_TASK_zgemm(
                        &options,
                        ChamConjTrans, ChamNoTrans,
                        tempmm, tempnn, A->mb, A->mb,
                        mzone, A(k, m),
                               A(k, n),
                        zone,  A(m, n));
                }
                RUNTIME_data_flush( sequence, A(k, m) );
            }

            RUNTIME_iteration_pop(chamctxt);
        }
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
