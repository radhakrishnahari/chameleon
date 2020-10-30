/**
 *
 * @file pzlauum.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlauum parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2020-03-03
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
/**
 *  Parallel UU' or L'L operation - dynamic scheduling
 */
void chameleon_pzlauum(cham_uplo_t uplo, CHAM_desc_t *A,
                          RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int k, m, n;
    int tempkm, tempkn;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);
    /*
     *  ChamLower
     */
    if (uplo == ChamLower) {
        for (k = 0; k < A->mt; k++) {
            tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
            for(n = 0; n < k; n++) {
                INSERT_TASK_zherk(
                    &options,
                    uplo, ChamConjTrans,
                    A->mb, tempkm, A->mb,
                    1.0, A(k, n),
                    1.0, A(n, n));

                for(m = n+1; m < k; m++) {
                    INSERT_TASK_zgemm(
                        &options,
                        ChamConjTrans, ChamNoTrans,
                        A->mb, A->nb, tempkm, A->mb,
                        1.0, A(k, m),
                             A(k, n),
                        1.0, A(m, n));
                }
            }
            for (n = 0; n < k; n++) {
                RUNTIME_data_flush( sequence, A(k, n) );
                INSERT_TASK_ztrmm(
                    &options,
                    ChamLeft, uplo, ChamConjTrans, ChamNonUnit,
                    tempkm, A->nb, A->mb,
                    1.0, A(k, k),
                         A(k, n));
            }
            RUNTIME_data_flush( sequence, A(k, k) );
            INSERT_TASK_zlauum(
                &options,
                uplo, tempkm, A->mb,
                A(k, k));
        }
    }
    /*
     *  ChamUpper
     */
    else {
        for (k = 0; k < A->mt; k++) {
            tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

            for (m = 0; m < k; m++) {
                INSERT_TASK_zherk(
                    &options,
                    uplo, ChamNoTrans,
                    A->mb, tempkn, A->mb,
                    1.0, A(m, k),
                    1.0, A(m, m));

                for (n = m+1; n < k; n++){
                    INSERT_TASK_zgemm(
                        &options,
                        ChamNoTrans, ChamConjTrans,
                        A->mb, A->nb, tempkn, A->mb,
                        1.0, A(m, k),
                             A(n, k),
                        1.0, A(m, n));
                }
            }
            for (m = 0; m < k; m++) {
                RUNTIME_data_flush( sequence, A(m, k) );
                INSERT_TASK_ztrmm(
                    &options,
                    ChamRight, uplo, ChamConjTrans, ChamNonUnit,
                    A->mb, tempkn, A->mb,
                    1.0, A(k, k),
                         A(m, k));
            }
            RUNTIME_data_flush( sequence, A(k, k) );
            INSERT_TASK_zlauum(
                &options,
                uplo, tempkn, A->mb,
                A(k, k));
        }
    }
    RUNTIME_options_finalize(&options, chamctxt);
}
