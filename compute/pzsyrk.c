/**
 *
 * @file pzsyrk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyrk parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
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
#define C(m,n) C,  m,  n
/**
 *  Parallel tile symmetric rank-k update - dynamic scheduling
 */
void chameleon_pzsyrk(cham_uplo_t uplo, cham_trans_t trans,
                         CHAMELEON_Complex64_t alpha, CHAM_desc_t *A,
                         CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int m, n, k;
    int tempnn, tempmm, tempkn, tempkm;

    CHAMELEON_Complex64_t zbeta;
    CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    for (n = 0; n < C->nt; n++) {
        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
        /*
         *  ChamNoTrans
         */
        if (trans == ChamNoTrans) {
            for (k = 0; k < A->nt; k++) {
                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                zbeta = k == 0 ? beta : zone;
                INSERT_TASK_zsyrk(
                    &options,
                    uplo, trans,
                    tempnn, tempkn, A->mb,
                    alpha, A(n, k), /* ldan * K */
                    zbeta, C(n, n)); /* ldc  * N */
            }
            /*
             *  ChamNoTrans / ChamLower
             */
            if (uplo == ChamLower) {
                for (m = n+1; m < C->mt; m++) {
                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                    for (k = 0; k < A->nt; k++) {
                        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                        zbeta = k == 0 ? beta : zone;
                        INSERT_TASK_zgemm(
                            &options,
                            trans, ChamTrans,
                            tempmm, tempnn, tempkn, A->mb,
                            alpha, A(m, k),  /* ldam * K */
                                   A(n, k),  /* ldan * K */
                            zbeta, C(m, n)); /* ldc  * N */
                    }
                }
            }
            /*
             *  ChamNoTrans / ChamUpper
             */
            else {
                for (m = n+1; m < C->mt; m++) {
                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                    for (k = 0; k < A->nt; k++) {
                        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                        zbeta = k == 0 ? beta : zone;
                        INSERT_TASK_zgemm(
                            &options,
                            trans, ChamTrans,
                            tempnn, tempmm, tempkn, A->mb,
                            alpha, A(n, k),  /* ldan * K */
                                   A(m, k),  /* ldam * M */
                            zbeta, C(n, m)); /* ldc  * M */
                    }
                }
            }
        }
        /*
         *  ChamTrans
         */
        else {
            for (k = 0; k < A->mt; k++) {
                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                zbeta = k == 0 ? beta : zone;
                INSERT_TASK_zsyrk(
                    &options,
                    uplo, trans,
                    tempnn, tempkm, A->mb,
                    alpha, A(k, n),  /* lda * N */
                    zbeta, C(n, n)); /* ldc * N */
            }
            /*
             *  ChamTrans / ChamLower
             */
            if (uplo == ChamLower) {
                for (m = n+1; m < C->mt; m++) {
                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                    for (k = 0; k < A->mt; k++) {
                        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                        zbeta = k == 0 ? beta : zone;
                        INSERT_TASK_zgemm(
                            &options,
                            trans, ChamNoTrans,
                            tempmm, tempnn, tempkm, A->mb,
                            alpha, A(k, m),  /* lda * M */
                                   A(k, n),  /* lda * N */
                            zbeta, C(m, n)); /* ldc * N */
                    }
                }
            }
            /*
             *  ChamTrans / ChamUpper
             */
            else {
                for (m = n+1; m < C->mt; m++) {
                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                    for (k = 0; k < A->mt; k++) {
                        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                        zbeta = k == 0 ? beta : zone;
                        INSERT_TASK_zgemm(
                            &options,
                            trans, ChamNoTrans,
                            tempnn, tempmm, tempkm, A->mb,
                            alpha, A(k, n),  /* lda * K */
                                   A(k, m),  /* lda * M */
                            zbeta, C(n, m)); /* ldc * M */
                    }
                }
            }
        }
    }
    RUNTIME_options_finalize(&options, chamctxt);
}
