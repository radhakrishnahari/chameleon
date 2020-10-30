/**
 *
 * @file pztrsmpl.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrsmpl parallel algorithm
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
#define B(m,n) B,  m,  n
#define L(m,n) L,  m,  n
#define IPIV(m,n) &(IPIV[(int64_t)A->nb*((int64_t)(m)+(int64_t)A->mt*(int64_t)(n))])
/**
 *  Parallel forward substitution for tile LU - dynamic scheduling
 */
void chameleon_pztrsmpl( CHAM_desc_t *A, CHAM_desc_t *B, CHAM_desc_t *L, int *IPIV,
                     RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int k, m, n;
    int tempkm, tempnn, tempkmin, tempmm, tempkn;
    int ib;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    ib = CHAMELEON_IB;
    for (k = 0; k < chameleon_min(A->mt, A->nt); k++) {
        tempkm   = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        tempkn   = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        tempkmin = k == chameleon_min(A->mt, A->nt)-1 ? chameleon_min(A->m, A->n)-k*A->mb : A->mb;
        for (n = 0; n < B->nt; n++) {
            tempnn = n == B->nt-1 ? B->n-n*B->nb : B->nb;
            INSERT_TASK_zgessm(
                &options,
                tempkm, tempnn, tempkmin, ib, L->nb,
                IPIV(k, k),
                L(k, k),
                A(k, k),
                B(k, n));
        }
        for (m = k+1; m < A->mt; m++) {
            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            for (n = 0; n < B->nt; n++) {
                tempnn  = n == B->nt-1 ? B->n-n*B->nb : B->nb;
                INSERT_TASK_zssssm(
                    &options,
                    A->nb, tempnn, tempmm, tempnn, tempkn, ib, L->nb,
                    B(k, n),
                    B(m, n),
                    L(m, k),
                    A(m, k),
                    IPIV(m, k));
            }
        }
    }
    RUNTIME_options_finalize(&options, chamctxt);
}
