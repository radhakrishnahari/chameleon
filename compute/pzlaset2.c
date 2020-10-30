/**
 *
 * @file pzlaset2.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset2 parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
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
 *  Parallel initializztion a 2-D array A to
 *  ALPHA on the offdiagonals.
 */
void chameleon_pzlaset2(cham_uplo_t uplo, CHAMELEON_Complex64_t alpha,
                           CHAM_desc_t *A,
                           RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int i, j;
    int tempim;
    int tempjm, tempjn;
    int minmn = chameleon_min(A->mt, A->nt);

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }

    RUNTIME_options_init(&options, chamctxt, sequence, request);

    if (uplo == ChamLower) {
       for (j = 0; j < minmn; j++){
           tempjm = j == A->mt-1 ? A->m-j*A->mb : A->mb;
           tempjn = j == A->nt-1 ? A->n-j*A->nb : A->nb;
           INSERT_TASK_zlaset2(
               &options,
               ChamLower, tempjm, tempjn, alpha,
               A(j, j));

           for (i = j+1; i < A->mt; i++){
               tempim = i == A->mt-1 ? A->m-i*A->mb : A->mb;
               INSERT_TASK_zlaset2(
                   &options,
                   ChamUpperLower, tempim, tempjn, alpha,
                   A(i, j));
           }
       }
    }
    else if (uplo == ChamUpper) {
       for (j = 1; j < A->nt; j++){
           tempjn = j == A->nt-1 ? A->n-j*A->nb : A->nb;
           for (i = 0; i < chameleon_min(j, A->mt); i++){
               tempim = i == A->mt-1 ? A->m-i*A->mb : A->mb;
               INSERT_TASK_zlaset2(
                   &options,
                   ChamUpperLower, tempim, tempjn, alpha,
                   A(i, j));
           }
       }
       for (j = 0; j < minmn; j++){
           tempjm = j == A->mt-1 ? A->m-j*A->mb : A->mb;
           tempjn = j == A->nt-1 ? A->n-j*A->nb : A->nb;
           INSERT_TASK_zlaset2(
               &options,
               ChamUpper, tempjm, tempjn, alpha,
               A(j, j));
       }
    }
    else {
       for (i = 0; i < A->mt; i++){
           tempim = i == A->mt-1 ? A->m-i*A->mb : A->mb;
           for (j = 0; j < A->nt; j++){
               tempjn = j == A->nt-1 ? A->n-j*A->nb : A->nb;
               INSERT_TASK_zlaset2(
                   &options,
                   ChamUpperLower, tempim, tempjn, alpha,
                   A(i, j));
           }
       }
    }
    RUNTIME_options_finalize(&options, chamctxt);
}
