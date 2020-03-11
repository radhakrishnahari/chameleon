/**
 *
 * @file pmap2.c
 *
 * @copyright 2018-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map2 parallel algorithm
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2020-03-03
 *
 */
#include "control/common.h"

#define A(m, n) A, m, n
#define B(m, n) B, m, n

/**
 *  chameleon_pmap2
 */
void chameleon_pmap2( cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B,
                      cham_binary_operator_t op_fct, void *op_args,
                      RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    int m, n;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS)
        return;
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    switch( uplo ) {
    case ChamUpper:
        for (n = 0; n < A->nt; n++) {
            for (m = 0; m < n; m++) {
                INSERT_TASK_map2(
                    &options, ChamUpperLower,
                    A(m, n), B(m, n),
                    op_fct, op_args );
            }
            INSERT_TASK_map2(
                &options, uplo,
                A(n, n), B(n, n),
                op_fct, op_args );
        }
        break;

    case ChamLower:
        for (n = 0; n < A->nt; n++) {
            INSERT_TASK_map2(
                &options, uplo,
                A(n, n), B(n, n),
                op_fct, op_args );
            for (m = n+1; m < A->mt; m++) {
                INSERT_TASK_map2(
                    &options, ChamUpperLower,
                    A(m, n), B(m, n),
                    op_fct, op_args );
            }
        }
        break;

    case ChamUpperLower:
    default:
        for (m = 0; m < A->mt; m++) {
            for (n = 0; n < A->nt; n++) {
                INSERT_TASK_map2(
                    &options, uplo,
                    A(m, n), B(m, n),
                    op_fct, op_args );
            }
        }
    }

    RUNTIME_options_finalize(&options, chamctxt);
}
