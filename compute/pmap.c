/**
 *
 * @file pmap.c
 *
 * @copyright 2018-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map parallel algorithm
 *
 * @version 1.2.0
 * @author Mathieu Faverge
 * @date 2022-02-22
 *
 */
#include "control/common.h"

#define A(m, n) A, m, n

/**
 *  chameleon_pmap
 */
void chameleon_pmap( cham_access_t access, cham_uplo_t uplo, CHAM_desc_t *A,
                     cham_unary_operator_t op_fct, void *op_args,
                     RUNTIME_sequence_t *sequence, RUNTIME_request_t *request, const char *name )
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
                INSERT_TASK_map(
                    &options, ChamUpperLower,
                    access, A(m, n),
                    op_fct, op_args, name );
            }
            INSERT_TASK_map(
                &options, uplo,
                access, A(n, n),
                op_fct, op_args, name );
        }
        break;

    case ChamLower:
        for (n = 0; n < A->nt; n++) {
            INSERT_TASK_map(
                &options, uplo,
                access, A(n, n),
                op_fct, op_args, name );
            for (m = n+1; m < A->mt; m++) {
                INSERT_TASK_map(
                    &options, ChamUpperLower,
                    access, A(m, n),
                    op_fct, op_args, name );
            }
        }
        break;

    case ChamUpperLower:
    default:
        for (m = 0; m < A->mt; m++) {
            for (n = 0; n < A->nt; n++) {
                INSERT_TASK_map(
                    &options, uplo,
                    access, A(m, n),
                    op_fct, op_args, name );
            }
        }
    }

    RUNTIME_options_finalize(&options, chamctxt);
}
