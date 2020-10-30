/**
 *
 * @file pzplghe.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplghe parallel algorithm
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Rade Mathis
 * @author Florent Pruvost
 * @date 2020-03-03
 * @precisions normal z -> c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
/**
 *  chameleon_pzplghe - Generate a random hermitian (positive definite if 'bump' is large enough) half-matrix by tiles.
 */
void chameleon_pzplghe( double bump, cham_uplo_t uplo, CHAM_desc_t *A,
                        unsigned long long int seed,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int m, n, minmn;
    int tempmm, tempnn;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    minmn = chameleon_min( A->mt, A->nt );
    switch ( uplo ) {
    case ChamLower:
        for (n = 0; n < minmn; n++) {
            tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

            for (m = n; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;

                options.priority = m + n;
                INSERT_TASK_zplghe(
                    &options,
                    bump, tempmm, tempnn, A(m, n),
                    A->m, m*A->mb, n*A->nb, seed );
            }
        }
        break;

    case ChamUpper:
        for (m = 0; m < minmn; m++) {
            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;

            for (n = m; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                options.priority = m + n;
                INSERT_TASK_zplghe(
                    &options,
                    bump, tempmm, tempnn, A(m, n),
                    A->m, m*A->mb, n*A->nb, seed );
            }
        }
        break;

    case ChamUpperLower:
    default:
        for (m = 0; m < A->mt; m++) {
            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;

            for (n = 0; n < A->nt; n++) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                options.priority = m + n;
                INSERT_TASK_zplghe(
                    &options,
                    bump, tempmm, tempnn, A(m, n),
                    A->m, m*A->mb, n*A->nb, seed );
            }
        }
    }
    RUNTIME_options_finalize(&options, chamctxt);
}
