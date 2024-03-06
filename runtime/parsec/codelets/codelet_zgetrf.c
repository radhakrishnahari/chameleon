/**
 *
 * @file parsec/codelet_zgetrf.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf PaRSEC codelet
 *
 * @version 1.2.0
 * @author Reazul Hoque
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2022-02-22
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zgetrf_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int m;
    int n;
    CHAMELEON_Complex64_t *A;
    int lda;
    int *IPIV;
    cham_bool_t *check_info;
    int iinfo;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    int info;

    parsec_dtd_unpack_args(
        this_task, &m, &n, &A, &lda, &IPIV, &check_info, &iinfo, &sequence, &request );

    CORE_zgetrf( m, n, A, lda, IPIV, &info );

    if ( (sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, sequence, request, iinfo+info );
    }

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zgetrf(const RUNTIME_option_t *options,
                       int m, int n, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       int *IPIV,
                       cham_bool_t check_info, int iinfo)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

    parsec_dtd_insert_task(
        PARSEC_dtd_taskpool, CORE_zgetrf_parsec, options->priority, PARSEC_DEV_CPU, "getrf",
        sizeof(int),                 &m,                          PARSEC_VALUE,
        sizeof(int),                 &n,                          PARSEC_VALUE,
        PASSED_BY_REF,               RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | PARSEC_INOUT | PARSEC_AFFINITY,
        sizeof(int), &(tileA->ld), PARSEC_VALUE,
        sizeof(int)*nb,              IPIV,                        PARSEC_SCRATCH,
        sizeof(cham_bool_t),         &check_info,                 PARSEC_VALUE,
        sizeof(int),                 &iinfo,                      PARSEC_VALUE,
        sizeof(RUNTIME_sequence_t*), &(options->sequence),        PARSEC_VALUE,
        sizeof(RUNTIME_request_t*),  &(options->request),         PARSEC_VALUE,
        PARSEC_DTD_ARG_END );
}
