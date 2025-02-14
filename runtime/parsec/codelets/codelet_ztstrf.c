/**
 *
 * @file parsec/codelet_ztstrf.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztstrf PaRSEC codelet
 *
 * @version 1.3.0
 * @author Reazul Hoque
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_ztstrf_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int m;
    int n;
    int ib;
    int nb;
    CHAMELEON_Complex64_t *U;
    int ldu;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex64_t *L;
    int ldl;
    int *IPIV;
    CHAMELEON_Complex64_t *WORK;
    int ldwork;
    cham_bool_t *check_info;
    int iinfo;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    int info;

    parsec_dtd_unpack_args(
        this_task, &m, &n, &ib, &nb, &U, &ldu, &A, &lda, &L, &ldl, &IPIV, &WORK, &ldwork, &check_info, &iinfo, &sequence, &request );

    CORE_ztstrf( m, n, ib, nb, U, ldu, A, lda, L, ldl, IPIV, WORK, ldwork, &info );

    if ( (sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, sequence, request, iinfo+info );
    }

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_ztstrf(const RUNTIME_option_t *options,
                       int m, int n, int ib, int nb,
                       const CHAM_desc_t *U, int Um, int Un,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *L, int Lm, int Ln,
                       int *IPIV,
                       cham_bool_t check_info, int iinfo)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileU = U->get_blktile( U, Um, Un );
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileL = L->get_blktile( L, Lm, Ln );

    parsec_dtd_insert_task(
        PARSEC_dtd_taskpool, CORE_ztstrf_parsec, options->priority, PARSEC_DEV_CPU, "tstrf",
        sizeof(int),                 &m,                                PARSEC_VALUE,
        sizeof(int),                 &n,                                PARSEC_VALUE,
        sizeof(int),                 &ib,                               PARSEC_VALUE,
        sizeof(int),                 &nb,                               PARSEC_VALUE,
        PASSED_BY_REF,               RTBLKADDR( U, CHAMELEON_Complex64_t, Um, Un ), chameleon_parsec_get_arena_index( U ) | PARSEC_INOUT,
        sizeof(int), &(tileU->ld), PARSEC_VALUE,
        PASSED_BY_REF,               RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | PARSEC_INOUT | PARSEC_AFFINITY,
        sizeof(int), &(tileA->ld), PARSEC_VALUE,
        PASSED_BY_REF,               RTBLKADDR( L, CHAMELEON_Complex64_t, Lm, Ln ), chameleon_parsec_get_arena_index( L ) | PARSEC_OUTPUT,
        sizeof(int), &(tileL->ld), PARSEC_VALUE,
        sizeof(int*),                &IPIV,                             PARSEC_VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb,    NULL,                 PARSEC_SCRATCH,
        sizeof(int),                 &nb,                               PARSEC_VALUE,
        sizeof(int),                 &check_info,                       PARSEC_VALUE,
        sizeof(int),                 &iinfo,                            PARSEC_VALUE,
        sizeof(RUNTIME_sequence_t*), &(options->sequence),              PARSEC_VALUE,
        sizeof(RUNTIME_request_t*),  &(options->request),               PARSEC_VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}
