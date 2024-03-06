/**
 *
 * @file parsec/codelet_zlange.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlange PaRSEC codelet
 *
 * @version 1.2.0
 * @author Reazul Hoque
 * @author Mathieu Faverge
 * @date 2022-02-22
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zlange_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_normtype_t norm;
    int M;
    int N;
    CHAMELEON_Complex64_t *A;
    int LDA;
    double *work;
    double *normA;

    parsec_dtd_unpack_args(
        this_task,   &norm,   &M,   &N, &A,   &LDA, &work, &normA );

    CORE_zlange( norm, M, N, A, LDA, work, normA );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zlange(const RUNTIME_option_t *options,
                       cham_normtype_t norm, int M, int N, int NB,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

    int szeW = chameleon_max( M, N );

    parsec_dtd_insert_task(
        PARSEC_dtd_taskpool, CORE_zlange_parsec, options->priority, PARSEC_DEV_CPU, "lange",
        sizeof(cham_normtype_t),            &norm,          PARSEC_VALUE,
        sizeof(int),                   &M,             PARSEC_VALUE,
        sizeof(int),                   &N,             PARSEC_VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | PARSEC_INPUT,
        sizeof(int), &(tileA->ld), PARSEC_VALUE,
        sizeof(double)*szeW,           NULL,           PARSEC_SCRATCH,
        PASSED_BY_REF,         RTBLKADDR( B, double, Bm, Bn ),            PARSEC_OUTPUT | PARSEC_AFFINITY,
        PARSEC_DTD_ARG_END );

    (void)NB;
}

static inline int
CORE_zlange_max_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    double *A;
    double *normA;

    parsec_dtd_unpack_args(
        this_task, &A, &normA );

    if ( *A > *normA )
        *normA = *A;

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zlange_max(const RUNTIME_option_t *options,
                           const CHAM_desc_t *A, int Am, int An,
                           const CHAM_desc_t *B, int Bm, int Bn)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_insert_task(
        PARSEC_dtd_taskpool, CORE_zlange_max_parsec, options->priority, PARSEC_DEV_CPU, "lange_max",
        PASSED_BY_REF,         RTBLKADDR( A, double, Am, An ), PARSEC_INPUT,
        PASSED_BY_REF,         RTBLKADDR( B, double, Bm, Bn ), PARSEC_OUTPUT | PARSEC_AFFINITY,
        PARSEC_DTD_ARG_END );
}
