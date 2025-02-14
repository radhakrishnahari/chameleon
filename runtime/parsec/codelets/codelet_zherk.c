/**
 *
 * @file parsec/codelet_zherk.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherk PaRSEC codelet
 *
 * @version 1.3.0
 * @author Reazul Hoque
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zherk_parsec( parsec_execution_stream_t *context,
                   parsec_task_t             *this_task )
{
    cham_uplo_t uplo;
    cham_trans_t trans;
    int n;
    int k;
    double alpha;
    CHAMELEON_Complex64_t *A;
    int lda;
    double beta;
    CHAMELEON_Complex64_t *C;
    int ldc;

    parsec_dtd_unpack_args(
        this_task, &uplo, &trans, &n, &k, &alpha, &A, &lda, &beta, &C, &ldc );

    CORE_zherk( uplo, trans, n, k,
                alpha, A, lda,
                beta,  C, ldc);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zherk(const RUNTIME_option_t *options,
                      cham_uplo_t uplo, cham_trans_t trans,
                      int n, int k, int nb,
                      double alpha, const CHAM_desc_t *A, int Am, int An,
                      double beta, const CHAM_desc_t *C, int Cm, int Cn)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileC = C->get_blktile( C, Cm, Cn );

    parsec_dtd_insert_task(
        PARSEC_dtd_taskpool, CORE_zherk_parsec, options->priority, PARSEC_DEV_CPU, "herk",
        sizeof(cham_uplo_t),    &uplo,                             PARSEC_VALUE,
        sizeof(cham_trans_t),    &trans,                            PARSEC_VALUE,
        sizeof(int),           &n,                                PARSEC_VALUE,
        sizeof(int),           &k,                                PARSEC_VALUE,
        sizeof(double),        &alpha,                            PARSEC_VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | PARSEC_INPUT,
        sizeof(int), &(tileA->ld), PARSEC_VALUE,
        sizeof(double),        &beta,                             PARSEC_VALUE,
        PASSED_BY_REF,         RTBLKADDR( C, CHAMELEON_Complex64_t, Cm, Cn ), chameleon_parsec_get_arena_index( C ) | PARSEC_INOUT | PARSEC_AFFINITY,
        sizeof(int), &(tileC->ld), PARSEC_VALUE,
        PARSEC_DTD_ARG_END );

    (void)nb;
}

