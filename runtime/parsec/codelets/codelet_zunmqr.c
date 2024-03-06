/**
 *
 * @file parsec/codelet_zunmqr.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zunmqr PaRSEC codelet
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
CORE_zunmqr_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    cham_side_t side;
    cham_trans_t trans;
    int m;
    int n;
    int k;
    int ib;
    CHAMELEON_Complex64_t *A;
    int lda;
    CHAMELEON_Complex64_t *T;
    int ldt;
    CHAMELEON_Complex64_t *C;
    int ldc;
    CHAMELEON_Complex64_t *WORK;
    int ldwork;

    parsec_dtd_unpack_args(
        this_task, &side, &trans, &m, &n, &k, &ib, &A, &lda, &T, &ldt, &C, &ldc, &WORK, &ldwork );

    CORE_zunmqr( side, trans, m, n, k, ib,
                 A, lda, T, ldt, C, ldc, WORK, ldwork);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zunmqr(const RUNTIME_option_t *options,
                       cham_side_t side, cham_trans_t trans,
                       int m, int n, int k, int ib, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *T, int Tm, int Tn,
                       const CHAM_desc_t *C, int Cm, int Cn)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileT = T->get_blktile( T, Tm, Tn );
    CHAM_tile_t *tileC = C->get_blktile( C, Cm, Cn );

    parsec_dtd_insert_task(
        PARSEC_dtd_taskpool, CORE_zunmqr_parsec, options->priority, PARSEC_DEV_CPU, "unmqr",
        sizeof(cham_side_t),    &side,                              PARSEC_VALUE,
        sizeof(cham_trans_t),    &trans,                             PARSEC_VALUE,
        sizeof(int),           &m,                                 PARSEC_VALUE,
        sizeof(int),           &n,                                 PARSEC_VALUE,
        sizeof(int),           &k,                                 PARSEC_VALUE,
        sizeof(int),           &ib,                                PARSEC_VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | PARSEC_INPUT,
        sizeof(int), &(tileA->ld), PARSEC_VALUE,
        PASSED_BY_REF,         RTBLKADDR( T, CHAMELEON_Complex64_t, Tm, Tn ), chameleon_parsec_get_arena_index( T ) | PARSEC_INPUT,
        sizeof(int), &(tileT->ld), PARSEC_VALUE,
        PASSED_BY_REF,         RTBLKADDR( C, CHAMELEON_Complex64_t, Cm, Cn ), chameleon_parsec_get_arena_index( C ) | PARSEC_INOUT | PARSEC_AFFINITY,
        sizeof(int), &(tileC->ld), PARSEC_VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb,   NULL,                          PARSEC_SCRATCH,
        sizeof(int),           &nb,                                PARSEC_VALUE,
        PARSEC_DTD_ARG_END );
}
