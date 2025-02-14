/**
 *
 * @file parsec/codelet_ztsmqr_hetra1.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztsmqr_hetra1 PaRSEC codelet
 *
 * @version 1.3.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @author Azzam Haidar
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_ztsmqr_hetra1_parsec( parsec_execution_stream_t *context,
                           parsec_task_t             *this_task )
{
    cham_side_t side;
    cham_trans_t trans;
    int m1;
    int n1;
    int m2;
    int n2;
    int k;
    int ib;
    CHAMELEON_Complex64_t *A1;
    int lda1;
    CHAMELEON_Complex64_t *A2;
    int lda2;
    CHAMELEON_Complex64_t *V;
    int ldv;
    CHAMELEON_Complex64_t *T;
    int ldt;
    CHAMELEON_Complex64_t *WORK;
    int ldwork;

    parsec_dtd_unpack_args(
        this_task, &side, &trans, &m1, &n1, &m2, &n2, &k, &ib, &A1, &lda1, &A2, &lda2, &V, &ldv, &T, &ldt, &WORK, &ldwork);

    CORE_ztsmqr_hetra1( side, trans, m1, n1, m2, n2, k, ib,
                        A1, lda1, A2, lda2,
                        V, ldv, T, ldt,
                        WORK, ldwork);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_ztsmqr_hetra1(const RUNTIME_option_t *options,
                              cham_side_t side, cham_trans_t trans,
                              int m1, int n1, int m2, int n2, int k, int ib, int nb,
                              const CHAM_desc_t *A1, int A1m, int A1n,
                              const CHAM_desc_t *A2, int A2m, int A2n,
                              const CHAM_desc_t *V, int Vm, int Vn,
                              const CHAM_desc_t *T, int Tm, int Tn)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA1 = A1->get_blktile( A1, A1m, A1n );
    CHAM_tile_t *tileA2 = A2->get_blktile( A2, A2m, A2n );
    CHAM_tile_t *tileV = V->get_blktile( V, Vm, Vn );
    CHAM_tile_t *tileT = T->get_blktile( T, Tm, Tn );
    int ldwork = side == ChamLeft ? ib : nb;

    parsec_dtd_insert_task(
        PARSEC_dtd_taskpool, CORE_ztsmqr_hetra1_parsec, options->priority, PARSEC_DEV_CPU, "tsmqr_hetra1",
        sizeof(cham_side_t), &side,   PARSEC_VALUE,
        sizeof(cham_trans_t), &trans,  PARSEC_VALUE,
        sizeof(int),        &m1,     PARSEC_VALUE,
        sizeof(int),        &n1,     PARSEC_VALUE,
        sizeof(int),        &m2,     PARSEC_VALUE,
        sizeof(int),        &n2,     PARSEC_VALUE,
        sizeof(int),        &k,      PARSEC_VALUE,
        sizeof(int),        &ib,     PARSEC_VALUE,
        PASSED_BY_REF,       RTBLKADDR(A1, CHAMELEON_Complex64_t, A1m, A1n), PARSEC_INOUT,
        sizeof(int), &(tileA1->ld), PARSEC_VALUE,
        PASSED_BY_REF,       RTBLKADDR(A2, CHAMELEON_Complex64_t, A2m, A2n), PARSEC_INOUT | PARSEC_AFFINITY,
        sizeof(int), &(tileA2->ld), PARSEC_VALUE,
        PASSED_BY_REF,       RTBLKADDR(V,  CHAMELEON_Complex64_t, Vm,  Vn),  PARSEC_INPUT,
        sizeof(int), &(tileV->ld), PARSEC_VALUE,
        PASSED_BY_REF,       RTBLKADDR(T,  CHAMELEON_Complex64_t, Tm,  Tn),  PARSEC_INPUT,
        sizeof(int), &(tileT->ld), PARSEC_VALUE,
        sizeof(CHAMELEON_Complex64_t)*ib*nb, NULL, PARSEC_SCRATCH,
        sizeof(int),        &ldwork, PARSEC_VALUE,
        PARSEC_DTD_ARG_END );
}
