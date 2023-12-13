/**
 *
 * @file parsec/codelet_map3.c
 *
 * @copyright 2018-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map3 PaRSEC codelet
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2020-03-03
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks.h"

static inline int
CORE_map3_parsec( parsec_execution_stream_t *context,
                  parsec_task_t             *this_task )
{
    const CHAM_desc_t *descA;
    const CHAM_desc_t *descB;
    const CHAM_desc_t *descC;
    cham_uplo_t uplo;
    int m;
    int n;
    void *dataA;
    void *dataB;
    void *dataC;
    cham_ternary_operator_t op_fct;
    void *op_args;

    parsec_dtd_unpack_args(
        this_task, &uplo, &m, &n,
        &descA, &dataA, &descB, &dataB, &descC, &dataC,
        &op_fct, &op_args );
    op_fct( uplo, m, n, descA, dataA, descB, dataB, descC, dataC, op_args );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_map3( const RUNTIME_option_t *options,
                       cham_uplo_t uplo,
                       cham_access_t accessA, cham_access_t accessB, cham_access_t accessC,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn,
                       const CHAM_desc_t *C, int Cm, int Cn,
                       cham_ternary_operator_t op_fct, void *op_args )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    int parsec_accessA = cham_to_parsec_access( accessA );
    int parsec_accessB = cham_to_parsec_access( accessB );
    int parsec_accessC = cham_to_parsec_access( accessC );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_map3_parsec, options->priority, "map3",
        sizeof(cham_uplo_t),              &uplo, VALUE,
        sizeof(int),                      &Am,   VALUE,
        sizeof(int),                      &An,   VALUE,
        sizeof(CHAM_desc_t*),             &A,    VALUE,
        PASSED_BY_REF, RTBLKADDR(A, void, Am, An), chameleon_parsec_get_arena_index( A ) | parsec_accessA,
        sizeof(CHAM_desc_t*),             &B,    VALUE,
        PASSED_BY_REF, RTBLKADDR(B, void, Bm, Bn), chameleon_parsec_get_arena_index( B ) | parsec_accessB,
        sizeof(CHAM_desc_t*),             &C,    VALUE,
        PASSED_BY_REF, RTBLKADDR(C, void, Cm, Cn), chameleon_parsec_get_arena_index( C ) | parsec_accessC,
        sizeof(cham_ternary_operator_t),  &op_fct,  VALUE,
        sizeof(void*),                    &op_args, VALUE,
        PARSEC_DTD_ARG_END );
}
