/**
 *
 * @file quark/codelet_map3.c
 *
 * @copyright 2018-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map3 Quark codelet
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2020-03-03
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks.h"

void CORE_map3_quark(Quark *quark)
{
    const CHAM_desc_t *descA;
    const CHAM_desc_t *descB;
    const CHAM_desc_t *descC;
    cham_uplo_t uplo;
    int m;
    int n;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileC;
    cham_ternary_operator_t op_fct;
    void *op_args;

    quark_unpack_args_11( quark, uplo, m, n,
                          descA, tileA, descB, tileB, descC, tileC,
                          op_fct, op_args );
    op_fct( uplo, m, n,
            descA, tileA, descB, tileB, descC, tileC,
            op_args );
}

void INSERT_TASK_map3( const RUNTIME_option_t *options,
                       cham_uplo_t uplo,
                       cham_access_t accessA, cham_access_t accessB, cham_access_t accessC,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn,
                       const CHAM_desc_t *C, int Cm, int Cn,
                       cham_ternary_operator_t op_fct, void *op_args )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);

    QUARK_Insert_Task(
        opt->quark, CORE_map3_quark, (Quark_Task_Flags*)opt,
        sizeof(cham_uplo_t),              &uplo, VALUE,
        sizeof(int),                      &Am,   VALUE,
        sizeof(int),                      &An,   VALUE,
        sizeof(CHAM_desc_t*),             &A,    VALUE,
        sizeof(void*), RTBLKADDR(A, void, Am, An), cham_to_quark_access( accessA ),
        sizeof(CHAM_desc_t*),             &B,    VALUE,
        sizeof(void*), RTBLKADDR(B, void, Bm, Bn), cham_to_quark_access( accessB ),
        sizeof(CHAM_desc_t*),             &C,    VALUE,
        sizeof(void*), RTBLKADDR(C, void, Cm, Cn), cham_to_quark_access( accessC ),
        sizeof(cham_ternary_operator_t),   &op_fct,  VALUE,
        sizeof(void*),                    &op_args, VALUE,
        0);
}
