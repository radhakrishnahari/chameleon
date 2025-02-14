/**
 *
 * @file quark/codelet_zbuild.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zbuild Quark codelet
 *
 * @version 1.3.0
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Guillaume Sylvand
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void CORE_zbuild_quark(Quark *quark)
{
    CHAM_tile_t *tileA;
    void *user_data;
    void (*user_build_callback)(int row_min, int row_max, int col_min, int col_max, void *buffer, int ld, void *user_data) ;
    int row_min, row_max, col_min, col_max;

    quark_unpack_args_7( quark, row_min, row_max, col_min, col_max, tileA, user_data, user_build_callback);

    user_build_callback(row_min, row_max, col_min, col_max, tileA->mat, tileA->ld, user_data);
}

void INSERT_TASK_zbuild( const RUNTIME_option_t *options,
                        const CHAM_desc_t *A, int Am, int An,
                        void *user_data, void* user_build_callback )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_BUILD;
    int row_min, row_max, col_min, col_max;
    row_min = Am*A->mb ;
    row_max = Am == A->mt-1 ? A->m-1 : row_min+A->mb-1 ;
    col_min = An*A->nb ;
    col_max = An == A->nt-1 ? A->n-1 : col_min+A->nb-1 ;

    QUARK_Insert_Task(opt->quark, CORE_zbuild_quark, (Quark_Task_Flags*)opt,
                      sizeof(int),                      &row_min,    VALUE,
                      sizeof(int),                      &row_max,    VALUE,
                      sizeof(int),                      &col_min,    VALUE,
                      sizeof(int),                      &col_max,    VALUE,
                      sizeof(void*), RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),         OUTPUT,
                      sizeof(void*),                    &user_data,  VALUE,
                      sizeof(void*),                    &user_build_callback,   VALUE,
                      0);
}
