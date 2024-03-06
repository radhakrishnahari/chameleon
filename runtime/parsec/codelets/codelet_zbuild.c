/**
 *
 * @file parsec/codelet_zbuild.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zbuild PaRSEC codelet
 *
 * @version 1.2.0
 * @author Reazul Hoque
 * @author Guillaume Sylvand
 * @author Mathieu Faverge
 * @date 2022-02-22
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zbuild_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    CHAMELEON_Complex64_t *A;
    int lda;
    void *user_data;
    void (*user_build_callback)( int row_min, int row_max, int col_min, int col_max,
                                 void *buffer, int ld, void *user_data );
    int row_min, row_max, col_min, col_max;

    parsec_dtd_unpack_args(
        this_task, &row_min, &row_max, &col_min, &col_max, &A, &lda, &user_data, &user_build_callback );

    user_build_callback(row_min, row_max, col_min, col_max, A, lda, user_data);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zbuild( const RUNTIME_option_t *options,
                        const CHAM_desc_t *A, int Am, int An,
                        void *user_data, void* user_build_callback )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    int row_min, row_max, col_min, col_max;
    row_min = Am*A->mb ;
    row_max = Am == A->mt-1 ? A->m-1 : row_min+A->mb-1 ;
    col_min = An*A->nb ;
    col_max = An == A->nt-1 ? A->n-1 : col_min+A->nb-1 ;

    parsec_dtd_insert_task(
        PARSEC_dtd_taskpool, CORE_zbuild_parsec, options->priority, PARSEC_DEV_CPU, "zbuild",
        sizeof(int),   &row_min,                          PARSEC_VALUE,
        sizeof(int),   &row_max,                          PARSEC_VALUE,
        sizeof(int),   &col_min,                          PARSEC_VALUE,
        sizeof(int),   &col_max,                          PARSEC_VALUE,
        PASSED_BY_REF,  RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | PARSEC_OUTPUT | PARSEC_AFFINITY,
        sizeof(int), &(tileA->ld), PARSEC_VALUE,
        sizeof(void*), &user_data,                        PARSEC_VALUE,
        sizeof(void*), &user_build_callback,              PARSEC_VALUE,
        PARSEC_DTD_ARG_END );
}
