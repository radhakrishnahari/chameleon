/**
 *
 * @file starpu/codelet_zbuild.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zbuild StarPU codelet
 *
 * @version 1.3.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Guillaume Sylvand
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zbuild_cpu_func(void *descr[], void *cl_arg)
{
    CHAM_tile_t *tileA;
    void *user_data;
    void (*user_build_callback)(int row_min, int row_max, int col_min, int col_max, void *buffer, int ld, void *user_data) ;
    int row_min, row_max, col_min, col_max;

    tileA = cti_interface_get(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &row_min, &row_max, &col_min, &col_max, &user_data, &user_build_callback );

    /* The callback 'user_build_callback' is expected to build the block of matrix [row_min, row_max] x [col_min, col_max]
     * (with both min and max values included in the intervals, index start at 0 like in C, NOT 1 like in Fortran)
     * and store it at the address 'buffer' with leading dimension 'ld'
     */
    user_build_callback(row_min, row_max, col_min, col_max, tileA->mat, tileA->ld, user_data);

}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zbuild, cl_zbuild_cpu_func)

void INSERT_TASK_zbuild( const RUNTIME_option_t *options,
                         const CHAM_desc_t *A, int Am, int An,
                         void *user_data, void* user_build_callback )
{

    struct starpu_codelet *codelet = &cl_zbuild;
    void (*callback)(void*) = options->profiling ? cl_zbuild_callback : NULL;
    int row_min, row_max, col_min, col_max;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_W(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    row_min = Am*A->mb ;
    row_max = Am == A->mt-1 ? A->m-1 : row_min+A->mb-1 ;
    col_min = An*A->nb ;
    col_max = An == A->nt-1 ? A->n-1 : col_min+A->nb-1 ;
    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,    &row_min,                      sizeof(int),
        STARPU_VALUE,    &row_max,                      sizeof(int),
        STARPU_VALUE,    &col_min,                      sizeof(int),
        STARPU_VALUE,    &col_max,                      sizeof(int),
        STARPU_W,         RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_VALUE,    &user_data,                    sizeof(void*),
        STARPU_VALUE,    &user_build_callback,          sizeof(void*),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        0);
}
