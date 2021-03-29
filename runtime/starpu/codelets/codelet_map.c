/**
 *
 * @file starpu/codelet_map.c
 *
 * @copyright 2018-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map StarPU codelet
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2023-07-06
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

CHAMELEON_CL_CB(map, cti_handle_get_m(task->handles[0]), cti_handle_get_n(task->handles[0]), 0, M*N)

#if !defined(CHAMELEON_SIMULATION)
static void cl_map_cpu_func(void *descr[], void *cl_arg)
{
    const CHAM_desc_t *desc;
    cham_uplo_t uplo;
    int m;
    int n;
    CHAM_tile_t *tile;
    cham_unary_operator_t op_fct;
    void *op_args;

    tile = cti_interface_get(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &desc, &uplo, &m, &n, &op_fct, &op_args );
    op_fct( desc, uplo, m, n, tile, op_args );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(map, cl_map_cpu_func)

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_access_t accessA, cham_uplo_t uplo, const CHAM_desc_t *A, int Am, int An,
                      cham_unary_operator_t op_fct, void *op_args,
                      const char *name )
{

    struct starpu_codelet *codelet = &cl_map;
    void (*callback)(void*) = options->profiling ? cl_map_callback : NULL;
    char                  *cl_name = (name == NULL) ? "map" : name;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    cl_name = chameleon_codelet_name( cl_name, 1,
                                      A->get_blktile( A, Am, An ) );

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,    &A,                      sizeof(CHAM_desc_t*),
        STARPU_VALUE,    &uplo,                   sizeof(cham_uplo_t),
        STARPU_VALUE,    &Am,                     sizeof(int),
        STARPU_VALUE,    &An,                     sizeof(int),
        cham_to_starpu_access(accessA), RTBLKADDR(A, ChamByte, Am, An),
        STARPU_VALUE,    &op_fct,                 sizeof(cham_unary_operator_t),
        STARPU_VALUE,    &op_args,                sizeof(void*),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, cl_name,
#endif
        0);
}
