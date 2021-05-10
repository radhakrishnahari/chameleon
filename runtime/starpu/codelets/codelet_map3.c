/**
 *
 * @file starpu/codelet_map3.c
 *
 * @copyright 2018-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map StarPU codelet
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2020-03-03
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

CHAMELEON_CL_CB(map3, cti_handle_get_m(task->handles[2]), cti_handle_get_n(task->handles[2]), 0, M*N)

#if !defined(CHAMELEON_SIMULATION)
static void cl_map3_cpu_func(void *descr[], void *cl_arg)
{
    const CHAM_desc_t *descA;
    const CHAM_desc_t *descB;
    const CHAM_desc_t *descC;
    cham_uplo_t uplo;
    int m, n;
    CHAM_tile_t *tileA, *tileB, *tileC;
    cham_ternary_operator_t op_fct;
    void *op_args;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);
    tileC = cti_interface_get(descr[2]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &m, &n, &descA, &descB, &descC, &op_fct, &op_args );
    op_fct( uplo, m, n,
            descA, tileA,
            descB, tileB,
            descC, tileC,
            op_args );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(map3, cl_map3_cpu_func)

void INSERT_TASK_map3( const RUNTIME_option_t *options,
                       cham_uplo_t uplo,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn,
                       const CHAM_desc_t *C, int Cm, int Cn,
                       cham_ternary_operator_t op_fct, void *op_args,
                       const char *name )
{

    struct starpu_codelet *codelet = &cl_map3;
    void (*callback)(void*) = options->profiling ? cl_map3_callback : NULL;
    starpu_option_request_t* schedopt = (starpu_option_request_t *)(options->request->schedopt);
    int workerid = (schedopt == NULL) ? -1 : schedopt->workerid;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R( A, Am, An);
    CHAMELEON_ACCESS_R( B, Bm, Bn);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,    &uplo,                   sizeof(cham_uplo_t),
        STARPU_VALUE,    &Am,                     sizeof(int),
        STARPU_VALUE,    &An,                     sizeof(int),
        STARPU_VALUE,    &A,                      sizeof(CHAM_desc_t*),
        STARPU_R,         RTBLKADDR(A, void, Am, An),
        STARPU_VALUE,    &B,                      sizeof(CHAM_desc_t*),
        STARPU_R,         RTBLKADDR(B, void, Bm, Bn),
        STARPU_VALUE,    &C,                      sizeof(CHAM_desc_t*),
        STARPU_RW,        RTBLKADDR(C, void, Cm, Cn),
        STARPU_VALUE,    &op_fct,                 sizeof(cham_ternary_operator_t),
        STARPU_VALUE,    &op_args,                sizeof(void*),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
        STARPU_EXECUTE_ON_WORKER, workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,      name == NULL ? "map3" : name,
#endif
        0);
}
