/**
 *
 * @file starpu/codelet_zgesum.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgesum StarPU codelet
 *
 * @version 1.3.0
 * @author Florent Pruvost
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zgesum_args_s {
    cham_store_t storev;
    int m;
    int n;
};

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgesum_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zgesum_args_s *clargs = (struct cl_zgesum_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileW;

    tileA = cti_interface_get(descr[0]);
    tileW = cti_interface_get(descr[1]);

    TCORE_zgesum( clargs->storev, clargs->m, clargs->n, tileA, tileW );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgesum, cl_zgesum_cpu_func)

void INSERT_TASK_zgesum( const RUNTIME_option_t *options,
                         cham_store_t storev, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SUMS, int SUMSm, int SUMSn )
{
    struct cl_zgesum_args_s *clargs = NULL;
    struct starpu_codelet *codelet = &cl_zgesum;
    void (*callback)(void*) = options->profiling ? cl_zgesum_callback : NULL;
    int exec = 0;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(SUMS, SUMSm, SUMSn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zgesum_args_s ) );
        clargs->storev = storev;
        clargs->m      = m;
        clargs->n      = n;
    }

    rt_starpu_insert_task(
        codelet,
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zgesum_args_s),
        STARPU_R,        RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_RW,       RTBLKADDR(SUMS, ChamComplexDouble, SUMSm, SUMSn),
        STARPU_PRIORITY, options->priority,
        STARPU_CALLBACK, callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        0 );
}
