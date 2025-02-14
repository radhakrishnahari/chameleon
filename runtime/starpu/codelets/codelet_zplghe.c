/**
 *
 * @file starpu/codelet_zplghe.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplghe StarPU codelet
 *
 * @version 1.3.0
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2024-10-18
 * @precisions normal z -> c
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zplghe_args_s {
    double bump;
    int m;
    int n;
    int bigM;
    int m0;
    int n0;
    unsigned long long int seed;
};

#if !defined(CHAMELEON_SIMULATION)
static void cl_zplghe_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zplghe_args_s *clargs = (struct cl_zplghe_args_s *)cl_arg;
    CHAM_tile_t *tileA;

    tileA = cti_interface_get(descr[0]);

    TCORE_zplghe( clargs->bump, clargs->m, clargs->n, tileA,
                  clargs->bigM, clargs->m0, clargs->n0, clargs->seed );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zplghe, cl_zplghe_cpu_func )

void INSERT_TASK_zplghe( const RUNTIME_option_t *options,
                         double bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed )
{
    struct cl_zplghe_args_s *clargs = NULL;
    void (*callback)(void*);
    int                      exec = 0;
    const char              *cl_name = "zplghe";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_W(A, Am, An);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zplghe_args_s ) );
        clargs->bump  = bump;
        clargs->m     = m;
        clargs->n     = n;
        clargs->bigM  = bigM;
        clargs->m0    = m0;
        clargs->n0    = n0;
        clargs->seed  = seed;
    }

    /* Callback fro profiling information */
    callback = options->profiling ? cl_zplghe_callback : NULL;

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zplghe,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zplghe_args_s),
        STARPU_W,      RTBLKADDR(A, ChamComplexDouble, Am, An),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_NAME,              cl_name,
        0 );

}
