/**
 *
 * @file starpu/codelet_zlauum.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlauum StarPU codelet
 *
 * @version 1.3.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zlauum_args_s {
    cham_uplo_t uplo;
    int         n;
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zlauum_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zlauum_args_s *clargs = (struct cl_zlauum_args_s *)cl_arg;
    CHAM_tile_t *tileA;

    tileA = cti_interface_get(descr[0]);

    TCORE_zlauum( clargs->uplo, clargs->n, tileA );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zlauum, cl_zlauum_cpu_func )

void INSERT_TASK_zlauum( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An )
{
    struct cl_zlauum_args_s *clargs = NULL;
    void (*callback)(void*);
    int                      exec = 0;
    const char              *cl_name = "zlauum";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zlauum_args_s ) );
        clargs->uplo  = uplo;
        clargs->n     = n;
    }

    /* Callback fro profiling information */
    callback = options->profiling ? cl_zlauum_callback : NULL;

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 1,
                                      A->get_blktile( A, Am, An ) );

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zlauum,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zlauum_args_s),
        STARPU_RW,     RTBLKADDR(A, ChamComplexDouble, Am, An),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_NAME,              cl_name,
        0 );

    (void)nb;
}
