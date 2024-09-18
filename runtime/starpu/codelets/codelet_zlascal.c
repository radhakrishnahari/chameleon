/**
 *
 * @file starpu/codelet_zlascal.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlascal StarPU codelet
 *
 * @version 1.3.0
 * @author Dalal Sukkari
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Samuel Thibault
 * @date 2024-03-11
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zlascal_args_s {
    cham_uplo_t uplo;
    int m;
    int n;
    CHAMELEON_Complex64_t alpha;
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zlascal_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_zlascal_args_s *clargs = (struct cl_zlascal_args_s *)cl_arg;
    CHAM_tile_t *tileA;

    tileA = cti_interface_get(descr[0]);

    TCORE_zlascal( clargs->uplo, clargs->m, clargs->n, clargs->alpha, tileA );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zlascal, cl_zlascal_cpu_func )

void INSERT_TASK_zlascal( const RUNTIME_option_t *options,
                          cham_uplo_t uplo,
                          int m, int n, int nb,
                          CHAMELEON_Complex64_t alpha,
                          const CHAM_desc_t *A, int Am, int An )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlaset( options, uplo, m, n,
                            alpha, alpha, A, Am, An );
        return;
    }
    else if ( alpha == 1. ) {
        return;
    }

    struct cl_zlascal_args_s *clargs = NULL;
    void (*callback)(void*);
    int                      exec = 0;
    const char              *cl_name = "zlascal";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zlascal_args_s ) );
        clargs->uplo  = uplo;
        clargs->m     = m;
        clargs->n     = n;
        clargs->alpha = alpha;
    }

    /* Callback fro profiling information */
    callback = options->profiling ? cl_zlascal_callback : NULL;

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zlascal,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zlascal_args_s),
        STARPU_RW,     RTBLKADDR(A, ChamComplexDouble, Am, An),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              cl_name,
#endif

        0 );

    (void)nb;
}
