/**
 *
 * @file starpu/codelet_zlascal.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zlascal_args_s {
    cham_uplo_t           uplo;
    int                   m;
    int                   n;
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

#if defined(CHAMELEON_STARPU_USE_INSERT)
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

    void (*callback)(void*);
    struct cl_zlascal_args_s *clargs  = NULL;
    int                       exec    = 0;
    const char               *cl_name = "zlascal";

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
        STARPU_NAME,              cl_name,
        0 );

    (void)nb;
}

#else

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

    INSERT_TASK_COMMON_PARAMETERS( zlascal, 1 );

    /*
     * Set the data handles and initialize exchanges if needed
     */
    starpu_cham_exchange_init_params( options, &params, A->get_rankof( A, Am, An ) );
    starpu_cham_exchange_data_before_execution( options, params, &nbdata, descrs, A, Am, An, STARPU_RW );

    /*
     * Not involved, let's return
     */
    if ( nbdata == 0 ) {
        return;
    }

    if ( params.do_execute )
    {
        int ret;
        struct starpu_task *task = starpu_task_create();
        task->cl = cl;

        /* Set codelet parameters */
        clargs = malloc( sizeof( struct cl_zlascal_args_s ) );
        clargs->uplo  = uplo;
        clargs->m     = m;
        clargs->n     = n;
        clargs->alpha = alpha;

        task->cl_arg      = clargs;
        task->cl_arg_size = sizeof( struct cl_zlascal_args_s );
        task->cl_arg_free = 1;

        /* Set common parameters */
        starpu_cham_task_set_options( options, task, nbdata, descrs, cl_zlascal_callback );

        /* Flops */
        //task->flops = flops_zlascal( uplo, m, n );

        /* Refine name */
        task->name = chameleon_codelet_name( cl_name, 1,
                                             A->get_blktile( A, Am, An ) );

        ret = starpu_task_submit( task );
        if ( ret == -ENODEV ) {
            task->destroy = 0;
            starpu_task_destroy( task );
            chameleon_error( "INSERT_TASK_zlascal", "Failed to submit the task to StarPU" );
            return;
        }
    }

    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    (void)nb;
}

#endif
