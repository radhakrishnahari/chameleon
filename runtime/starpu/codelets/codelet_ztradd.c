/**
 *
 * @file starpu/codelet_ztradd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztradd StarPU codelet
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_ztradd_args_s {
    cham_uplo_t           uplo;
    cham_trans_t          trans;
    int                   m;
    int                   n;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t beta;
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_ztradd_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_ztradd_args_s *clargs = (struct cl_ztradd_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    TCORE_ztradd( clargs->uplo, clargs->trans, clargs->m, clargs->n,
                  clargs->alpha, tileA, clargs->beta, tileB );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( ztradd, cl_ztradd_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)
void INSERT_TASK_ztradd( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, uplo, m, n, nb,
                             beta, B, Bm, Bn );
        return;
    }

    void (*callback)(void*);
    struct cl_ztradd_args_s *clargs  = NULL;
    int                      exec    = 0;
    const char              *cl_name = "ztradd";
    int                      accessB;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(B, Bm, Bn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_ztradd_args_s ) );
        clargs->uplo  = uplo;
        clargs->trans = trans;
        clargs->m     = m;
        clargs->n     = n;
        clargs->alpha = alpha;
        clargs->beta  = beta;
    }

    /* Callback fro profiling information */
    callback = options->profiling ? cl_ztradd_callback : NULL;

    /* Reduce the B access if needed */
    accessB = ( beta == 0. ) ? STARPU_W : STARPU_RW;

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_ztradd,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_ztradd_args_s),
        STARPU_R,      RTBLKADDR(A, ChamComplexDouble, Am, An),
        accessB,       RTBLKADDR(B, ChamComplexDouble, Bm, Bn),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_NAME,              cl_name,
        0 );

    (void)nb;
}

#else

void INSERT_TASK_ztradd( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, uplo, m, n, nb,
                             beta, B, Bm, Bn );
        return;
    }

    INSERT_TASK_COMMON_PARAMETERS( ztradd, 2 );
    int accessB;


    /* Reduce the B access if needed */
    accessB = ( beta == 0. ) ? STARPU_W : STARPU_RW;

    /*
     * Set the data handles and initialize exchanges if needed
     */
    starpu_cham_exchange_init_params( options, &params, B->get_rankof( B, Bm, Bn ) );
    starpu_cham_exchange_data_before_execution( options, params, &nbdata, descrs, A, Am, An, STARPU_R );
    starpu_cham_exchange_data_before_execution( options, params, &nbdata, descrs, B, Bm, Bn, accessB  );

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
        clargs = malloc( sizeof( struct cl_ztradd_args_s ) );
        clargs->uplo  = uplo;
        clargs->trans = trans;
        clargs->m     = m;
        clargs->n     = n;
        clargs->alpha = alpha;
        clargs->beta  = beta;

        task->cl_arg      = clargs;
        task->cl_arg_size = sizeof( struct cl_ztradd_args_s );
        task->cl_arg_free = 1;

        /* Set common parameters */
        starpu_cham_task_set_options( options, task, nbdata, descrs, cl_ztradd_callback );

        /* Flops */
        //task->flops = flops_ztradd( m, n );

        /* Refine name */
        task->name = chameleon_codelet_name( cl_name, 2,
                                             A->get_blktile( A, Am, An ),
                                             B->get_blktile( B, Bm, Bn ) );

        ret = starpu_task_submit( task );
        if ( ret == -ENODEV ) {
            task->destroy = 0;
            starpu_task_destroy( task );
            chameleon_error( "INSERT_TASK_ztradd", "Failed to submit the task to StarPU" );
            return;
        }
    }

    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    (void)nb;
}

#endif
