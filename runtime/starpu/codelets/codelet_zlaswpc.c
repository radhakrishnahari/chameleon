/**
 *
 * @file starpu/codelet_zlaswpc.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets to apply zlaswpc on a panel
 *
 * @version 1.3.0
 * @author Matteo Marcos
 * @date 2025-06-16
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zlaswpc_args_s {
    int n0;
    int m;
    int n;
    int k;
};

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlaswpc_get_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_zlaswpc_args_s *clargs = (struct cl_zlaswpc_args_s *)cl_arg;
    int         *perm;
    CHAM_tile_t *A, *B;

    perm = (int *)STARPU_VECTOR_GET_PTR( descr[0] );
    A    = (CHAM_tile_t *) cti_interface_get( descr[1] );
    B    = (CHAM_tile_t *) cti_interface_get( descr[2] );

    TCORE_zlaswpc_get( clargs->n0, clargs->m, clargs->n, clargs->k, A, B, perm );
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU( zlaswpc_get, cl_zlaswpc_get_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zlaswpc_get( const RUNTIME_option_t *options,
                              cham_dir_t dir, int n0, int m, int n, int k,
                              const CHAM_ipiv_t *ipiv, int ipivk,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *U, int Um, int Un )
{
    void                  *ipiv_handle;
    struct starpu_codelet *codelet = &cl_zlaswpc_get;
    if ( A->get_rankof( A, Am, An) != A->myrank ) {
        return;
    }

    struct cl_zlaswpc_args_s *clargs;
    clargs = malloc( sizeof( struct cl_zlaswpc_args_s ) );
    clargs->n0 = n0;
    clargs->m  = m;
    clargs->n  = n;
    clargs->k  = k;

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_perm_getaddr( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_invp_getaddr( ipiv, ipivk );
    }
    //void (*callback)(void*) = options->profiling ? cl_zlaswp_get_callback : NULL;

    rt_starpu_insert_task(
        codelet,
        STARPU_CL_ARGS,             clargs, sizeof(struct cl_zlaswpc_args_s),
        STARPU_R,                   ipiv_handle,
        STARPU_R,                   RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_RW | STARPU_COMMUTE, RTBLKADDR(U, ChamComplexDouble, Um, Un),
        STARPU_PRIORITY,            options->priority,
        //STARPU_CALLBACK,            callback,
        STARPU_EXECUTE_ON_WORKER,   options->workerid,
        0 );
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zlaswpc_get( const RUNTIME_option_t *options,
                              cham_dir_t dir, int n0, int m, int n, int k,
                              const CHAM_ipiv_t *ipiv, int ipivk,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *U, int Um, int Un )
{
    int                 ret;
    struct starpu_task *task;
    void               *ipiv_handle;

    if ( A->get_rankof( A, Am, An) != A->myrank ) {
        return;
    }

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_perm_getaddr( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_invp_getaddr( ipiv, ipivk );
    }

    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zlaswpc_get, zlaswpc_get, zlaswpc, 3);

    /*
     * Register the data handles, might need to receive perm and invp
     */
    starpu_cham_exchange_init_params( options, &params, U->get_rankof( U, Um, Un ) );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  ipiv_handle, STARPU_R );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( A, ChamComplexDouble, Am, An ), STARPU_R );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( U, ChamComplexDouble, Um, Un ),
                                STARPU_RW | STARPU_COMMUTE );

    task = starpu_task_create();
    task->cl = cl;

    clargs = malloc( sizeof( struct cl_zlaswpc_args_s ) );
    clargs->n0 = n0;
    clargs->m  = m;
    clargs->n  = n;
    clargs->k  = k;

    task->cl_arg      = clargs;
    task->cl_arg_size = sizeof( struct cl_zlaswpc_args_s );
    task->cl_arg_free = 1;

    starpu_cham_task_set_options( options, task, nbdata, descrs, NULL );

    /* Flops */
    task->flops = 0.;

    /* Refine name */
    task->name = cl_name;

    ret = starpu_task_submit( task );
    if ( ret == -ENODEV ) {
        task->destroy = 0;
        starpu_task_destroy( task );
        chameleon_error( "INSERT_TASK_zlaswp_get", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlaswpc_set_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_zlaswpc_args_s *clargs = (struct cl_zlaswpc_args_s *)cl_arg;
    int         *invp;
    CHAM_tile_t *A, *B;

    invp = (int *)STARPU_VECTOR_GET_PTR( descr[0] );
    A    = (CHAM_tile_t *) cti_interface_get( descr[1] );
    B    = (CHAM_tile_t *) cti_interface_get( descr[2] );

    TCORE_zlaswpc_set( clargs->n0, clargs->m, clargs->n, clargs->k, A, B, invp );
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU( zlaswpc_set, cl_zlaswpc_set_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zlaswpc_set( const RUNTIME_option_t *options,
                              cham_dir_t dir, int n0, int m, int n, int k,
                              const CHAM_ipiv_t *ipiv, int ipivk,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *B, int Bm, int Bn )
{
    void                  *ipiv_handle;
    struct starpu_codelet *codelet = &cl_zlaswpc_set;
    if ( B->get_rankof( B, Bm, Bn) != A->myrank ) {
        return;
    }

    struct cl_zlaswpc_args_s *clargs;
    clargs = malloc( sizeof( struct cl_zlaswpc_args_s ) );
    clargs->n0 = n0;
    clargs->m  = m;
    clargs->n  = n;
    clargs->k  = k;

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_invp_getaddr( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_perm_getaddr( ipiv, ipivk );
    }

    //void (*callback)(void*) = options->profiling ? cl_zlaswp_set_callback : NULL;

    rt_starpu_insert_task(
        codelet,
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_zlaswpc_args_s),
        STARPU_R,                 ipiv_handle,
        STARPU_R,                 RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_RW,                RTBLKADDR(B, ChamComplexDouble, Bm, Bn),
        STARPU_PRIORITY,          options->priority,
        //STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        0 );
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zlaswpc_set( const RUNTIME_option_t *options,
                              cham_dir_t dir, int n0, int m, int n, int k,
                              const CHAM_ipiv_t *ipiv, int ipivk,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *B, int Bm, int Bn )
{
    int                 ret;
    struct starpu_task *task;
    void               *ipiv_handle;

    if ( B->get_rankof( B, Bm, Bn) != A->myrank ) {
        return;
    }

    if( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_invp_getaddr( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_perm_getaddr( ipiv, ipivk );
    }

    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zlaswpc_set, zlaswpc_set, zlaswpc, 3);

    /*
     * Register the data handles, might need to receive perm and invp
     */
    starpu_cham_exchange_init_params( options, &params, B->get_rankof( B, Bm, Bn ) );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  ipiv_handle, STARPU_R );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( A, ChamComplexDouble, Am, An ), STARPU_R );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( B, ChamComplexDouble, Bm, Bn ), STARPU_RW );

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    clargs = malloc( sizeof( struct cl_zlaswpc_args_s ) );
    clargs->n0 = n0;
    clargs->m  = m;
    clargs->n  = n;
    clargs->k  = k;

    task->cl_arg      = clargs;
    task->cl_arg_size = sizeof( struct cl_zlaswpc_args_s );
    task->cl_arg_free = 1;

    /* Set common parameters */
    starpu_cham_task_set_options( options, task, nbdata, descrs, NULL );

    /* Flops */
    task->flops = 0.;

    /* Refine name */
    task->name = cl_name;

    ret = starpu_task_submit( task );
    if ( ret == -ENODEV ) {
        task->destroy = 0;
        starpu_task_destroy( task );
        chameleon_error( "INSERT_TASK_zlaswp_set", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
}
#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

