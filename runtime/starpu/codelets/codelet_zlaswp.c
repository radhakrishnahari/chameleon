/**
 *
 * @file starpu/codelet_zlaswp.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets to apply zlaswp on a panel
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-10-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zlaswp_args_s {
    cham_side_t side;
    int         m0;
    int         m;
    int         n;
    int         k;
};

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlaswp_get_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_zlaswp_args_s *clargs = (struct cl_zlaswp_args_s *)cl_arg;
    int         *perm;
    CHAM_tile_t *A, *B;

    perm = (int *)STARPU_VECTOR_GET_PTR( descr[0] );
    A    = (CHAM_tile_t *) cti_interface_get( descr[1] );
    B    = (CHAM_tile_t *) cti_interface_get( descr[2] );

    TCORE_zlaswp_get( clargs->side, clargs->m0, clargs->m, clargs->n, clargs->k, A, B, perm );
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU( zlaswp_get, cl_zlaswp_get_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zlaswp_get( const RUNTIME_option_t *options,
                             cham_side_t side, cham_dir_t dir,
                             int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *U, int Um, int Un )
{
    void                  *ipiv_handle;
    struct starpu_codelet *codelet = &cl_zlaswp_get;
    if ( A->get_rankof( A, Am, An ) != A->myrank ) {
        return;
    }

    struct cl_zlaswp_args_s *clargs;
    clargs = malloc( sizeof( struct cl_zlaswp_args_s ) );
    clargs->side = side;
    clargs->m0   = m0;
    clargs->m    = m;
    clargs->n    = n;
    clargs->k    = k;

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_perm_getaddr( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_invp_getaddr( ipiv, ipivk );
    }
    //void (*callback)(void*) = options->profiling ? cl_zlaswp_get_callback : NULL;

    rt_starpu_insert_task(
        codelet,
        STARPU_CL_ARGS,             clargs, sizeof(struct cl_zlaswp_args_s),
        STARPU_R,                   ipiv_handle,
        STARPU_R,                   RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_RW | STARPU_COMMUTE, RTBLKADDR(U, ChamComplexDouble, Um, Un),
        STARPU_PRIORITY,            options->priority,
        //STARPU_CALLBACK,            callback,
        STARPU_EXECUTE_ON_WORKER,   options->workerid,
        0 );
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zlaswp_get( const RUNTIME_option_t *options,
                             cham_side_t side, cham_dir_t dir,
                             int m0, int m, int n, int k,
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

    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zlaswp_get, zlaswp_get, zlaswp, 3);

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

    clargs = malloc( sizeof( struct cl_zlaswp_args_s ) );
    clargs->side = side;
    clargs->m0   = m0;
    clargs->m    = m;
    clargs->n    = n;
    clargs->k    = k;

    task->cl_arg      = clargs;
    task->cl_arg_size = sizeof( struct cl_zlaswp_args_s );
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
static void cl_zlaswp_set_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_zlaswp_args_s *clargs = (struct cl_zlaswp_args_s *)cl_arg;
    int         *invp;
    CHAM_tile_t *A, *B;

    invp = (int *)STARPU_VECTOR_GET_PTR( descr[0] );
    A    = (CHAM_tile_t *) cti_interface_get( descr[1] );
    B    = (CHAM_tile_t *) cti_interface_get( descr[2] );

    TCORE_zlaswp_set( clargs->side, clargs->m0, clargs->m, clargs->n, clargs->k, A, B, invp );
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU( zlaswp_set, cl_zlaswp_set_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zlaswp_set( const RUNTIME_option_t *options,
                             cham_side_t             side,
                             cham_dir_t dir, int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn )
{
    void                  *ipiv_handle;
    struct starpu_codelet *codelet = &cl_zlaswp_set;
    if ( B->get_rankof( B, Bm, Bn) != A->myrank ) {
        return;
    }

    struct cl_zlaswp_args_s *clargs;
    clargs = malloc( sizeof( struct cl_zlaswp_args_s ) );
    clargs->side = side;
    clargs->m0   = m0;
    clargs->m    = m;
    clargs->n    = n;
    clargs->k    = k;

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_invp_getaddr( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_perm_getaddr( ipiv, ipivk );
    }

    //void (*callback)(void*) = options->profiling ? cl_zlaswp_set_callback : NULL;

    rt_starpu_insert_task(
        codelet,
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_zlaswp_args_s),
        STARPU_R,                 ipiv_handle,
        STARPU_R,                 RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_RW,                RTBLKADDR(B, ChamComplexDouble, Bm, Bn),
        STARPU_PRIORITY,          options->priority,
        //STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        0 );
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zlaswp_set( const RUNTIME_option_t *options,
                             cham_side_t             side,
                             cham_dir_t dir, int m0, int m, int n, int k,
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

    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zlaswp_set, zlaswp_set, zlaswp, 3);

    /*
     * Register the data handles, might need to receive perm and invp
     */
    starpu_cham_exchange_init_params( options, &params, B->get_rankof( B, Bm, Bn ) );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  ipiv_handle, STARPU_R );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( A, ChamComplexDouble, Am, An), STARPU_R );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( B, ChamComplexDouble, Bm, Bn ), STARPU_RW );

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    clargs = malloc( sizeof( struct cl_zlaswp_args_s ) );
    clargs->side = side;
    clargs->m0   = m0;
    clargs->m    = m;
    clargs->n    = n;
    clargs->k    = k;

    task->cl_arg      = clargs;
    task->cl_arg_size = sizeof( struct cl_zlaswp_args_s );
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

#if defined(CHAMELEON_USE_MPI)

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlaswp_ret_cpu_func( void *descr[], void *cl_arg )
{
    CHAM_tile_t           *A_tile;
    cpui_interface_t      *ws     = (cpui_interface_t*) descr[1];
    cham_side_t            side   = ws->side;
    int                    ldb    = ws->n;
    int                   *index  = ws->ws.index;
    CHAMELEON_Complex64_t *rows   = ws->ws.rows;
    CHAMELEON_Complex64_t *A;
    int                    i, lda, A_inc;

    A_tile = (CHAM_tile_t *) cti_interface_get( descr[0] );
    A      = CHAM_tile_get_ptr( A_tile );

    lda   = ( side == ChamLeft ) ? A_tile->ld : 1;
    A_inc = ( side == ChamLeft ) ? 1          : A_tile->ld;

    for ( i = 0; i < ws->m; i++ ) {
        if ( index[i] != -1 ) {
            cblas_zcopy( ws->n, rows + index[i] * ldb,   1,
                                A    + i        * A_inc, lda );
        }
    }

    (void)cl_arg;
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU( zlaswp_ret, cl_zlaswp_ret_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zlaswp_ret( const RUNTIME_option_t *options,
                             CHAM_perm_t       *ws, int Wm, int Wn,
                             const CHAM_desc_t *A,  int Am, int An)
{
    struct starpu_codelet *codelet = &cl_zlaswp_ret;
    if ( A->get_rankof( A, Am, An) != A->myrank ) {
        return;
    }
    //void (*callback)(void*) = options->profiling ? cl_zlaswp_get_callback : NULL;

    rt_starpu_insert_task(
        codelet,
        STARPU_W,                   RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_R,                   RUNTIME_cpui_getaddr( ws, Wm, Wn ),
        STARPU_PRIORITY,            options->priority,
        //STARPU_CALLBACK,            callback,
        STARPU_EXECUTE_ON_WORKER,   options->workerid,
        0 );
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zlaswp_ret( const RUNTIME_option_t *options,
                             CHAM_perm_t       *ws, int Wm, int Wn,
                             const CHAM_desc_t *A,  int Am, int An)
{
    int                 rank = A->get_rankof( A, Am, An );
    int                 ret;
    struct starpu_task *task;

    if ( A->get_rankof( A, Am, An) != A->myrank ) {
        return;
    }

    INSERT_TASK_COMMON_PARAMETERS( zlaswp_ret, 2);

    //void (*callback)(void*) = options->profiling ? cl_zlaswp_get_callback : NULL;
    starpu_cham_exchange_init_params( options, &params, rank );

    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( A, ChamComplexDouble, Am, An ), STARPU_W );
    starpu_cham_register_descr( &nbdata, descrs, RUNTIME_cpui_getaddr( ws, Wm, Wn ),
                                STARPU_R );

    task = starpu_task_create();
    task->cl = cl;

    starpu_cham_task_set_options( options, task, nbdata, descrs, NULL );

    /* Flops */
    task->flops = 0.;

    /* Refine name */
    task->name = cl_name;

    ret = starpu_task_submit( task );
    if ( ret == -ENODEV ) {
        task->destroy = 0;
        starpu_task_destroy( task );
        chameleon_error( "INSERT_TASK_zlaswp_ret", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
    (void)clargs;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

#endif /* defined(CHAMELEON_USE_MPI) */
