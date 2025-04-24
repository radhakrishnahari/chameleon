/**
 *
 * @file starpu/codelet_zperm_reduce.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets to do the reduction
 *
 * @version 1.3.0
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-04-11
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"
#include <coreblas/cblas_wrapper.h>

#if defined(CHAMELEON_USE_MPI)

struct cl_redux_args_s {
    int tempmm;
    int mb;
    int n;
    int p;
    int q;
    int me;
};

static void
cl_zperm_reduce_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_redux_args_s      *clargs     = (struct cl_redux_args_s *) cl_arg;
    const CHAM_tile_t           *tileUinout = cti_interface_get( descr[0] );
    const CHAM_tile_t           *tileUin    = cti_interface_get( descr[1] );
    const int                   *perm       = (int *)STARPU_VECTOR_GET_PTR( descr[2] );
    CHAMELEON_Complex64_t       *Uinout     = CHAM_tile_get_ptr( tileUinout );
    const CHAMELEON_Complex64_t *Uin        = CHAM_tile_get_ptr( tileUin );

    int tempmm  = clargs->tempmm;
    int mb      = clargs->mb;
    int n       = clargs->n;
    int p       = clargs->p;
    int q       = clargs->q;
    int me      = clargs->me;
    int nb      = tileUinout->n;
    int i, m, ownerp;

    for ( i = 0; i < tempmm; i++ ) {
        m      = perm[ i ] / mb;
        ownerp = ( (m % p) * q + (n % q) );

        if ( me != ownerp )
        {
            cblas_zcopy( nb, Uin    + i, tileUin->ld,
                             Uinout + i, tileUinout->ld );
        }
    }
}

CODELETS_CPU( zperm_reduce, cl_zperm_reduce_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)

static void
INSERT_TASK_zperm_reduce_send( const RUNTIME_option_t *options,
                               CHAM_desc_t            *U,
                               int                     me,
                               int                     dst,
                               int                     n )
{
    rt_starpu_insert_task(
        NULL,
        STARPU_EXECUTE_ON_NODE, dst,
        STARPU_R,               RTBLKADDR(U, CHAMELEON_Complex64_t, me, n),
        STARPU_PRIORITY,        options->priority,
        0 );
}

static void
INSERT_TASK_zperm_reduce_recv( const RUNTIME_option_t *options,
                               cham_dir_t              dir,
                               CHAM_desc_t            *U,
                               CHAM_ipiv_t            *ipiv,
                               int                     ipivk,
                               int                     me,
                               int                     src,
                               int                     n,
                               int                     tempmm,
                               int                     p,
                               int                     q )
{
    struct cl_redux_args_s *clargs;
    void                   *ipiv_handle;

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_perm_getaddr( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_invp_getaddr( ipiv, ipivk );
    }

    clargs = malloc( sizeof( struct cl_redux_args_s ) );
    clargs->tempmm  = tempmm;
    clargs->mb      = U->mb;
    clargs->n       = n;
    clargs->p       = p;
    clargs->q       = q;
    clargs->me      = me;

    rt_starpu_insert_task(
        &cl_zperm_reduce,
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_redux_args_s),
        STARPU_RW,                RTBLKADDR(U, CHAMELEON_Complex64_t, me,  n),
        STARPU_R,                 RTBLKADDR(U, CHAMELEON_Complex64_t, src, n),
        STARPU_R,                 ipiv_handle,
        STARPU_EXECUTE_ON_NODE,   me,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_PRIORITY,          options->priority,
        0 );
    starpu_mpi_cache_flush( options->sequence->comm, RTBLKADDR(U, CHAMELEON_Complex64_t, src, n) );
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

static void
INSERT_TASK_zperm_reduce_send( const RUNTIME_option_t *options,
                               CHAM_desc_t            *U,
                               int                     me,
                               int                     dst,
                               int                     n )
{
    INSERT_TASK_COMMON_PARAMETERS_CLNULL( zperm_reduce_send, 1 );

    starpu_cham_exchange_init_params( options, &params, dst );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( U, ChamComplexDouble, me, n ),
                                                  STARPU_R );
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
    (void)cl;
    (void)cl_name;
}

static void
INSERT_TASK_zperm_reduce_recv( const RUNTIME_option_t *options,
                               cham_dir_t              dir,
                               CHAM_desc_t            *U,
                               CHAM_ipiv_t            *ipiv,
                               int                     ipivk,
                               int                     me,
                               int                     src,
                               int                     n,
                               int                     tempmm,
                               int                     p,
                               int                     q )
{
    int                 ret;
    struct starpu_task *task;
    void               *ipiv_handle;

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_perm_getaddr( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_invp_getaddr( ipiv, ipivk );
    }

    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zperm_reduce_send, zperm_reduce, redux, 3 );

    starpu_cham_exchange_init_params( options, &params, me );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( U, ChamComplexDouble, me, n ),
                                                  STARPU_RW );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( U, ChamComplexDouble, src, n ),
                                                  STARPU_R );
    starpu_cham_register_descr( &nbdata, descrs, ipiv_handle, STARPU_R );

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    clargs = malloc( sizeof( struct cl_redux_args_s ) );
    clargs->tempmm  = tempmm;
    clargs->mb      = U->mb;
    clargs->n       = n;
    clargs->p       = p;
    clargs->q       = q;
    clargs->me      = me;

    task->cl_arg      = clargs;
    task->cl_arg_size = sizeof( struct cl_redux_args_s );
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
        chameleon_error( "INSERT_TASK_zperm_reduce", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
    starpu_mpi_cache_flush( options->sequence->comm, RTBLKADDR(U, CHAMELEON_Complex64_t, src, n) );
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

static void
zperm_reduce_chameleon_starpu_task( const RUNTIME_option_t     *options,
                                    cham_dir_t                  dir,
                                    const CHAM_desc_t          *A,
                                    CHAM_desc_t                *U,
                                    int                         Um,
                                    int                         Un,
                                    CHAM_ipiv_t                *ipiv,
                                    int                         ipivk,
                                    int                         k,
                                    int                         n,
                                    CHAM_reduce_t              *reduce )
{
    int *proc_involved = reduce->proc_involved;
    int  np_involved   = chameleon_min( chameleon_desc_datadist_get_iparam(A, 0), A->mt - k );
    int  p_recv, p_send, p/*, me*/;

    if( np_involved == 1 ){
        assert( proc_involved[0] == A->myrank );
        return;
    }

    if ( A->myrank == A->get_rankof( A, k, n ) ) {
        for ( p = 0; p < np_involved; p ++ ) {
            p_recv = proc_involved[p];

            if ( p_recv == A->myrank ) {
                continue;
            }

            INSERT_TASK_zperm_reduce_recv( options, dir, U, ipiv, ipivk, A->myrank, p_recv,
                                           n, k == (A->mt-1) ? A->m - k * A->mb : A->mb,
                                           chameleon_desc_datadist_get_iparam(A, 0),
                                           chameleon_desc_datadist_get_iparam(A, 1) );

        }
    }
    else {
        p_send = A->get_rankof( A, k, n );

        INSERT_TASK_zperm_reduce_send( options, U, A->myrank, p_send, n );
    }

}

void
INSERT_TASK_zperm_reduce_row( const RUNTIME_option_t *options,
                              cham_dir_t              dir,
                              const CHAM_desc_t      *A,
                              CHAM_desc_t            *U,
                              int                     Um,
                              int                     Un,
                              CHAM_ipiv_t            *ipiv,
                              int                     ipivk,
                              int                     k,
                              int                     n,
                              void                   *ws )
{
    struct chameleon_pzlaswp_s *tmp = (struct chameleon_pzlaswp_s *)ws;
    cham_getrf_allreduce_t      alg = tmp->reduce.alg_allreduce;
    switch( alg ) {
    case ChamStarPUTasks:
    default:
        zperm_reduce_chameleon_starpu_task( options, dir, A, U, Um, Un, ipiv, ipivk, k, n, &(tmp->reduce) );
    }
}

#else

void
INSERT_TASK_zperm_reduce_row( const RUNTIME_option_t *options,
                              cham_dir_t              dir,
                              const CHAM_desc_t      *A,
                              CHAM_desc_t            *U,
                              int                     Um,
                              int                     Un,
                              CHAM_ipiv_t            *ipiv,
                              int                     ipivk,
                              int                     k,
                              int                     n,
                              void                   *ws )
{
    (void)options;
    (void)A;
    (void)U;
    (void)Um;
    (void)Un;
    (void)ipiv;
    (void)ipivk;
    (void)k;
    (void)n;
    (void)ws;
}

#endif
