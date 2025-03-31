/**
 *
 * @file starpu/codelet_zperm_allreduce.c
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
 * @author Pierre Esterie
 * @date 2025-03-24
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
    int p_first;
    int me;
    int shift;
    int np_inv;
};

static void
cl_zperm_allreduce_cpu_func( void *descr[], void *cl_arg )
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
    int p_first = clargs->p_first / q;
    int shift   = clargs->shift;
    int np      = clargs->np_inv;
    int me      = ( p <= np ) ? clargs->me / q : ( ( clargs->me / q ) - p_first + p ) % p;
    int nb      = tileUinout->n;
    int first   = me - 2 * shift + 1;
    int last    = me -     shift;
    int i, m, ownerp;

    for ( i = 0; i < tempmm; i++ ) {
        m      = perm[ i ] / mb;
        ownerp = ( p <= np ) ? ( (m % p) * q + (n % q) ) / q : ( ( (m % p) * q + (n % q) ) / q - p_first + p ) % p;

        if ( ( (first    <= ownerp) && (ownerp <= last   ) ) ||
             ( (first+np <= ownerp) && (ownerp <= last+np) ) )
        {
            cblas_zcopy( nb, Uin    + i, tileUin->ld,
                             Uinout + i, tileUinout->ld );
        }
    }
}

CODELETS_CPU( zperm_allreduce, cl_zperm_allreduce_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)

static void
INSERT_TASK_zperm_allreduce_send( const RUNTIME_option_t *options,
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
INSERT_TASK_zperm_allreduce_recv( const RUNTIME_option_t *options,
                                  cham_dir_t              dir,
                                  CHAM_desc_t            *U,
                                  CHAM_ipiv_t            *ipiv,
                                  int                     ipivk,
                                  int                     me,
                                  int                     src,
                                  int                     n,
                                  int                     tempmm,
                                  int                     p,
                                  int                     q,
                                  int                     shift,
                                  int                     np,
                                  int                     p_first )
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
    clargs->p_first = p_first;
    clargs->me      = me;
    clargs->shift   = shift;
    clargs->np_inv  = np;

    rt_starpu_insert_task(
        &cl_zperm_allreduce,
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
INSERT_TASK_zperm_allreduce_send( const RUNTIME_option_t *options,
                                  CHAM_desc_t            *U,
                                  int                     me,
                                  int                     dst,
                                  int                     n )
{
    INSERT_TASK_COMMON_PARAMETERS_CLNULL( zperm_allreduce_send, 1 );

    starpu_cham_exchange_init_params( options, &params, dst );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( U, ChamComplexDouble, me, n ),
                                                  STARPU_R );
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
    (void)cl;
    (void)cl_name;
}

static void
INSERT_TASK_zperm_allreduce_recv( const RUNTIME_option_t *options,
                                  cham_dir_t              dir,
                                  CHAM_desc_t            *U,
                                  CHAM_ipiv_t            *ipiv,
                                  int                     ipivk,
                                  int                     me,
                                  int                     src,
                                  int                     n,
                                  int                     tempmm,
                                  int                     p,
                                  int                     q,
                                  int                     shift,
                                  int                     np,
                                  int                     p_first )
{
    int                 ret;
    struct starpu_task *task;
    void               *ipiv_handle

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_perm_getaddr( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_invp_getaddr( ipiv, ipivk );
    }

    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zperm_allreduce_send, zperm_allreduce, redux, 3 );

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
    clargs->p_first = p_first;
    clargs->me      = me;
    clargs->shift   = shift;
    clargs->np_inv  = np;

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
        chameleon_error( "INSERT_TASK_zperm_allreduce", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
    starpu_mpi_cache_flush( options->sequence->comm, RTBLKADDR(U, CHAMELEON_Complex64_t, src, n) );
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

static void
zperm_allreduce_chameleon_starpu_task( const RUNTIME_option_t     *options,
                                       cham_dir_t                  dir,
                                       const CHAM_desc_t          *A,
                                       CHAM_desc_t                *U,
                                       int                         Um,
                                       int                         Un,
                                       CHAM_ipiv_t                *ipiv,
                                       int                         ipivk,
                                       int                         k,
                                       int                         n,
                                       struct chameleon_pzgetrf_s *ws )
{
    int *proc_involved = ws->proc_involved;
    int  np_involved   = chameleon_min( chameleon_desc_datadist_get_iparam(A, 0), A->mt - k );
    int  np_iter       = np_involved;
    int  p_recv, p_send, me, p_first;
    int  shift = 1;

    if ( np_involved == 1 ) {
        assert( proc_involved[0] == A->myrank );
    }
    else {
        p_first = proc_involved[0];
        for( me = 0; me < np_involved; me++ ) {
            if ( proc_involved[me] == A->myrank ) {
                break;
            }
        }
        assert( me < np_involved );
        while ( np_iter > 1 ) {
            p_send = proc_involved[ ( me + shift               ) % np_involved ];
            p_recv = proc_involved[ ( me - shift + np_involved ) % np_involved ];

            INSERT_TASK_zperm_allreduce_send( options, U, A->myrank, p_send, n );
            INSERT_TASK_zperm_allreduce_recv( options, dir, U, ipiv, ipivk, A->myrank, p_recv,
                                              n, k == (A->mt-1) ? A->m - k * A->mb : A->mb,
                                              chameleon_desc_datadist_get_iparam(A, 0),
                                              chameleon_desc_datadist_get_iparam(A, 1),
                                              shift, np_involved, p_first );

            shift   = shift << 1;
            np_iter = chameleon_ceil( np_iter, 2 );
        }
    }
}

void
INSERT_TASK_zperm_allreduce_row( const RUNTIME_option_t *options,
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
    struct chameleon_pzgetrf_s *tmp = (struct chameleon_pzgetrf_s *)ws;
    cham_getrf_allreduce_t alg = tmp->alg_allreduce;
    switch( alg ) {
    case ChamStarPUTasks:
    default:
        zperm_allreduce_chameleon_starpu_task( options, dir, A, U, Um, Un, ipiv, ipivk, k, n, tmp );
    }
}

void
INSERT_TASK_zperm_allreduce_send_A( const RUNTIME_option_t *options,
                                    CHAM_desc_t            *A,
                                    int                     Am,
                                    int                     An,
                                    int                     myrank,
                                    int                     np,
                                    int                    *proc_involved )
{
    int p;

    for ( p = 0; p < np; p ++ ) {
        if ( proc_involved[ p ] == myrank ) {
            continue;
        }
        starpu_mpi_get_data_on_node_detached( options->sequence->comm,
                                              RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
                                              proc_involved[ p ], NULL, NULL );
    }
}

void
INSERT_TASK_zperm_allreduce_send_perm( const RUNTIME_option_t *options,
                                       cham_dir_t              dir,
                                       CHAM_ipiv_t            *ipiv,
                                       int                     ipivk,
                                       int                     myrank,
                                       int                     np,
                                       int                    *proc_involved )
{
    int   p;
    void *ipiv_handle;

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_perm_getaddr( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_invp_getaddr( ipiv, ipivk );
    }

    for ( p = 0; p < np; p++ ) {
        if ( proc_involved[ p ] == myrank ) {
            continue;
        }
        starpu_mpi_get_data_on_node_detached( options->sequence->comm,
                                              ipiv_handle,
                                              proc_involved[ p ], NULL, NULL );
    }
}

void
INSERT_TASK_zperm_allreduce_send_invp_row( const RUNTIME_option_t *options,
                                           cham_dir_t              dir,
                                           CHAM_ipiv_t            *ipiv,
                                           int                     ipivk,
                                           const CHAM_desc_t      *A,
                                           int                     k,
                                           int                     n )
{
    int   b, rank;
    void *ipiv_handle;

    if ( dir == ChamDirForward ) {
        ipiv_handle = RUNTIME_invp_getaddr( ipiv, ipivk );
    }
    else {
        ipiv_handle = RUNTIME_perm_getaddr( ipiv, ipivk );
    }

    for ( b = k+1; (b < A->mt) && ((b-(k+1)) < chameleon_desc_datadist_get_iparam(A, 0)); b ++ ) {
        rank = A->get_rankof( A, b, n );
        if ( rank == A->myrank ) {
            continue;
        }
        starpu_mpi_get_data_on_node_detached( options->sequence->comm,
                                              ipiv_handle,
                                              rank, NULL, NULL );
    }
}

#else
void
INSERT_TASK_zperm_allreduce_send_A( const RUNTIME_option_t *options,
                                    CHAM_desc_t            *A,
                                    int                     Am,
                                    int                     An,
                                    int                     myrank,
                                    int                     np,
                                    int                    *proc_involved )
{
    (void)options;
    (void)A;
    (void)Am;
    (void)An;
    (void)myrank;
    (void)np;
    (void)proc_involved;
}

void
INSERT_TASK_zperm_allreduce_send_perm( const RUNTIME_option_t *options,
                                       cham_dir_t              dir,
                                       CHAM_ipiv_t            *ipiv,
                                       int                     ipivk,
                                       int                     myrank,
                                       int                     np,
                                       int                    *proc_involved )
{
    (void)options;
    (void)ipiv;
    (void)ipivk;
    (void)myrank;
    (void)np;
    (void)proc_involved;
}

void
INSERT_TASK_zperm_allreduce_send_invp_row( const RUNTIME_option_t *options,
                                           cham_dir_t              dir,
                                           CHAM_ipiv_t            *ipiv,
                                           int                     ipivk,
                                           const CHAM_desc_t      *A,
                                           int                     k,
                                           int                     n )
{
    (void)options;
    (void)ipiv;
    (void)ipivk;
    (void)A;
    (void)k;
    (void)n;
}

void
INSERT_TASK_zperm_allreduce_row( const RUNTIME_option_t *options,
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
