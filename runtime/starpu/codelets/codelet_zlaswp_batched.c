/**
 *
 * @file starpu/codelet_zlaswp_batched.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets to apply zlaswp on a panel
 *
 * @version 1.3.0
 * @author Alycia Lisito
 * @date 2024-11-12
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zlaswp_batched_args_s {
    int                      tasks_nbr;
    int                      minmn;
    int                      m0[CHAMELEON_BATCH_SIZE];
    struct starpu_data_descr handle_mode[CHAMELEON_BATCH_SIZE];
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zlaswp_batched_cpu_func( void *descr[],
                            void *cl_arg )
{
    int          i, m0, minmn, *perm, *invp;
    CHAM_tile_t *A, *U, *B;
    struct cl_zlaswp_batched_args_s *clargs = ( struct cl_zlaswp_batched_args_s * ) cl_arg;

    minmn = clargs->minmn;
    perm = (int *)STARPU_VECTOR_GET_PTR( descr[0] );
    invp = (int *)STARPU_VECTOR_GET_PTR( descr[1] );
    U    = (CHAM_tile_t *) cti_interface_get( descr[2] );
    B    = (CHAM_tile_t *) cti_interface_get( descr[3] );

    for ( i = 0; i < clargs->tasks_nbr; i++ ) {
        A  = (CHAM_tile_t *) cti_interface_get( descr[ i + 4 ] );
        m0 = clargs->m0[ i ];
        TCORE_zlaswp_get( m0, A->m, A->n, minmn, A, U, perm );
        TCORE_zlaswp_set( m0, A->m, A->n, minmn, B, A, invp );
    }
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU( zlaswp_batched, cl_zlaswp_batched_cpu_func )

void INSERT_TASK_zlaswp_batched( const RUNTIME_option_t *options,
                                 int                     m0,
                                 int                     minmn,
                                 void                   *ws,
                                 const CHAM_ipiv_t      *ipiv,
                                 int                     ipivk,
                                 const CHAM_desc_t      *Am,
                                 int                     Amm,
                                 int                     Amn,
                                 const CHAM_desc_t      *Ak,
                                 int                     Akm,
                                 int                     Akn,
                                 const CHAM_desc_t      *U,
                                 int                     Um,
                                 int                     Un,
                                 void                  **clargs_ptr )
{
    int task_num   = 0;
    int batch_size = ((struct chameleon_pzgetrf_s *)ws)->batch_size_swap;
    struct cl_zlaswp_batched_args_s *clargs = *clargs_ptr;
    if ( Am->get_rankof( Am, Amm, Amn) != Am->myrank ) {
        return;
    }

    if( clargs == NULL ) {
        clargs = malloc( sizeof( struct cl_zlaswp_batched_args_s ) ) ;
        clargs->tasks_nbr = 0;
        clargs->minmn     = minmn;
        *clargs_ptr       = clargs;
    }

    task_num               = clargs->tasks_nbr;
    clargs->m0[ task_num ] = m0;
    clargs->handle_mode[ task_num ].handle = RTBLKADDR(Am, CHAMELEON_Complex64_t, Amm, Amn);
    clargs->handle_mode[ task_num ].mode   = STARPU_RW;
    clargs->tasks_nbr ++;

    if ( clargs->tasks_nbr == batch_size ) {
        INSERT_TASK_zlaswp_batched_flush( options, ipiv, ipivk, Ak, Akm, Akn, U, Um, Un, clargs_ptr );
    }
}

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zlaswp_batched_flush( const RUNTIME_option_t *options,
                                       const CHAM_ipiv_t      *ipiv,
                                       int                     ipivk,
                                       const CHAM_desc_t      *Ak,
                                       int                     Akm,
                                       int                     Akn,
                                       const CHAM_desc_t      *U,
                                       int                     Um,
                                       int                     Un,
                                       void                  **clargs_ptr )
{
    struct cl_zlaswp_batched_args_s *clargs   = *clargs_ptr;
    int                             nhandles;

    if( clargs == NULL ) {
        return;
    }

    nhandles = clargs->tasks_nbr;
    rt_starpu_insert_task(
        &cl_zlaswp_batched,
        STARPU_CL_ARGS,             clargs, sizeof(struct cl_zlaswp_batched_args_s),
        STARPU_R,                   RUNTIME_perm_getaddr( ipiv, ipivk ),
        STARPU_R,                   RUNTIME_invp_getaddr( ipiv, ipivk ),
        STARPU_RW | STARPU_COMMUTE, RTBLKADDR(U, ChamComplexDouble, Um, Un),
        STARPU_R,                   RTBLKADDR(Ak, ChamComplexDouble, Akm, Akn),
        STARPU_DATA_MODE_ARRAY,     clargs->handle_mode, nhandles,
        STARPU_PRIORITY,            options->priority,
        STARPU_EXECUTE_ON_WORKER,   options->workerid,
        0 );

    /* clargs is freed by starpu. */
    *clargs_ptr = NULL;
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zlaswp_batched_flush( const RUNTIME_option_t *options,
                                       const CHAM_ipiv_t      *ipiv,
                                       int                     ipivk,
                                       const CHAM_desc_t      *Ak,
                                       int                     Akm,
                                       int                     Akn,
                                       const CHAM_desc_t      *U,
                                       int                     Um,
                                       int                     Un,
                                       void                  **clargs_ptr )
{
    int ret, k;
    struct starpu_task *task;
    struct cl_zlaswp_batched_args_s *myclargs = *clargs_ptr;

    if( myclargs == NULL ) {
        return;
    }

    INSERT_TASK_COMMON_PARAMETERS( zlaswp_batched, myclargs->tasks_nbr + 4 );

    /*
     * Register the data handles, might need to receive perm and invp
     */
    starpu_cham_exchange_init_params( options, &params, Ak->myrank );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RUNTIME_perm_getaddr( ipiv, ipivk ),
                                                  STARPU_R );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RUNTIME_invp_getaddr( ipiv, ipivk ),
                                                  STARPU_R );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( U, ChamComplexDouble, Um, Un ),
                                STARPU_RW | STARPU_COMMUTE );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( Ak, ChamComplexDouble, Akm, Akn ), STARPU_R );
    for ( k = 0; k < myclargs->tasks_nbr; k++ ) {
        starpu_cham_register_descr( &nbdata, descrs, myclargs->handle_mode[ k ].handle, STARPU_RW );
    }

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    task->cl_arg      = myclargs;
    task->cl_arg_size = sizeof( struct cl_zlaswp_batched_args_s );
    task->cl_arg_free = 1;

    /* Set common parameters */
    starpu_cham_task_set_options( options, task, nbdata, descrs, NULL );

    /* Flops */
    task->flops = 0.;

    ret = starpu_task_submit( task );
    if ( ret == -ENODEV ) {
        task->destroy = 0;
        starpu_task_destroy( task );
        chameleon_error( "INSERT_TASK_zlaswp_batched", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );

    /* clargs is freed by starpu. */
    *clargs_ptr = NULL;
    (void)clargs;
    (void)cl_name;
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */
