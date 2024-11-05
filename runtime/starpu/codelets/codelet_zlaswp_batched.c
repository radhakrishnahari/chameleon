/**
 *
 * @file starpu/codelet_zlaswp_batched.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets to apply zlaswp on a panel
 *
 * @version 1.3.0
 * @author Alycia Lisito
 * @date 2024-10-21
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_laswp_batched_args_t {
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
    struct cl_laswp_batched_args_t *clargs = ( struct cl_laswp_batched_args_t * ) cl_arg;

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
                                 int                     k,
                                 int                     m,
                                 int                     n,
                                 void                   *ws,
                                 const CHAM_ipiv_t      *ipiv,
                                 int                     ipivk,
                                 const CHAM_desc_t      *A,
                                 const CHAM_desc_t      *Wu,
                                 void                  **clargs_ptr )
{
    int task_num   = 0;
    int batch_size = ((struct chameleon_pzgetrf_s *)ws)->batch_size;
    int nhandles;
    struct cl_laswp_batched_args_t *clargs = *clargs_ptr;
    if ( A->get_rankof( A, m, n) != A->myrank ) {
        return;
    }

    if( clargs == NULL ) {
        clargs = malloc( sizeof( struct cl_laswp_batched_args_t ) ) ;
        clargs->tasks_nbr = 0;
        clargs->minmn     = minmn;
        *clargs_ptr       = clargs;
    }

    task_num               = clargs->tasks_nbr;
    clargs->m0[ task_num ] = m0;
    clargs->handle_mode[ task_num ].handle = RTBLKADDR(A, CHAMELEON_Complex64_t, m, n);
    clargs->handle_mode[ task_num ].mode   = STARPU_RW;
    clargs->tasks_nbr ++;

    if ( clargs->tasks_nbr == batch_size ) {
        nhandles = clargs->tasks_nbr;
        rt_starpu_insert_task(
            &cl_zlaswp_batched,
            STARPU_CL_ARGS,             clargs, sizeof(struct cl_laswp_batched_args_t),
            STARPU_R,                   RUNTIME_perm_getaddr( ipiv, ipivk ),
            STARPU_R,                   RUNTIME_invp_getaddr( ipiv, ipivk ),
            STARPU_RW | STARPU_COMMUTE, RTBLKADDR(Wu, ChamComplexDouble, A->myrank, n),
            STARPU_R,                   RTBLKADDR(A, ChamComplexDouble, k, n),
            STARPU_DATA_MODE_ARRAY,     clargs->handle_mode, nhandles,
            STARPU_PRIORITY,            options->priority,
            STARPU_EXECUTE_ON_WORKER,   options->workerid,
            0 );

        /* clargs is freed by starpu. */
        *clargs_ptr = NULL;
    }
}

void INSERT_TASK_zlaswp_batched_flush( const RUNTIME_option_t *options,
                                       int                     k,
                                       int                     n,
                                       const CHAM_ipiv_t      *ipiv,
                                       int                     ipivk,
                                       const CHAM_desc_t      *A,
                                       const CHAM_desc_t      *U,
                                       void                  **clargs_ptr )
{
    struct cl_laswp_batched_args_t *clargs   = *clargs_ptr;
    int                             nhandles;

    if( clargs == NULL ) {
        return;
    }

    nhandles = clargs->tasks_nbr;
    rt_starpu_insert_task(
        &cl_zlaswp_batched,
        STARPU_CL_ARGS,             clargs, sizeof(struct cl_laswp_batched_args_t),
        STARPU_R,                   RUNTIME_perm_getaddr( ipiv, ipivk ),
        STARPU_R,                   RUNTIME_invp_getaddr( ipiv, ipivk ),
        STARPU_RW | STARPU_COMMUTE, RTBLKADDR(U, ChamComplexDouble, k, n),
        STARPU_R,                   RTBLKADDR(A, ChamComplexDouble, k, n),
        STARPU_DATA_MODE_ARRAY,     clargs->handle_mode, nhandles,
        STARPU_PRIORITY,            options->priority,
        STARPU_EXECUTE_ON_WORKER,   options->workerid,
        0 );

    /* clargs is freed by starpu. */
    *clargs_ptr = NULL;
}
