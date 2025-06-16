/**
 *
 * @file starpu/codelet_zgetrf_blocked.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpanel StarPU codelets
 *
 * @version 1.3.0
 * @comment Codelets to perform panel factorization with partial pivoting
 *
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-06-16
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zgetrf_blocked_args_s {
        int                 m;
        int                 n;
        int                 h;
        int                 m0;
        int                 ib;
        RUNTIME_sequence_t *sequence;
        RUNTIME_request_t  *request;
};

CHAMELEON_CL_CB( zgetrf_blocked_diag,    cti_handle_get_m(task->handles[0]), 0, 0, M )
CHAMELEON_CL_CB( zgetrf_blocked_offdiag, cti_handle_get_m(task->handles[0]), 0, 0, M )
CHAMELEON_CL_CB( zgetrf_blocked_trsm,    cti_handle_get_m(task->handles[0]), 0, 0, M )

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_blocked_diag_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zgetrf_blocked_args_s *clargs = (struct cl_zgetrf_blocked_args_s *)cl_arg;
    CHAM_tile_t           *tileA;
    CHAM_tile_t           *tileU;
    int                   *ipiv;
    cppi_interface_t      *nextpiv;
    cppi_interface_t      *prevpiv;
    CHAMELEON_Complex64_t *U   = NULL;
    int                    ldu = -1;;

    tileA   = cti_interface_get(descr[0]);
    ipiv    = (int *)STARPU_VECTOR_GET_PTR(descr[1]);
    nextpiv = (cppi_interface_t*) descr[2];
    prevpiv = (cppi_interface_t*) descr[3];
    if ( descr[4] != NULL ) {
        tileU = cti_interface_get(descr[4]);
        U     = CHAM_tile_get_ptr( tileU );
        ldu   = tileU->ld;
    }

    if ( clargs->h > 0 ) {
        cppi_display_dbg( prevpiv, stderr, "Prevpiv before call: " );
    }
    if ( clargs->h < tileA->n ) {
        cppi_display_dbg( nextpiv, stderr, "Nextpiv before call: " );
    }

    /*
     * Make sure the nextpiv interface store the right information about the
     * column and diagonal row for the reduction
     */
    nextpiv->h        = clargs->h;
    nextpiv->has_diag = 1;

    coreblas_kernel_trace( tileA );
    CORE_zgetrf_panel_diag( clargs->m, clargs->n, clargs->h, clargs->m0, clargs->ib,
                            CHAM_tile_get_ptr( tileA ), tileA->ld,
                            U, ldu,
                            ipiv, &(nextpiv->pivot), &(prevpiv->pivot) );

    if ( clargs->h > 0 ) {
        cppi_display_dbg( prevpiv, stderr, "Prevpiv after call: " );
    }
    if ( clargs->h < tileA->n ) {
        cppi_display_dbg( nextpiv, stderr, "Nextpiv after call: " );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zgetrf_blocked_diag, cl_zgetrf_blocked_diag_cpu_func )

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zgetrf_blocked_diag( const RUNTIME_option_t *options,
                                      int m, int n, int h, int m0, int ib,
                                      CHAM_desc_t *A, int Am, int An,
                                      CHAM_desc_t *U, int Um, int Un,
                                      CHAM_ipiv_t *ipiv,
                                      CHAM_desc_pivot_t *pivot )
{
#if !defined(HAVE_STARPU_NONE_NONZERO)
    /* STARPU_NONE can't be equal to 0 */
    fprintf( stderr, "INSERT_TASK_zgetrf_blocked_diag: STARPU_NONE can not be equal to 0\n" );
    assert( 0 );
#endif
    void (*callback)(void*) = options->profiling ? cl_zgetrf_blocked_diag_callback : NULL;
    const char *cl_name = "zgetrf_blocked_diag";
    int rankA           = A->get_rankof(A, Am, An);

#if defined ( CHAMELEON_USE_MPI )
    if ( ( h % ib == 0 ) && ( h > 0 ) ) {
        starpu_mpi_cache_flush( options->sequence->comm, RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un) );
    }

    if ( rankA != A->myrank ) {
        return;
    }
#endif

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW( A, Am, An );
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Set codelet parameters */
    struct cl_zgetrf_blocked_args_s *clargs;
    clargs = malloc( sizeof( struct cl_zgetrf_blocked_args_s ) );
    clargs->m        = m;
    clargs->n        = n;
    clargs->h        = h;
    clargs->m0       = m0;
    clargs->ib       = ib;
    clargs->sequence = options->sequence;
    clargs->request  = options->request;

    int access_ipiv = ( h == 0 )        ? STARPU_W    : STARPU_RW;
    int access_npiv = ( h == pivot->n ) ? STARPU_R    : STARPU_REDUX;
    int access_ppiv = ( h == 0 )        ? STARPU_NONE : STARPU_R;
    int accessU     = STARPU_RW;

    if ( h == 0 ) {
        accessU = STARPU_NONE;
        /* U can be set after ppiv because they are both none together, so it won't shift the buffers */
    }
    else if ( h%ib == 0 ) {
        accessU = STARPU_R;
    }
    else if ( ( h%ib == 1 ) || ( ib == 1 ) ) {
        accessU = STARPU_W;
    }

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 1,
                                      A->get_blktile( A, Am, An ) );

    rt_starpu_insert_task(
        &cl_zgetrf_blocked_diag,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zgetrf_blocked_args_s),

        /* Task handles */
        STARPU_RW,                RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        access_ipiv,              RUNTIME_ipiv_getaddr( ipiv, An ),
        access_npiv,              RUNTIME_pivot_getaddr( pivot, rankA, An, h ),
        access_ppiv,              RUNTIME_pivot_getaddr( pivot, rankA, An, h-1 ),
        accessU,                  RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_NAME,              cl_name,
        0 );
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zgetrf_blocked_diag( const RUNTIME_option_t *options,
                                      int m, int n, int h, int m0, int ib,
                                      CHAM_desc_t *A, int Am, int An,
                                      CHAM_desc_t *U, int Um, int Un,
                                      CHAM_ipiv_t *ipiv,
                                      CHAM_desc_pivot_t *pivot )
{
    int ret, access_ipiv, access_npiv, access_ppiv, accessU;
    struct starpu_task *task;
    int rankA = A->get_rankof(A, Am, An);

#if defined ( CHAMELEON_USE_MPI )
    if ( ( h % ib == 0 ) && ( h > 0 ) ) {
        starpu_mpi_cache_flush( options->sequence->comm, RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un) );
    }

    if ( rankA != A->myrank ) {
        return;
    }
#endif

    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zgetrf_blocked_diag, zgetrf_blocked_diag, zgetrf_blocked, 5 );

    access_ipiv = ( h == 0 )        ? STARPU_W    : STARPU_RW;
    access_npiv = ( h == pivot->n ) ? STARPU_R    : STARPU_REDUX;
    access_ppiv = ( h == 0 )        ? STARPU_NONE : STARPU_R;
    accessU     = STARPU_RW;
    if ( h == 0 ) {
        accessU = STARPU_NONE;
        /* U can be set after ppiv because they are both none together, so it won't shift the buffers */
    }
    else if ( h%ib == 0 ) {
        accessU = STARPU_R;
    }
    else if ( ( h%ib == 1 ) || ( ib == 1 ) ) {
        accessU = STARPU_W;
    }

    /*
     * Register the data handles, no exchange needed
     */
    starpu_cham_exchange_init_params( options, &params, rankA );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( A, ChamComplexDouble, Am, An ),     STARPU_RW );
    starpu_cham_register_descr( &nbdata, descrs, RUNTIME_ipiv_getaddr( ipiv, An),               access_ipiv );
    starpu_cham_register_descr( &nbdata, descrs, RUNTIME_pivot_getaddr( pivot, rankA, An, h ),   access_npiv );
    starpu_cham_register_descr( &nbdata, descrs, RUNTIME_pivot_getaddr( pivot, rankA, An, h-1 ), access_ppiv );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( U, ChamComplexDouble, Um, Un ),     accessU );

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    clargs = malloc( sizeof( struct cl_zgetrf_blocked_args_s ) );
    clargs->m        = m;
    clargs->n        = n;
    clargs->h        = h;
    clargs->m0       = m0;
    clargs->ib       = ib;
    clargs->sequence = options->sequence;
    clargs->request  = options->request;

    task->cl_arg      = clargs;
    task->cl_arg_size = sizeof( struct cl_zgetrf_blocked_args_s );
    task->cl_arg_free = 1;

    /* Set common parameters */
    starpu_cham_task_set_options( options, task, nbdata, descrs, cl_zgetrf_blocked_diag_callback );

    /* Flops */
    task->flops = flops_zgetrf_blocked_diag( m, n, h, ib );

    /* Refine name */
    task->name = chameleon_codelet_name( cl_name, 1, A->get_blktile( A, Am, An ) );

    ret = starpu_task_submit( task );
    if ( ret == -ENODEV ) {
        task->destroy = 0;
        starpu_task_destroy( task );
        chameleon_error( "INSERT_TASK_zgetrf_blocked_diag", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_blocked_offdiag_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zgetrf_blocked_args_s *clargs = (struct cl_zgetrf_blocked_args_s *)cl_arg;
    CHAM_tile_t           *tileA;
    CHAM_tile_t           *tileU;
    cppi_interface_t      *nextpiv;
    cppi_interface_t      *prevpiv;
    CHAMELEON_Complex64_t *U   = NULL;
    int                    ldu = -1;

    tileA   = cti_interface_get(descr[0]);
    nextpiv = (cppi_interface_t*) descr[1];
    prevpiv = (cppi_interface_t*) descr[2];
    if ( descr[3] != NULL ) {
        tileU = cti_interface_get(descr[3]);
        U     = CHAM_tile_get_ptr( tileU );
        ldu   = tileU->ld;
    }

    if ( clargs->h > 0 ) {
        cppi_display_dbg( prevpiv, stderr, "Prevpiv offdiag before call: " );
    }
    if ( clargs->h < tileA->n ) {
        cppi_display_dbg( nextpiv, stderr, "Nextpiv offdiag before call: " );
    }

    nextpiv->h = clargs->h; /* Initialize in case it uses a copy */
    nextpiv->has_diag = chameleon_max( -1, nextpiv->has_diag);

    coreblas_kernel_trace( tileA );
    CORE_zgetrf_panel_offdiag( clargs->m, clargs->n, clargs->h, clargs->m0, clargs->ib,
                               CHAM_tile_get_ptr(tileA), tileA->ld,
                               U, ldu,
                               &(nextpiv->pivot), &(prevpiv->pivot) );

    if ( clargs->h > 0 ) {
        cppi_display_dbg( prevpiv, stderr, "Prevpiv offdiag after call: " );
    }
    if ( clargs->h < tileA->n ) {
        cppi_display_dbg( nextpiv, stderr, "Nextpiv offdiag after call: " );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgetrf_blocked_offdiag, cl_zgetrf_blocked_offdiag_cpu_func)

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zgetrf_blocked_offdiag( const RUNTIME_option_t *options,
                                         int m, int n, int h, int m0, int ib,
                                         CHAM_desc_t *A, int Am, int An,
                                         CHAM_desc_t *U, int Um, int Un,
                                         CHAM_desc_pivot_t *pivot )
{
#if !defined(HAVE_STARPU_NONE_NONZERO)
    /* STARPU_NONE can't be equal to 0 */
    fprintf( stderr, "INSERT_TASK_zgetrf_blocked_diag: STARPU_NONE can not be equal to 0\n" );
    assert( 0 );
#endif
    int access_npiv = ( h == pivot->n )        ? STARPU_R    : STARPU_REDUX;
    int access_ppiv = ( h == 0 )               ? STARPU_NONE : STARPU_R;
    int accessU     = ((h%ib == 0) && (h > 0)) ? STARPU_R : STARPU_NONE;
    int rankA       = A->get_rankof(A, Am, An);

#if defined ( CHAMELEON_USE_MPI )
    if ( rankA != A->myrank ) {
        if ( ( accessU != STARPU_NONE ) &&
             ( A->myrank == A->get_rankof( A, An, An ) ) )
        {
            starpu_mpi_get_data_on_node_detached( options->sequence->comm,
                                                  RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un),
                                                  rankA, NULL, NULL );
        }
        return;
    }
#endif

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW( A, Am, An );
    if ((h%ib == 0) && (h > 0)) {
        CHAMELEON_ACCESS_R( U, Um, Un );
    }
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Set codelet parameters */
    struct cl_zgetrf_blocked_args_s *clargs;
    clargs = malloc( sizeof( struct cl_zgetrf_blocked_args_s ) );
    clargs->m        = m;
    clargs->n        = n;
    clargs->h        = h;
    clargs->m0       = m0;
    clargs->ib       = ib;
    clargs->sequence = options->sequence;
    clargs->request  = options->request;

    void (*callback)(void*) = options->profiling ? cl_zgetrf_blocked_offdiag_callback : NULL;
    const char *cl_name = "zgetrf_blocked_offdiag";

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 1,
                                      A->get_blktile( A, Am, An ) );

    rt_starpu_insert_task(
        &cl_zgetrf_blocked_offdiag,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zgetrf_blocked_args_s),

        /* Task handles */
        STARPU_RW,                RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        access_npiv,              RUNTIME_pivot_getaddr( pivot, rankA, An, h ),
        access_ppiv,              RUNTIME_pivot_getaddr( pivot, rankA, An, h-1 ),
        accessU,                  RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_NAME,              cl_name,
        0 );
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zgetrf_blocked_offdiag( const RUNTIME_option_t *options,
                                         int m, int n, int h, int m0, int ib,
                                         CHAM_desc_t *A, int Am, int An,
                                         CHAM_desc_t *U, int Um, int Un,
                                         CHAM_desc_pivot_t *pivot )
{
    int ret;
    struct starpu_task *task;
    int rankA       = A->get_rankof(A, Am, An);
    int access_npiv = ( h == pivot->n ) ? STARPU_R    : STARPU_REDUX;
    int access_ppiv = ( h == 0 )       ? STARPU_NONE : STARPU_R;
    int accessU     = ((h%ib == 0) && (h > 0)) ? STARPU_R : STARPU_NONE;

#if defined ( CHAMELEON_USE_MPI )
    if ( rankA != A->myrank ) {
        if ( ( accessU != STARPU_NONE ) &&
             ( A->myrank == A->get_rankof( A, An, An ) ) )
        {
            starpu_mpi_get_data_on_node_detached( options->sequence->comm,
                                                  RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un),
                                                  rankA, NULL, NULL );
        }
        return;
    }
#endif

    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zgetrf_blocked_offdiag, zgetrf_blocked_offdiag, zgetrf_blocked, 4 );

    /*
     * Register the data handles, exchange needed only for U
     */
    starpu_cham_exchange_init_params( options, &params, rankA );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( A, ChamComplexDouble, Am, An ),     STARPU_RW );
    starpu_cham_register_descr( &nbdata, descrs, RUNTIME_pivot_getaddr( pivot, rankA, An, h ),   access_npiv );
    starpu_cham_register_descr( &nbdata, descrs, RUNTIME_pivot_getaddr( pivot, rankA, An, h-1 ), access_ppiv );
    starpu_cham_exchange_handle_before_execution( options, &params, &nbdata, descrs,
                                                  RTBLKADDR( U, ChamComplexDouble, Um, Un ),
                                                  accessU );

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    clargs = malloc( sizeof( struct cl_zgetrf_blocked_args_s ) );
    clargs->m        = m;
    clargs->n        = n;
    clargs->h        = h;
    clargs->m0       = m0;
    clargs->ib       = ib;
    clargs->sequence = options->sequence;
    clargs->request  = options->request;

    task->cl_arg      = clargs;
    task->cl_arg_size = sizeof( struct cl_zgetrf_blocked_args_s );
    task->cl_arg_free = 1;

    /* Set common parameters */
    starpu_cham_task_set_options( options, task, nbdata, descrs, cl_zgetrf_blocked_offdiag_callback );

    /* Flops */
    task->flops = flops_zgetrf_blocked_offdiag( m, n, h, ib );

    /* Refine name */
    task->name = chameleon_codelet_name( cl_name, 1, A->get_blktile( A, Am, An ) );

    ret = starpu_task_submit( task );
    if ( ret == -ENODEV ) {
        task->destroy = 0;
        starpu_task_destroy( task );
        chameleon_error( "INSERT_TASK_zgetrf_blocked_offdiag", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */

#if !defined(CHAMELEON_SIMULATION)
static const CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;

static void cl_zgetrf_blocked_trsm_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zgetrf_blocked_args_s *clargs = (struct cl_zgetrf_blocked_args_s *)cl_arg;
    CHAM_tile_t           *tileU;
    cppi_interface_t      *prevpiv;
    CHAMELEON_Complex64_t *U;
    int                    ldu;

    tileU   = cti_interface_get(descr[0]);
    prevpiv = (cppi_interface_t*) descr[1];
    U       = CHAM_tile_get_ptr( tileU );
    ldu     = tileU->ld;

    coreblas_kernel_trace( tileU );

    /* Copy the final max line of the block and solve */
    cblas_zcopy( clargs->n, prevpiv->pivot.pivrow, 1,
                            U +  clargs->m - 1, ldu );

    if ( ( clargs->n - clargs->h ) > 0 ) {
        cblas_ztrsm( CblasColMajor,
                     CblasLeft, CblasLower,
                     CblasNoTrans, CblasUnit,
                     clargs->ib, clargs->n - clargs->h,
                     CBLAS_SADDR(zone), U + (clargs->h-clargs->ib) * ldu, ldu,
                                        U +  clargs->h             * ldu, ldu );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgetrf_blocked_trsm, cl_zgetrf_blocked_trsm_cpu_func)

#if defined(CHAMELEON_STARPU_USE_INSERT)

void INSERT_TASK_zgetrf_blocked_trsm( const RUNTIME_option_t *options,
                                      int m, int n, int h, int ib,
                                      CHAM_desc_t *U, int Um, int Un,
                                      CHAM_desc_pivot_t *pivot )
{
    void (*callback)(void*) = options->profiling ? cl_zgetrf_blocked_trsm_callback : NULL;
    const char *cl_name = "zgetrf_blocked_trsm";
    int rankU = U->get_rankof(U, Um, Un);

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW( U, Um, Un );
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 1,
                                      U->get_blktile( U, Um, Un ) );

    if ( U->myrank != U->get_rankof(U, Um, Un) ) {
        return;
    }

    /* Set codelet parameters */
    struct cl_zgetrf_blocked_args_s *clargs;
    clargs = malloc( sizeof( struct cl_zgetrf_blocked_args_s ) );
    clargs->m  = m;
    clargs->n  = n;
    clargs->h  = h;
    clargs->ib = ib;

    rt_starpu_insert_task(
        &cl_zgetrf_blocked_trsm,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zgetrf_blocked_args_s),

        /* Task handles */
        STARPU_RW,                RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un),
        STARPU_R,                 RUNTIME_pivot_getaddr( pivot, rankU, Un, h-1 ),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_NAME,              cl_name,
        0 );
}

#else /* defined(CHAMELEON_STARPU_USE_INSERT) */

void INSERT_TASK_zgetrf_blocked_trsm( const RUNTIME_option_t *options,
                                      int m, int n, int h, int ib,
                                      CHAM_desc_t *U, int Um, int Un,
                                      CHAM_desc_pivot_t *pivot )
{
    int ret;
    struct starpu_task *task;
    int rankU = U->get_rankof(U, Um, Un);

    if ( U->myrank != rankU ) {
        return;
    }

    INSERT_TASK_COMMON_PARAMETERS_EXTENDED( zgetrf_blocked_trsm, zgetrf_blocked_trsm, zgetrf_blocked, 2 );

    /*
     * Register the data handles, no exchange needed
     */
    starpu_cham_exchange_init_params( options, &params, rankU );
    starpu_cham_register_descr( &nbdata, descrs, RTBLKADDR( U, ChamComplexDouble, Um, Un ),     STARPU_RW );
    starpu_cham_register_descr( &nbdata, descrs, RUNTIME_pivot_getaddr( pivot, rankU, Un, h-1 ), STARPU_R  );

    task = starpu_task_create();
    task->cl = cl;

    /* Set codelet parameters */
    clargs = malloc( sizeof( struct cl_zgetrf_blocked_args_s ) );
    clargs->m  = m;
    clargs->n  = n;
    clargs->h  = h;
    clargs->ib = ib;

    task->cl_arg      = clargs;
    task->cl_arg_size = sizeof( struct cl_zgetrf_blocked_args_s );
    task->cl_arg_free = 1;

    /* Set common parameters */
    starpu_cham_task_set_options( options, task, nbdata, descrs, cl_zgetrf_blocked_trsm_callback );

    /* Flops */
    task->flops = flops_zgetrf_trsm( m, n, h, ib );

    /* Refine name */
    task->name = chameleon_codelet_name( cl_name, 1, U->get_blktile( U, Um, Un ) );

    ret = starpu_task_submit( task );
    if ( ret == -ENODEV ) {
        task->destroy = 0;
        starpu_task_destroy( task );
        chameleon_error( "INSERT_TASK_zgetrf_blocked_diag", "Failed to submit the task to StarPU" );
        return;
    }
    starpu_cham_task_exchange_data_after_execution( options, params, nbdata, descrs );
}

#endif /* defined(CHAMELEON_STARPU_USE_INSERT) */
