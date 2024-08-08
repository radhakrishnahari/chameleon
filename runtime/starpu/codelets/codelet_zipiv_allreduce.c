/**
 *
 * @file starpu/codelet_zipiv_allreduce.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets to do the reduction
 *
 * @version 1.3.0
 * @author Alycia Lisito
 * @date 2024-06-11
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if defined(CHAMELEON_USE_MPI)
struct cl_redux_args_t {
    int h;
    int n;
};

struct cl_redux_mpi_args_t {
    int       h;
    int       n;
    int       nb;
    int       k;
    int       p;
    int       q;
    int       rankA;
    int       np_involved;
    MPI_Comm  comm_panel;
    int       proc_involved[1];
};

struct zipiv_interf_s
{
    int                   has_diag;
    int                   h;
    int                   n;
    int                   blkm0;
    int                   blkidx;
    CHAMELEON_Complex64_t pivrow[1];
};

static void
zipiv_allreduce_copy( cppi_interface_t      *cppi_int,
                      struct zipiv_interf_s *cppi_exc,
                      int                    n )
{
    cppi_exc->has_diag = cppi_int->has_diag;
    cppi_exc->h        = cppi_int->h;
    assert( cppi_int->n >= n );
    cppi_exc->n        = n;
    cppi_exc->blkm0    = cppi_int->pivot.blkm0;
    cppi_exc->blkidx   = cppi_int->pivot.blkidx;

    CHAMELEON_Complex64_t  *pivrow  = (CHAMELEON_Complex64_t *)(cppi_int->pivot.pivrow);
    CHAMELEON_Complex64_t  *diagrow = (CHAMELEON_Complex64_t *)(cppi_int->pivot.diagrow);
    for ( int k = 0; k < cppi_exc->n; k++ ) {
        cppi_exc->pivrow[ k ]               = pivrow[ k ];
        cppi_exc->pivrow[ k + cppi_exc->n ] = diagrow[ k ];
    }
}

static void
zipiv_allreduce_cpu_func( cppi_interface_t      *cppi_me,
                          CHAMELEON_Complex64_t *pivrow_src,
                          CHAMELEON_Complex64_t *diagrow_src,
                          int                    has_diag_src,
                          int                    h_src,
                          int                    n_src,
                          int                    blkm0_src,
                          int                    blkidx_src,
                          int                    h,
                          int                    n )
{
    CHAM_pivot_t          *nextpiv_me = &(cppi_me->pivot);
    CHAMELEON_Complex64_t *pivrow_me  = (CHAMELEON_Complex64_t *)(nextpiv_me->pivrow);

    cppi_display_dbg( cppi_me,  stderr, "Global redux Inout: ");

    assert( cppi_me->n >= n_src );
    assert( cppi_me->h == h_src );

    if ( cabs( pivrow_src[ h ] ) > cabs( pivrow_me[ h ] ) ) {
        nextpiv_me->blkm0  = blkm0_src;
        nextpiv_me->blkidx = blkidx_src;
        cblas_zcopy( n, pivrow_src, 1, pivrow_me, 1 );
    }

    /* Let's copy the diagonal row if needed */
    if ( ( has_diag_src == 1 ) &&
         ( cppi_me->has_diag  == -1 ) )
    {
        cblas_zcopy( n, diagrow_src, 1, nextpiv_me->diagrow, 1 );
        cppi_me->has_diag = 1;
    }

    cppi_display_dbg( cppi_me,  stderr, "Global redux Inout(After): ");
}

static void
zipiv_allreduce_mpi_cpu_func( cppi_interface_t      *cppi_me,
                              struct zipiv_interf_s *cppi_src,
                              int                    h,
                              int                    n )
{
    CHAMELEON_Complex64_t *pivrow_src   = cppi_src->pivrow;
    CHAMELEON_Complex64_t *diagrow_src  = cppi_src->pivrow + cppi_src->n;
    int                    has_diag_src = cppi_src->has_diag;
    int                    h_src        = cppi_src->h;
    int                    n_src        = cppi_src->n;
    int                    blkm0_src    = cppi_src->blkm0;
    int                    blkidx_src   = cppi_src->blkidx;

    assert( n_src == n );
    zipiv_allreduce_cpu_func( cppi_me, pivrow_src, diagrow_src, has_diag_src,
                              h_src, n_src, blkm0_src, blkidx_src, h, n );
}

static void
zipiv_allreduce_spu_cpu_func( cppi_interface_t *cppi_me,
                              cppi_interface_t *cppi_src,
                              int               h,
                              int               n )
{
    CHAM_pivot_t          *nextpiv_src  = &(cppi_src->pivot);
    CHAMELEON_Complex64_t *pivrow_src   = (CHAMELEON_Complex64_t *)(nextpiv_src->pivrow);
    CHAMELEON_Complex64_t *diagrow_src  = (CHAMELEON_Complex64_t *)(nextpiv_src->diagrow);
    int                    has_diag_src = cppi_src->has_diag;
    int                    h_src        = cppi_src->h;
    int                    n_src        = cppi_src->n;
    int                    blkm0_src    = nextpiv_src->blkm0;
    int                    blkidx_src   = nextpiv_src->blkidx;

    zipiv_allreduce_cpu_func( cppi_me, pivrow_src, diagrow_src, has_diag_src,
                              h_src, n_src, blkm0_src, blkidx_src, h, n );
}

static void
cl_zipiv_allreduce_mpi_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_redux_mpi_args_t *clargs = (struct cl_redux_mpi_args_t *) cl_arg;
    cppi_interface_t      *cppi_send_int = ((cppi_interface_t *)descr[0]);
    int                    h             = clargs->h;
    int                    n             = clargs->n;
    int                    nb            = clargs->nb;
    int                    k             = clargs->k;
    int                    q             = clargs->q;
    int                    p             = clargs->p;
    int                    rankA         = clargs->rankA;
    int                    np_involved   = clargs->np_involved;
    int                    np_iter       = clargs->np_involved;
    int                   *proc_involved = clargs->proc_involved;
    MPI_Comm               comm_panel    = clargs->comm_panel;
    int                    shift         = 1;
    int                    p_recv, p_send, me, size;
    int                    tag_send, tag_recv, where_i_am;
    struct zipiv_interf_s *cppi_send, *cppi_recv;
    MPI_Request            request;
    MPI_Status             status;

    for( me = 0; me < np_involved; me++ ) {
        if ( proc_involved[me] == rankA ) {
            break;
        }
    }
    assert( me < np_involved );

    size      = sizeof( struct zipiv_interf_s ) + (2*n-1) * sizeof( CHAMELEON_Complex64_t );
    cppi_send = malloc( size );
    cppi_recv = malloc( size );
    tag_send  = rankA / q + k * nb * p + h * p;

    while ( np_iter > 1 ) {
        p_send = proc_involved[ ( me + shift               ) % np_involved ] / q;
        p_recv = proc_involved[ ( me - shift + np_involved ) % np_involved ] / q;
        tag_recv = p_recv + k * nb * p + h * p;

        zipiv_allreduce_copy( cppi_send_int, cppi_send, n );
        MPI_Isend( cppi_send, size, MPI_BYTE, p_send, tag_send, comm_panel, &request );
        MPI_Recv( cppi_recv, size, MPI_BYTE, p_recv, tag_recv, comm_panel, MPI_STATUS_IGNORE );
        zipiv_allreduce_mpi_cpu_func( cppi_send_int, cppi_recv, h, n );

        shift   = shift << 1;
        np_iter = chameleon_ceil( np_iter, 2 );
        MPI_Wait( &request, &status );
    }

    free( cppi_recv );
    free( cppi_send );
}

CODELETS_CPU( zipiv_allreduce_mpi, cl_zipiv_allreduce_mpi_cpu_func )

static void
cl_zipiv_allreduce_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_redux_args_t *clargs   = (struct cl_redux_args_t *) cl_arg;
    cppi_interface_t       *cppi_me  = ((cppi_interface_t *) descr[0]);
    cppi_interface_t       *cppi_src = ((cppi_interface_t *) descr[1]);

    zipiv_allreduce_spu_cpu_func( cppi_me, cppi_src, clargs->h, clargs->n );
}

CODELETS_CPU( zipiv_allreduce, cl_zipiv_allreduce_cpu_func )

static void
INSERT_TASK_zipiv_allreduce_send( const RUNTIME_option_t *options,
                                  CHAM_ipiv_t            *ipiv,
                                  int                     me,
                                  int                     dst,
                                  int                     k,
                                  int                     h )
{
    rt_starpu_insert_task(
        NULL,
        STARPU_EXECUTE_ON_NODE, dst,
        STARPU_R,               RUNTIME_pivot_getaddr( ipiv, me, k, h ),
        STARPU_PRIORITY,        options->priority,
        0 );
}

static void
INSERT_TASK_zipiv_allreduce_recv( const RUNTIME_option_t *options,
                                  CHAM_ipiv_t            *ipiv,
                                  int                     me,
                                  int                     src,
                                  int                     k,
                                  int                     h,
                                  int                     n )
{
    struct cl_redux_args_t *clargs;
    clargs    = malloc( sizeof( struct cl_redux_args_t ) );
    clargs->h = h;
    clargs->n = n;

    rt_starpu_insert_task(
        &cl_zipiv_allreduce,
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_redux_args_t),
        STARPU_RW,                RUNTIME_pivot_getaddr( ipiv, me,  k, h ),
        STARPU_R,                 RUNTIME_pivot_getaddr( ipiv, src, k, h ),
        STARPU_EXECUTE_ON_NODE,   me,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_PRIORITY,          options->priority,
        0 );
    starpu_mpi_cache_flush( options->sequence->comm, RUNTIME_pivot_getaddr( ipiv, src, k, h ) );
}

static void
zipiv_allreduce_chameleon_starpu_task( const RUNTIME_option_t *options,
                                       CHAM_desc_t            *A,
                                       CHAM_ipiv_t            *ipiv,
                                       int                    *proc_involved,
                                       int                     k,
                                       int                     h,
                                       int                     n )
{
    int np_involved = chameleon_min( chameleon_desc_datadist_get_iparam(A, 0), A->mt - k);
    int np_iter     = np_involved;
    int p_recv, p_send, me;
    int shift = 1;

    if ( h > 0 ) {
        starpu_data_invalidate_submit( RUNTIME_pivot_getaddr( ipiv, A->myrank, k, h-1 ) );
    }
    if ( h >= ipiv->n ) {
        return;
    }

    if ( np_involved == 1 ) {
        assert( proc_involved[0] == A->myrank );
    }
    else {
        for( me = 0; me < np_involved; me++ ) {
            if ( proc_involved[me] == A->myrank ) {
                break;
            }
        }
        assert( me < np_involved );
        while ( np_iter > 1 ) {
            p_send = proc_involved[ ( me + shift               ) % np_involved ];
            p_recv = proc_involved[ ( me - shift + np_involved ) % np_involved ];

            INSERT_TASK_zipiv_allreduce_send( options, ipiv, A->myrank, p_send, k, h    );
            INSERT_TASK_zipiv_allreduce_recv( options, ipiv, A->myrank, p_recv, k, h, n );

            shift   = shift << 1;
            np_iter = chameleon_ceil( np_iter, 2 );
        }
    }
}

static void
zipiv_allreduce_chameleon_starpu_mpi_task( CHAM_desc_t            *A,
                                           const RUNTIME_option_t *options,
                                           CHAM_ipiv_t            *ipiv,
                                           int                    *proc_involved,
                                           int                     k,
                                           int                     h,
                                           int                     n,
                                           MPI_Comm                comm_panel )
{
    int    np_involved = chameleon_min( chameleon_desc_datadist_get_iparam(A, 0), A->mt - k);
    size_t size        = ( np_involved - 1 ) * sizeof( int );
    struct cl_redux_mpi_args_t *clargs;
    int i;
    if ( h > 0 ) {
        starpu_data_invalidate_submit( RUNTIME_pivot_getaddr( ipiv, A->myrank, k, h-1 ) );
    }
    if ( h >= ipiv->n ) {
        return;
    }
    if ( np_involved == 1 ) {
        assert( proc_involved[0] == A->myrank );
        return;
    }

    clargs = malloc( sizeof( struct cl_redux_mpi_args_t ) + size );
    clargs->h             = h;
    clargs->n             = n;
    clargs->nb            = A->nb;
    clargs->k             = k;
    clargs->p             = chameleon_desc_datadist_get_iparam(A, 0);
    clargs->q             = chameleon_desc_datadist_get_iparam(A, 1);
    clargs->rankA         = A->myrank;
    clargs->np_involved   = np_involved;
    clargs->comm_panel    = comm_panel;
    for ( i = 0; i < np_involved; i ++ ) {
        clargs->proc_involved[i] = proc_involved[i];
    }

    rt_starpu_insert_task(
        &cl_zipiv_allreduce_mpi,
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_redux_mpi_args_t) + size,
        STARPU_RW,                RUNTIME_pivot_getaddr( ipiv, A->myrank, k, h ),
        STARPU_EXECUTE_ON_NODE,   A->myrank,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_PRIORITY,          options->priority,
        0 );
}

void
INSERT_TASK_zipiv_allreduce( const RUNTIME_option_t *options,
                             CHAM_desc_t            *A,
                             CHAM_ipiv_t            *ipiv,
                             int                     k,
                             int                     h,
                             int                     n,
                             void                   *ws )
{
    struct chameleon_pzgetrf_s *tmp = (struct chameleon_pzgetrf_s *)ws;
    cham_getrf_allreduce_t alg = tmp->alg_allreduce;
    switch( alg ) {
    case ChamStarPUMPITasks:
        zipiv_allreduce_chameleon_starpu_mpi_task( A, options, ipiv, tmp->proc_involved, k, h, n, tmp->comm_panel );
    break;
    case ChamStarPUTasks:
    default:
        zipiv_allreduce_chameleon_starpu_task( options, A, ipiv, tmp->proc_involved, k, h, n );
    }
}
#else
void
INSERT_TASK_zipiv_allreduce( const RUNTIME_option_t *options,
                             CHAM_desc_t            *A,
                             CHAM_ipiv_t            *ipiv,
                             int                     k,
                             int                     h,
                             int                     n,
                             void                   *ws )
{
    if ( h > 0 ) {
        starpu_data_invalidate_submit( RUNTIME_pivot_getaddr( ipiv, A->myrank, k, h-1 ) );
    }

    (void)options;
    (void)n;
    (void)ws;
}
#endif
