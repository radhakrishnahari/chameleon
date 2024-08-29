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
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"
#include <coreblas/cblas_wrapper.h>

#if defined(CHAMELEON_USE_MPI)
struct cl_redux_args_t {
    int h;
    int n;
    int k;
};

static void cl_zipiv_allreduce_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_redux_args_t *clargs      = (struct cl_redux_args_t *) cl_arg;
    cppi_interface_t       *cppi_me     = ((cppi_interface_t *) descr[0]);
    cppi_interface_t       *cppi_src    = ((cppi_interface_t *) descr[1]);
    CHAM_pivot_t           *nextpiv_me  = &(cppi_me->pivot);
    CHAM_pivot_t           *nextpiv_src = &(cppi_src->pivot);
    CHAMELEON_Complex64_t  *pivrow_me   = (CHAMELEON_Complex64_t *)(nextpiv_me->pivrow);
    CHAMELEON_Complex64_t  *pivrow_src  = (CHAMELEON_Complex64_t *)(nextpiv_src->pivrow);

    cppi_display_dbg( cppi_me,  stderr, "Global redux Inout: ");
    cppi_display_dbg( cppi_src, stderr, "Global redux Input: ");

    assert( cppi_me->n         == cppi_src->n         );
    assert( cppi_me->h         == cppi_src->h         );
    assert( cppi_me->flttype   == cppi_src->flttype   );
    assert( cppi_me->arraysize == cppi_src->arraysize );

    if ( cabs( pivrow_src[ clargs->h ] ) > cabs( pivrow_me[ clargs->h ] ) ) {
        nextpiv_me->blkm0  = nextpiv_src->blkm0;
        nextpiv_me->blkidx = nextpiv_src->blkidx;
        cblas_zcopy( clargs->n, pivrow_src, 1, pivrow_me, 1 );
    }

    /* Let's copy the diagonal row if needed */
    if ( ( cppi_src->has_diag == 1 ) &&
         ( cppi_me->has_diag  == -1 ) )
    {
        cblas_zcopy( clargs->n, nextpiv_src->diagrow, 1, nextpiv_me->diagrow, 1 );
        assert( cppi_src->arraysize == clargs->n * sizeof(CHAMELEON_Complex64_t) );
        cppi_me->has_diag = 1;
    }

    cppi_display_dbg( cppi_me,  stderr, "Global redux Inout(After): ");
}

CODELETS_CPU( zipiv_allreduce, cl_zipiv_allreduce_cpu_func )

void
INSERT_TASK_zipiv_allreduce_send( CHAM_ipiv_t *ipiv,
                                  int          me,
                                  int          dst,
                                  int          k,
                                  int          h,
                                  const RUNTIME_option_t *options )
{
    rt_starpu_insert_task(
        NULL,
        STARPU_EXECUTE_ON_NODE, dst,
        STARPU_R,               RUNTIME_pivot_getaddr( ipiv, me, k, h ),
        STARPU_PRIORITY,        options->priority,
        0 );
}

void
INSERT_TASK_zipiv_allreduce_recv( CHAM_ipiv_t *ipiv,
                                  int          me,
                                  int          src,
                                  int          k,
                                  int          h,
                                  int          n,
                                  const RUNTIME_option_t *options )
{
    struct cl_redux_args_t *clargs;
    clargs = malloc( sizeof( struct cl_redux_args_t ) );
    clargs->h = h;
    clargs->n = n;
    clargs->k = k;

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

void INSERT_TASK_zipiv_allreduce( CHAM_desc_t            *A,
                                  const RUNTIME_option_t *options,
                                  CHAM_ipiv_t            *ipiv,
                                  int                    *proc_involved,
                                  int                     k,
                                  int                     h,
                                  int                     n )
{
    int np_involved   = chameleon_min( A->p, A->mt - k);
    int np_iter       = np_involved;
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

            INSERT_TASK_zipiv_allreduce_send( ipiv, A->myrank, p_send, k, h,    options );
            INSERT_TASK_zipiv_allreduce_recv( ipiv, A->myrank, p_recv, k, h, n, options );

            shift   = shift << 1;
            np_iter = chameleon_ceil( np_iter, 2 );
        }
    }
}
#else
void INSERT_TASK_zipiv_allreduce( CHAM_desc_t            *A,
                                  const RUNTIME_option_t *options,
                                  CHAM_ipiv_t            *ipiv,
                                  int                    *proc_involved,
                                  int                     k,
                                  int                     h,
                                  int                     n )
{
    if ( h > 0 ) {
        starpu_data_invalidate_submit( RUNTIME_pivot_getaddr( ipiv, A->myrank, k, h-1 ) );
    }

    (void)options;
    (void)proc_involved;
    (void)n;
}
#endif
