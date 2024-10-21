/**
 *
 * @file starpu/codelet_zperm_allreduce.c
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
#include <coreblas/cblas_wrapper.h>

#if defined(CHAMELEON_USE_MPI)
struct cl_redux_args_t {
    int tempmm;
    int n;
    int p;
    int q;
    int p_first;
    int me;
    int shift;
    int np_inv;
};

struct cl_redux_mpi_args_t {
    int       mb;
    int       nt;
    int       k;
    int       n;
    int       p;
    int       q;
    int       tag;
    int       myrank;
    int       np_involved;
    MPI_Comm  comm_panel;
    int       proc_involved[1];
};

static void
zperm_allreduce_buffer_send( int          *perm,
                             CHAM_tile_t  *tileU,
                             int           myrank,
                             int           mb_full,
                             int           n,
                             int           p,
                             int           q,
                             int          *idx,
                             char         *send_data )
{
    int64_t  i, m;
    int64_t  count = 0;
    char    *send  = send_data;
    int      mb    = tileU->m;
    int      nb    = tileU->n;
    CHAMELEON_Complex64_t *ptr_U = CHAM_tile_get_ptr( tileU );

    for ( i = 0; i < mb; i++ ) {
        m = perm[ i ] / mb_full;
        if ( (m % p) * q + (n % q) == myrank ) {
            count++;
        }
    }

    memcpy( send_data, &count, sizeof(int64_t) );
    send_data += sizeof(int64_t);
    for ( i = 0; i < mb; i++ ) {
        m = perm[ i ] / mb_full;
        if ( (m % p) * q + (n % q) == myrank ) {
            memcpy( send_data, &i, sizeof(int64_t) );
            send_data += sizeof(int64_t);
            cblas_zcopy( nb, ptr_U + i, tileU->ld, (void *)send_data, 1 );
            send_data += nb * sizeof(CHAMELEON_Complex64_t);
            idx[ i ] = 1;
        }
        else {
            idx[ i ] = 0;
        }
    }
    send += ( count + 1 ) * sizeof(int64_t) + count * nb * sizeof(CHAMELEON_Complex64_t);
    assert( (char *)send == send_data );
}

static void
zperm_allreduce_copy( char        *recv_rows,
                      char        *send_data,
                      CHAM_tile_t *tileU,
                      int         *idx )
{
    int64_t  i, j;
    char *recv          = recv_rows;
    int64_t  count_recv = *(int64_t *)recv_rows;
    int64_t  count_send = *(int64_t *)send_data;
    int64_t  new_count  = count_recv + count_send;
    int      nb         = tileU->n;
    CHAMELEON_Complex64_t *ptr_U = CHAM_tile_get_ptr( tileU );

    recv_rows += sizeof(int64_t) ;
    memcpy( send_data, &new_count, sizeof(int64_t) );
    send_data += ( count_send + 1 ) * sizeof(int64_t) + count_send * nb * sizeof(CHAMELEON_Complex64_t);
    for ( j = 0; j < count_recv; j++ ) {
        i = *(int64_t *)recv_rows;
        recv_rows += sizeof(int64_t);
        if ( idx[ i ] == 0 ) {
            cblas_zcopy( nb, (void *)recv_rows, 1, ptr_U + i, tileU->ld );
            memcpy( send_data, &i, sizeof(int64_t) );
            send_data += sizeof(int64_t);
            memcpy( send_data, recv_rows, nb * sizeof(CHAMELEON_Complex64_t) );
            send_data += nb * sizeof(CHAMELEON_Complex64_t);
            idx[ i ] = 1;
        }
        recv_rows += nb * sizeof(CHAMELEON_Complex64_t);
    }
    recv += ( count_recv + 1 ) * sizeof(int64_t) + count_recv * nb * sizeof(CHAMELEON_Complex64_t);
    assert( (char *)recv == recv_rows );
}

static void
cl_zperm_allreduce_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_redux_args_t      *clargs     = (struct cl_redux_args_t *) cl_arg;
    const CHAM_tile_t           *tileUinout = cti_interface_get( descr[0] );
    const CHAM_tile_t           *tileUin    = cti_interface_get( descr[1] );
    const int                   *perm       = (int *)STARPU_VECTOR_GET_PTR( descr[2] );
    CHAMELEON_Complex64_t       *Uinout     = CHAM_tile_get_ptr( tileUinout );
    const CHAMELEON_Complex64_t *Uin        = CHAM_tile_get_ptr( tileUin );

    int tempmm  = clargs->tempmm;
    int n       = clargs->n;
    int p       = clargs->p;
    int q       = clargs->q;
    int p_first = clargs->p_first / q;
    int shift   = clargs->shift;
    int np      = clargs->np_inv;
    int me      = ( p <= np ) ? clargs->me / q : ( ( clargs->me / q ) - p_first + p ) % p;
    int nb      = tileUinout->n;
    int mb      = tileUinout->m;
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

static void
cl_zperm_allreduce_mpi_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_redux_mpi_args_t *clargs = (struct cl_redux_mpi_args_t *) cl_arg;
    CHAM_tile_t  *tileU         = cti_interface_get( descr[0] );
    int          *perm          = (int *)STARPU_VECTOR_GET_PTR( descr[1] );
    int           mb_full       = clargs->mb;
    int           nb            = tileU->n;
    int           mb            = tileU->m;
    int           nt            = clargs->nt;
    int           k             = clargs->k;
    int           n             = clargs->n;
    int           p             = clargs->p;
    int           q             = clargs->q;
    int           tag           = clargs->tag;
    int           myrank        = clargs->myrank;
    int           np_involved   = clargs->np_involved;
    int           np_iter       = clargs->np_involved;
    int          *proc_involved = clargs->proc_involved;
    MPI_Comm      comm_panel    = clargs->comm_panel;
    int           shift         = 1;
    int           p_recv, p_send, me, size, size_max, count;
    int           tag_send, tag_recv, where_i_am;
    int           idx[ mb ];
    char         *send_data;
    char         *recv_data;
    MPI_Request   request;
    MPI_Status    status;

    for( me = 0; me < np_involved; me++ ) {
        if ( proc_involved[me] == myrank ) {
            break;
        }
    }
    assert( me < np_involved );

    size_max  = sizeof( int64_t ) * ( mb + 1 ) + mb * nb * sizeof( CHAMELEON_Complex64_t );
    send_data = malloc( size_max );
    recv_data = malloc( size_max );
    memset( send_data, 0, size_max );
    memset( recv_data, 0, size_max );
    zperm_allreduce_buffer_send( perm, tileU, myrank, mb_full, n, p, q, idx, send_data );
    tag_send  = myrank / q + tag + k * nt * p + n * p;

    while ( np_iter > 1 ) {
        p_send   = proc_involved[ ( me + shift               ) % np_involved ] / q;
        p_recv   = proc_involved[ ( me - shift + np_involved ) % np_involved ] / q;
        tag_recv = p_recv + tag + k * nt * p + n * p;
        count    = *((int64_t *)send_data);
        size     = sizeof( int64_t ) * ( count + 1 ) + count * nb * sizeof( CHAMELEON_Complex64_t );

        MPI_Isend( send_data, size,    MPI_BYTE, p_send, tag_send, comm_panel, &request );
        MPI_Recv( recv_data, size_max, MPI_BYTE, p_recv, tag_recv, comm_panel, MPI_STATUS_IGNORE );
        MPI_Wait( &request, &status );
        zperm_allreduce_copy( recv_data, send_data, tileU, idx );

        shift   = shift << 1;
        np_iter = chameleon_ceil( np_iter, 2 );
    }

    free( send_data );
    free( recv_data );
}

CODELETS_CPU( zperm_allreduce_mpi, cl_zperm_allreduce_mpi_cpu_func )

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
    struct cl_redux_args_t *clargs;
    clargs = malloc( sizeof( struct cl_redux_args_t ) );
    clargs->tempmm  = tempmm;
    clargs->n       = n;
    clargs->p       = p;
    clargs->q       = q;
    clargs->p_first = p_first;
    clargs->me      = me;
    clargs->shift   = shift;
    clargs->np_inv  = np;

    rt_starpu_insert_task(
        &cl_zperm_allreduce,
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_redux_args_t),
        STARPU_RW,                RTBLKADDR(U, CHAMELEON_Complex64_t, me,  n),
        STARPU_R,                 RTBLKADDR(U, CHAMELEON_Complex64_t, src, n),
        STARPU_R,                 RUNTIME_perm_getaddr( ipiv, ipivk ),
        STARPU_EXECUTE_ON_NODE,   me,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_PRIORITY,          options->priority,
        0 );
    starpu_mpi_cache_flush( options->sequence->comm, RTBLKADDR(U, CHAMELEON_Complex64_t, src, n) );
}

static void
zperm_allreduce_chameleon_starpu_task( const RUNTIME_option_t     *options,
                                       const CHAM_desc_t          *A,
                                       CHAM_desc_t                *U,
                                       int                         Um,
                                       int                         Un,
                                       CHAM_ipiv_t                *ipiv,
                                       int                         ipivk,
                                       int                         k,
                                       int                         n,
                                       struct chameleon_pzgetrf_s *ws)
{
    int *proc_involved = ws->proc_involved;
    int  np_involved   = chameleon_min( chameleon_desc_datadist_get_iparam(A, 0), A->mt - k);
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
            INSERT_TASK_zperm_allreduce_recv( options, U, ipiv, ipivk, A->myrank, p_recv,
                                              n, k == (A->mt-1) ? A->m - k * A->mb : A->mb,
                                              chameleon_desc_datadist_get_iparam(A, 0), chameleon_desc_datadist_get_iparam(A, 1), shift, np_involved, p_first );

            shift   = shift << 1;
            np_iter = chameleon_ceil( np_iter, 2 );
        }
    }
}

void
zperm_allreduce_chameleon_starpu_mpi_task( const RUNTIME_option_t     *options,
                                           const CHAM_desc_t          *A,
                                           CHAM_ipiv_t                *ipiv,
                                           int                         ipivk,
                                           int                         k,
                                           int                         n,
                                           CHAM_desc_t                *U,
                                           int                         Um,
                                           int                         Un,
                                           struct chameleon_pzgetrf_s *ws )
{
    int    np_involved = chameleon_min( chameleon_desc_datadist_get_iparam(A, 0), A->mt - k);
    size_t size        = ( np_involved - 1 ) * sizeof( int );
    struct cl_redux_mpi_args_t *clargs;
    int i, size_task;

    if ( np_involved == 1 ) {
        assert( ws->proc_involved[0] == A->myrank );
        return;
    }

    clargs = malloc( sizeof( struct cl_redux_mpi_args_t ) + size );
    clargs->mb            = A->mb;
    clargs->nt            = A->nt;
    clargs->k             = k;
    clargs->n             = n;
    clargs->p             = chameleon_desc_datadist_get_iparam(A, 0);
    clargs->q             = chameleon_desc_datadist_get_iparam(A, 1);
    clargs->tag           = ( ( chameleon_min( A->mt, A->nt ) + 1 ) * A->nb + 1 ) * chameleon_desc_datadist_get_iparam(A, 0);
    clargs->myrank        = A->myrank;
    clargs->np_involved   = np_involved;
    clargs->comm_panel    = ws->comm_panel;
    for ( i = 0; i < np_involved; i ++ ) {
        clargs->proc_involved[i] = ws->proc_involved[i];
    }

    uint64_t tag = k * A->nt * chameleon_desc_datadist_get_iparam(A, 0) + n * chameleon_desc_datadist_get_iparam(A, 0) + A->myrank / chameleon_desc_datadist_get_iparam(A, 1);
    if ( ws->tag != -1 ) {
        starpu_tag_declare_deps( (starpu_tag_t)tag, 1, (starpu_tag_t)ws->tag );
    }
    ws->tag = tag;

    rt_starpu_insert_task(
        &cl_zperm_allreduce_mpi,
        STARPU_CL_ARGS,           clargs, sizeof(struct cl_redux_mpi_args_t) + size,
        STARPU_RW,                RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un),
        STARPU_R,                 RUNTIME_perm_getaddr( ipiv, ipivk ),
        STARPU_TAG,               (starpu_tag_t)tag,
        STARPU_EXECUTE_ON_NODE,   A->myrank,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_PRIORITY,          options->priority,
        0 );
}

void
INSERT_TASK_zperm_allreduce( const RUNTIME_option_t *options,
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
    case ChamStarPUMPITasks:
        zperm_allreduce_chameleon_starpu_mpi_task( options, A, ipiv, ipivk, k, n, U, Um, Un, tmp );
    break;
    case ChamStarPUTasks:
    default:
        zperm_allreduce_chameleon_starpu_task( options, A, U, Um, Un, ipiv, ipivk, k, n, tmp );
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
                                       CHAM_ipiv_t            *ipiv,
                                       int                     ipivk,
                                       int                     myrank,
                                       int                     np,
                                       int                    *proc_involved )
{
    int p;

    for ( p = 0; p < np; p++ ) {
        if ( proc_involved[ p ] == myrank ) {
            continue;
        }
        starpu_mpi_get_data_on_node_detached( options->sequence->comm,
                                              RUNTIME_perm_getaddr( ipiv, ipivk ),
                                              proc_involved[ p ], NULL, NULL );
    }
}

void
INSERT_TASK_zperm_allreduce_send_invp( const RUNTIME_option_t *options,
                                       CHAM_ipiv_t            *ipiv,
                                       int                     ipivk,
                                       const CHAM_desc_t      *A,
                                       int                     k,
                                       int                     n )
{
    int b, rank;

    for ( b = k+1; (b < A->mt) && ((b-(k+1)) < chameleon_desc_datadist_get_iparam(A, 0)); b ++ ) {
        rank = A->get_rankof( A, b, n );
        if ( rank == A->myrank ) {
            continue;
        }
        starpu_mpi_get_data_on_node_detached( options->sequence->comm,
                                              RUNTIME_invp_getaddr( ipiv, ipivk ),
                                              rank, NULL, NULL );
    }
}

void
INSERT_TASK_zperm_allreduce_tag_free( )
{
    starpu_tag_clear();
}
#else
void
INSERT_TASK_zperm_allreduce_send_A( const RUNTIME_option_t *options,
                                    CHAM_desc_t            *A,
                                    int                     Am,
                                    int                     An,
                                    int                     myrank,
                                    int                     np,
                                    int                    *proc_involved  )
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
                                       CHAM_ipiv_t            *ipiv,
                                       int                     ipivk,
                                       int                     myrank,
                                       int                     np,
                                       int                    *proc_involved  )
{
    (void)options;
    (void)ipiv;
    (void)ipivk;
    (void)myrank;
    (void)np;
    (void)proc_involved;
}

void
INSERT_TASK_zperm_allreduce_send_invp( const RUNTIME_option_t *options,
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
INSERT_TASK_zperm_allreduce( const RUNTIME_option_t *options,
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

void
INSERT_TASK_zperm_allreduce_tag_free( )
{}
#endif
