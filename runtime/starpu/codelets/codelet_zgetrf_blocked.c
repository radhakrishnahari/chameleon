/**
 *
 * @file starpu/codelet_zgetrf_blocked.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

CHAMELEON_CL_CB( zgetrf_blocked_diag,    cti_handle_get_m(task->handles[0]), 0, 0, M )
CHAMELEON_CL_CB( zgetrf_blocked_offdiag, cti_handle_get_m(task->handles[0]), 0, 0, M )
CHAMELEON_CL_CB( zgetrf_blocked_trsm,    cti_handle_get_m(task->handles[0]), 0, 0, M )

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_blocked_diag_cpu_func(void *descr[], void *cl_arg)
{
    int                    m, n, h, m0, ib;
    RUNTIME_sequence_t    *sequence;
    RUNTIME_request_t     *request;
    CHAM_tile_t           *tileA;
    CHAM_tile_t           *tileU;
    int                   *ipiv;
    cppi_interface_t      *nextpiv;
    cppi_interface_t      *prevpiv;
    CHAMELEON_Complex64_t *U   = NULL;
    int                    ldu = -1;;

    starpu_codelet_unpack_args( cl_arg, &m, &n, &h, &m0, &ib,
                                &sequence, &request );

    tileA   = cti_interface_get(descr[0]);
    ipiv    = (int *)STARPU_VECTOR_GET_PTR(descr[1]);
    nextpiv = (cppi_interface_t*) descr[2];
    prevpiv = (cppi_interface_t*) descr[3];
    if ( descr[4] != NULL ) {
        tileU = cti_interface_get(descr[4]);
        U     = CHAM_tile_get_ptr( tileU );
        ldu   = tileU->ld;
    }

    if ( h > 0 ) {
        cppi_display_dbg( prevpiv, stderr, "Prevpiv before call: " );
    }
    if ( h < tileA->n ) {
        cppi_display_dbg( nextpiv, stderr, "Nextpiv before call: " );
    }

    /*
     * Make sure the nextpiv interface store the right information about the
     * column and diagonal row for the reduction
     */
    nextpiv->h        = h;
    nextpiv->has_diag = 1;

    coreblas_kernel_trace( tileA );
    CORE_zgetrf_panel_diag( m, n, h, m0, ib,
                            CHAM_tile_get_ptr( tileA ), tileA->ld,
                            U, ldu,
                            ipiv, &(nextpiv->pivot), &(prevpiv->pivot) );

    if ( h > 0 ) {
        cppi_display_dbg( prevpiv, stderr, "Prevpiv after call: " );
    }
    if ( h < tileA->n ) {
        cppi_display_dbg( nextpiv, stderr, "Nextpiv after call: " );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zgetrf_blocked_diag, cl_zgetrf_blocked_diag_cpu_func )

void INSERT_TASK_zgetrf_blocked_diag( const RUNTIME_option_t *options,
                                      int m, int n, int h, int m0, int ib,
                                      CHAM_desc_t *A, int Am, int An,
                                      CHAM_desc_t *U, int Um, int Un,
                                      CHAM_ipiv_t *ipiv )
{
    struct starpu_codelet *codelet = &cl_zgetrf_blocked_diag;
    void (*callback)(void*) = options->profiling ? cl_zgetrf_blocked_diag_callback : NULL;
    const char *cl_name = "zgetrf_blocked_diag";
    int rankA           = A->get_rankof(A, Am, An);
#if !defined(HAVE_STARPU_NONE_NONZERO)
    /* STARPU_NONE can't be equal to 0 */
    fprintf( stderr, "INSERT_TASK_zgetrf_percol_diag: STARPU_NONE can not be equal to 0\n" );
    assert( 0 );
#endif

#if defined ( CHAMELEON_USE_MPI )
    if ( ( h % ib == 0 ) && ( h > 0 ) ) {
        starpu_mpi_cache_flush( options->sequence->comm, RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un) );
    }

    if ( rankA != A->myrank ) {
        return;
    }
#endif

    int access_ipiv = ( h == 0 )       ? STARPU_W    : STARPU_RW;
    int access_npiv = ( h == ipiv->n ) ? STARPU_R    : STARPU_REDUX;
    int access_ppiv = ( h == 0 )       ? STARPU_NONE : STARPU_R;
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

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW( A, Am, An );
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 1,
                                      A->get_blktile( A, Am, An ) );

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,             &m,                   sizeof(int),
        STARPU_VALUE,             &n,                   sizeof(int),
        STARPU_VALUE,             &h,                   sizeof(int),
        STARPU_VALUE,             &m0,                  sizeof(int),
        STARPU_VALUE,             &ib,                  sizeof(int),
        STARPU_VALUE,             &(options->sequence), sizeof(RUNTIME_sequence_t*),
        STARPU_VALUE,             &(options->request),  sizeof(RUNTIME_request_t*),
        STARPU_RW,                RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        access_ipiv,              RUNTIME_ipiv_getaddr( ipiv, An ),
        access_npiv,              RUNTIME_pivot_getaddr( ipiv, rankA, An, h ),
        access_ppiv,              RUNTIME_pivot_getaddr( ipiv, rankA, An, h-1 ),
        accessU,                  RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un),
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_NAME,              cl_name,
        0 );
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_blocked_offdiag_cpu_func(void *descr[], void *cl_arg)
{
    int                    m, n, h, k, m0, ib;
    RUNTIME_sequence_t    *sequence;
    RUNTIME_request_t     *request;
    CHAM_tile_t           *tileA;
    CHAM_tile_t           *tileU;
    cppi_interface_t      *nextpiv;
    cppi_interface_t      *prevpiv;
    CHAMELEON_Complex64_t *U   = NULL;
    int                    ldu = -1;

    starpu_codelet_unpack_args( cl_arg, &m, &n, &h, &k, &m0, &ib, &sequence, &request );

    tileA   = cti_interface_get(descr[0]);
    nextpiv = (cppi_interface_t*) descr[1];
    prevpiv = (cppi_interface_t*) descr[2];
    if ( descr[3] != NULL ) {
        tileU = cti_interface_get(descr[3]);
        U     = CHAM_tile_get_ptr( tileU );
        ldu   = tileU->ld;
    }

    if ( h > 0 ) {
        cppi_display_dbg( prevpiv, stderr, "Prevpiv offdiag before call: " );
    }
    if ( h < tileA->n ) {
        cppi_display_dbg( nextpiv, stderr, "Nextpiv offdiag before call: " );
    }

    nextpiv->h = h; /* Initialize in case it uses a copy */
    nextpiv->has_diag = chameleon_max( -1, nextpiv->has_diag);

    coreblas_kernel_trace( tileA );
    CORE_zgetrf_panel_offdiag( m, n, h, m0, ib,
                               CHAM_tile_get_ptr(tileA), tileA->ld,
                               U, ldu,
                               &(nextpiv->pivot), &(prevpiv->pivot) );

    if ( h > 0 ) {
        cppi_display_dbg( prevpiv, stderr, "Prevpiv offdiag after call: " );
    }
    if ( h < tileA->n ) {
        cppi_display_dbg( nextpiv, stderr, "Nextpiv offdiag after call: " );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgetrf_blocked_offdiag, cl_zgetrf_blocked_offdiag_cpu_func)

void INSERT_TASK_zgetrf_blocked_offdiag( const RUNTIME_option_t *options,
                                         int m, int n, int h, int m0, int ib,
                                         CHAM_desc_t *A, int Am, int An,
                                         CHAM_desc_t *U, int Um, int Un,
                                         CHAM_ipiv_t *ipiv )
{
    struct starpu_codelet *codelet = &cl_zgetrf_blocked_offdiag;

    int access_npiv = ( h == ipiv->n ) ? STARPU_R    : STARPU_REDUX;
    int access_ppiv = ( h == 0 )       ? STARPU_NONE : STARPU_R;
    int accessU     = ((h%ib == 0) && (h > 0)) ? STARPU_R : STARPU_NONE;
    int rankA       = A->get_rankof(A, Am, An);
#if !defined(HAVE_STARPU_NONE_NONZERO)
    /* STARPU_NONE can't be equal to 0 */
    fprintf( stderr, "INSERT_TASK_zgetrf_percol_diag: STARPU_NONE can not be equal to 0\n" );
    assert( 0 );
#endif

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

    void (*callback)(void*) = options->profiling ? cl_zgetrf_blocked_offdiag_callback : NULL;
    const char *cl_name = "zgetrf_blocked_offdiag";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW( A, Am, An );
    if ((h%ib == 0) && (h > 0)) {
        CHAMELEON_ACCESS_R( U, Um, Un );
    }
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 1,
                                      A->get_blktile( A, Am, An ) );

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,             &m,                   sizeof(int),
        STARPU_VALUE,             &n,                   sizeof(int),
        STARPU_VALUE,             &h,                   sizeof(int),
        STARPU_VALUE,             &An,                  sizeof(int),
        STARPU_VALUE,             &m0,                  sizeof(int),
        STARPU_VALUE,             &ib,                  sizeof(int),
        STARPU_VALUE,             &(options->sequence), sizeof(RUNTIME_sequence_t *),
        STARPU_VALUE,             &(options->request),  sizeof(RUNTIME_request_t *),
        STARPU_RW,                RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        access_npiv,              RUNTIME_pivot_getaddr( ipiv, rankA, An, h ),
        access_ppiv,              RUNTIME_pivot_getaddr( ipiv, rankA, An, h-1 ),
        accessU,                  RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un),
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_NAME,              cl_name,
        0 );
}

static const CHAMELEON_Complex64_t zone  = (CHAMELEON_Complex64_t)1.0;

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_blocked_trsm_cpu_func(void *descr[], void *cl_arg)
{
    int                    m, n, h, ib;
    CHAM_tile_t           *tileU;
    cppi_interface_t      *prevpiv;
    CHAMELEON_Complex64_t *U;
    int                    ldu;

    starpu_codelet_unpack_args( cl_arg, &m, &n, &h, &ib );

    tileU   = cti_interface_get(descr[0]);
    prevpiv = (cppi_interface_t*) descr[1];
    U       = CHAM_tile_get_ptr( tileU );
    ldu     = tileU->ld;

    coreblas_kernel_trace( tileU );

    /* Copy the final max line of the block and solve */
    cblas_zcopy( n, prevpiv->pivot.pivrow, 1,
                    U + m - 1, ldu );

    if ( (n-h) > 0 ) {
        cblas_ztrsm( CblasColMajor,
                     CblasLeft, CblasLower,
                     CblasNoTrans, CblasUnit,
                     ib, n - h,
                     CBLAS_SADDR(zone), U + (h-ib) * ldu, ldu,
                                        U +  h     * ldu, ldu );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgetrf_blocked_trsm, cl_zgetrf_blocked_trsm_cpu_func)

void INSERT_TASK_zgetrf_blocked_trsm( const RUNTIME_option_t *options,
                                      int m, int n, int h, int ib,
                                      CHAM_desc_t *U, int Um, int Un,
                                      CHAM_ipiv_t *ipiv )
{
    struct starpu_codelet *codelet = &cl_zgetrf_blocked_trsm;

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

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,             &m,                   sizeof(int),
        STARPU_VALUE,             &n,                   sizeof(int),
        STARPU_VALUE,             &h,                   sizeof(int),
        STARPU_VALUE,             &ib,                  sizeof(int),
        STARPU_RW,                RTBLKADDR(U, CHAMELEON_Complex64_t, Um, Un),
        STARPU_R,                 RUNTIME_pivot_getaddr( ipiv, rankU, Un, h-1 ),
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_NAME,              cl_name,
        0 );
}
