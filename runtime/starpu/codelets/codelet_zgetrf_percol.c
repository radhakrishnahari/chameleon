/**
 *
 * @file starpu/codelet_zgetrf_percol.c
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
 * @date 2024-03-11
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"
#include <coreblas/cblas_wrapper.h>

CHAMELEON_CL_CB( zgetrf_percol_diag,    cti_handle_get_m(task->handles[0]), 0, 0, M );
CHAMELEON_CL_CB( zgetrf_percol_offdiag, cti_handle_get_m(task->handles[0]), 0, 0, M );

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_percol_diag_cpu_func(void *descr[], void *cl_arg)
{
    int                 m, n, h, m0;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t  *request;
    CHAM_tile_t        *tileA;
    int                *ipiv;
    cppi_interface_t   *nextpiv;
    cppi_interface_t   *prevpiv;

    starpu_codelet_unpack_args( cl_arg, &m, &n, &h, &m0, &sequence, &request );

    tileA   = cti_interface_get(descr[0]);
    ipiv    = (int *)STARPU_VECTOR_GET_PTR(descr[1]);
    nextpiv = (cppi_interface_t*) descr[2];
    prevpiv = (cppi_interface_t*) descr[3];

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

    CORE_zgetrf_panel_diag( m, n, h, m0, tileA->n,
                            CHAM_tile_get_ptr( tileA ), tileA->ld,
                            NULL, -1,
                            ipiv, &(nextpiv->pivot), &(prevpiv->pivot) );

    if ( h > 0 ) {
        cppi_display_dbg( prevpiv, stderr, "Prevpiv after call: " );
    }
    if ( h < n ) {
        cppi_display_dbg( nextpiv, stderr, "Nextpiv after call: " );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zgetrf_percol_diag, cl_zgetrf_percol_diag_cpu_func );

void INSERT_TASK_zgetrf_percol_diag( const RUNTIME_option_t *options,
                                     int m, int n, int h, int m0,
                                     CHAM_desc_t *A, int Am, int An,
                                     CHAM_ipiv_t *ipiv )
{
    struct starpu_codelet *codelet = &cl_zgetrf_percol_diag;
    void (*callback)(void*) = options->profiling ? cl_zgetrf_percol_diag_callback : NULL;
    const char *cl_name = "zgetrf_percol_diag";
    int rankA           = A->get_rankof(A, Am, An);
#if !defined ( HAVE_STARPU_NONE_NONZERO )
    /* STARPU_NONE can't be equal to 0 */
    fprintf( stderr, "INSERT_TASK_zgetrf_percol_diag: STARPU_NONE can not be equal to 0\n" );
    assert( 0 );
#endif

    if ( rankA != A->myrank ) {
        return;
    }

    int access_ipiv = ( h == 0 )       ? STARPU_W    : STARPU_RW;
    int access_npiv = ( h == ipiv->n ) ? STARPU_R    : STARPU_REDUX;
    int access_ppiv = ( h == 0 )       ? STARPU_NONE : STARPU_R;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW( A, Am, An );
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 1, A->get_blktile( A, Am, An ) );

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,             &m,                   sizeof(int),
        STARPU_VALUE,             &n,                   sizeof(int),
        STARPU_VALUE,             &h,                   sizeof(int),
        STARPU_VALUE,             &m0,                  sizeof(int),
        STARPU_VALUE,             &(options->sequence), sizeof(RUNTIME_sequence_t*),
        STARPU_VALUE,             &(options->request),  sizeof(RUNTIME_request_t*),
        STARPU_RW,                RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        access_ipiv,              RUNTIME_ipiv_getaddr( ipiv, An ),
        access_npiv,              RUNTIME_pivot_getaddr( ipiv, rankA, An, h   ),
        access_ppiv,              RUNTIME_pivot_getaddr( ipiv, rankA, An, h-1 ),
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              cl_name,
#endif
        0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_percol_offdiag_cpu_func(void *descr[], void *cl_arg)
{
    int                 m, n, h, m0;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t  *request;
    CHAM_tile_t        *tileA;
    cppi_interface_t   *nextpiv;
    cppi_interface_t   *prevpiv;

    starpu_codelet_unpack_args( cl_arg, &m, &n, &h, &m0, &sequence, &request );

    tileA   = cti_interface_get(descr[0]);
    nextpiv = (cppi_interface_t*) descr[1];
    prevpiv = (cppi_interface_t*) descr[2];

    nextpiv->h = h; /* Initialize in case it uses a copy */
    nextpiv->has_diag = chameleon_max( -1, nextpiv->has_diag);

    CORE_zgetrf_panel_offdiag( m, n, h, m0, tileA->n,
                               CHAM_tile_get_ptr(tileA), tileA->ld,
                               NULL, -1,
                               &(nextpiv->pivot), &(prevpiv->pivot) );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgetrf_percol_offdiag, cl_zgetrf_percol_offdiag_cpu_func)

void INSERT_TASK_zgetrf_percol_offdiag( const RUNTIME_option_t *options,
                                        int m, int n, int h, int m0,
                                        CHAM_desc_t *A, int Am, int An,
                                        CHAM_ipiv_t *ipiv )
{
    struct starpu_codelet *codelet = &cl_zgetrf_percol_offdiag;

    void (*callback)(void*) = options->profiling ? cl_zgetrf_percol_offdiag_callback : NULL;
    const char *cl_name = "zgetrf_percol_offdiag";
    int access_npiv = ( h == ipiv->n ) ? STARPU_R    : STARPU_REDUX;
    int access_ppiv = ( h == 0 )       ? STARPU_NONE : STARPU_R;
    int rankA       = A->get_rankof(A, Am, An);
#if !defined ( HAVE_STARPU_NONE_NONZERO )
    /* STARPU_NONE can't be equal to 0 */
    fprintf( stderr, "INSERT_TASK_zgetrf_percol_diag: STARPU_NONE can not be equal to 0\n" );
    assert( 0 );
#endif

    if ( rankA != A->myrank ) {
        return;
    }

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW( A, Am, An );
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 1, A->get_blktile( A, Am, An ) );

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,             &m,                   sizeof(int),
        STARPU_VALUE,             &n,                   sizeof(int),
        STARPU_VALUE,             &h,                   sizeof(int),
        STARPU_VALUE,             &m0,                  sizeof(int),
        STARPU_VALUE,             &(options->sequence), sizeof(RUNTIME_sequence_t *),
        STARPU_VALUE,             &(options->request),  sizeof(RUNTIME_request_t *),
        STARPU_RW,                RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        access_npiv,              RUNTIME_pivot_getaddr( ipiv, rankA, An, h   ),
        access_ppiv,              RUNTIME_pivot_getaddr( ipiv, rankA, An, h-1 ),
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              cl_name,
#endif
        0);
}
