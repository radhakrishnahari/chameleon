/**
 *
 * @file starpu/codelet_dlag2h.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon dlag2h StarPU codelet
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2024-10-18
 * @precisions normal d -> d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_d.h"

#if !defined(CHAMELEON_SIMULATION)
#if defined(CHAMELEON_USE_CUDA)
static void
cl_dlag2h_cuda_func( void *descr[], void *cl_arg )
{
    cublasHandle_t handle = starpu_cublas_get_local_handle();
    CHAM_tile_t   *tileA;
    CHAM_tile_t   *tileB;
    int            m, n;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );

    starpu_codelet_unpack_args( cl_arg, &m, &n );

    int rc = CUDA_dlag2h(
        m, n,
        tileA->mat, tileA->ld,
        tileB->mat, tileB->ld,
        handle );

    if ( rc != 0 ) {
        fprintf( stderr, "core_dlag2h failed with info(%d)\n", rc );
    }
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS( dlag2h, NULL, cl_dlag2h_cuda_func, STARPU_CUDA_ASYNC )

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_dlag2h( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_dlag2h;
    void (*callback)(void*) = options->profiling ? cl_dlag2h_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_R,         RTBLKADDR(A, ChamRealDouble, Am, An),
        STARPU_W,         RTBLKADDR(B, ChamRealHalf, Bm, Bn),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        0);
}

#if !defined(CHAMELEON_SIMULATION)
#if defined(CHAMELEON_USE_CUDA)
static void
cl_hlag2d_cuda_func( void *descr[], void *cl_arg )
{
    cublasHandle_t handle = starpu_cublas_get_local_handle();
    CHAM_tile_t   *tileA;
    CHAM_tile_t   *tileB;
    int            m, n;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );

    starpu_codelet_unpack_args( cl_arg, &m, &n );

    int rc = CUDA_hlag2d(
        m, n,
        tileA->mat, tileA->ld,
        tileB->mat, tileB->ld,
        handle );

    if ( rc != 0 ) {
        fprintf( stderr, "core_hlag2d failed with info(%d)\n", rc );
    }
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS( hlag2d, NULL, cl_hlag2d_cuda_func, STARPU_CUDA_ASYNC )

void INSERT_TASK_hlag2d( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_hlag2d;
    void (*callback)(void*) = options->profiling ? cl_hlag2d_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R( A, Am, An );
    CHAMELEON_ACCESS_W( B, Bm, Bn );
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_R,         RTBLKADDR(A, ChamComplexFloat, Am, An),
        STARPU_W,         RTBLKADDR(B, ChamComplexDouble, Bm, Bn),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        0);
}
