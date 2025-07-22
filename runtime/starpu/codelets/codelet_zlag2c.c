/**
 *
 * @file starpu/codelet_zlag2c.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlag2c StarPU codelet
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2025-07-22
 * @precisions mixed zc -> ds
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_zc.h"

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zlag2c_cpu_func( void *descr[], void *cl_arg )
{
    int info = 0;
    int m;
    int n;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    starpu_codelet_unpack_args( cl_arg, &m, &n );
    TCORE_zlag2c( m, n, tileA, tileB, &info );
}

#if defined(CHAMELEON_USE_CUDA) && defined(GPUCUBLAS_HAVE_CUDA_TOOLKIT)
static void
cl_zlag2c_cuda_func( void *descr[], void *cl_arg )
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

    int rc = CUDA_zlag2c(
        m, n,
        tileA->mat, tileA->ld,
        tileB->mat, tileB->ld,
        handle );

    if ( rc != 0 ) {
        fprintf( stderr, "core_zlag2c failed with info(%d)\n", rc );
    }
}
#endif /* defined(CHAMELEON_USE_CUDA) && defined(GPUCUBLAS_HAVE_CUDA_TOOLKIT) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
#if defined(CHAMELEON_USE_CUDA) && defined(GPUCUBLAS_HAVE_CUDA_TOOLKIT)
CODELETS( zlag2c, cl_zlag2c_cpu_func, cl_zlag2c_cuda_func, STARPU_CUDA_ASYNC )
#else
CODELETS_CPU( zlag2c, cl_zlag2c_cpu_func )
#endif

/**
 *
 * @ingroup INSERT_TASK_Complex64_t
 *
 */
void INSERT_TASK_zlag2c(const RUNTIME_option_t *options,
                       int m, int n, int nb,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zlag2c;
    void (*callback)(void*) = options->profiling ? cl_zlag2c_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,    &m,                 sizeof(int),
        STARPU_VALUE,    &n,                 sizeof(int),
        STARPU_R,         RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_W,         RTBLKADDR(B, ChamComplexFloat, Bm, Bn),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_clag2z_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &m, &n);
    TCORE_clag2z( m, n, tileA, tileB);
}

#if defined(CHAMELEON_USE_CUDA) && defined(GPUCUBLAS_HAVE_CUDA_TOOLKIT)
static void
cl_clag2z_cuda_func( void *descr[], void *cl_arg )
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

    int rc = CUDA_clag2z(
        m, n,
        tileA->mat, tileA->ld,
        tileB->mat, tileB->ld,
        handle );
    if ( rc != 0 ) {
        fprintf( stderr, "core_clag2z failed with info(%d)\n", rc );
    }
}
#endif /* defined(CHAMELEON_USE_CUDA) && defined(GPUCUBLAS_HAVE_CUDA_TOOLKIT) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
#if defined(CHAMELEON_USE_CUDA) && defined(GPUCUBLAS_HAVE_CUDA_TOOLKIT)
CODELETS( clag2z, cl_clag2z_cpu_func, cl_clag2z_cuda_func, STARPU_CUDA_ASYNC )
#else
CODELETS_CPU( clag2z, cl_clag2z_cpu_func )
#endif

void INSERT_TASK_clag2z( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_clag2z;
    void (*callback)(void*) = options->profiling ? cl_clag2z_callback : NULL;

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
