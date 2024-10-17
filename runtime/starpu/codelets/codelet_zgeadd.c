/**
 *
 * @file starpu/codelet_zgeadd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgeadd StarPU codelet
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zgeadd_cpu_func( void *descr[], void *cl_arg )
{
    cham_trans_t          trans;
    int                   M;
    int                   N;
    CHAMELEON_Complex64_t alpha;
    CHAM_tile_t          *tileA;
    CHAMELEON_Complex64_t beta;
    CHAM_tile_t          *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    starpu_codelet_unpack_args( cl_arg, &trans, &M, &N, &alpha, &beta );
    TCORE_zgeadd( trans, M, N, alpha, tileA, beta, tileB );

    return;
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_zgeadd_cuda_func( void *descr[], void *cl_arg )
{
    cublasHandle_t  handle = starpu_cublas_get_local_handle();
    cham_trans_t    trans;
    int             M;
    int             N;
    cuDoubleComplex alpha;
    CHAM_tile_t    *tileA;
    cuDoubleComplex beta;
    CHAM_tile_t    *tileB;

    tileA = cti_interface_get( descr[0] );
    tileB = cti_interface_get( descr[1] );
    starpu_codelet_unpack_args( cl_arg, &trans, &M, &N, &alpha, &beta );

    CUDA_zgeadd( trans, M, N,
                 &alpha, tileA->mat, tileA->ld,
                 &beta,  tileB->mat, tileB->ld,
                 handle );

    return;
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS( zgeadd, cl_zgeadd_cpu_func, cl_zgeadd_cuda_func, STARPU_CUDA_ASYNC );

void INSERT_TASK_zgeadd( const RUNTIME_option_t *options,
                         cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, ChamUpperLower, m, n, nb,
                             beta, B, Bm, Bn );
        return;
    }

    struct starpu_codelet *codelet = &cl_zgeadd;
    void (*callback)(void*) = options->profiling ? cl_zgeadd_callback : NULL;
    int accessB = ( beta == 0. ) ? STARPU_W : STARPU_RW;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,    &trans,              sizeof(int),
        STARPU_VALUE,    &m,                  sizeof(int),
        STARPU_VALUE,    &n,                  sizeof(int),
        STARPU_VALUE,    &alpha,              sizeof(CHAMELEON_Complex64_t),
        STARPU_R,         RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_VALUE,    &beta,               sizeof(CHAMELEON_Complex64_t),
        accessB,          RTBLKADDR(B, ChamComplexDouble, Bm, Bn),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        0 );

    (void)nb;
}
