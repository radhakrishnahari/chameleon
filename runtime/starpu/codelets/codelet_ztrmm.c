/**
 *
 * @file starpu/codelet_ztrmm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrmm StarPU codelet
 *
 * @version 1.3.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Loris Lucido
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_ztrmm_args_s {
    cham_side_t side;
    cham_uplo_t uplo;
    cham_trans_t transA;
    cham_diag_t diag;
    int m;
    int n;
    CHAMELEON_Complex64_t alpha;
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_ztrmm_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_ztrmm_args_s *clargs = (struct cl_ztrmm_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    TCORE_ztrmm( clargs->side, clargs->uplo, clargs->transA, clargs->diag,
                 clargs->m, clargs->n, clargs->alpha, tileA, tileB );
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_ztrmm_cuda_func(void *descr[], void *cl_arg)
{
    cublasHandle_t handle = starpu_cublas_get_local_handle();
    struct cl_ztrmm_args_s *clargs = (struct cl_ztrmm_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    CUDA_ztrmm(
        clargs->side, clargs->uplo, clargs->transA, clargs->diag,
        clargs->m, clargs->n,
        (cuDoubleComplex*)&(clargs->alpha),
        tileA->mat, tileA->ld,
        tileB->mat, tileB->ld,
        handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#if defined(CHAMELEON_USE_HIP)
static void
cl_ztrmm_hip_func(void *descr[], void *cl_arg)
{
    struct cl_ztrmm_args_s *clargs = (struct cl_ztrmm_args_s *)cl_arg;
    hipblasHandle_t         handle = starpu_hipblas_get_local_handle();
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    HIP_ztrmm(
        clargs->side, clargs->uplo, clargs->transA, clargs->diag,
        clargs->m, clargs->n,
        (hipblasDoubleComplex*)&(clargs->alpha),
        tileA->mat, tileA->ld,
        tileB->mat, tileB->ld,
        handle );
}
#endif /* defined(CHAMELEON_USE_HIP) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
#if defined(CHAMELEON_USE_HIP)
CODELETS_GPU( ztrmm, cl_ztrmm_cpu_func, cl_ztrmm_hip_func, STARPU_HIP_ASYNC )
#else
CODELETS( ztrmm, cl_ztrmm_cpu_func, cl_ztrmm_cuda_func, STARPU_CUDA_ASYNC )
#endif

void INSERT_TASK_ztrmm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn )
{
    struct cl_ztrmm_args_s *clargs = NULL;
    void (*callback)(void*);
    int                      exec = 0;
    const char              *cl_name = "ztrmm";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(B, Bm, Bn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_ztrmm_args_s ) );
        clargs->side   = side;
        clargs->uplo   = uplo;
        clargs->transA = transA;
        clargs->diag   = diag;
        clargs->m      = m;
        clargs->n      = n;
        clargs->alpha  = alpha;
    }

    /* Callback fro profiling information */
    callback = options->profiling ? cl_ztrmm_callback : NULL;

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 2,
                                      A->get_blktile( A, Am, An ),
                                      B->get_blktile( B, Bm, Bn ) );

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_ztrmm,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_ztrmm_args_s),
        STARPU_R,      RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_RW,     RTBLKADDR(B, ChamComplexDouble, Bm, Bn),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_NAME,              cl_name,
        0 );

    (void)nb;
}
