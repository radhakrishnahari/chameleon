/**
 *
 * @file starpu/codelet_zherk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherk StarPU codelet
 *
 * @version 1.3.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Gwenole Lucas
 * @author Loris Lucido
 * @date 2024-10-18
 * @precisions normal z -> c
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zherk_args_s {
    cham_uplo_t uplo;
    cham_trans_t trans;
    int n;
    int k;
    double alpha;
    double beta;
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zherk_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zherk_args_s *clargs = (struct cl_zherk_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileC = cti_interface_get(descr[1]);

    assert( tileA->flttype == ChamComplexDouble );
    assert( tileC->flttype == ChamComplexDouble );

    TCORE_zherk( clargs->uplo, clargs->trans, clargs->n, clargs->k,
                 clargs->alpha, tileA, clargs->beta, tileC );
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_zherk_cuda_func(void *descr[], void *cl_arg)
{
    cublasHandle_t handle = starpu_cublas_get_local_handle();
    struct cl_zherk_args_s *clargs = (struct cl_zherk_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileC = cti_interface_get(descr[1]);

    CUDA_zherk(
        clargs->uplo, clargs->trans, clargs->n, clargs->k,
        &(clargs->alpha),
        tileA->mat, tileA->ld,
        &(clargs->beta),
        tileC->mat, tileC->ld,
        handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#if defined(CHAMELEON_USE_HIP)
static void
cl_zherk_hip_func(void *descr[], void *cl_arg)
{
    hipblasHandle_t handle = starpu_hipblas_get_local_handle();
    struct cl_zherk_args_s *clargs = (struct cl_zherk_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileC = cti_interface_get(descr[1]);

    HIP_zherk(
        clargs->uplo, clargs->trans, clargs->n, clargs->k,
        &(clargs->alpha),
        tileA->mat, tileA->ld,
        &(clargs->beta),
        tileC->mat, tileC->ld,
        handle );
}
#endif /* defined(CHAMELEON_USE_HIP) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
#if defined(CHAMELEON_USE_HIP)
CODELETS_GPU( zherk, cl_zherk_cpu_func, cl_zherk_hip_func, STARPU_HIP_ASYNC )
#else
CODELETS( zherk, cl_zherk_cpu_func, cl_zherk_cuda_func, STARPU_CUDA_ASYNC )
#endif

void INSERT_TASK_zherk( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_trans_t trans,
                        int n, int k, int nb,
                        double alpha, const CHAM_desc_t *A, int Am, int An,
                        double beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, uplo, n, n, nb,
                             beta, C, Cm, Cn );
        return;
    }

    struct cl_zherk_args_s *clargs = NULL;
    void (*callback)(void*);
    int                      accessC;
    int                      exec = 0;
    const char              *cl_name = "zherk";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zherk_args_s ) );
        clargs->uplo  = uplo;
        clargs->trans = trans;
        clargs->n     = n;
        clargs->k     = k;
        clargs->alpha = alpha;
        clargs->beta  = beta;
    }

    /* Callback fro profiling information */
    callback = options->profiling ? cl_zherk_callback : NULL;

    /* Reduce the C access if needed */
    accessC = ( beta == 0. ) ? STARPU_W : STARPU_RW;

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 2,
                                      A->get_blktile( A, Am, An ),
                                      C->get_blktile( C, Cm, Cn ) );

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zherk,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zherk_args_s),
        STARPU_R,      RTBLKADDR(A, ChamComplexDouble, Am, An),
        accessC,       RTBLKADDR(C, ChamComplexDouble, Cm, Cn),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_POSSIBLY_PARALLEL, options->parallel,
        STARPU_NAME,              cl_name,
        0 );

    (void)nb;
}
