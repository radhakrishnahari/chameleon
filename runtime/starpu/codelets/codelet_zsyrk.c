/**
 *
 * @file starpu/codelet_zsyrk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyrk StarPU codelet
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
 * @author Terry Cojean
 * @author Loris Lucido
 * @date 2023-07-06
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

struct cl_zsyrk_args_s {
    cham_uplo_t uplo;
    cham_trans_t trans;
    int n;
    int k;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t beta;
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zsyrk_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zsyrk_args_s *clargs = (struct cl_zsyrk_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileC = cti_interface_get(descr[1]);

    TCORE_zsyrk( clargs->uplo, clargs->trans, clargs->n, clargs->k,
                 clargs->alpha, tileA, clargs->beta, tileC );
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_zsyrk_cuda_func(void *descr[], void *cl_arg)
{
    cublasHandle_t handle = starpu_cublas_get_local_handle();
    struct cl_zsyrk_args_s *clargs = (struct cl_zsyrk_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileC = cti_interface_get(descr[1]);

    CUDA_zsyrk(
        clargs->uplo, clargs->trans, clargs->n, clargs->k,
        (cuDoubleComplex*)&(clargs->alpha),
        tileA->mat, tileA->ld,
        (cuDoubleComplex*)&(clargs->beta),
        tileC->mat, tileC->ld,
        handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#if defined(CHAMELEON_USE_HIP)
static void
cl_zsyrk_hip_func(void *descr[], void *cl_arg)
{
    hipblasHandle_t handle = starpu_hipblas_get_local_handle();
    struct cl_zsyrk_args_s *clargs = (struct cl_zsyrk_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileC;

    tileA = cti_interface_get(descr[0]);
    tileC = cti_interface_get(descr[1]);

    HIP_zsyrk(
        clargs->uplo, clargs->trans, clargs->n, clargs->k,
        (hipblasDoubleComplex*)&(clargs->alpha),
        tileA->mat, tileA->ld,
        (hipblasDoubleComplex*)&(clargs->beta),
        tileC->mat, tileC->ld,
        handle );
}
#endif /* defined(CHAMELEON_USE_HIP) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
#if defined(CHAMELEON_USE_HIP)
CODELETS_GPU( zsyrk, cl_zsyrk_cpu_func, cl_zsyrk_hip_func, STARPU_HIP_ASYNC )
#else
CODELETS( zsyrk, cl_zsyrk_cpu_func, cl_zsyrk_cuda_func, STARPU_CUDA_ASYNC )
#endif

/*
 * Model per code. For now the number of cores is fixed to match the machine I
 * used (with 32 cores). For each worker, we read the power consumption from a file.
 */
static cham_fixdbl_t zsyrk_energy_array [32];

cham_fixdbl_t syrk_worker_cost_function(struct starpu_task *t, unsigned workerid, unsigned i)
{
    (void)t; (void)workerid; (void)i;
    return zsyrk_energy_array[workerid];
}

static struct starpu_perfmodel energy_model =
{
    .type = STARPU_PER_WORKER,
    .worker_cost_function = syrk_worker_cost_function,
    .symbol = "zsyrk"
};

__attribute__((constructor))
static void init_cl_zsyrk() {
    /* copy energy values from env variables, only if variables are
     defined. Otherwise no energy model will be used*/
    if(getenv("STARPU_SCHED_GAMMA"))
    {
        FILE *fp;
        char* line = NULL;
        size_t len = 0;
        ssize_t read;
        int i = 0;
        fp = fopen("/root/syrk_energy_file", "r");
        if(fp == NULL){
            printf("failed to open file syrk_energy_file");
            exit(EXIT_FAILURE);
        }
        while((read = getline(&line, &len, fp)) != -1) {
            zsyrk_energy_array[i] = atof(line);
            i ++;
        }
        fclose(fp);
        if(line)
            free(line);
        cl_zsyrk.energy_model = &energy_model ;
    }
}

void INSERT_TASK_zsyrk( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_trans_t trans,
                        int n, int k, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        CHAMELEON_Complex64_t beta,  const CHAM_desc_t *C, int Cm, int Cn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, uplo, n, n, nb,
                             beta, C, Cm, Cn );
        return;
    }

    struct cl_zsyrk_args_s *clargs = NULL;
    void (*callback)(void*);
    int                      accessC;
    int                      exec = 0;
    char                    *cl_name = "zsyrk";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(C, Cm, Cn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zsyrk_args_s ) );
        clargs->uplo  = uplo;
        clargs->trans = trans;
        clargs->n     = n;
        clargs->k     = k;
        clargs->alpha = alpha;
        clargs->beta  = beta;
    }

    /* Callback fro profiling information */
    callback = options->profiling ? cl_zsyrk_callback : NULL;

    /* Reduce the C access if needed */
    accessC = ( beta == 0. ) ? STARPU_W : STARPU_RW;

    /* Refine name */
    cl_name = chameleon_codelet_name( cl_name, 2,
                                      A->get_blktile( A, Am, An ),
                                      C->get_blktile( C, Cm, Cn ) );

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zsyrk,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zsyrk_args_s),
        STARPU_R,      RTBLKADDR(A, ChamComplexDouble, Am, An),
        accessC,       RTBLKADDR(C, ChamComplexDouble, Cm, Cn),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_POSSIBLY_PARALLEL, options->parallel,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              cl_name,
#endif

        0 );

    (void)nb;
}
