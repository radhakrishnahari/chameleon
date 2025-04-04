/**
 *
 * @file starpu/codelet_zge2ps.c
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zge2ps StarPU codelet
 *
 * @version 1.3.0
 * @author Abel Calluaud
 * @date 2024-03-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

struct cl_zge2ps_args_s {
    int m;
    int n;
    int lda;
    int ldb;
};

#if !defined(CHAMELEON_SIMULATION)
static void cl_zge2ps_starpu_func(void *descr[], void *cl_arg)
{
    static const struct starpu_data_interface_ops *interface_ops = &starpu_interface_cham_tile_ops;
    const struct starpu_data_copy_methods         *copy_methods  = interface_ops->copy_methods;
    struct cl_zge2ps_args_s                       *clargs        = (struct cl_zge2ps_args_s *)cl_arg;

    int      workerid    = starpu_worker_get_id_check();
    unsigned memory_node = starpu_worker_get_memory_node( workerid );

    void *src_interface = descr[0];
    void *dst_interface = descr[1];

    int rc;

    assert( clargs->displA == 0 );
    assert( clargs->displB == 0 );

    rc = copy_methods->any_to_any( src_interface, memory_node,
                                   dst_interface, memory_node, NULL );
    assert( rc == 0 );
}

static void
cl_zge2ps_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zge2ps_args_s *clargs = (struct cl_zge2ps_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileR;

    tileA = cti_interface_get(descr[0]);
    tileR = cti_interface_get(descr[1]);

    int *rank = CHAM_tile_get_ptr(tileR);
    *rank = TCORE_zlrgetrk( tileA );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zge2ps,  cl_zge2ps_cpu_func  )
CODELETS( zge2ps_starpu, cl_zge2ps_starpu_func, cl_zge2ps_starpu_func, STARPU_CUDA_ASYNC )

static inline void
insert_task_zge2ps_on_local_node( const RUNTIME_option_t *options,
                                  starpu_data_handle_t handleA,
                                  starpu_data_handle_t handleB )
{
    void (*callback)(void*) = options->profiling ? cl_zge2ps_callback : NULL;
#if defined(CHAMELEON_RUNTIME_SYNC)
    starpu_data_cpy_priority( handleB, handleA, 0, callback, NULL, options->priority );
#else
    starpu_data_cpy_priority( handleB, handleA, 1, callback, NULL, options->priority );
#endif
}

#if defined(CHAMELEON_USE_MPI)
static inline void
insert_task_zge2ps_on_remote_node( const RUNTIME_option_t *options,
                                   starpu_data_handle_t handleA,
                                   starpu_data_handle_t handleB )
{
    void (*callback)(void*) = options->profiling ? cl_zge2ps_callback : NULL;
#if defined(CHAMELEON_RUNTIME_SYNC)
    starpu_mpi_data_cpy_priority( handleB, handleA, MPI_COMM_WORLD, 0, callback, NULL, options->priority );
#else
    starpu_mpi_data_cpy_priority( handleB, handleA, MPI_COMM_WORLD, 1, callback, NULL, options->priority );
#endif
}
#endif

void INSERT_TASK_zge2ps( const RUNTIME_option_t *options,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *R, int Rm, int Rn )
{
    int          exec    = 0;
    char        *cl_name = "zge2ps";
    CHAM_tile_t *tileA   = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileR   = R->get_blktile( R, Rm, Rn );

        /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(R, Rm, Rn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Insert the task */
    struct cl_zge2ps_args_s *clargs = NULL;
    void (*callback)(void*);

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zge2ps_args_s ) );
        clargs->lda    = tileA->ld;
        clargs->ldb    = tileB->ld;
    }

    /* Callback for profiling information */
    callback = options->profiling ? cl_zge2ps_callback : NULL;

    rt_starpu_insert_task(
        &cl_zge2ps,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zge2ps_args_s),
        STARPU_R,      RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_W,      RTBLKADDR(B, ChamComplexDouble, Bm, Bn),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME,              cl_name,
#endif

        0 );
}
