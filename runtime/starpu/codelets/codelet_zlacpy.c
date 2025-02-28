/**
 *
 * @file starpu/codelet_zlacpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy StarPU codelet
 *
 * @version 1.3.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @author Alycia Lisito
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zlacpy_args_s {
    cham_uplo_t uplo;
    int m;
    int n;
    int displA;
    int displB;
    int lda;
    int ldb;
};

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlacpy_starpu_func(void *descr[], void *cl_arg)
{
    static const struct starpu_data_interface_ops *interface_ops = &starpu_interface_cham_tile_ops;
    const struct starpu_data_copy_methods         *copy_methods  = interface_ops->copy_methods;
    struct cl_zlacpy_args_s                       *clargs        = (struct cl_zlacpy_args_s *)cl_arg;
    (void)clargs;

    int      workerid    = starpu_worker_get_id_check();
    unsigned memory_node = starpu_worker_get_memory_node( workerid );

    void *src_interface = descr[0];
    void *dst_interface = descr[1];

    int rc;
    (void)rc;

    assert( clargs->displA == 0 );
    assert( clargs->displB == 0 );

    rc = copy_methods->any_to_any( src_interface, memory_node,
                                   dst_interface, memory_node, NULL );
    assert( rc == 0 );
}

static void
cl_zlacpy_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zlacpy_args_s *clargs = (struct cl_zlacpy_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    assert( clargs->displA == 0 );
    assert( clargs->displB == 0 );

    TCORE_zlacpy( clargs->uplo, clargs->m, clargs->n, tileA, tileB );
}

static void
cl_zlacpyx_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_zlacpy_args_s *clargs = (struct cl_zlacpy_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    TCORE_zlacpyx( clargs->uplo, clargs->m, clargs->n, clargs->displA,
                   tileA, clargs->lda, clargs->displB, tileB, clargs->ldb );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zlacpy,  cl_zlacpy_cpu_func  )
CODELETS_CPU( zlacpyx, cl_zlacpyx_cpu_func )
CODELETS( zlacpy_starpu, cl_zlacpy_starpu_func, cl_zlacpy_starpu_func, STARPU_CUDA_ASYNC )

static inline void
insert_task_zlacpy_on_local_node( const RUNTIME_option_t *options,
                                  starpu_data_handle_t handleA,
                                  starpu_data_handle_t handleB )
{
    void (*callback)(void*) = NULL; // options->profiling ? cl_zlacpy_callback : NULL;
#if defined(CHAMELEON_RUNTIME_SYNC)
    starpu_data_cpy_priority( handleB, handleA, 0, callback, NULL, options->priority );
#else
    starpu_data_cpy_priority( handleB, handleA, 1, callback, NULL, options->priority );
#endif
}

#if defined(CHAMELEON_USE_MPI)
static inline void
insert_task_zlacpy_on_remote_node( const RUNTIME_option_t *options,
                                   starpu_data_handle_t handleA,
                                   starpu_data_handle_t handleB )
{
    void (*callback)(void*) = NULL; // options->profiling ? cl_zlacpy_callback : NULL;
#if defined(CHAMELEON_RUNTIME_SYNC)
    starpu_mpi_data_cpy_priority( handleB, handleA, options->sequence->comm, 0, callback, NULL, options->priority );
#else
    starpu_mpi_data_cpy_priority( handleB, handleA, options->sequence->comm, 1, callback, NULL, options->priority );
#endif
}
#endif

void INSERT_TASK_zlacpyx( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n,
                          int displA, const CHAM_desc_t *A, int Am, int An, int lda,
                          int displB, const CHAM_desc_t *B, int Bm, int Bn, int ldb )
{
    int          exec    = 0;
    char        *cl_name = "zlacpyx";
    CHAM_tile_t *tileA   = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB   = B->get_blktile( B, Bm, Bn );

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R( A, Am, An );
    CHAMELEON_ACCESS_W( B, Bm, Bn );
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

#if !defined(CHAMELEON_USE_MPI) || defined(HAVE_STARPU_MPI_DATA_CPY_PRIORITY)
    /* Insert the task */
    if ( (!options->withlacpy) &&
         (uplo == ChamUpperLower) &&
         (tileA->m == m) && (tileA->n == n) &&
         (tileB->m == m) && (tileB->n == n) &&
         (displA == 0) && (displB == 0) )
    {
#if defined(CHAMELEON_USE_MPI)
        insert_task_zlacpy_on_remote_node( options,
                                           RTBLKADDR(A, ChamComplexDouble, Am, An),
                                           RTBLKADDR(B, ChamComplexDouble, Bm, Bn) );
#else
        insert_task_zlacpy_on_local_node( options,
                                          RTBLKADDR(A, ChamComplexDouble, Am, An),
                                          RTBLKADDR(B, ChamComplexDouble, Bm, Bn) );
#endif
    }
    else
#endif
    {
        struct cl_zlacpy_args_s *clargs = NULL;
        void (*callback)(void*);

        if ( exec ) {
            clargs = malloc( sizeof( struct cl_zlacpy_args_s ) );
            clargs->uplo   = uplo;
            clargs->m      = m;
            clargs->n      = n;
            clargs->displA = displA;
            clargs->displB = displB;
            clargs->lda    = lda;
            clargs->ldb    = ldb;
        }

        /* Callback fro profiling information */
        callback = options->profiling ? cl_zlacpyx_callback : NULL;

        /* Insert the task */
        rt_starpu_insert_task(
            &cl_zlacpyx,
            /* Task codelet arguments */
            STARPU_CL_ARGS, clargs, sizeof(struct cl_zlacpy_args_s),
            STARPU_R,      RTBLKADDR(A, ChamComplexDouble, Am, An),
            STARPU_W,      RTBLKADDR(B, ChamComplexDouble, Bm, Bn),

            /* Common task arguments */
            STARPU_PRIORITY,          options->priority,
            STARPU_CALLBACK,          callback,
            STARPU_EXECUTE_ON_WORKER, options->workerid,
            STARPU_NAME,              cl_name,
            0 );
    }

    (void)tileA;
    (void)tileB;
}

void INSERT_TASK_zlacpy( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn )
{
    int          exec    = 0;
    char        *cl_name = "zlacpy";
    CHAM_tile_t *tileA   = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB   = B->get_blktile( B, Bm, Bn );

        /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_W(B, Bm, Bn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

#if !defined(CHAMELEON_USE_MPI) || defined(HAVE_STARPU_MPI_DATA_CPY_PRIORITY)
    /* Insert the task */
    if ( (!options->withlacpy) &&
         (uplo == ChamUpperLower) &&
         (tileA->m == m) && (tileA->n == n) &&
         (tileB->m == m) && (tileB->n == n) )
    {
#if defined(CHAMELEON_USE_MPI)
        insert_task_zlacpy_on_remote_node( options,
                                           RTBLKADDR(A, ChamComplexDouble, Am, An),
                                           RTBLKADDR(B, ChamComplexDouble, Bm, Bn) );
#else
        insert_task_zlacpy_on_local_node( options,
                                          RTBLKADDR(A, ChamComplexDouble, Am, An),
                                          RTBLKADDR(B, ChamComplexDouble, Bm, Bn) );
#endif
    }
    else
#endif
    {
        struct cl_zlacpy_args_s *clargs = NULL;
        void (*callback)(void*);

        if ( exec ) {
            clargs = malloc( sizeof( struct cl_zlacpy_args_s ) );
            clargs->uplo   = uplo;
            clargs->m      = m;
            clargs->n      = n;
            clargs->displA = 0;
            clargs->displB = 0;
            clargs->lda    = tileA->ld;
            clargs->ldb    = tileB->ld;
        }

        /* Callback for profiling information */
        callback = options->profiling ? cl_zlacpy_callback : NULL;

        rt_starpu_insert_task(
            &cl_zlacpy,
            /* Task codelet arguments */
            STARPU_CL_ARGS, clargs, sizeof(struct cl_zlacpy_args_s),
            STARPU_R,      RTBLKADDR(A, ChamComplexDouble, Am, An),
            STARPU_W,      RTBLKADDR(B, ChamComplexDouble, Bm, Bn),

            /* Common task arguments */
            STARPU_PRIORITY,          options->priority,
            STARPU_CALLBACK,          callback,
            STARPU_EXECUTE_ON_WORKER, options->workerid,
            STARPU_NAME,              cl_name,
            0 );
    }
}
