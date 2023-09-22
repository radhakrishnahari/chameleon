/**
 *
 * @file starpu/codelet_zlaswp.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets to apply zlaswp on a panel
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @date 2023-08-31
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlaswp_get_cpu_func( void *descr[], void *cl_arg )
{
    int          m0, k, *perm;
    CHAM_tile_t *A, *B;

    starpu_codelet_unpack_args( cl_arg, &m0, &k );

    perm = (int *)STARPU_VECTOR_GET_PTR( descr[0] );
    A    = (CHAM_tile_t *) cti_interface_get( descr[1] );
    B    = (CHAM_tile_t *) cti_interface_get( descr[2] );

    TCORE_zlaswp_get( m0, A->m, A->n, k, A, B, perm );
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU( zlaswp_get, cl_zlaswp_get_cpu_func )

void INSERT_TASK_zlaswp_get( const RUNTIME_option_t *options,
                             int m0, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *U, int Um, int Un )
{
    struct starpu_codelet *codelet = &cl_zlaswp_get;

    //void (*callback)(void*) = options->profiling ? cl_zlaswp_get_callback : NULL;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,               &m0, sizeof(int),
        STARPU_VALUE,               &k,  sizeof(int),
        STARPU_R,                   RUNTIME_perm_getaddr( ipiv, ipivk ),
        STARPU_R,                   RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_MPI_REDUX, RTBLKADDR(U, ChamComplexDouble, Um, Un),
        STARPU_PRIORITY,            options->priority,
        //STARPU_CALLBACK,            callback,
        STARPU_EXECUTE_ON_WORKER,   options->workerid,
        0 );
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zlaswp_set_cpu_func( void *descr[], void *cl_arg )
{
    int          m0, k, *invp;
    CHAM_tile_t *A, *B;

    starpu_codelet_unpack_args( cl_arg, &m0, &k );

    invp = (int *)STARPU_VECTOR_GET_PTR( descr[0] );
    A    = (CHAM_tile_t *) cti_interface_get( descr[1] );
    B    = (CHAM_tile_t *) cti_interface_get( descr[2] );

    TCORE_zlaswp_set( m0, B->m, B->n, k, A, B, invp );
}
#endif

/*
 * Codelet definition
 */
CODELETS_CPU( zlaswp_set, cl_zlaswp_set_cpu_func )

void INSERT_TASK_zlaswp_set( const RUNTIME_option_t *options,
                             int m0, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn )
{
    struct starpu_codelet *codelet = &cl_zlaswp_set;

    //void (*callback)(void*) = options->profiling ? cl_zlaswp_set_callback : NULL;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,             &m0, sizeof(int),
        STARPU_VALUE,             &k,  sizeof(int),
        STARPU_R,                 RUNTIME_invp_getaddr( ipiv, ipivk ),
        STARPU_R,                 RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_RW,                RTBLKADDR(B, ChamComplexDouble, Bm, Bn),
        STARPU_PRIORITY,          options->priority,
        //STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        0 );
}
