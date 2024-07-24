/**
 *
 * @file starpu/codelet_ipiv.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets to work with ipiv array
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @date 2024-03-16
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelets.h"

static void cl_ipiv_init_cpu_func(void *descr[], void *cl_arg)
{
    int *ipiv = (int *)STARPU_VECTOR_GET_PTR(descr[0]);

#if !defined(CHAMELEON_SIMULATION)
    {
        int i, m0, n;
        starpu_codelet_unpack_args( cl_arg, &m0, &n );

        for( i=0; i<n; i++ ) {
            ipiv[i] = m0 + i + 1;
        }
    }
#endif
}

struct starpu_codelet cl_ipiv_init = {
    .where     = STARPU_CPU,
    .cpu_func  = cl_ipiv_init_cpu_func,
    .nbuffers  = 1,
};

void INSERT_TASK_ipiv_init( const RUNTIME_option_t *options,
                            CHAM_ipiv_t            *ipiv )
{
    int64_t mt = ipiv->mt;
    int64_t mb = ipiv->mb;
    int     m;

    for (m = 0; m < mt; m++) {
        starpu_data_handle_t ipiv_src = RUNTIME_ipiv_getaddr( ipiv, m );
        int m0 = m * mb;
        int n  = (m == (mt-1)) ? ipiv->m - m0 : mb;

        rt_starpu_insert_task(
            &cl_ipiv_init,
            STARPU_VALUE, &m0, sizeof(int),
            STARPU_VALUE, &n,  sizeof(int),
            STARPU_W, ipiv_src,
            0);
    }
}

void INSERT_TASK_ipiv_reducek( const RUNTIME_option_t *options,
                               CHAM_ipiv_t *ipiv, int k, int h, int rank )
{
    starpu_data_handle_t prevpiv = RUNTIME_pivot_getaddr( ipiv, rank, k, h-1 );

#if defined(HAVE_STARPU_MPI_REDUX) && defined(CHAMELEON_USE_MPI)
#if !defined(HAVE_STARPU_MPI_REDUX_WRAPUP)
    starpu_data_handle_t nextpiv = RUNTIME_pivot_getaddr( ipiv, rank, k, h   );
    if ( h < ipiv->n ) {
        starpu_mpi_redux_data_prio_tree( options->sequence->comm, nextpiv,
                                         options->priority, 2 /* Binary tree */ );
    }
#endif
#endif

    /* Invalidate the previous pivot structure for correct initialization in later reuse */
    if ( h > 0 ) {
        starpu_data_invalidate_submit( prevpiv );
    }

    (void)options;
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_ipiv_to_perm_cpu_func( void *descr[], void *cl_arg )
{
    int m0, m, k;
    int *ipiv, *perm, *invp;

    starpu_codelet_unpack_args( cl_arg, &m0, &m, &k );

    ipiv = (int*)STARPU_VECTOR_GET_PTR(descr[0]);
    perm = (int*)STARPU_VECTOR_GET_PTR(descr[1]);
    invp = (int*)STARPU_VECTOR_GET_PTR(descr[2]);

    CORE_ipiv_to_perm( m0, m, k, ipiv, perm, invp );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
* Codelet definition
*/
static struct starpu_codelet cl_ipiv_to_perm = {
    .where        = STARPU_CPU,
#if defined(CHAMELEON_SIMULATION)
    .cpu_funcs[0] = (starpu_cpu_func_t)1,
#else
    .cpu_funcs[0] = cl_ipiv_to_perm_cpu_func,
#endif
    .nbuffers     = 3,
    .model        = NULL,
    .name         = "ipiv_to_perm"
};

void INSERT_TASK_ipiv_to_perm( const RUNTIME_option_t *options,
                               int m0, int m, int k,
                               const CHAM_ipiv_t *ipivdesc, int ipivk )
{
    struct starpu_codelet *codelet = &cl_ipiv_to_perm;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,             &m0,  sizeof(int),
        STARPU_VALUE,             &m,   sizeof(int),
        STARPU_VALUE,             &k,   sizeof(int),
        STARPU_R,                 RUNTIME_ipiv_getaddr( ipivdesc, ipivk ),
        STARPU_W,                 RUNTIME_perm_getaddr( ipivdesc, ipivk ),
        STARPU_W,                 RUNTIME_invp_getaddr( ipivdesc, ipivk ),
        STARPU_PRIORITY,          options->priority,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        0 );
}
