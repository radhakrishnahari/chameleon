/**
 *
 * @file starpu/codelet_ztradd.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztradd StarPU codelet
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_ztradd_args_s {
    cham_uplo_t uplo;
    cham_trans_t trans;
    int m;
    int n;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t beta;
};

#if !defined(CHAMELEON_SIMULATION)
static void
cl_ztradd_cpu_func(void *descr[], void *cl_arg)
{
    struct cl_ztradd_args_s *clargs = (struct cl_ztradd_args_s *)cl_arg;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;

    tileA = cti_interface_get(descr[0]);
    tileB = cti_interface_get(descr[1]);

    TCORE_ztradd( clargs->uplo, clargs->trans, clargs->m, clargs->n,
                  clargs->alpha, tileA, clargs->beta, tileB );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( ztradd, cl_ztradd_cpu_func )

void INSERT_TASK_ztradd( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn )
{
    if ( alpha == 0. ) {
        INSERT_TASK_zlascal( options, uplo, m, n, nb,
                             beta, B, Bm, Bn );
        return;
    }

    struct cl_ztradd_args_s *clargs = NULL;
    void (*callback)(void*);
    int                      accessB;
    int                      exec = 0;
    const char              *cl_name = "ztradd";

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(A, Am, An);
    CHAMELEON_ACCESS_RW(B, Bm, Bn);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_ztradd_args_s ) );
        clargs->uplo  = uplo;
        clargs->trans = trans;
        clargs->m     = m;
        clargs->n     = n;
        clargs->alpha = alpha;
        clargs->beta  = beta;
    }

    /* Callback fro profiling information */
    callback = options->profiling ? cl_ztradd_callback : NULL;

    /* Reduce the B access if needed */
    accessB = ( beta == 0. ) ? STARPU_W : STARPU_RW;

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_ztradd,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_ztradd_args_s),
        STARPU_R,      RTBLKADDR(A, ChamComplexDouble, Am, An),
        accessB,       RTBLKADDR(B, ChamComplexDouble, Bm, Bn),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_NAME,              cl_name,
        0 );

    (void)nb;
}
