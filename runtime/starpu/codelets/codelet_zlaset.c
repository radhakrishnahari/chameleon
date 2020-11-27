/**
 *
 * @file starpu/codelet_zlaset.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset StarPU codelet
 *
 * @version 1.3.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @author Gwenole Lucas
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_zlaset_args_s {
    cham_uplo_t uplo;
    int m;
    int n;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t beta;
    CHAM_tile_t *tileA;
};

#if defined(CHAMELEON_USE_BUBBLE)
static inline int
cl_zlaset_is_bubble( struct starpu_task *t, void *_args )
{
    struct cl_zlaset_args_s *clargs = (struct cl_zlaset_args_s *)(t->cl_arg);
    (void)_args;

    return( clargs->tileA->format & CHAMELEON_TILE_DESC );
}

static void
cl_zlaset_bubble_func( struct starpu_task *t, void *_args )
{
    struct cl_zlaset_args_s *clargs  = (struct cl_zlaset_args_s *)(t->cl_arg);
    bubble_args_t           *b_args  = (bubble_args_t *)_args;
    RUNTIME_request_t        request = RUNTIME_REQUEST_INITIALIZER;

    /* Register the task parent */
    request.parent = t;

#if defined(CHAMELEON_BUBBLE_PARALLEL_INSERT)
    request.dependency = t;
    starpu_task_end_dep_add( t, 1 );
#endif

    chameleon_pzlaset( clargs->uplo, clargs->alpha, clargs->beta, clargs->tileA->mat,
                       b_args->sequence, &request );

    free( _args );
}
#endif /* defined(CHAMELEON_USE_BUBBLE) */

#if !defined(CHAMELEON_SIMULATION)
static void
cl_zlaset_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_zlaset_args_s *clargs = (struct cl_zlaset_args_s *)cl_arg;
    CHAM_tile_t *tileA;

    tileA = cti_interface_get(descr[0]);

    TCORE_zlaset( clargs->uplo, clargs->m, clargs->n, clargs->alpha, clargs->beta, tileA );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU( zlaset, cl_zlaset_cpu_func )

void INSERT_TASK_zlaset( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta,
                         const CHAM_desc_t *A, int Am, int An )
{
    struct cl_zlaset_args_s *clargs = NULL;
    void (*callback)(void*);
    RUNTIME_request_t       *request  = options->request;
    int                      is_bubble;
    int                      exec = 0;
    const char              *cl_name = "zlaset";
    bubble_args_t           *b_args = NULL;

    /* Handle cache */
    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_W(A, Am, An);
    exec = __chameleon_need_exec;
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs = malloc( sizeof( struct cl_zlaset_args_s ) );
        clargs->uplo  = uplo;
        clargs->m     = m;
        clargs->n     = n;
        clargs->alpha = alpha;
        clargs->beta  = beta;
        clargs->tileA = A->get_blktile( A, Am, An );
    }

    /* Callback fro profiling information */
    callback = options->profiling ? cl_zlaset_callback : NULL;

    /* Check if this is a bubble */
    is_bubble = ( clargs->tileA->format & CHAMELEON_TILE_DESC );
    if ( is_bubble ) {
        b_args = malloc( sizeof(bubble_args_t) + sizeof(struct cl_zlaset_args_s) );
        b_args->sequence = options->sequence;
        b_args->parent   = request->parent;
        memcpy( &(b_args->clargs), clargs, sizeof(struct cl_zlaset_args_s) );
        cl_name = "zlaset_bubble";
    }

    /* Insert the task */
    rt_starpu_insert_task(
        &cl_zlaset,
        /* Task codelet arguments */
        STARPU_CL_ARGS, clargs, sizeof(struct cl_zlaset_args_s),
        STARPU_W,      RTBLKADDR(A, ChamComplexDouble, Am, An),

        /* Common task arguments */
        STARPU_PRIORITY,          options->priority,
        STARPU_CALLBACK,          callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        STARPU_NAME,              cl_name,

        /* Bubble management */
#if defined(CHAMELEON_USE_BUBBLE)
        STARPU_BUBBLE_FUNC,             is_bubble_func,
        STARPU_BUBBLE_FUNC_ARG,         b_args,
        STARPU_BUBBLE_GEN_DAG_FUNC,     cl_zlaset_bubble_func,
        STARPU_BUBBLE_GEN_DAG_FUNC_ARG, b_args,

#if defined(CHAMELEON_BUBBLE_PROFILE)
        STARPU_BUBBLE_PARENT, request->parent,
#endif

#if defined(CHAMELEON_BUBBLE_PARALLEL_INSERT)
        STARPU_CALLBACK_WITH_ARG_NFREE, callback_end_dep_release, request->dependency,
#endif
#endif
        0 );

    /* Dependency is used only by the first submitted task and should not be reused */
    request->dependency = NULL;
}
