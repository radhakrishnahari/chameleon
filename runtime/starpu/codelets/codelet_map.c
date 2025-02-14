/**
 *
 * @file starpu/codelet_map.c
 *
 * @copyright 2018-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map StarPU codelet
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2025-01-29
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

struct cl_map_args_s {
    cham_uplo_t          uplo;
    int                  m, n;
    cham_map_operator_t *op_fcts;
    void                *op_args;
    const CHAM_desc_t   *desc[1];
};

/*
 * Map with a single tile as parameter
 */
#if !defined(CHAMELEON_SIMULATION)
static void cl_map_one_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    CHAM_tile_t          *tileA;

    tileA = cti_interface_get( descr[0] );
    clargs->op_fcts->cpufunc( clargs->op_args, clargs->uplo, clargs->m, clargs->n, 1,
                              clargs->desc[0], tileA );
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_map_one_cuda_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    cublasHandle_t        handle = starpu_cublas_get_local_handle();
    CHAM_tile_t          *tileA;

    tileA = cti_interface_get( descr[0] );
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    clargs->op_fcts->cudafunc( handle, clargs->op_args, clargs->uplo, clargs->m, clargs->n, 1,
                               clargs->desc[0], tileA );
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#if defined(CHAMELEON_USE_HIP)
static void
cl_map_one_hip_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    hipblasHandle_t       handle = starpu_hipblas_get_local_handle();
    CHAM_tile_t          *tileA;

    tileA = cti_interface_get( descr[0] );
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    clargs->op_fcts->hipfunc( handle, clargs->op_args, clargs->uplo, clargs->m, clargs->n, 1,
                              clargs->desc[0], tileA );
}
#endif /* defined(CHAMELEON_USE_HIP) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CHAMELEON_CL_CB( map_one, cti_handle_get_m( task->handles[0] ), cti_handle_get_n( task->handles[0] ), 0, M * N )
#if defined(CHAMELEON_USE_HIP)
    CODELETS_GPU( map_one, cl_map_one_cpu_func, cl_map_one_hip_func, STARPU_HIP_ASYNC )
#else
CODELETS( map_one, cl_map_one_cpu_func, cl_map_one_cuda_func, STARPU_CUDA_ASYNC )
#endif

/*
 * Map with two tiles as parameter
 */
#if !defined(CHAMELEON_SIMULATION)
static void cl_map_two_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    CHAM_tile_t          *tileA;
    CHAM_tile_t          *tileB;

    tileA = cti_interface_get( descr[0] );
    tileB = cti_interface_get( descr[1] );
    clargs->op_fcts->cpufunc( clargs->op_args, clargs->uplo, clargs->m, clargs->n, 2,
                              clargs->desc[0], tileA, clargs->desc[1], tileB );
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_map_two_cuda_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    cublasHandle_t        handle = starpu_cublas_get_local_handle();
    CHAM_tile_t          *tileA;
    CHAM_tile_t          *tileB;

    tileA = cti_interface_get( descr[0] );
    tileB = cti_interface_get( descr[1] );
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );
    clargs->op_fcts->cudafunc( handle, clargs->op_args, clargs->uplo, clargs->m, clargs->n, 2,
                               clargs->desc[0], tileA, clargs->desc[1], tileB );
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#if defined(CHAMELEON_USE_HIP)
static void
cl_map_two_hip_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    hipblasHandle_t       handle = starpu_hipblas_get_local_handle();
    CHAM_tile_t          *tileA;
    CHAM_tile_t          *tileB;

    tileA = cti_interface_get( descr[0] );
    tileB = cti_interface_get( descr[1] );
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );
    clargs->op_fcts->hipfunc( handle, clargs->op_args, clargs->uplo, clargs->m, clargs->n, 2,
                              clargs->desc[0], tileA, clargs->desc[1], tileB );
}
#endif /* defined(CHAMELEON_USE_HIP) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CHAMELEON_CL_CB( map_two, cti_handle_get_m( task->handles[0] ), cti_handle_get_n( task->handles[0] ), 0, M * N )
#if defined(CHAMELEON_USE_HIP)
    CODELETS_GPU( map_two, cl_map_two_cpu_func, cl_map_two_hip_func, STARPU_HIP_ASYNC )
#else
CODELETS( map_two, cl_map_two_cpu_func, cl_map_two_cuda_func, STARPU_CUDA_ASYNC )
#endif

/*
 * Map with three tiles as parameter
 */
#if !defined(CHAMELEON_SIMULATION)
static void cl_map_three_cpu_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    CHAM_tile_t          *tileA;
    CHAM_tile_t          *tileB;
    CHAM_tile_t          *tileC;

    tileA = cti_interface_get( descr[0] );
    tileB = cti_interface_get( descr[1] );
    tileC = cti_interface_get( descr[2] );
    clargs->op_fcts->cpufunc( clargs->op_args, clargs->uplo, clargs->m, clargs->n, 3,
                              clargs->desc[0], tileA, clargs->desc[1], tileB,
                              clargs->desc[2], tileC );
}

#if defined(CHAMELEON_USE_CUDA)
static void
cl_map_three_cuda_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    cublasHandle_t        handle = starpu_cublas_get_local_handle();
    CHAM_tile_t          *tileA;
    CHAM_tile_t          *tileB;
    CHAM_tile_t          *tileC;

    tileA = cti_interface_get( descr[0] );
    tileB = cti_interface_get( descr[1] );
    tileC = cti_interface_get( descr[2] );
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );
    assert( tileC->format & CHAMELEON_TILE_FULLRANK );
    clargs->op_fcts->cudafunc( handle, clargs->op_args, clargs->uplo, clargs->m, clargs->n, 3,
                               clargs->desc[0], tileA, clargs->desc[1], tileB,
                               clargs->desc[2], tileC );
}
#endif /* defined(CHAMELEON_USE_CUDA) */

#if defined(CHAMELEON_USE_HIP)
static void
cl_map_three_hip_func( void *descr[], void *cl_arg )
{
    struct cl_map_args_s *clargs = (struct cl_map_args_s*)cl_arg;
    hipblasHandle_t       handle = starpu_hipblas_get_local_handle();
    CHAM_tile_t          *tileA;
    CHAM_tile_t          *tileB;
    CHAM_tile_t          *tileC;

    tileA = cti_interface_get( descr[0] );
    tileB = cti_interface_get( descr[1] );
    tileC = cti_interface_get( descr[2] );
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );
    assert( tileB->format & CHAMELEON_TILE_FULLRANK );
    assert( tileC->format & CHAMELEON_TILE_FULLRANK );
    clargs->op_fcts->hipfunc( handle, clargs->op_args, clargs->uplo, clargs->m, clargs->n, 3,
                              clargs->desc[0], tileA, clargs->desc[1], tileB,
                              clargs->desc[2], tileC );
}
#endif /* defined(CHAMELEON_USE_HIP) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CHAMELEON_CL_CB( map_three, cti_handle_get_m( task->handles[0] ), cti_handle_get_n( task->handles[0] ), 0, M * N )
#if defined(CHAMELEON_USE_HIP)
    CODELETS_GPU( map_three, cl_map_three_cpu_func, cl_map_three_hip_func, STARPU_HIP_ASYNC )
#else
CODELETS( map_three, cl_map_three_cpu_func, cl_map_three_cuda_func, STARPU_CUDA_ASYNC )
#endif

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, int m, int n,
                      int ndata, cham_map_data_t *data,
                      cham_map_operator_t *op_fcts, void *op_args )
{
    struct cl_map_args_s *clargs  = NULL;
    const char           *cl_name = (op_fcts->name == NULL) ? "map" : op_fcts->name;
    int                   exec    = 0;
    int                   i, readonly = 1;
    size_t                clargs_size = 0;
    uint32_t              where       = 0;
    void (*callback)(void*);

    if ( ( ndata < 0 ) || ( ndata > 3 ) ) {
        fprintf( stderr, "INSERT_TASK_map() can handle only 1 to 3 parameters\n" );
        return;
    }

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    for( i=0; i<ndata; i++ ) {
        if ( data[i].access == ChamRW ) {
            CHAMELEON_ACCESS_RW( data[i].desc, m, n );
            readonly = 0;
        }
        else if ( data[i].access == ChamW ) {
            CHAMELEON_ACCESS_W( data[i].desc, m, n );
            readonly = 0;
        }
        else {
            CHAMELEON_ACCESS_R( data[i].desc, m, n );
        }
    }
    exec = __chameleon_need_exec;
    /* Force execution for read-only functions */
    if ( readonly && __chameleon_need_submit ) {
        exec = 1;
    }
    CHAMELEON_END_ACCESS_DECLARATION;

    if ( exec ) {
        clargs_size = sizeof( struct cl_map_args_s ) + sizeof( CHAM_desc_t * ) * (ndata - 1);
        clargs = malloc( clargs_size );
        clargs->uplo    = uplo;
        clargs->m       = m;
        clargs->n       = n;
        clargs->op_fcts = op_fcts;
        clargs->op_args = op_args;
        for( i=0; i<ndata; i++ ) {
            clargs->desc[i] = data[i].desc;
        }
    }

    /* Refine name */
    for( i=0; i<ndata; i++ ) {
        cl_name = chameleon_codelet_name( cl_name, 1,
                                          (data[i].desc)->get_blktile( data[i].desc, m, n ) );
    }

    /* Where to execute */
    if ( op_fcts->cpufunc ) {
        where |= STARPU_CPU;
    }
    if ( op_fcts->cudafunc ) {
        where |= STARPU_CUDA;
    }
    if ( op_fcts->hipfunc ) {
        where |= STARPU_HIP;
    }

    /* Insert the task */
    switch( ndata ) {
    case 1:
        callback = options->profiling ? cl_map_one_callback : NULL;
        rt_starpu_insert_task(
            &cl_map_one,
            /* Task codelet arguments */
            STARPU_CL_ARGS, clargs, clargs_size,

            /* Task handles */
            cham_to_starpu_access( data[0].access ), RTBLKADDR( data[0].desc, ChamByte, m, n ),

            /* Common task arguments */
            STARPU_PRIORITY,          options->priority,
            STARPU_CALLBACK,          callback,
            STARPU_EXECUTE_ON_WORKER, options->workerid,
            STARPU_EXECUTE_WHERE,     where,
            STARPU_NAME,              cl_name,
            0 );
        break;

    case 2:
        callback = options->profiling ? cl_map_two_callback : NULL;
        rt_starpu_insert_task(
            &cl_map_two,
            /* Task codelet arguments */
            STARPU_CL_ARGS, clargs, clargs_size,

            /* Task handles */
            cham_to_starpu_access( data[0].access ), RTBLKADDR( data[0].desc, ChamByte, m, n ),
            cham_to_starpu_access( data[1].access ), RTBLKADDR( data[1].desc, ChamByte, m, n ),

            /* Common task arguments */
            STARPU_PRIORITY,          options->priority,
            STARPU_CALLBACK,          callback,
            STARPU_EXECUTE_ON_WORKER, options->workerid,
            STARPU_EXECUTE_WHERE,     where,
            STARPU_NAME,              cl_name,
            0 );
        break;

    case 3:
        callback = options->profiling ? cl_map_three_callback : NULL;
        rt_starpu_insert_task(
            &cl_map_three,
            /* Task codelet arguments */
            STARPU_CL_ARGS, clargs, clargs_size,

            /* Task handles */
            cham_to_starpu_access( data[0].access ), RTBLKADDR( data[0].desc, ChamByte, m, n ),
            cham_to_starpu_access( data[1].access ), RTBLKADDR( data[1].desc, ChamByte, m, n ),
            cham_to_starpu_access( data[2].access ), RTBLKADDR( data[2].desc, ChamByte, m, n ),

            /* Common task arguments */
            STARPU_PRIORITY,          options->priority,
            STARPU_CALLBACK,          callback,
            STARPU_EXECUTE_ON_WORKER, options->workerid,
            STARPU_EXECUTE_WHERE,     where,
            STARPU_NAME,              cl_name,
            0 );
        break;
    }
}
