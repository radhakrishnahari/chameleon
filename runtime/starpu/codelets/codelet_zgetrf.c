/**
 *
 * @file starpu/codelet_zgetrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf StarPU codelet
 *
 * @version 1.3.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2024-10-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    CHAM_tile_t *tileA;
    int *IPIV;
    cham_bool_t check_info;
    int iinfo;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    int info = 0;

    tileA = cti_interface_get(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &m, &n, &IPIV, &check_info, &iinfo, &sequence, &request);
    TCORE_zgetrf( m, n, tileA, IPIV, &info );

    if ( (sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, sequence, request, iinfo+info );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgetrf, cl_zgetrf_cpu_func)

void INSERT_TASK_zgetrf( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo )
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgetrf;
    void (*callback)(void*) = options->profiling ? cl_zgetrf_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,             &m,                        sizeof(int),
        STARPU_VALUE,             &n,                        sizeof(int),
        STARPU_RW,                     RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_VALUE,                  &IPIV,                      sizeof(int*),
        STARPU_VALUE,    &check_info,                sizeof(cham_bool_t),
        STARPU_VALUE,         &iinfo,                        sizeof(int),
        STARPU_VALUE,    &(options->sequence),       sizeof(RUNTIME_sequence_t*),
        STARPU_VALUE,    &(options->request),        sizeof(RUNTIME_request_t*),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        0 );
}
