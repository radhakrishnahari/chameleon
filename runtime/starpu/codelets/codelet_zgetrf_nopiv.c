/**
 *
 * @file starpu/codelet_zgetrf_nopiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_nopiv StarPU codelet
 *
 * @version 1.3.0
 * @author Omar Zenati
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

/*
 * Codelet CPU
 */
#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_nopiv_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    int ib;
    CHAM_tile_t *tileA;
    int iinfo;
    RUNTIME_sequence_t *sequence;
    RUNTIME_request_t *request;
    int info = 0;

    tileA = cti_interface_get(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &m, &n, &ib, &iinfo, &sequence, &request);
    TCORE_zgetrf_nopiv(m, n, ib, tileA, &info);

    if ( (sequence->status == CHAMELEON_SUCCESS) && (info != 0) ) {
        RUNTIME_sequence_flush( NULL, sequence, request, iinfo+info );
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgetrf_nopiv, cl_zgetrf_nopiv_cpu_func)

void INSERT_TASK_zgetrf_nopiv(const RUNTIME_option_t *options,
                              int m, int n, int ib, int nb,
                              const CHAM_desc_t *A, int Am, int An,
                              int iinfo)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgetrf_nopiv;
    void (*callback)(void*) = options->profiling ? cl_zgetrf_nopiv_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE,    &m,                         sizeof(int),
        STARPU_VALUE,    &n,                         sizeof(int),
        STARPU_VALUE,    &ib,                        sizeof(int),
        STARPU_RW,        RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_VALUE,    &iinfo,                     sizeof(int),
        STARPU_VALUE,    &(options->sequence),       sizeof(RUNTIME_sequence_t*),
        STARPU_VALUE,    &(options->request),        sizeof(RUNTIME_request_t*),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
        0 );
}
