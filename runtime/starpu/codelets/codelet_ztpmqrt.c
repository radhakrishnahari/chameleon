/**
 *
 * @file starpu/codelet_ztpmqrt.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon ztpmqrt StarPU codelet
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Lucas Barros de Assis
 * @author Florent Pruvost
 * @date 2023-07-06
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"

#if !defined(CHAMELEON_SIMULATION)
static void cl_ztpmqrt_cpu_func(void *descr[], void *cl_arg)
{
    cham_side_t side;
    cham_trans_t trans;
    int M;
    int N;
    int K;
    int L;
    int ib;
    size_t lwork;
    CHAM_tile_t *tileV;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileW;

    tileV = cti_interface_get(descr[0]);
    tileT = cti_interface_get(descr[1]);
    tileA = cti_interface_get(descr[2]);
    tileB = cti_interface_get(descr[3]);
    tileW = cti_interface_get(descr[4]); /* ib * nb */
    starpu_codelet_unpack_args( cl_arg, &side, &trans, &M, &N, &K, &L, &ib, &lwork );

    TCORE_ztpmqrt( side, trans, M, N, K, L, ib,
                   tileV, tileT, tileA, tileB, tileW->mat );
}

#if defined(CHAMELEON_USE_CUDA)
static void cl_ztpmqrt_cuda_func(void *descr[], void *cl_arg)
{
    cublasHandle_t handle = starpu_cublas_get_local_handle();
    cham_side_t side;
    cham_trans_t trans;
    int M;
    int N;
    int K;
    int L;
    int ib;
    size_t lwork;
    CHAM_tile_t *tileV;
    CHAM_tile_t *tileT;
    CHAM_tile_t *tileA;
    CHAM_tile_t *tileB;
    CHAM_tile_t *tileW;

    tileV = cti_interface_get(descr[0]);
    tileT = cti_interface_get(descr[1]);
    tileA = cti_interface_get(descr[2]);
    tileB = cti_interface_get(descr[3]);
    tileW = cti_interface_get(descr[4]); /* 3*ib*nb */

    starpu_codelet_unpack_args( cl_arg, &side, &trans, &M, &N, &K, &L, &ib, &lwork );

    CUDA_ztpmqrt(
            side, trans, M, N, K, L, ib,
            tileV->mat, tileV->ld,
            tileT->mat, tileT->ld,
            tileA->mat, tileA->ld,
            tileB->mat, tileB->ld,
            tileW->mat, lwork, handle );
}
#endif /* defined(CHAMELEON_USE_CUDA) */
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS(ztpmqrt, cl_ztpmqrt_cpu_func, cl_ztpmqrt_cuda_func, STARPU_CUDA_ASYNC)

void INSERT_TASK_ztpmqrt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int M, int N, int K, int L, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn )
{
    struct starpu_codelet *codelet = &cl_ztpmqrt;
    void (*callback)(void*) = options->profiling ? cl_ztpmqrt_callback : NULL;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_R(V, Vm, Vn);
    CHAMELEON_ACCESS_R(T, Tm, Tn);
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_ACCESS_RW(B, Bm, Bn);
    CHAMELEON_END_ACCESS_DECLARATION;

    rt_starpu_insert_task(
        codelet,
        STARPU_VALUE, &side,  sizeof(int),
        STARPU_VALUE, &trans, sizeof(int),
        STARPU_VALUE, &M,     sizeof(int),
        STARPU_VALUE, &N,     sizeof(int),
        STARPU_VALUE, &K,     sizeof(int),
        STARPU_VALUE, &L,     sizeof(int),
        STARPU_VALUE, &ib,     sizeof(int),
        STARPU_VALUE, &(options->ws_wsize), sizeof(size_t),
        STARPU_R,      RTBLKADDR(V, ChamComplexDouble, Vm, Vn),
        STARPU_R,      RTBLKADDR(T, ChamComplexDouble, Tm, Tn),
        STARPU_RW,     RTBLKADDR(A, ChamComplexDouble, Am, An),
        STARPU_RW,     RTBLKADDR(B, ChamComplexDouble, Bm, Bn),
        /* Other options */
        STARPU_SCRATCH,   options->ws_worker,
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
        STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_USE_MPI)
        STARPU_EXECUTE_ON_NODE, B->get_rankof(B, Bm, Bn),
#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, (( L == 0 ) ? "ztsmqr" : "ztpmqrt"),
#endif
        0);

    (void)ib; (void)nb;
}
