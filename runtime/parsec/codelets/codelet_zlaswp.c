/**
 *
 * @file parsec/codelet_zlaswp.c
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC codelets to apply zlaswp on a panel
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zlaswp_get_parsec( parsec_execution_stream_t *context,
                        parsec_task_t             *this_task )
{
    int          m0, m, n, k, lda, ldb, *perm;
    CHAMELEON_Complex64_t *A, *B;

    parsec_dtd_unpack_args( this_task, &m0, &m, &n, &k, &A, &lda, &B, &ldb, &perm );

    CORE_zlaswp_get( m0, m, n, k, A, lda, B, ldb, perm );
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zlaswp_get( const RUNTIME_option_t *options,
                             int m0, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *U, int Um, int Un )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileU = U->get_blktile( U, Um, Un );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zlaswp_get_parsec, options->priority, "laswp_get",
        sizeof(int),          &m0,         VALUE,
        sizeof(int),          &(tileA->m), VALUE,
        sizeof(int),          &(tileA->n), VALUE,
        sizeof(int),          &k,          VALUE,
        PASSED_BY_REF, RTBLKADDR(A, ChamComplexDouble, Am, An), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),         &(tileA->ld), VALUE,
        PASSED_BY_REF, RTBLKADDR(U, ChamComplexDouble, Um, Un), chameleon_parsec_get_arena_index( U ) | INOUT,
        sizeof(int),         &(tileU->ld), VALUE,
        PASSED_BY_REF, RUNTIME_perm_getaddr( ipiv, ipivk ),     chameleon_parsec_get_arena_index_perm( ipiv ) | INPUT,
        PARSEC_DTD_ARG_END );
}

static inline int
CORE_zlaswp_set_parsec( parsec_execution_stream_t *context,
                        parsec_task_t             *this_task )
{
    int          m0, m, n, k, lda, ldb, *invp;
    CHAMELEON_Complex64_t *A, *B;

    parsec_dtd_unpack_args( this_task, &m0, &m, &n, &k, &A, &lda, &B, &ldb, &invp );

    CORE_zlaswp_set( m0, m, n, k, A, lda, B, ldb, invp );
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_zlaswp_set( const RUNTIME_option_t *options,
                             int m0, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zlaswp_set_parsec, options->priority, "laswp_set",
        sizeof(int),          &m0,         VALUE,
        sizeof(int),          &(tileB->m), VALUE,
        sizeof(int),          &(tileB->n), VALUE,
        sizeof(int),          &k,          VALUE,
        PASSED_BY_REF, RTBLKADDR(A, ChamComplexDouble, Am, An), chameleon_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),         &(tileA->ld), VALUE,
        PASSED_BY_REF, RTBLKADDR(B, ChamComplexDouble, Bm, Bn), chameleon_parsec_get_arena_index( B ) | INOUT,
        sizeof(int),         &(tileB->ld), VALUE,
        PASSED_BY_REF, RUNTIME_invp_getaddr( ipiv, ipivk ),     chameleon_parsec_get_arena_index_invp( ipiv ) | INPUT,
        PARSEC_DTD_ARG_END );
}
