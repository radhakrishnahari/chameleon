/**
 *
 * @file zgemm_batch.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon batch zgemm wrappers
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2020-03-11
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#if !defined(CHAMELEON_SIMULATION)
#include "coreblas/coreblas_ztile.h"
#endif

struct zgemm_batch_args_s {
    cham_trans_t          transA;
    cham_trans_t          transB;
    CHAMELEON_Complex64_t alpha;
    CHAMELEON_Complex64_t beta;
};
typedef struct zgemm_batch_args_s zgemm_batch_args_t;

static int
chameleon_zgemm_batch_operator( cham_uplo_t uplo, int m, int n,
                                const CHAM_desc_t *descA, CHAM_tile_t *tileA,
                                const CHAM_desc_t *descB, CHAM_tile_t *tileB,
                                const CHAM_desc_t *descC, CHAM_tile_t *tileC,
                                void *op_args )
{
    int tempkk;
    zgemm_batch_args_t *args = (zgemm_batch_args_t*)op_args;

    tempkk = ( args->transA == ChamNoTrans ) ? descA->nb : descA->mb;

#if !defined(CHAMELEON_SIMULATION)
    TCORE_zgemm(
        args->transA, args->transB, descC->mb, descC->nb, tempkk,
        args->alpha, tileA, tileB, args->beta, tileC );
#endif

    (void)descB;
    (void)uplo;
    (void)m;
    (void)n;
    return 0;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zgemm_batch_Tile - Performs multiple matrix multiplication in parallel.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Specifies whether the tiles from A are transposed, not transposed or conjugate transposed:
 *          = ChamNoTrans:   tiles from A are not transposed;
 *          = ChamTrans:     tiles from A are transposed;
 *          = ChamConjTrans: tiles from A are conjugate transposed.
 *
 * @param[in] transB
 *          Specifies whether the tiles from B are transposed, not transposed or conjugate transposed:
 *          = ChamNoTrans:   tiles from B are not transposed;
 *          = ChamTrans:     tiles from B are transposed;
 *          = ChamConjTrans: tiles from B are conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a collection of mt-by-nt tiles of size A->mb by A->nb
 *
 * @param[in] B
 *          B is a collection of mt-by-nt tiles of size B->mb by B->nb
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a collection of mt-by-nt tiles of size C->mb by C->nb
 *          On exit, each tile Cij is overwritten by the matrix:
 *          \f[ alpha * op( A[i,j] )*op( B[i,j] ) * C[i,j] \f]
 *
 *******************************************************************************
 *
 * @return CHAMELEON_SUCCESS on successful exit
 * @return CHAMELEON_ERR_... on error
 *
 */
int CHAMELEON_zgemm_batch_Tile( cham_trans_t transA, cham_trans_t transB,
                                CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                                CHAMELEON_Complex64_t beta,  CHAM_desc_t *C )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    zgemm_batch_args_t params = { transA, transB, alpha, beta };
    int status;
    int accessC;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgemm_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    /* Reduce the C access if needed */
    accessC = ( beta == 0. ) ? ChamW : ChamRW;

    chameleon_pmap3( ChamR, ChamR, accessC,
                     ChamUpperLower, A, B, C,
                     chameleon_zgemm_batch_operator, &params,
                     sequence, &request, "zgemm" );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( B, sequence );
    CHAMELEON_Desc_Flush( C, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}
