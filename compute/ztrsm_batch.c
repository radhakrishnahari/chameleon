/**
 *
 * @file ztrsm_batch.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon batch ztrsm wrappers
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

struct ztrsm_batch_args_s {
    cham_side_t           side;
    cham_uplo_t           uplo;
    cham_trans_t          transA;
    cham_diag_t           diag;
    CHAMELEON_Complex64_t alpha;
};
typedef struct ztrsm_batch_args_s ztrsm_batch_args_t;

static int
chameleon_ztrsm_batch_operator( cham_uplo_t uplo, int m, int n,
                                const CHAM_desc_t *descA, CHAM_tile_t *tileA,
                                const CHAM_desc_t *descB, CHAM_tile_t *tileB,
                                void *op_args )
{
    ztrsm_batch_args_t *args = (ztrsm_batch_args_t*)op_args;

#if !defined(CHAMELEON_SIMULATION)
    TCORE_ztrsm(
        args->side, args->uplo, args->transA, args->diag,
        descB->mb, descB->nb, args->alpha, tileA, tileB );
#endif

    (void)descA;
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
 *  CHAMELEON_ztrsm_batch_Tile - Performs multiple triangular solves in parallel.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether tiles of A appears on the left or on the right of tiles of X:
 *          = ChamLeft:  A[i,j] * X[i,j]  = B[i,j]
 *          = ChamRight: X[i,j] * A[i,j]  = B[i,j]
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = ChamUpper: Upper triangle of tiles of A are stored;
 *          = ChamLower: Lower triangle of tiles of A are stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = ChamNoTrans:   tiles of A are transposed;
 *          = ChamTrans:     tiles of A are not transposed;
 *          = ChamConjTrans: tiles of A are conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = ChamNonUnit: tiles of A are non unit;
 *          = ChamUnit:    tiles of A are unit.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a collection of mt-by-nt tiles of size A->mb by A->nb
 *
 * @param[in,out] B
 *          B is a collection of mt-by-nt tiles of size B->mb by B->nb
 *          On exit, each tile B[i,j] is overwritten by X[i,j]
 *
 *******************************************************************************
 *
 * @return CHAMELEON_SUCCESS on successful exit
 * @return CHAMELEON_ERR_... on error
 *
 */
int CHAMELEON_ztrsm_batch_Tile( cham_side_t side, cham_uplo_t uplo,
                                cham_trans_t trans, cham_diag_t diag,
                                CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    ztrsm_batch_args_t params = { side, uplo, trans, diag, alpha };
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_ztrsm_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    chameleon_pmap2( ChamUpperLower, A, B,
                     chameleon_ztrsm_batch_operator, &params,
                     sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( B, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}
