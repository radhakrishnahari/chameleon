/**
 *
 * @file zherk_batch.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon batch zherk wrappers
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

struct zherk_batch_args_s {
    cham_uplo_t  uplo;
    cham_trans_t trans;
    double       alpha;
    double       beta;
};
typedef struct zherk_batch_args_s zherk_batch_args_t;

static int
chameleon_zherk_batch_operator( cham_uplo_t uplo, int m, int n,
                                const CHAM_desc_t *descA, CHAM_tile_t *tileA,
                                const CHAM_desc_t *descB, CHAM_tile_t *tileB,
                                void *op_args )
{
    int tempkk;
    zherk_batch_args_t *args = (zherk_batch_args_t*)op_args;

    tempkk = ( args->trans == ChamNoTrans ) ? descA->nb : descA->mb;

#if !defined(CHAMELEON_SIMULATION)
    TCORE_zherk(
        args->uplo, args->trans, descB->nb, tempkk,
        args->alpha, tileA, args->beta, tileB );
#endif

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
 *  CHAMELEON_zherk_batch_Tile - Performs multiple rank-k update in parallel.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = ChamUpper: Upper triangle of C tiles are stored;
 *          = ChamLower: Lower triangle of C tiles are stored.
 *
 * @param[in] trans
 *          Specifies whether the tiles of A are transposed or conjugate transposed:
 *          = ChamNoTrans:   tiles of A are not transposed;
 *          = ChamConjTrans: tiles of A are conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a collection of mt-by-nt tiles of size A->mb by A->nb
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] B
 *          B is a collection of mt-by-nt tiles of size B->mb by B->nb
 *          On exit, each tile B[i,j] is overwritten by
 *          \f[ B = \alpha [ op( A ) \times conjg( op( A )' )] + \beta B \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = conjg( X' )
 *
 *******************************************************************************
 *
 * @return CHAMELEON_SUCCESS on successful exit
 * @return CHAMELEON_ERR_... on error
 *
 */
int CHAMELEON_zherk_batch_Tile( cham_uplo_t uplo, cham_trans_t trans,
                                double alpha, CHAM_desc_t *A,
                                double beta,  CHAM_desc_t *B )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    zherk_batch_args_t params = { uplo, trans, alpha, beta };
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zherk_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    chameleon_pmap2( ChamUpperLower, A, B,
                     chameleon_zherk_batch_operator, &params,
                     sequence, &request, "zherk" );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( B, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}
