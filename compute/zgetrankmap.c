/**
 *
 * @file zgetrankmap.c
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrankmap wrappers
 *
 * @version 1.2.0
 * @author Abel Calluaud
 * @date 2024-03-18
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

int chameleon_zgetrankmap_op( const CHAM_desc_t *desc,
                              cham_uplo_t uplo, int m, int n,
                              CHAM_tile_t *data, void *op_args ) {
    CHAM_tile_t *tile = tile->get_blktile( desc, m, n );
    CHAM_desc_t *rkmap = (CHAM_desc_t*) op_args;
    assert( desc->format & CHAMELEON_TILE_DEC & CHAMELEON_TILE_LOWRANK );
    assert( rkmap->format & CHAMELEON_TILE_DEC );
    assert( rkmat->dtype == ChamInteger );
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 * @brief Tile equivalent of CHAMELEON_zgetrankmap().
 *
 * Operates on matrices stored by tiles.
 * All matrices are passed through descriptors.
 * All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] norm
 *          = ChamMaxNorm: Max norm
 *          = ChamOneNorm: One norm
 *          = ChamInfNorm: Infinity norm
 *          = ChamFrobeniusNorm: Frobenius norm
 *
 * @param[in] A
 *          On entry, the triangular factor U or L.
 *          On exit, if UPLO = 'U', the upper triangle of A is
 *          overwritten with the upper triangle of the product U * U';
 *          if UPLO = 'L', the lower triangle of A is overwritten with
 *          the lower triangle of the product L' * L.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrankmap
 * @sa CHAMELEON_zgetrankmap_Tile_Async
 * @sa CHAMELEON_cgetrankmap_Tile
 * @sa CHAMELEON_dgetrankmap_Tile
 * @sa CHAMELEON_sgetrankmap_Tile
 *
 */
int CHAMELEON_zgetrankmap_Tile( CHAM_desc_t *A, CHAM_desc_t *rkmap )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgetrankmap_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zgetrankmap_Tile_Async(A, &rkmap, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return CHAMELEON_SUCCESS;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 * @brief Non-blocking equivalent of CHAMELEON_zgetrankmap_Tile().
 *
 * @Warning Note that this algorithm includes a RUNTIME_Sequence_wait to cleanup
 * workspaces and thus, will return only when the result is computed. It does
 * not allow to do pipelining.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrankmap
 * @sa CHAMELEON_zgetrankmap_Tile
 * @sa CHAMELEON_cgetrankmap_Tile_Async
 * @sa CHAMELEON_dgetrankmap_Tile_Async
 * @sa CHAMELEON_sgetrankmap_Tile_Async
 *
 */
int CHAMELEON_zgetrankmap_Tile_Async( CHAM_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *rkmap,
                                 RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgetrankmap_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zgetrankmap_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zgetrankmap_Tile", "NULL request");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == CHAMELEON_SUCCESS) {
        request->status = CHAMELEON_SUCCESS;
    }
    else {
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED);
    }

    /* Check descriptors for correctness */
    if (chameleon_desc_check(A) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgetrankmap_Tile", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    if (chameleon_desc_check(rkmap) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zgetrankmap_Tile", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    chameleon_map_Tile_Async( ChamR, ChamUpperLower, A, chameleon_zgetrankmap_op, (void*)rkmap, sequence, request );

    return CHAMELEON_SUCCESS;
}
