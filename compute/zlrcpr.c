/**
 *
 * @file zlrcpr.c
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon low rank matrix compression ratio calculation.
 *
 * @version 1.2.0
 * @author Abel Calluaud
 * @date 2024-05-02
 * @precisions normal z -> s d c
 *
 */
#include "coreblas.h"
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zlrcpr_Tile compute the low rank matrix 
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          The M-by-N matrix A. If uplo = ChamUpper, only the upper trapezium
 *          is accessed; if UPLO = ChamLower, only the lower trapezium is
 *          accessed.
 *
 * @param[out] B
 *          The M-by-N low rank matrix B.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlrcpr
 * @sa CHAMELEON_zlrcpr_Tile_Async
 * @sa CHAMELEON_clrcpr_Tile
 * @sa CHAMELEON_dlrcpr_Tile
 * @sa CHAMELEON_slrcpr_Tile
 *
 */
#include "chameleon.h"
#include "runtime_rpk.h"
#include "coreblas.h"

int zlrcpr_op_cpu( void *args, cham_uplo_t uplo, int m, int n, int ndata,
                   const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... ) {
    rpk_matrix_t *Ara;
    double tile_footprint;

    assert( tileA->format == CHAMELEON_TILE_LOWRANK );

    Ara = (rpk_matrix_t *) CHAM_tile_get_ptr( tileA );
    if (Ara->rk == -1) {
    	tile_footprint = tileA->m * tileA->n;
    }
    else {
	tile_footprint = Ara->rkmax * ((double) tileA->m + tileA->n);
    }

    return CHAMELEON_SUCCESS;
}

int CHAMELEON_zlrcpr_Tile( cham_uplo_t uplo, CHAM_desc_t *descAt ) {
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status = CHAMELEON_SUCCESS;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlrcpr_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    cham_map_operator_t zlrcpr_op;
    zlrcpr_op.name = "zlrcpr";
    zlrcpr_op.cpufunc = zlrcpr_op_cpu;
    zlrcpr_op.cudafunc = NULL;
    zlrcpr_op.hipfunc = NULL;

    cham_map_data_t zlrcpr_data[2];
    zlrcpr_data[0].access = ChamR;
    zlrcpr_data[0].desc = descAt;

    CHAMELEON_mapv_Tile_Async( uplo, 1, zlrcpr_data, &zlrcpr_op, NULL, sequence, &request );

    CHAMELEON_Desc_Flush( descAt, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}
