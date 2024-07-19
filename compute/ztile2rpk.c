/**
 *
 * @file ztile2rpk.c
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon tile to rapack converion routine
 *
 * @version 1.2.0
 * @author Abel Calluaud
 * @date 2024-02-19
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include "coreblas.h"
#include "runtime_lrmem.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_ztile2rpk_Tile constructs a low rank matrix B from a full rank
 *  matrix A.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          The M-by-N full rank matrix A. If uplo = ChamUpper, only the upper trapezium
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
 * @sa CHAMELEON_ztile2rpk
 * @sa CHAMELEON_ztile2rpk_Tile_Async
 * @sa CHAMELEON_ctile2rpk_Tile
 * @sa CHAMELEON_dtile2rpk_Tile
 * @sa CHAMELEON_stile2rpk_Tile
 *
 */
#include "chameleon.h"
#include "runtime_rpk.h"
#include "coreblas.h"
#include "rapack.h"

static int ztile2rpk_op_cpu( void *args, cham_uplo_t uplo, int m, int n, int ndata,
                   const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... ) {
    /*
	 * Update the tile data structure to rapack format
	 */
	runtime_ztile2rpk( tileA );

    return CHAMELEON_SUCCESS;
}

int CHAMELEON_ztile2rpk_Tile( cham_uplo_t uplo, CHAM_desc_t *descAt ) {
    int status = CHAMELEON_SUCCESS;

    int i;
    int j;

	cham_map_operator_t ztile2rpk_op = {0};
    ztile2rpk_op.name = "ztile2rpk";
    ztile2rpk_op.cpufunc = ztile2rpk_op_cpu;

    cham_map_data_t ztile2rpk_data;
    ztile2rpk_data.access = ChamRW;
    ztile2rpk_data.desc = descAt;

    /*
     * Convert each tile inplace to rapack format.
     * The original approximate rank structure is kept.
     */
	CHAMELEON_mapv_Tile(
		uplo, 1, &ztile2rpk_data, &ztile2rpk_op, NULL
	);

	for (i = 0; i < descAt->mt; i++) {
		for (j = 0; j < descAt->nt; j++) {
			CHAM_tile_t *tile = descAt->get_blktile( descAt, i, j );
	        assert( tile->format != CHAMELEON_TILE_LOWRANK );
            assert( tile->mat == NULL );
            tile->format = CHAMELEON_TILE_LOWRANK;
		}
	}

    descAt->format = ChamTLR;

	return status;
}
