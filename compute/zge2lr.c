/**
 *
 * @file zge2lr.c
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon low rank assembly routines.
 *
 * @version 1.2.0
 * @author Abel Calluaud
 * @date 2024-02-05
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
 *  CHAMELEON_zge2lr_Tile constructs a low rank matrix B from a full rank
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
 * @sa CHAMELEON_zge2lr
 * @sa CHAMELEON_zge2lr_Tile_Async
 * @sa CHAMELEON_cge2lr_Tile
 * @sa CHAMELEON_dge2lr_Tile
 * @sa CHAMELEON_sge2lr_Tile
 *
 */
#include "chameleon.h"
#include "runtime_rpk.h"
#include "coreblas.h"

int zge2lr_op_cpu( void *args, cham_uplo_t uplo, int m, int n, int ndata,
                   const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... ) {
	const rpk_ctx_t *ctx = runtime_rpk_zctx_get();
    int compressible = (m != n);

    rpk_matrix_t *Ara = (rpk_matrix_t*)CHAM_tile_get_ptr( tileA );

    (void)args;

	/*
	 * Compress the tile if needed
	 */
    if (compressible) {
        rpk_complex64_t *Ala = malloc( tileA->m * tileA->n * sizeof(rpk_complex64_t) );
        rpkx_zlr2ge( ctx, RapackNoTrans, tileA->m, tileA->n, Ara, Ala, tileA->m );

        rpkx_zlrfree( ctx, tileA->m, tileA->n, Ara );
        CORE_zge2lr(
            chameleon_min(tileA->m, tileA->n),
            tileA->m,
            tileA->n,
            Ala,
            tileA->m,
            Ara
        );

        free(Ala);
    }
    return CHAMELEON_SUCCESS;
}

int CHAMELEON_zge2lr_Tile( cham_uplo_t uplo, CHAM_desc_t *descAt ) {
    CHAM_context_t *chamctxt;
    int status = CHAMELEON_SUCCESS;
	int i;

    cham_map_operator_t zge2lr_op;
    zge2lr_op.name = "zge2lr";
    zge2lr_op.cpufunc = zge2lr_op_cpu;
    zge2lr_op.cudafunc = NULL;
    zge2lr_op.hipfunc = NULL;

    cham_map_data_t zge2lr_data[1];
    zge2lr_data[0].access = ChamRW;
    zge2lr_data[0].desc = descAt;

    assert( descAt->format == ChamRapack );

    CHAMELEON_mapv_Tile( uplo, 1, zge2lr_data, &zge2lr_op, NULL);

    return status;
}

