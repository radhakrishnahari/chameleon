/**
 *
 * @file zrpk2fmt.c
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon RAPACK descriptor format conversion routines
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
 *  CHAMELEON_zrpk2fmt_Tile constructs a low rank matrix B from a full rank
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
 * @sa CHAMELEON_zrpk2fmt
 * @sa CHAMELEON_zrpk2fmt_Tile_Async
 * @sa CHAMELEON_crpk2fmt_Tile
 * @sa CHAMELEON_drpk2fmt_Tile
 * @sa CHAMELEON_srpk2fmt_Tile
 *
 */
#include "chameleon.h"
#include "runtime_rpk.h"
#include "coreblas.h"

typedef struct zrpk2fmt_args_s {
    cham_mtxfmt_t fmtA;
} zrpk2fmt_args_t;

static int zrpk2fmt_op_cpu( void *args, cham_uplo_t uplo, int m, int n, int ndata,
                   const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... ) {
    zrpk2fmt_args_t *zrpk2fmt_args = (zrpk2fmt_args_t*)args;
    cham_mtxfmt_t fmtA = zrpk2fmt_args->fmtA;
    int compressible;
    CHAMELEON_Complex64_t *Ala;
	rpk_matrix_t *Ara;
    int status = CHAMELEON_SUCCESS;
	
	const rpk_ctx_t *ctx = runtime_rpk_zctx_get();
    
	assert( tileA->m == tileA->ld );
    
    /* Hierarchical data format is not supported yet */
    assert( tileA->format == CHAMELEON_TILE_LOWRANK );

    switch( fmtA ) {
        case ChamTFR:
            {
                Ara = (rpk_matrix_t*)CHAM_tile_get_ptr( tileA );

                if (Ara->rk == -1) {
                    return status;
                }

                Ala = (CHAMELEON_Complex64_t*)runtime_rpk_malloc( ctx, tileA->m * tileA->n * sizeof(CHAMELEON_Complex64_t) );

                rpkx_zlr2ge( ctx, RapackNoTrans, tileA->m, tileA->n, Ara, Ala, tileA->m );

                rpkx_zlrsze( ctx, 0, tileA->m, tileA->n, Ara, -1, -1, -1 );
                
                assert( Ara->rk == -1 );

                CORE_zlacpy(
                    uplo,
                    tileA->m, tileA->n,
                    Ala, tileA->m,
                    Ara->u, tileA->m
                );

                runtime_rpk_free( ctx, Ala, tileA->m * tileA->n * sizeof(CHAMELEON_Complex64_t) );
            }
            break;
        case ChamTLR:
            break;
        case ChamTile:
            fprintf(stderr, "error: could not convert to tile format\n");
            exit(1);
            break;
#if defined(CHAMELEON_USE_HMATOSS)
        case ChamHmatOSS:
            fprintf(stderr, "error: Could not convert to Hmat-OSS format\n");
            exit(1);
#endif 
            break;
    }

    return status;
}

int CHAMELEON_zrpk2fmt_Tile( cham_mtxfmt_t fmtA, cham_uplo_t uplo, CHAM_desc_t *descAt ) {
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status = CHAMELEON_SUCCESS;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlacpy_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

	CHAMELEON_zrpk2fmt_Tile_Async( fmtA, uplo, descAt, sequence, &request );

	CHAMELEON_Desc_Flush( descAt, sequence );

	chameleon_sequence_wait( chamctxt, sequence );
	status = sequence->status;
	chameleon_sequence_destroy( chamctxt, sequence );
	return status;
}

int CHAMELEON_zrpk2fmt_Tile_Async( cham_mtxfmt_t fmtA, cham_uplo_t uplo, CHAM_desc_t *descAt, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request ) {
    CHAM_context_t *chamctxt;
    int status = CHAMELEON_SUCCESS;
	int i;

	chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zgemm_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zgemm_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zgemm_Tile_Async", "NULL request");
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
    if (chameleon_desc_check(descAt) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zrpk2fmt_Tile_Async", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

	cham_map_operator_t zrpk2fmt_op;
    zrpk2fmt_op.name = "zrpk2fmt";
    zrpk2fmt_op.cpufunc = zrpk2fmt_op_cpu;
    zrpk2fmt_op.cudafunc = NULL;
    zrpk2fmt_op.hipfunc = NULL;

    cham_map_data_t zrpk2fmt_data[1];
    zrpk2fmt_data[0].access = ChamRW;
    zrpk2fmt_data[0].desc = descAt;

    zrpk2fmt_args_t zrpk2fmt_args;
    zrpk2fmt_args.fmtA = fmtA;
    
    assert( (descAt->format == ChamTLR) || (descAt->format == ChamTFR) );

	CHAMELEON_mapv_Tile(
		uplo, 1, zrpk2fmt_data, &zrpk2fmt_op, (void*)&zrpk2fmt_args
	);

	/* TODO: Find out why the async version is not working */
	
/*
 *     CHAMELEON_mapv_Tile_Async(
 *         uplo, 1, zrpk2fmt_data, &zrpk2fmt_op, (void*)&zrpk2fmt_args,
 *         sequence, request
 *     );
 */

	return status;
}
