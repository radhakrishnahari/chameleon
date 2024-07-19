/**
 *
 * @file zfmt2fmt.c
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon RAPACK tile format conversion routine
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
 *  CHAMELEON_zfmt2fmt_Tile constructs a low rank matrix B from a full rank
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
 * @sa CHAMELEON_zfmt2fmt
 * @sa CHAMELEON_zfmt2fmt_Tile_Async
 * @sa CHAMELEON_cfmt2fmt_Tile
 * @sa CHAMELEON_dfmt2fmt_Tile
 * @sa CHAMELEON_sfmt2fmt_Tile
 *
 */
#include "chameleon.h"
#include "runtime_rpk.h"
#include "coreblas.h"

static int zfmt2fmt_op_cpu( void *args, cham_uplo_t uplo, int m, int n, int ndata,
                   const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... ) {
    int compressible;
    CHAMELEON_Complex64_t *Ala;
	rpk_matrix_t *Ara;
    int status = CHAMELEON_SUCCESS;
    const CHAM_desc_t *descB;
    CHAM_tile_t *tileB;
    cham_mtxfmt_t fmtA;
    cham_mtxfmt_t fmtB;
   
    va_list varg_list;
    va_start(varg_list, tileA);
    descB = va_arg(varg_list, const CHAM_desc_t*);
    tileB = va_arg(varg_list, CHAM_tile_t*);
    va_end(varg_list);

    fmtA = descA->format;
    fmtB = descB->format;

    if( 0 ) {
    }
#if defined(CHAMELEON_USE_HMATOSS) && defined(CHAMELEON_USE_RAPACK)
    else if( (fmtA == ChamHmatOSS) && (fmtB == ChamRapack) ) {
        TCORE_zhmat2rpk( tileA, tileB );
    }
#endif
    else if ( (fmtA == ChamTile) && (fmtB == ChamRapack) ) {
        TCORE_ztile2rpk( tileA, tileB );
    }
    else {
        fprintf(stderr, "Unsupported format conversion from %d to %d\n", fmtA, fmtB);
        chameleon_error("CHAMELEON_zfmt2fmt_Tile", "unsupported format conversion");
        return CHAMELEON_ERR_NOT_SUPPORTED;
    }

    return status;
}

int CHAMELEON_zfmt2fmt_Tile( cham_uplo_t uplo, CHAM_desc_t *descAt, CHAM_desc_t *descBt) {
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

	CHAMELEON_zfmt2fmt_Tile_Async( uplo, descAt, descBt, sequence, &request );

	CHAMELEON_Desc_Flush( descAt, sequence );

	chameleon_sequence_wait( chamctxt, sequence );
	status = sequence->status;
	chameleon_sequence_destroy( chamctxt, sequence );
	return status;
}

int CHAMELEON_zfmt2fmt_Tile_Async( cham_uplo_t uplo, CHAM_desc_t *descAt, CHAM_desc_t *descBt, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request ) {
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

    assert(descAt->format != descBt->format);

	/* Check descriptors for correctness */
    if (chameleon_desc_check(descAt) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zfmt2fmt_Tile_Async", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

	cham_map_operator_t zfmt2fmt_op = {0};
    zfmt2fmt_op.name = "zfmt2fmt";
    zfmt2fmt_op.cpufunc = zfmt2fmt_op_cpu;

    cham_map_data_t zfmt2fmt_data[2];
    zfmt2fmt_data[0].access = ChamR;
    zfmt2fmt_data[0].desc = descAt;
    zfmt2fmt_data[1].access = ChamW;
    zfmt2fmt_data[1].desc = descBt;

	CHAMELEON_mapv_Tile(
		uplo, 2, zfmt2fmt_data, &zfmt2fmt_op, NULL
	);

	return status;
}
