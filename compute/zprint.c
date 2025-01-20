/**
 *
 * @file zprint.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zprint wrappers
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @date 2024-03-14
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#if !defined(CHAMELEON_SIMULATION)
#include <coreblas/coreblas_z.h>
#endif

struct zprint_args_s {
    FILE       *file;
    const char *header;
};

static inline int
zprint_cpu( void *op_args,
            cham_uplo_t uplo, int m, int n, int ndata,
            const CHAM_desc_t *descA, CHAM_tile_t *tileA, ... )
{
    struct zprint_args_s  *options = (struct zprint_args_s *)op_args;
    CHAMELEON_Complex64_t *A = CHAM_tile_get_ptr( tileA );

    int tempmm = descA->get_blkdim( descA, m, DIM_m, descA->m );
    int tempnn = descA->get_blkdim( descA, n, DIM_n, descA->n );
    int lda    = tileA->ld;

    if ( ndata > 1 ) {
        fprintf( stderr, "zprint_cpu: supports only one piece of data and %d have been given\n", ndata );
    }
    assert( tileA->format & CHAMELEON_TILE_FULLRANK );

#if !defined(CHAMELEON_SIMULATION)
    CORE_zprint( options->file, options->header, uplo,
                 tempmm, tempnn, m, n, A, lda );
#endif

    return 0;
}

static cham_map_operator_t zprint_map = {
    .name     = "zprint",
    .cpufunc  = zprint_cpu,
    .cudafunc = NULL,
    .hipfunc  = NULL,
};

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 * @brief Print a matrix tile by tile.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in] A
 *          The M-by-N matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @retval the two-norm estimate.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlatms_Tile
 * @sa CHAMELEON_zlatms_Tile_Async
 * @sa CHAMELEON_clange
 * @sa CHAMELEON_dlange
 * @sa CHAMELEON_slange
 *
 */
int CHAMELEON_zprint( FILE *file, const char *header,
                      cham_uplo_t uplo, int M, int N,
                      CHAMELEON_Complex64_t *A, int LDA )
{
    int                  NB, status;
    CHAM_context_t      *chamctxt;
    RUNTIME_sequence_t  *sequence = NULL;
    RUNTIME_request_t    request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t          descAl, descAt;
    cham_map_data_t      data;
    struct zprint_args_s options = {
        .file   = file,
        .header = header,
    };

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlatms", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        chameleon_error("CHAMELEON_zlatms", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zlatms", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zlatms", "illegal value of LDA");
        return -4;
    }

    /* Quick return */
    if (chameleon_min(N, M) == 0) {
        return (double)0.0;
    }

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zlatms", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, &descAl, &descAt, ChamDescInput, uplo,
                         A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    data.access = ChamR;
    data.desc   = &descAt;
    chameleon_pmap( uplo, 1, &data, &zprint_map, &options, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInput, uplo, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );

    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 * @brief Tile equivalent of CHAMELEON_zprint().
 *
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          On entry, the input matrix A.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlatms
 * @sa CHAMELEON_zlatms_Tile_Async
 * @sa CHAMELEON_clange_Tile
 * @sa CHAMELEON_dlange_Tile
 * @sa CHAMELEON_slange_Tile
 *
 */
int CHAMELEON_zprint_Tile( FILE *file, const char *header,
                           cham_uplo_t uplo, CHAM_desc_t *A )
{
    CHAM_context_t      *chamctxt;
    RUNTIME_sequence_t  *sequence = NULL;
    RUNTIME_request_t    request = RUNTIME_REQUEST_INITIALIZER;
    cham_map_data_t      data;
    struct zprint_args_s options = {
        .file   = file,
        .header = header,
    };
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zlatms_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    data.access = ChamR;
    data.desc   = A;

    chameleon_pmap( uplo, 1, &data, &zprint_map, &options, sequence, &request );
    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}
