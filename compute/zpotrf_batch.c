/**
 *
 * @file zpotrf_batch.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon batch zpotrf wrappers
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

struct zpotrf_batch_args_s {
    cham_uplo_t uplo;
};
typedef struct zpotrf_batch_args_s zpotrf_batch_args_t;

static int
chameleon_zpotrf_batch_operator( cham_uplo_t uplo, int m, int n,
                                 const CHAM_desc_t *descA, CHAM_tile_t *tileA,
                                 void *op_args )
{
    int info = 0;
    zpotrf_batch_args_t *args = (zpotrf_batch_args_t*)op_args;

#if !defined(CHAMELEON_SIMULATION)
    TCORE_zpotrf(
        args->uplo, descA->nb, tileA, &info );
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
 *  CHAMELEON_zpotrf_batch_Tile - Performs multiple Cholesky factorization in parallel.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = ChamUpper: Upper triangle of A tiles are stored;
 *          = ChamLower: Lower triangle of A tiles are stored.
 *
 * @param[in,out] A
 *          A is a collection of mt-by-nt tiles of size A->mb by A->nb
 *          On exit, each tile is factorized with Cholesky factorization.
 *
 *******************************************************************************
 *
 * @return CHAMELEON_SUCCESS on successful exit
 * @return CHAMELEON_ERR_... on error
 *
 */
int CHAMELEON_zpotrf_batch_Tile( cham_uplo_t uplo, CHAM_desc_t *A )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    zpotrf_batch_args_t params = { uplo };
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zpotrf_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    chameleon_pmap( ChamRW, ChamUpperLower, A,
                    chameleon_zpotrf_batch_operator, &params,
                    sequence, &request, "zpotrf" );

    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}
