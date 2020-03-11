/**
 *
 * @file zplghe_batch.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon batch zplghe wrappers
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

struct zplghe_batch_args_s {
    double                 bump;
    unsigned long long int seed;
};
typedef struct zplghe_batch_args_s zplghe_batch_args_t;

static int
chameleon_zplghe_batch_operator( cham_uplo_t uplo, int m, int n,
                                 const CHAM_desc_t *descA, CHAM_tile_t *tileA,
                                 void *op_args )
{
    zplghe_batch_args_t *args = (zplghe_batch_args_t*)op_args;

#if !defined(CHAMELEON_SIMULATION)
    TCORE_zplghe( args->bump, descA->mb, descA->nb, tileA,
                  descA->m, m * descA->mb, n * descA->nb, args->seed );
#endif

    (void)uplo;
    return 0;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zplghe_batch_Tile - Performs multiple hermitian matrix generation
 *  in parallel.
 *
 *******************************************************************************
 *
 * @param[in] bump
 *          The value to add to the diagonal of each tile to be sure
 *          they are positive definite matrices.
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 * @param[in,out] A
 *          A is a collection of mt-by-nt tiles of size A->mb by A->nb
 *          On exit, each tile is initialized as an hermitian matrix.
 *
 *******************************************************************************
 *
 * @return CHAMELEON_SUCCESS on successful exit
 * @return CHAMELEON_ERR_... on error
 *
 */
int CHAMELEON_zplghe_batch_Tile( double bump, CHAM_desc_t *A,
                                 unsigned long long int seed )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    zplghe_batch_args_t params = { bump, seed };
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zplghe_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    chameleon_pmap( ChamUpperLower, A,
                    chameleon_zplghe_batch_operator, &params,
                    sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}
