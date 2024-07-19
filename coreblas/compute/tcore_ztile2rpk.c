/**
 *
 * @file ztile2rpk.cc
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon tile to RAPACK format conversion
 *
 * @version 1.2.0
 * @author Abel Calluaud
 * @date 2024-06-28
 * @precisions normal z -> c d s
 *
 */

#include "runtime_rpk.h"
#include "runtime_lrmem.h"
#include "chameleon/config.h"
#include "chameleon/constants.h"
#include "chameleon/types.h"
#include "chameleon/struct.h"
#include "coreblas.h"
#include <rapack.h>
#include <starpu.h>
#include "cham_tile_interface.h"
#include "runtime_rpk.h"

void
TCORE_ztile2rpk( CHAM_tile_t *tileA, CHAM_tile_t *tileB )
{
    assert( tileA->m == tileB->m );
    assert( tileA->n == tileB->n );
    assert( tileA->format == CHAMELEON_TILE_FULLRANK );
    assert( tileB->format == CHAMELEON_TILE_LOWRANK );

    rpk_matrix_t *Bra = (rpk_matrix_t *) CHAM_tile_get_ptr( tileB );
    CHAMELEON_Complex64_t *Ala = (CHAMELEON_Complex64_t *) CHAM_tile_get_ptr( tileA );
    int compressible = (tileB->m != tileB->n);
    const rpk_ctx_t *ctx = runtime_rpk_zctx_get();

    /* TODO: Should we try to compress ? */
    compressible = 0;

    if (compressible) {
        rpkx_zlrfree( ctx, tileB->m, tileB->n, Bra );
        CORE_zge2lr(
            chameleon_min(tileB->m, tileB->n),
            tileB->m,
            tileB->n,
            Ala,
            tileB->m,
            Bra
        );
    }
    else {
        rpkx_zlrsze( ctx, 0, tileB->m, tileB->n, Bra, -1, -1, -1 );
        CORE_zlacpy(
            ChamUpperLower,
            tileB->m, tileB->n,
            Ala,
            tileA->ld,
            Bra->u,
            tileB->m
        );
    }
}
