/**
 *
 * @file core_zlrcrk.c
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon CPU kernel interface from CHAM_tile_t layout to the real one.
 *
 * @version 1.3.0
 * @author Abel Calluaud
 * @date 2023-08-31
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"
#include "coreblas/coreblas_ztile.h"
#if defined(CHAMELEON_USE_HMATOSS)
#include "coreblas/hmat.h"
#endif

int
TCORE_zlrcrk( const CHAM_tile_t *tileA )
{
    int rank = -1;
    if ( tileA->format & CHAMELEON_TILE_LOWRANK ) {
        rpk_matrix_t *Ara = CHAM_tile_get_ptr( tileA );
        rank = Ara->rk;
    }
#if defined(CHAMELEON_USE_HMATOSS)
    else if ( tileA->format & CHAMELEON_TILE_HMAT ) {
        hmat_matrix_t *hmat = CHAM_tile_get_ptr( tileA );
        rank = hmat_zlrcrk( hmat );
    }
#endif
    else if ( tileA->format & CHAMELEON_TILE_FULLRANK ) {
        rank = -1;
    }
    else {
        fprintf(stderr, "error: unknown tile format\n");
        exit(2);
    }
    return rank;
}
