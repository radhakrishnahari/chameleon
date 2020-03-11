/**
 *
 * @file openmp/codelet_map2.c
 *
 * @copyright 2018-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map2 OpenMP codelet
 *
 * @version 1.0.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @date 2020-03-03
 *
 */
#include "chameleon_openmp.h"

void INSERT_TASK_map2( const RUNTIME_option_t *options,
                       cham_uplo_t uplo,
                       const CHAM_desc_t *A, int Am, int An,
                       const CHAM_desc_t *B, int Bm, int Bn,
                       cham_binary_operator_t op_fct, void *op_args )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    CHAM_tile_t *tileB = B->get_blktile( B, Bm, Bn );

#pragma omp task depend( in:tileA[0] ) depend( inout:tileB[0] )
    {
        op_fct( uplo, Am, An, A, tileA, B, tileB, op_args );
    }
}
