/**
 *
 * @file openmp/codelet_map.c
 *
 * @copyright 2018-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon map OpenMP codelet
 *
 * @version 1.2.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @date 2022-02-22
 *
 */
#include "chameleon_openmp.h"

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_access_t accessA, cham_uplo_t uplo, const CHAM_desc_t *A, int Am, int An,
                      cham_unary_operator_t op_fct, void *op_args )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

    switch( accessA ) {
    case ChamW:
#pragma omp task depend( out: tileA[0] )
    {
        op_fct( A, uplo, Am, An, tileA, op_args );
    }
    break;

    case ChamR:
#pragma omp task depend( in: tileA[0] )
    {
        op_fct( uplo, Am, An, A, tileA, op_args );
    }

    break;

    case ChamRW:
    default:
#pragma omp task depend( inout: tileA[0] )
    {
        op_fct( A, uplo, Am, An, tileA, op_args );
    }
    }

    (void)options;
}
