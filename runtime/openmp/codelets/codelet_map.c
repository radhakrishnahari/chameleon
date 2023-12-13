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

// Define the access function interface
typedef void (*access_fct_t)(cham_uplo_t uplo,
                                 CHAM_tile_t *tileA, const CHAM_desc_t *A, int Am, int An,
                                 cham_unary_operator_t op_fct, void *op_args);

void access_fct_R(cham_uplo_t uplo,
                 CHAM_tile_t *tileA, const CHAM_desc_t *A, int Am, int An,
                 cham_unary_operator_t op_fct, void *op_args) {
#pragma omp task depend( in:tileA[0])
    {
        op_fct( A, uplo, Am, An, tileA, op_args );
    }
}

void access_fct_W(cham_uplo_t uplo,
               CHAM_tile_t *tileA, const CHAM_desc_t *A, int Am, int An,
               cham_unary_operator_t op_fct, void *op_args) {
#pragma omp task depend( out:tileA[0])
    {
        op_fct( A, uplo, Am, An, tileA, op_args );
    }
}

void access_fct_RW(cham_uplo_t uplo,
               CHAM_tile_t *tileA, const CHAM_desc_t *A, int Am, int An,
               cham_unary_operator_t op_fct, void *op_args) {
#pragma omp task depend( inout:tileA[0])
    {
        op_fct( A, uplo, Am, An, tileA, op_args );
    }
}

static access_fct_t array_access_fct_t[3] = { access_fct_R, access_fct_W, access_fct_RW };

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_access_t accessA, cham_uplo_t uplo, const CHAM_desc_t *A, int Am, int An,
                      cham_unary_operator_t op_fct, void *op_args, const char *name )
{
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );

    array_access_fct_t[accessA - 1](uplo, tileA, A, Am, An, op_fct, op_args);

    (void)options;
}
