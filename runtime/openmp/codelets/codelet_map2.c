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

// Define the access function interface
typedef void (*access_fct_t)(cham_uplo_t uplo,
                                 CHAM_tile_t *tileA, const CHAM_desc_t *A, int Am, int An,
                                 CHAM_tile_t *tileB, const CHAM_desc_t *B,
                                 cham_binary_operator_t op_fct, void *op_args);

void access_fct_R_R(cham_uplo_t uplo,
                 CHAM_tile_t *tileA, const CHAM_desc_t *A, int Am, int An,
                 CHAM_tile_t *tileB, const CHAM_desc_t *B,
                 cham_binary_operator_t op_fct, void *op_args) {
#pragma omp task depend( in:tileA[0], tileB[0])
    {
        op_fct(uplo, Am, An, A, tileA, B, tileB, op_args);
    }
}

void access_fct_R_W(cham_uplo_t uplo,
                 CHAM_tile_t *tileA, const CHAM_desc_t *A, int Am, int An,
                 CHAM_tile_t *tileB, const CHAM_desc_t *B,
                 cham_binary_operator_t op_fct, void *op_args) {
#pragma omp task depend( in:tileA[0] ) depend (out: tileB[0] )
    {
        op_fct(uplo, Am, An, A, tileA, B, tileB, op_args);
    }
}

void access_fct_R_RW(cham_uplo_t uplo,
                  CHAM_tile_t *tileA, const CHAM_desc_t *A, int Am, int An,
                  CHAM_tile_t *tileB, const CHAM_desc_t *B,
                  cham_binary_operator_t op_fct, void *op_args) {
#pragma omp task depend( in:tileA[0] ) depend (inout: tileB[0] )
    {
        op_fct(uplo, Am, An, A, tileA, B, tileB, op_args);
    }
}

void access_fct_W_R(cham_uplo_t uplo,
                 CHAM_tile_t *tileA, const CHAM_desc_t *A, int Am, int An,
                 CHAM_tile_t *tileB, const CHAM_desc_t *B,
                 cham_binary_operator_t op_fct, void *op_args) {
#pragma omp task depend( out:tileA[0]) depend( in:tileB[0] )
    {
        op_fct(uplo, Am, An, A, tileA, B, tileB, op_args);
    }
}

void access_fct_W_W(cham_uplo_t uplo,
                 CHAM_tile_t *tileA, const CHAM_desc_t *A, int Am, int An,
                 CHAM_tile_t *tileB, const CHAM_desc_t *B,
                 cham_binary_operator_t op_fct, void *op_args) {
#pragma omp task depend( out: tileA[0], tileB[0])
    {
        op_fct(uplo, Am, An, A, tileA, B, tileB, op_args);
    }
}

void access_fct_W_RW(cham_uplo_t uplo,
                  CHAM_tile_t *tileA, const CHAM_desc_t *A, int Am, int An,
                  CHAM_tile_t *tileB, const CHAM_desc_t *B,
                  cham_binary_operator_t op_fct, void *op_args) {
#pragma omp task depend( out:tileA[0] ) depend (inout: tileB[0] )
    {
        op_fct(uplo, Am, An, A, tileA, B, tileB, op_args);
    }
}

void access_fct_RW_R(cham_uplo_t uplo,
                  CHAM_tile_t *tileA, const CHAM_desc_t *A, int Am, int An,
                  CHAM_tile_t *tileB, const CHAM_desc_t *B,
                  cham_binary_operator_t op_fct, void *op_args) {
#pragma omp task depend( inout:tileA[0]) depend( in:tileB[0] )
    {
        op_fct(uplo, Am, An, A, tileA, B, tileB, op_args);
    }
}


void access_fct_RW_W(cham_uplo_t uplo,
                  CHAM_tile_t *tileA, const CHAM_desc_t *A, int Am, int An,
                  CHAM_tile_t *tileB, const CHAM_desc_t *B,
                  cham_binary_operator_t op_fct, void *op_args) {
#pragma omp task depend( inout:tileA[0]) depend( out:tileB[0] )
    {
        op_fct(uplo, Am, An, A, tileA, B, tileB, op_args);
    }
}


void access_fct_RW_RW(cham_uplo_t uplo,
                   CHAM_tile_t *tileA, const CHAM_desc_t *A, int Am, int An,
                   CHAM_tile_t *tileB, const CHAM_desc_t *B,
                   cham_binary_operator_t op_fct, void *op_args) {
#pragma omp task depend( inout:tileA[0], tileB[0] )
    {
        op_fct(uplo, Am, An, A, tileA, B, tileB, op_args);
    }
}


static access_fct_t array_access_fct_t[3][3] = {
        {
                access_fct_R_R, access_fct_R_W, access_fct_R_RW
        },
        {
                access_fct_W_R, access_fct_W_W, access_fct_W_RW
        },
        {
                access_fct_RW_R, access_fct_RW_W, access_fct_RW_RW
        }
};

void INSERT_TASK_map2(const RUNTIME_option_t *options,
                      cham_access_t accessA, cham_access_t accessB,
                      cham_uplo_t uplo,
                      const CHAM_desc_t *A, int Am, int An,
                      const CHAM_desc_t *B, int Bm, int Bn,
                      cham_binary_operator_t op_fct, void *op_args,
                      const char *name) {
    CHAM_tile_t *tileA = A->get_blktile(A, Am, An);
    CHAM_tile_t *tileB = B->get_blktile(B, Bm, Bn);

    array_access_fct_t[accessA - 1][accessB - 1](uplo, tileA, A, Am, An, tileB, B, op_fct, op_args);

    (void) options;
}
