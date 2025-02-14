/**
 *
 * @file openmp/codelet_zbuild.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zbuild OpenMP codelet
 *
 * @version 1.3.0
 * @author Philippe Virouleau
 * @author Mathieu Faverge
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

void INSERT_TASK_zbuild( const RUNTIME_option_t *options,
                         const CHAM_desc_t *A, int Am, int An,
                         void *user_data, void* user_build_callback )
{
    int row_min, row_max, col_min, col_max;
    row_min = Am*A->mb ;
    row_max = Am == A->mt-1 ? A->m-1 : row_min+A->mb-1 ;
    col_min = An*A->nb ;
    col_max = An == A->nt-1 ? A->n-1 : col_min+A->nb-1 ;
    CHAM_tile_t *tileA = A->get_blktile( A, Am, An );
    void ( *callback )( int row_min, int row_max, int col_min, int col_max, void *buffer, int ld, void *user_data ) ;
    callback = user_build_callback;

#pragma omp task firstprivate( row_min, row_max, col_min, col_max, tileA, user_data ) depend( inout:tileA[0] )
    callback( row_min, row_max, col_min, col_max, tileA->mat, tileA->ld, user_data );

    (void)options;
}
