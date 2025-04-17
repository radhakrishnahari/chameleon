/**
 *
 * @file parsec/codelet_zlaswpc.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC codelets to apply zlaswpc on a panel
 *
 * @version 1.3.0
 * @author Matteo Marcos
 * @date 2025-03-27
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

void INSERT_TASK_zlaswpc_get( const RUNTIME_option_t *options,
                              cham_dir_t dir, int k, int n0,
                              const CHAM_ipiv_t *ipiv, int ipivk,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *U, int Um, int Un )
{
    (void)options;
    (void)dir;
    (void)k;
    (void)n0;
    (void)ipiv;
    (void)ipivk;
    (void)A;
    (void)Am;
    (void)An;
    (void)U;
    (void)Um;
    (void)Un;
}

void INSERT_TASK_zlaswpc_set( const RUNTIME_option_t *options,
                              cham_dir_t dir, int k, int n0,
                              const CHAM_ipiv_t *ipiv, int ipivk,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *B, int Bm, int Bn )
{
    (void)options;
    (void)dir;
    (void)k;
    (void)n0;
    (void)ipiv;
    (void)ipivk;
    (void)A;
    (void)Am;
    (void)An;
    (void)B;
    (void)Bm;
    (void)Bn;
}

