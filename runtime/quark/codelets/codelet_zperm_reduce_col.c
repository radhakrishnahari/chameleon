/**
 *
 * @file quark/codelet_zperm_reduce_col.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark codelets to do the reduction
 *
 * @version 1.4.0
 * @author Matteo Marcos
 * @date 2025-04-11
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"

void
INSERT_TASK_zperm_reduce_col( const RUNTIME_option_t *options,
                              cham_dir_t              dir,
                              const CHAM_desc_t      *A,
                              CHAM_desc_t            *U,
                              int                     Um,
                              int                     Un,
                              CHAM_ipiv_t            *ipiv,
                              int                     ipivk,
                              int                     m,
                              int                     k,
                              void                   *ws )
{
    (void)options;
    (void)A;
    (void)U;
    (void)Um;
    (void)Un;
    (void)ipiv;
    (void)ipivk;
    (void)m;
    (void)k;
    (void)ws;
}
