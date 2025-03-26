/**
 *
 * @file parsec/codelet_zperm_allreduce_col.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC codelets to do the reduction
 *
 * @version 1.3.0
 * @author Matteo Marcos
 * @date 2025-03-27
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"

void
INSERT_TASK_zperm_allreduce_send_invp_col( const RUNTIME_option_t *options,
                                           cham_dir_t              dir,
                                           CHAM_ipiv_t            *ipiv,
                                           int                     ipivk,
                                           const CHAM_desc_t      *A,
                                           int                     m,
                                           int                     k )
{
    (void)options;
    (void)ipiv;
    (void)ipivk;
    (void)A;
    (void)m;
    (void)k;
}

void
INSERT_TASK_zperm_allreduce_col( const RUNTIME_option_t *options,
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
