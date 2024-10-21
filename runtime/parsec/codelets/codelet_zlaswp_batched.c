/**
 *
 * @file parsec/codelet_zlaswp_batched.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Parsec codelets to apply zlaswp on a panel
 *
 * @version 1.3.0
 * @author Alycia Lisito
 * @date 2024-10-21
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"

void INSERT_TASK_zlaswp_batched( const RUNTIME_option_t *options,
                                 int                     m0,
                                 int                     minmn,
                                 int                     k,
                                 int                     m,
                                 int                     n,
                                 void                   *ws,
                                 const CHAM_ipiv_t      *ipiv,
                                 int                     ipivk,
                                 const CHAM_desc_t      *A,
                                 const CHAM_desc_t      *Wu,
                                 void                  **clargs_ptr )
{
    (void)options;
    (void)m0;
    (void)minmn;
    (void)k;
    (void)m;
    (void)n;
    (void)ws;
    (void)ipiv;
    (void)ipivk;
    (void)A;
    (void)Wu;
    (void)clargs_ptr;
}

void INSERT_TASK_zlaswp_batched_flush( const RUNTIME_option_t *options,
                                       int                     k,
                                       int                     n,
                                       const CHAM_ipiv_t      *ipiv,
                                       int                     ipivk,
                                       const CHAM_desc_t      *A,
                                       const CHAM_desc_t      *U,
                                       void                  **clargs_ptr )
{
    (void)options;
    (void)k;
    (void)n;
    (void)ipiv;
    (void)ipivk;
    (void)A;
    (void)U;
    (void)clargs_ptr;
}
