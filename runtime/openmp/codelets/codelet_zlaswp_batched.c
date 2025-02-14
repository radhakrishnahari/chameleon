/**
 *
 * @file openmp/codelet_zlaswp_batched.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon OpenMP codelets to apply zlaswp on a panel
 *
 * @version 1.3.0
 * @author Alycia Lisito
 * @date 2024-11-12
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks_z.h"

void INSERT_TASK_zlaswp_batched( const RUNTIME_option_t *options,
                                 int                     m0,
                                 int                     minmn,
                                 void                   *ws,
                                 const CHAM_ipiv_t      *ipiv,
                                 int                     ipivk,
                                 const CHAM_desc_t      *Am,
                                 int                     Amm,
                                 int                     Amn,
                                 const CHAM_desc_t      *Ak,
                                 int                     Akm,
                                 int                     Akn,
                                 const CHAM_desc_t      *U,
                                 int                     Um,
                                 int                     Un,
                                 void                  **clargs_ptr )
{
    (void)options;
    (void)m0;
    (void)minmn;
    (void)ws;
    (void)ipiv;
    (void)ipivk;
    (void)Am;
    (void)Amm;
    (void)Amn;
    (void)Ak;
    (void)Akm;
    (void)Akn;
    (void)U;
    (void)Um;
    (void)Un;
    (void)clargs_ptr;
}

void INSERT_TASK_zlaswp_batched_flush( const RUNTIME_option_t *options,
                                       const CHAM_ipiv_t      *ipiv,
                                       int                     ipivk,
                                       const CHAM_desc_t      *Ak,
                                       int                     Akm,
                                       int                     Akn,
                                       const CHAM_desc_t      *U,
                                       int                     Um,
                                       int                     Un,
                                       void                  **clargs_ptr )
{
    (void)options;
    (void)ipiv;
    (void)ipivk;
    (void)Ak;
    (void)Akm;
    (void)Akn;
    (void)U;
    (void)Um;
    (void)Un;
    (void)clargs_ptr;
}
