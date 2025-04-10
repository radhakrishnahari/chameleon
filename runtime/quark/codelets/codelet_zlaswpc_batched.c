/**
 *
 * @file quark/codelet_zlaswpc_batched.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon quark codelets to apply zlaswp on a panel
 *
 * @version 1.3.0
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2024-11-12
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"

void INSERT_TASK_zlaswpc_batched( const RUNTIME_option_t *options,
                                  cham_dir_t              dir,
                                  int                     n0,
                                  int                     n,
                                  void                   *ws,
                                  const CHAM_ipiv_t      *ipiv,
                                  int                     ipivk,
                                  const CHAM_desc_t      *An,
                                  int                     Anm,
                                  int                     Ann,
                                  const CHAM_desc_t      *Ak,
                                  int                     Akm,
                                  int                     Akn,
                                  const CHAM_desc_t      *U,
                                  int                     Um,
                                  int                     Un,
                                  void                  **clargs_ptr )
{
    (void)options;
    (void)dir;
    (void)n0;
    (void)n;
    (void)ws;
    (void)ipiv;
    (void)ipivk;
    (void)An;
    (void)Anm;
    (void)Ann;
    (void)Ak;
    (void)Akm;
    (void)Akn;
    (void)U;
    (void)Um;
    (void)Un;
    (void)clargs_ptr;
}

void INSERT_TASK_zlaswpc_batched_flush( const RUNTIME_option_t *options,
                                        cham_dir_t              dir,
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
    (void)dir;
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
