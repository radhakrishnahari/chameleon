/**
 *
 * @file parsec/codelet_zlaswp_batched.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Parsec codelets to apply zlaswp on a panel
 *
 * @version 1.3.0
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-07-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"

void INSERT_TASK_zlaswp_batched( const RUNTIME_option_t *options,
                                 cham_side_t             side,
                                 cham_dir_t              dir,
                                 int                     m0,
                                 int                     m,
                                 int                     n,
                                 int                     k,
                                 void                   *ws,
                                 const CHAM_ipiv_t      *ipiv, int ipivk,
                                 const CHAM_desc_t      *Am,   int Amm, int Amn,
                                 const CHAM_desc_t      *Ak,   int Akm, int Akn,
                                 const CHAM_desc_t      *U,    int Um,  int Un,
                                 void                  **clargs_ptr )
{
    assert( 0 );
    (void)options;
    (void)dir;
    (void)m0;
    (void)m;
    (void)n;
    (void)k;
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
                                       cham_dir_t              dir,
                                       const CHAM_ipiv_t      *ipiv, int ipivk,
                                       const CHAM_desc_t      *Ak,   int Akm, int Akn,
                                       const CHAM_desc_t      *U,    int Um,  int Un,
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
