/**
 *
 * @file core_zge2lr.c
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zge2lr CPU kernel
 *
 * @version 1.2.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Abel Calluaud
 * @date 2024-03-20
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"
#include "runtime_rpk.h"

void CORE_zge2lr( rpk_int_t        rklimit,
                  rpk_int_t        m,
                  rpk_int_t        n,
                  const CHAMELEON_Complex64_t* A,
                  rpk_int_t        lda,
                  rpk_matrix_t    *Alr ) {

    const rpk_ctx_t * ctx = runtime_rpk_zctx_get();
    ctx->rpk_ge2lr(
        ctx,
        rklimit,
        m,
        n,
        (void*)A,
        lda,
        Alr
    );
}

/**/








