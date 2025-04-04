/**
 *
 * @file runtime_lrmem.h
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon rapack memory management functions
 *
 * @version 1.3.0
 * @author Abel Calluaud
 * @date 2024-04-24
 *
 */
#ifndef _rpk_lrmem_h_
#define _rpk_lrmem_h_

#include "chameleon/config.h"
#include <stdint.h>
#include <stdlib.h>

BEGIN_C_DECLS

struct rpk_ctx_s;
struct rpk_matrix_s;

typedef struct rpk_ctx_s rpk_ctx_t;
typedef struct rpk_matrix_s rpk_matrix_t;
typedef int rpk_int_t;

#if 0
int runtime_rpk_lrmem_alloc(
    const rpk_ctx_t *ctx,
    void           **new_u,
    void           **new_v,
    size_t          u_size,
    size_t          v_size,
    size_t          elem_size
);

int runtime_rpk_lrmem_free(
    const rpk_ctx_t *ctx,
    void            *original_u,
    void            *original_v,
    size_t           u_size,
    size_t           v_size,
    size_t           elem_size
);

int runtime_rpk_lrmem_set(
    const rpk_ctx_t *ctx,
    rpk_matrix_t    *Ara,
    void            *new_u,
    void            *new_v
); 
#endif

void *
runtime_rpk_malloc( const rpk_ctx_t *ctx, size_t size );

void
runtime_rpk_free( const rpk_ctx_t *ctx, void *ptr, size_t size );

END_C_DECLS

#endif
