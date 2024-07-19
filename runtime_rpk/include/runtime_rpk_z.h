/**
 *
 * @file runtime_rpk_z.h
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon rapack context management routines
 *
 * @version 1.3.0
 * @author Abel Calluaud
 * @date 2024-04-26
 * @precisions normal z -> c d s
 *
 */

#ifndef _runtime_zrpk_h_
#define _runtime_zrpk_h_

struct rpk_ctx_s;
typedef struct rpk_ctx_s rpk_ctx_t;

#include "chameleon/config.h"
#include "chameleon/struct.h"
#include <starpu.h>

BEGIN_C_DECLS

extern rpk_ctx_t *cham_rpk_zctx;

/* Get the current double complex rapack ctx */
const rpk_ctx_t *
runtime_rpk_zctx_get( void );

/* Set the current runtime handle ctx */
void
runtime_rpk_zhandle_set( starpu_data_handle_t handle );	

/* Convert tile from lowrank to rapack format */
void
runtime_ztile2rpk( CHAM_tile_t *tile );

END_C_DECLS

#endif
