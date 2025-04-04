/**
 *
 * @file runtime_rpk_zctx.c
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon rapack context configuration
 *
 * @version 1.2.0
 * @author Abel Calluaud
 * @date 2024-03-27
 * @precisions normal z -> c d s
 *
 */

#include "runtime_rpk.h"
#include "runtime_lrmem.h"
#include "chameleon/constants.h"
#include "chameleon/types.h"
#include "chameleon/struct.h"
#include <rapack.h>
#include <starpu.h>
#include "cham_tile_interface.h"

const rpk_ctx_t *
runtime_rpk_zctx_get( void ) {
	return runtime_rpk_ctx_get( ChamComplexDouble );
}

/*
void
runtime_rpk_zhandle_set( starpu_data_handle_t handle ) {
	int wid = starpu_worker_get_id_check();
	runtime_rpk_handle_set( ChamComplexDouble, wid, handle );
}
*/
