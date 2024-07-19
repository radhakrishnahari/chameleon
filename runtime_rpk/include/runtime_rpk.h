/**
 *
 * @file runtime_rpk.h
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
 *
 */

#ifndef _runtime_rpk_h_
#define _runtime_rpk_h_

#include "chameleon/config.h"
#include "chameleon/constants.h"
#include "runtime_rpk/runtime_rpk_z.h"
#include "runtime_rpk/runtime_rpk_c.h"
#include "runtime_rpk/runtime_rpk_d.h"
#include "runtime_rpk/runtime_rpk_s.h"

BEGIN_C_DECLS

struct cham_starpu_rpk_ctx_s {
	int node; /*< memory node of a runtime-managed Rapack matrix */
	starpu_data_handle_t handle; /*< handle number of the matrix in task parameters */
};

typedef struct cham_starpu_rpk_ctx_s cham_starpu_rpk_ctx_t;

void runtime_rpk_ctx_init( void );

void runtime_rpk_ctx_exit( void );

/*
void runtime_rpk_handle_set( cham_flttype_t flttype, int wid, starpu_data_handle_t handle );

void runtime_rpk_node_set( cham_flttype_t flttype, int wid, unsigned node );

void runtime_rpk_handle_get( cham_flttype_t flttype, int wid, int node );

const rpk_ctx_t *runtime_rpk_ctx_get( cham_flttype_t type, int wid );
*/

const rpk_ctx_t *runtime_rpk_ctx_get( cham_flttype_t type );

END_C_DECLS

#endif
