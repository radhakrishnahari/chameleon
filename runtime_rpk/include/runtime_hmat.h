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

#ifndef _runtime_hmat_h_
#define _runtime_hmat_h_

#include "hmat/hmat.h"
#include "chameleon/config.h"
#include "runtime_rpk/runtime_hmat_z.h"
#include "runtime_rpk/runtime_hmat_c.h"
#include "runtime_rpk/runtime_hmat_d.h"
#include "runtime_rpk/runtime_hmat_s.h"

BEGIN_C_DECLS

void runtime_hmat_interface_init( void );

void runtime_hmat_interface_exit( void );

hmat_interface_t *runtime_hmat_interface_get( hmat_value_t ftype );

END_C_DECLS

#endif
