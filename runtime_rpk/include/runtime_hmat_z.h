/**
 *
 * @file runtime_hmat_z.h
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Hmat-oss management routines
 *
 * @version 1.3.0
 * @author Abel Calluaud
 * @date 2024-05-20
 * @precisions normal z -> c d s
 *
 */

#ifndef _runtime_zhmat_h_
#define _runtime_zhmat_h_

#include "chameleon/config.h"
#include "chameleon/struct.h"
#include <starpu.h>

BEGIN_C_DECLS

/* Get the current double complex rapack interface */
hmat_interface_t *
runtime_hmat_zinterface_get( void );

END_C_DECLS

#endif
