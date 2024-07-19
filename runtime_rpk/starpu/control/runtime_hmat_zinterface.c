/**
 *
 * @file runtime_hmat_zinterface.c
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
 * @date 2024-05-19
 * @precisions normal z -> c d s
 *
 */

#include "runtime_hmat.h"

hmat_interface_t *
runtime_hmat_zinterface_get( void ) {
	return runtime_hmat_interface_get( HMAT_DOUBLE_COMPLEX );
}
