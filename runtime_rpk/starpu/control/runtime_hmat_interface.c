/**
 *
 * @file runtime_hmat_interface.c
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
 * @date 2024-05-20
 *
 */

#include "hmat/hmat.h"
#include <stdio.h>

static hmat_interface_t cham_hmat_zinterface;
static hmat_interface_t cham_hmat_cinterface;
static hmat_interface_t cham_hmat_dinterface;
static hmat_interface_t cham_hmat_sinterface;

void
runtime_hmat_interface_init( void ) {
    hmat_init_default_interface( &cham_hmat_zinterface, HMAT_DOUBLE_COMPLEX );
    hmat_init_default_interface( &cham_hmat_cinterface, HMAT_SIMPLE_COMPLEX );
    hmat_init_default_interface( &cham_hmat_dinterface, HMAT_DOUBLE_PRECISION );
    hmat_init_default_interface( &cham_hmat_sinterface, HMAT_SIMPLE_PRECISION );
}

hmat_interface_t *
runtime_hmat_interface_get( hmat_value_t ftype ) {
    switch( ftype ) {
    case HMAT_DOUBLE_COMPLEX:
        return &cham_hmat_zinterface;
    case HMAT_SIMPLE_COMPLEX:
        return &cham_hmat_cinterface;
    case HMAT_DOUBLE_PRECISION:
        return &cham_hmat_dinterface;
    case HMAT_SIMPLE_PRECISION:
        return &cham_hmat_sinterface;
    default:
        fprintf(stderr, "Error: unknown hmat value type\n");
        exit(2);
    }
}

void runtime_hmat_interface_exit( void ) {

}
