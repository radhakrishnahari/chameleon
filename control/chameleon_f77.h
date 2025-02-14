/**
 *
 * @file chameleon_f77.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Fortran77 naming macros
 *
 * @version 1.3.0
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Philippe Virouleau
 * @date 2024-02-18
 *
 */
#ifndef _chameleon_f77_h_
#define _chameleon_f77_h_

#include "chameleon/mangling.h"

/**
 *  Determine FORTRAN names
 */
#define CHAMELEON_FNAME(lcname, UCNAME) CHAMELEON_GLOBAL(chameleon_##lcname, CHAMELEON_##UCNAME)
#define CHAMELEON_TILE_FNAME(lcname, UCNAME) CHAMELEON_GLOBAL(chameleon_##lcname##_tile, CHAMELEON_##UCNAME##_TILE)
#define CHAMELEON_ASYNC_FNAME(lcname, UCNAME) CHAMELEON_GLOBAL(chameleon_##lcname##_tile_async, CHAMELEON_##UCNAME##_TILE_ASYNC)
#define CHAMELEON_WS_FNAME(lcname, UCNAME) CHAMELEON_GLOBAL(chameleon_alloc_workspace_##lcname, CHAMELEON_ALLOC_WORKSPACE_##UCNAME)
#define CHAMELEON_WST_FNAME(lcname, UCNAME) CHAMELEON_GLOBAL(chameleon_alloc_workspace_##lcname##_tile, CHAMELEON_ALLOC_WORKSPACE_##UCNAME##_TILE)

#define __CHAMELEON_INIT CHAMELEON_GLOBAL(__chameleon_init, __CHAMELEON_INIT)
#define __CHAMELEON_FINALIZE CHAMELEON_GLOBAL(__chameleon_finalize, __CHAMELEON_FINALIZE)
#define CHAMELEON_ENABLE CHAMELEON_GLOBAL(chameleon_enable, CHAMELEON_ENABLE)
#define CHAMELEON_DISABLE CHAMELEON_GLOBAL(chameleon_disable, CHAMELEON_DISABLE)
#define CHAMELEON_SET CHAMELEON_GLOBAL(chameleon_set, CHAMELEON_SET)
#define CHAMELEON_GET CHAMELEON_GLOBAL(chameleon_get, CHAMELEON_GET)
#define CHAMELEON_DEALLOC_HANDLE CHAMELEON_GLOBAL(chameleon_dealloc_handle, CHAMELEON_DEALLOC_HANDLE)
#define CHAMELEON_VERSION CHAMELEON_GLOBAL(chameleon_version, CHAMELEON_VERSION)
#define CHAMELEON_DESC_CREATE CHAMELEON_GLOBAL(chameleon_desc_create, CHAMELEON_DESC_CREATE)
#define CHAMELEON_DESC_CREATE_OOC CHAMELEON_GLOBAL(chameleon_desc_create_ooc, CHAMELEON_DESC_CREATE_OOC)
#define CHAMELEON_DESC_CREATE_USER CHAMELEON_GLOBAL(chameleon_desc_create_user, CHAMELEON_DESC_CREATE_USER)
#define CHAMELEON_DESC_CREATE_OOC_USER CHAMELEON_GLOBAL(chameleon_desc_create_ooc_user, CHAMELEON_DESC_CREATE_OOC_USER)
#define CHAMELEON_DESC_DESTROY CHAMELEON_GLOBAL(chameleon_desc_destroy, CHAMELEON_DESC_DESTROY)
#define CHAMELEON_LAPACK_TO_TILE CHAMELEON_GLOBAL(chameleon_lapack_to_tile, CHAMELEON_LAPACK_TO_TILE)
#define CHAMELEON_TILE_TO_LAPACK CHAMELEON_GLOBAL(chameleon_tile_to_lapack, CHAMELEON_TILE_TO_LAPACK)

#if defined(CHAMELEON_SCHED_OPENMP)
#define CHAMELEON_INIT(nworkers, ncudas)\
    CALL __CHAMELEON_INIT(nworkers, ncudas)\
    !$omp parallel\
    !$omp master
#define CHAMELEON_FINALIZE()\
    !$omp end master\
    !$omp end parallel\
    CALL __CHAMELEON_FINALIZE()
#else
#define CHAMELEON_INIT(nworkers, ncudas)\
    CALL __CHAMELEON_INIT(nworkers, ncudas)
#define CHAMELEON_FINALIZE()\
    CALL __CHAMELEON_FINALIZE()
#endif

#endif /* _chameleon_f77_h_ */
