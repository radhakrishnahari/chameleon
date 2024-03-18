/**
 *
 * @file starpu/chameleon_starpu_internal.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU runtime header
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Philippe Swartvagher
 * @author Samuel Thibault
 * @author Loris Lucido
 * @author Terry Cojean
 * @author Matthieu Kuhn
 * @date 2024-03-16
 *
 */
#ifndef _chameleon_starpu_internal_h_
#define _chameleon_starpu_internal_h_

#include "control/common.h"
#include "chameleon_starpu.h"

/* Chameleon interfaces for StarPU */
#include "cham_tile_interface.h"
#include "cppi_interface.h"

/**
 * @brief StarPU specific Chameleon structure stored in the schedopt field
 */
typedef struct CHAM_context_starpu_s
{
    struct starpu_conf                    starpu_conf; /**< StarPU main configuration structure   */
    struct starpu_parallel_worker_config *pw_config;   /**< StarPU parallel workers configuration */
} CHAM_context_starpu_t;

/**
 * @brief StarPU specific request field stored in the schedopt field
 */
typedef struct RUNTIME_request_starpu_s {
    int workerid; // to force task execution on a specific workerid
} RUNTIME_request_starpu_t;

/**
 * @brief Convert the Chameleon access enum to the StarPU one
 */
static inline int cham_to_starpu_access( cham_access_t accessA ) {
    assert( (enum starpu_data_access_mode)ChamR  == STARPU_R  );
    assert( (enum starpu_data_access_mode)ChamW  == STARPU_W  );
    assert( (enum starpu_data_access_mode)ChamRW == STARPU_RW );
    return accessA;
}

#endif /* _chameleon_starpu_internal_h_ */
