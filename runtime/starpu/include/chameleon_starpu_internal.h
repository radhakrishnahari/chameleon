/**
 *
 * @file starpu/chameleon_starpu_internal.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
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
 * @date 2024-09-17
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

#include "runtime_codelets.h"
#include "runtime_profiling.h"
#include "runtime_codelet_profile.h"
#include "runtime_workspace.h"

void *RUNTIME_data_getaddr_withconversion( const RUNTIME_option_t *options,
                                           cham_access_t access, cham_flttype_t flttype,
                                           const CHAM_desc_t *A, int m, int n );

/*
 * MPI Redefinitions
 */
#if defined(CHAMELEON_USE_MPI)

#if defined(CHAMELEON_RUNTIME_SYNC)
#define rt_starpu_insert_task( _codelet_, ... )                         \
    starpu_mpi_insert_task( options->sequence->comm, (_codelet_), STARPU_TASK_SYNCHRONOUS, 1, ##__VA_ARGS__ )
#else
#define rt_starpu_insert_task( _codelet_, ... )                         \
    starpu_mpi_insert_task( options->sequence->comm, (_codelet_), ##__VA_ARGS__ )
#endif

#else

#if defined(CHAMELEON_RUNTIME_SYNC)
#define rt_starpu_insert_task( _codelet_, ... )                         \
    starpu_insert_task( (_codelet_), STARPU_TASK_SYNCHRONOUS, 1, ##__VA_ARGS__ )
#else
#define rt_starpu_insert_task( _codelet_, ... )                         \
    starpu_insert_task( (_codelet_), ##__VA_ARGS__ )
#endif

#endif

#if defined(CHAMELEON_RUNTIME_SYNC)
#define rt_shm_starpu_insert_task( _codelet_, ... )                         \
    starpu_insert_task( (_codelet_), STARPU_TASK_SYNCHRONOUS, 1, ##__VA_ARGS__ )
#else
#define rt_shm_starpu_insert_task( _codelet_, ... )                         \
    starpu_insert_task( (_codelet_), ##__VA_ARGS__ )
#endif

/*
 * Enable codelets names
 */
#if (STARPU_MAJOR_VERSION > 1) || ((STARPU_MAJOR_VERSION == 1) && (STARPU_MINOR_VERSION > 1))
#define CHAMELEON_CODELETS_HAVE_NAME
#endif

/**
 * MPI tag management
 */
void    chameleon_starpu_tag_init( );
int64_t chameleon_starpu_tag_book( int64_t nbtags );
void    chameleon_starpu_tag_release( int64_t min );

/**
 * Access to block pointer and leading dimension
 */
#define RTBLKADDR( desc, type, m, n ) ( (starpu_data_handle_t)RUNTIME_data_getaddr( desc, m, n ) )

void RUNTIME_set_reduction_methods(starpu_data_handle_t handle, cham_flttype_t dtyp);

#include "runtime_mpi.h"
#include "runtime_wontuse.h"

#if defined(CHAMELEON_USE_MPI) && defined(HAVE_STARPU_MPI_CACHED_RECEIVE)
static inline int
chameleon_starpu_data_iscached(const CHAM_desc_t *A, int m, int n)
{
    int64_t mm = m + (A->i / A->mb);
    int64_t nn = n + (A->j / A->nb);

    starpu_data_handle_t *ptrtile = A->schedopt;
    ptrtile += ((int64_t)A->lmt) * nn + mm;

    if (!(*ptrtile)) {
        return 0;
    }

    return starpu_mpi_cached_receive(*ptrtile);
}

#define RUNTIME_ACCESS_WRITE_CACHED(A, Am, An) do {                 \
        if (chameleon_starpu_data_iscached(A, Am, An)) __chameleon_need_submit = 1; } while(0)

#else

#if defined(CHAMELEON_USE_MPI)
#warning "WAR dependencies need starpu_mpi_cached_receive support from StarPU 1.2.1 or greater"
#endif
#define RUNTIME_ACCESS_WRITE_CACHED(A, Am, An) do {} while (0)

#endif

#ifdef CHAMELEON_ENABLE_PRUNING_STATS

#define RUNTIME_PRUNING_STATS_BEGIN_ACCESS_DECLARATION \
    int __chameleon_exec = 0; \
    int __chameleon_changed = 0;

#define RUNTIME_PRUNING_STATS_ACCESS_W(A, Am, An) \
    if (chameleon_desc_islocal(A, Am, An)) \
        __chameleon_exec = 1;

#define RUNTIME_PRUNING_STATS_END_ACCESS_DECLARATION \
    RUNTIME_total_tasks++; \
    if (__chameleon_exec) \
        RUNTIME_exec_tasks++; \
    else if (__chameleon_need_submit) \
        RUNTIME_comm_tasks++; \
    else if (__chameleon_changed) \
        RUNTIME_changed_tasks++;

#define RUNTIME_PRUNING_STATS_RANK_CHANGED(rank) \
    int __chameleon_myrank; \
    RUNTIME_comm_rank(&__chameleon_myrank); \
    __chameleon_exec = (rank) == __chameleon_myrank; \
    __chameleon_changed = 1; \

#else
#define RUNTIME_PRUNING_STATS_BEGIN_ACCESS_DECLARATION
#define RUNTIME_PRUNING_STATS_ACCESS_W(A, Am, An)
#define RUNTIME_PRUNING_STATS_END_ACCESS_DECLARATION
#define RUNTIME_PRUNING_STATS_RANK_CHANGED(rank)
#endif

#define RUNTIME_BEGIN_ACCESS_DECLARATION        \
    RUNTIME_PRUNING_STATS_BEGIN_ACCESS_DECLARATION

#define RUNTIME_ACCESS_R(A, Am, An)

#define RUNTIME_ACCESS_W(A, Am, An)             \
    RUNTIME_PRUNING_STATS_ACCESS_W(A, Am, An);  \
    RUNTIME_ACCESS_WRITE_CACHED(A, Am, An)

#define RUNTIME_ACCESS_RW(A, Am, An)            \
    RUNTIME_PRUNING_STATS_ACCESS_W(A, Am, An);  \
    RUNTIME_ACCESS_WRITE_CACHED(A, Am, An)

#define RUNTIME_RANK_CHANGED(rank)              \
    RUNTIME_PRUNING_STATS_RANK_CHANGED(rank)

#define RUNTIME_END_ACCESS_DECLARATION          \
    RUNTIME_PRUNING_STATS_END_ACCESS_DECLARATION;

#endif /* _chameleon_starpu_internal_h_ */
