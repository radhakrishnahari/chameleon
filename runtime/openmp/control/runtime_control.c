/**
 *
 * @file openmp/runtime_control.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon OpenMP control routines
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Philippe Virouleau
 * @author Philippe Swartvagher
 * @author Matthieu Kuhn
 * @date 2024-02-18
 *
 */
#include "chameleon_openmp.h"

/**
 *
 */
int RUNTIME_init( CHAM_context_t *chamctxt,
                  int ncpus,
                  int ncudas,
                  int nthreads_per_worker )
{
    if ( ncudas > 0 ) {
        chameleon_warning( "RUNTIME_init_scheduler(OpenMP)", "GPUs are not supported for now");
    }

    if ( nthreads_per_worker > 0 ) {
        chameleon_warning( "RUNTIME_init_scheduler(OpenMP)", "Multi-threaded kernels are not supported for now");
    }

    chamctxt->nworkers = ncpus;
    chamctxt->nthreads_per_worker = nthreads_per_worker;

    return CHAMELEON_SUCCESS;
}

/**
 *
 */
void RUNTIME_finalize( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return;
}

/**
 *  To suspend the processing of new tasks by workers
 */
void RUNTIME_pause( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return;
}

/**
 *  This is the symmetrical call to RUNTIME_pause,
 *  used to resume the workers polling for new tasks.
 */
void RUNTIME_resume( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return;
}

/**
 *  Busy-waiting barrier
 */
void RUNTIME_barrier( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
#pragma omp barrier
}

/**
 *  Display a progress information when executing the tasks
 */
void RUNTIME_progress( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return;
}


/**
 * Thread rank.
 */
int RUNTIME_thread_rank( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return omp_get_thread_num();
}

/**
 * Number of threads.
 */
int RUNTIME_thread_size( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return omp_get_num_threads();
}

/**
 *  The process rank
 */
int RUNTIME_comm_rank( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return 0;
}

/**
 *  This returns the size of the distributed computation
 */
int RUNTIME_comm_size( CHAM_context_t *chamctxt )
{
    (void)chamctxt;
    return 1;
}

void RUNTIME_set_minmax_submitted_tasks( int min, int max ) {
    (void)min;
    (void)max;
}
