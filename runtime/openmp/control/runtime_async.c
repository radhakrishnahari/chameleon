/**
 *
 * @file openmp/runtime_async.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon OpenMP asynchronous routines
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Philippe Virouleau
 * @date 2024-02-18
 *
 */
#include "chameleon_openmp.h"

/**
 *  Create a sequence
 */
int RUNTIME_sequence_create( CHAM_context_t  *chamctxt,
                             RUNTIME_sequence_t *sequence )
{
    (void)chamctxt;
    (void)sequence;
    return CHAMELEON_SUCCESS;
}

/**
 *  Destroy a sequence
 */
int RUNTIME_sequence_destroy( CHAM_context_t  *chamctxt,
                              RUNTIME_sequence_t *sequence )
{
    (void)chamctxt;
    (void)sequence;
    return CHAMELEON_SUCCESS;
}

/**
 *  Wait for the completion of a sequence
 */
int RUNTIME_sequence_wait( CHAM_context_t  *chamctxt,
                           RUNTIME_sequence_t *sequence )
{
    (void)chamctxt;
    (void)sequence;

#pragma omp taskwait
    return CHAMELEON_SUCCESS;
}

/**
 *  Terminate a sequence
 */
void RUNTIME_sequence_flush( CHAM_context_t  *chamctxt,
                             RUNTIME_sequence_t *sequence,
                             RUNTIME_request_t  *request,
                             int status )
{
    (void)chamctxt;
    assert( status == CHAMELEON_SUCCESS );
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    return;
}

/**
 *  Create a request
 */
int RUNTIME_request_create( CHAM_context_t  *chamctxt,
                            RUNTIME_request_t *request )
{
    (void)chamctxt;
    request->schedopt = NULL;
    request->status = CHAMELEON_SUCCESS;
    return CHAMELEON_SUCCESS;
}

/**
 *  Destroy a request
 */
int RUNTIME_request_destroy( CHAM_context_t  *chamctxt,
                             RUNTIME_request_t *request )
{
    (void)chamctxt;
    (void)request;
    return CHAMELEON_SUCCESS;
}

/**
 *  Set runtime parameter for a request
 */
int RUNTIME_request_set( CHAM_context_t  *chamctxt,
                         RUNTIME_request_t *request,
                         int param, int value )
{
    (void)chamctxt;
    (void)request;
    (void)param;
    (void)value;
    return CHAMELEON_SUCCESS;
}
