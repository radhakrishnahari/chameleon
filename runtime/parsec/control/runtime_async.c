/**
 *
 * @file parsec/runtime_async.c
 *
 * @copyright 2012-2017 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon PaRSEC asynchronous routines
 *
 * @version 1.2.0
 * @author Reazul Hoque
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2022-02-22
 *
 */
#include "chameleon_parsec.h"

/**
 *  Create a sequence
 */
int RUNTIME_sequence_create( CHAM_context_t     *chamctxt,
                             RUNTIME_sequence_t *sequence )
{
    parsec_context_t  *parsec        = (parsec_context_t *)(chamctxt->schedopt);
    parsec_taskpool_t *parsec_dtd_tp = parsec_dtd_taskpool_new();

    parsec_context_add_taskpool( parsec, (parsec_taskpool_t *)parsec_dtd_tp );
    sequence->schedopt = parsec_dtd_tp;

    parsec_context_start(parsec);

    return CHAMELEON_SUCCESS;
}

/**
 *  Destroy a sequence
 */
int RUNTIME_sequence_destroy( CHAM_context_t     *chamctxt,
                              RUNTIME_sequence_t *sequence )
{
    parsec_taskpool_t *parsec_dtd_tp = (parsec_taskpool_t *)(sequence->schedopt);

    assert( parsec_dtd_tp );
    parsec_taskpool_free( parsec_dtd_tp );

    sequence->schedopt = NULL;

    (void)chamctxt;
    return CHAMELEON_SUCCESS;
}

/**
 *  Wait for the completion of a sequence
 */
int RUNTIME_sequence_wait( CHAM_context_t  *chamctxt,
                           RUNTIME_sequence_t *sequence )
{
    parsec_context_t  *parsec = (parsec_context_t *)chamctxt->schedopt;
    parsec_taskpool_t *parsec_dtd_tp = (parsec_taskpool_t *) sequence->schedopt;

    assert( parsec_dtd_tp );
    parsec_taskpool_wait( parsec_dtd_tp );

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
    assert( status == CHAMELEON_SUCCESS );
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    (void)chamctxt;
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
