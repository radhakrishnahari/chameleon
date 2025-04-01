/**
 *
 * @file descriptor_ipiv.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon descriptors routines
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Florent Pruvost
 * @date 2025-03-24
 *
 ***
 *
 * @defgroup Descriptor
 * @brief Group descriptor routines exposed to users to manipulate IPIV data structures
 *
 */
#define _GNU_SOURCE 1
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "control/common.h"
#include "control/descriptor.h"
#include "chameleon/runtime.h"

/**
 ******************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Internal function to create tiled descriptor associated to a pivot array.
 *
 ******************************************************************************
 *
 * @param[in,out] ipiv
 *          The pointer to the ipiv descriptor to initialize.
 *
 * @param[in] desc
 *          The tile descriptor for which an associated ipiv descriptor must be generated.
 *
 * @param[in] m
 *          The size of the pivot array.
 *
 * @param[in] data
 *          The pointer to the original vector where to store the pivot values.
 *
 ******************************************************************************
 *
 * @return CHAMELEON_SUCCESS on success, CHAMELEON_ERR_NOT_INITIALIZED otherwise.
 *
 */
int chameleon_ipiv_init( CHAM_ipiv_t *ipiv, const CHAM_desc_t *desc, int m, void *data )
{
    CHAM_context_t *chamctxt;
    int rc = CHAMELEON_SUCCESS;

    memset( ipiv, 0, sizeof(CHAM_ipiv_t) );

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Desc_Create", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    ipiv->desc = desc;
    ipiv->data = data;
    ipiv->i    = 0;
    ipiv->m    = m;
    ipiv->mb   = desc->mb;
    ipiv->mt   = chameleon_ceil( ipiv->m, ipiv->mb );

    /* Create runtime specific structure like registering data */
    RUNTIME_ipiv_create( ipiv );

    return rc;
}

/**
 ******************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Internal function to create tiled descriptor associated to a pivot.
 *
 ******************************************************************************
 *
 * @param[in,out] pivot
 *          The pointer to the pivot descriptor to initialize.
 *
 * @param[in] desc
 *          The tile descriptor for which an associated pivot descriptor must be generated.
 *
 ******************************************************************************
 *
 * @return CHAMELEON_SUCCESS on success, CHAMELEON_ERR_NOT_INITIALIZED otherwise.
 *
 */
int chameleon_pivot_init( CHAM_desc_pivot_t *pivot, const CHAM_desc_t *desc )
{
    CHAM_context_t *chamctxt;
    int rc = CHAMELEON_SUCCESS;

    memset( pivot, 0, sizeof(CHAM_desc_pivot_t) );

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Desc_Create", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    pivot->P       = chameleon_desc_datadist_get_iparam( desc, 0 );
    pivot->Q       = chameleon_desc_datadist_get_iparam( desc, 1 );
    pivot->n       = chameleon_min(desc->mb, desc->nb);
    pivot->nb      = desc->mb;
    pivot->dtyp    = desc->dtyp;

    /* Create runtime specific structure like registering data */
    RUNTIME_pivot_create( pivot );

    return rc;
}

/**
 ******************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Internal function to destroy a tiled descriptor associated to a pivot array.
 *
 ******************************************************************************
 *
 * @param[in,out] ipiv
 *          The pointer to the ipiv descriptor to destroy.
 *
 */
void chameleon_ipiv_destroy( CHAM_ipiv_t       *ipiv,
                             const CHAM_desc_t *desc )
{
    RUNTIME_ipiv_destroy( ipiv, desc );
}

/**
 ******************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Internal function to destroy a tiled descriptor associated to a pivot array.
 *
 ******************************************************************************
 *
 * @param[in,out] pivot
 *          The pointer to the pivot descriptor to destroy.
 *
 * @param[in] desc
 *          The tile descriptor for which an associated pivot descriptor must be generated.
 *
 */
void chameleon_pivot_destroy( CHAM_desc_pivot_t *pivot )
{
    RUNTIME_pivot_destroy( pivot );
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Create a tiled ipiv descriptor associated to a given matrix.
 *
 ******************************************************************************
 *
 * @param[in,out] ipiv
 *          The pointer to the ipiv descriptor to initialize.
 *
 * @param[in] desc
 *          The tile descriptor for which an associated ipiv descriptor must be generated.
 *
 * @param[in] m
 *          The size of the pivot array.
 *
 * @param[in] data
 *          The pointer to the original vector where to store the pivot values.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS on successful exit
 * @retval CHAMELEON_ERR_NOT_INITIALIZED if failed to initialize the descriptor.
 * @retval CHAMELEON_ERR_OUT_OF_RESOURCES if failed to allocated some ressources.
 *
 */
int CHAMELEON_Ipiv_Create( CHAM_ipiv_t **ipivptr, const CHAM_desc_t *desc, int m, void *data )
{
    CHAM_context_t *chamctxt;
    CHAM_ipiv_t *ipiv;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Ipiv_Create", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    /* Allocate memory and initialize the ipivriptor */
    ipiv = (CHAM_ipiv_t*)malloc(sizeof(CHAM_ipiv_t));
    if (ipiv == NULL) {
        chameleon_error("CHAMELEON_Ipiv_Create", "malloc() failed");
        return CHAMELEON_ERR_OUT_OF_RESOURCES;
    }

    chameleon_ipiv_init( ipiv, desc, m, data );

    *ipivptr = ipiv;
    return CHAMELEON_SUCCESS;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  @brief initialize the IPIV descriptor.
 *
 *******************************************************************************
 *
 * @param[in] descA
 *          Descriptor of the matrix A.
 *
 * @param[in,out] descIPIV
 *          Descriptor of the pivot array. Should be initialized using
 *          CHAMELEON_Ipiv_Create() with data filled with the vector of pivot.
 *
 *******************************************************************************
 *
 *
 */
void CHAMELEON_Ipiv_Init( const CHAM_desc_t *descA,
                          CHAM_ipiv_t       *descIPIV )
{

    RUNTIME_option_t    options;
    RUNTIME_request_t   request  = RUNTIME_REQUEST_INITIALIZER;
    RUNTIME_sequence_t *sequence = NULL;
    CHAM_context_t     *chamctxt;

    chamctxt = chameleon_context_self();
    chameleon_sequence_create( chamctxt, &sequence );
    RUNTIME_options_init( &options, chamctxt, sequence, &request );

    INSERT_TASK_ipiv_init_data( &options, descIPIV );

    chameleon_sequence_wait( chamctxt, sequence );
    chameleon_sequence_destroy( chamctxt, sequence );
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Destroys an ipiv tile descriptor.
 *
 ******************************************************************************
 *
 * @param[in] ipivptr
 *          The Ipiv tile descriptor to destroy.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Ipiv_Destroy( CHAM_ipiv_t **ipivptr,
                            const CHAM_desc_t *desc )
{
    CHAM_context_t *chamctxt;
    CHAM_ipiv_t *ipiv;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Ipiv_Destroy", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    if ((ipivptr == NULL) || (*ipivptr == NULL)) {
        chameleon_error("CHAMELEON_Ipiv_Destroy", "attempting to destroy a NULL descriptor");
        return CHAMELEON_ERR_UNALLOCATED;
    }

    ipiv = *ipivptr;
    chameleon_ipiv_destroy( ipiv, desc );
    free(ipiv);
    *ipivptr = NULL;
    return CHAMELEON_SUCCESS;
}

 /**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Flushes the data in the sequence when they won't be reused. This calls
 * cleans up the distributed communication caches, and transfer the data back to
 * the CPU.
 *
 ******************************************************************************
 *
 * @param[in] ipiv
 *          ipiv descriptor.
 *
 * @param[in] sequence
 *          The sequence in which to submit the calls to flush the data.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Ipiv_Flush( const CHAM_ipiv_t        *ipiv,
                          const RUNTIME_sequence_t *sequence )
{
    RUNTIME_ipiv_flush( sequence, ipiv );
    return CHAMELEON_SUCCESS;
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 * @brief Gathers an IPIV tile descriptor in a single vector on the given root node.
 *
 ******************************************************************************
 *
 * @param[in] ipivdesc
 *          the ipiv vector descriptor to gather.
 *
 * @param[in] ipiv
 *          The ipiv vector where to store the result. Allocated vector of size
 *          ipivdesc->m on root, not referenced on other nodes.
 *
 * @param[in] root
 *          root node on which to gather the data.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Ipiv_Gather( CHAM_ipiv_t *ipivdesc, int *ipiv, int root )
{
    CHAM_context_t     *chamctxt = chameleon_context_self();
    RUNTIME_sequence_t *sequence = NULL;
    int                 status;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the tasks to collect the ipiv array on root */
    RUNTIME_ipiv_gather( sequence, ipivdesc, ipiv, root );

    /* Wait for the end */
    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}
