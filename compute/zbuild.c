/**
 *
 * @file zbuild.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @brief Chameleon zbuild wrappers
 *
 * @version 1.3.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Guillaume Sylvand
 * @author Florent Pruvost
 * @date 2025-10-16
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zbuild - Generate a matrix by calling user provided function.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *            = ChamUpperLower: All the matrix A
 *            = ChamUpper: Upper triangular part
 *            = ChamLower: Lower triangular part
 *
 * @param[in] M
 *          The number of rows of A.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[out] A
 *          On exit, The matrix A generated.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] user_data
 *          The user data used in the matrix generation, it will be passed by chameleon
 *          to the build callback function (see below).
 *
 * @param[in] user_build_callback
 *          The user function to call to generate tiles.
 *          The prototype of the callback is :
 *          void myFcn(int row_min, int row_max, int col_min, int col_max, void *buffer, int ld, void *user_data)
 *          It is expected to build the block of matrix [row_min, row_max] x [col_min, col_max]
 *          (with both min and max values included in the intervals,
 *          index start at 0 like in C, NOT 1 like in Fortran)
 *          and store it at the adresse 'buffer' with leading dimension 'ld'
 *          The argument 'user_data' is an opaque pointer on any user data, it is passed by
 *          the user to Cham_zbuild (see above) and transmitted by chameleon to the callback.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zbuild_Tile
 * @sa CHAMELEON_zbuild_Tile_Async
 * @sa CHAMELEON_cbuild
 * @sa CHAMELEON_dbuild
 * @sa CHAMELEON_sbuild
 *
 */
int CHAMELEON_zbuild( cham_uplo_t uplo, int M, int N,
                  CHAMELEON_Complex64_t *A, int LDA,
                  void *user_data, void* user_build_callback )
{
    int NB;
    int status;
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t descAl, descAt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zbuild", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        chameleon_error("CHAMELEON_zbuild", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        chameleon_error("CHAMELEON_zbuild", "illegal value of N");
        return -2;
    }
    if (LDA < chameleon_max(1, M)) {
        chameleon_error("CHAMELEON_zbuild", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (chameleon_min(M, N) == 0)
        return CHAMELEON_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, M, N, 0);
    if (status != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zbuild", "chameleon_tune() failed");
        return status;
    }

    /* Set NT */
    NB = CHAMELEON_NB;
    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, "A", &descAl, &descAt, ChamDescOutput, uplo,
                     A, NB, NB, LDA, N, M, N, sequence, &request );

    /* Call the tile interface */
    chameleon_pzbuild( uplo, &descAt, user_data, user_build_callback, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                     ChamDescOutput, uplo, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );

    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  CHAMELEON_zbuild_Tile - Generate a matrix by tiles by calling user provided function.
 *  Tile equivalent of CHAMELEON_zbuild().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *            = ChamUpperLower: All the matrix A
 *            = ChamUpper: Upper triangular part
 *            = ChamLower: Lower triangular part
 *
 * @param[in] A
 *          On exit, The matrix A generated.
 *
 * @param[in] user_data
 *          The data used in the matrix generation.
 *
 * @param[in] user_build_callback
 *          The function called by the codelet to fill the tiles (see CHAMELEON_zbuild)
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zbuild
 * @sa CHAMELEON_zbuild_Tile_Async
 * @sa CHAMELEON_cbuild_Tile
 * @sa CHAMELEON_dbuild_Tile
 * @sa CHAMELEON_sbuild_Tile
 *
 */
int CHAMELEON_zbuild_Tile( cham_uplo_t uplo, CHAM_desc_t *A,
                       void *user_data, void* user_build_callback )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zbuild_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    chameleon_pzbuild( uplo, A, user_data, user_build_callback, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 *  CHAMELEON_zbuild_Tile_Async - Generate a matrix by tiles by calling user provided function.
 *  Non-blocking equivalent of CHAMELEON_zbuild_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *            = ChamUpperLower: All the matrix A
 *            = ChamUpper: Upper triangular part
 *            = ChamLower: Lower triangular part
 *
 * @param[in] A
 *          On exit, The matrix A generated.
 *
 * @param[in] user_data
 *          The data used in the matrix generation.
 *
 * @param[in] user_build_callback
 *          The function called by the codelet to fill the tiles (see CHAMELEON_zbuild)
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zbuild
 * @sa CHAMELEON_zbuild_Tile
 * @sa CHAMELEON_cbuild_Tile_Async
 * @sa CHAMELEON_dbuild_Tile_Async
 * @sa CHAMELEON_sbuild_Tile_Async
 *
 */
int CHAMELEON_zbuild_Tile_Async( cham_uplo_t uplo, CHAM_desc_t     *A,
                             void *user_data, void* user_build_callback,
                             RUNTIME_sequence_t *sequence,
                             RUNTIME_request_t  *request )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zbuild_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        chameleon_fatal_error("CHAMELEON_zbuild_Tile", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        chameleon_fatal_error("CHAMELEON_zbuild_Tile", "NULL request");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == CHAMELEON_SUCCESS) {
        request->status = CHAMELEON_SUCCESS;
    }
    else {
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED);
    }

    /* Check descriptors for correctness */
    if (chameleon_desc_check(A) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zbuild_Tile", "invalid descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (chameleon_min( A->m, A->n ) == 0) {
        return CHAMELEON_SUCCESS;
    }

    chameleon_pzbuild( uplo, A, user_data, user_build_callback, sequence, request );

    return CHAMELEON_SUCCESS;
}
