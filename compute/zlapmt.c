/**
 *
 * @file zlapmt.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlapmt wrappers
 *
 * @version 1.3.0
 * @author Matteo Marcos
 * @date 2025-04-28
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  @brief Computes the permutation P*op(A) or op(A)*P where P is the permutation
 *         matrix generated from PERM.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the permutation is done on the rows or the columns.
 *          = ChamLeft:  op(A) = A
 *          = ChamRight: op(A) = A^T
 *
 * @param[in] dir
 *          Specifies the order of the permutation.
 *          = ChamDirForward:  Natural order. P*op(A)
 *          = ChamDirBackward: Reverse order. op(A)*P
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          The M-by-N matrix A.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] PERM
 *          Vector of pivot indices.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlapmt_Tile
 * @sa CHAMELEON_zlapmt_Tile_Async
 * @sa CHAMELEON_clapmt
 * @sa CHAMELEON_dlapmt
 * @sa CHAMELEON_slapmt
 *
 */
int CHAMELEON_zlapmt( cham_side_t            side,
                      cham_dir_t             dir,
                      int                    M,
                      int                    N,
                      CHAMELEON_Complex64_t *A,
                      int                    LDA,
                      int                   *PERM )
{
    int                 status;
    int                 NB;
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request  = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t         descAl, descAt;
    CHAM_ipiv_t        *descIPIV;
    int                 K = ( side == ChamLeft ) ? M : N;
    int                 P, Q;
    void               *ws;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error("CHAMELEON_zlapmt", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( M < 0 ) {
        chameleon_error("CHAMELEON_zlapmt", "illegal value of M");
        return -2;
    }
    if ( N < 0 ) {
        chameleon_error("CHAMELEON_zlapmt", "illegal value of N");
        return -3;
    }
    if ( LDA < chameleon_max( 1, M ) ) {
        chameleon_error("CHAMELEON_zlapmt", "illegal value of LDA");
        return -5;
    }

    /* Quick return */
    if ( chameleon_min( N, M ) == 0 ) {
        return (double)0.0;
    }

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, M, N, 0);
    if ( status != CHAMELEON_SUCCESS ) {
        chameleon_error("CHAMELEON_zlapmt", "chameleon_tune() failed");
        return status;
    }

    /* Set NB */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, &descAl, &descAt, ChamDescInput, ChamUpperLower,
                         A, NB, NB, LDA, N, M, N, sequence, &request );

    P = chameleon_desc_datadist_get_iparam( &descAt, 0 );
    Q = chameleon_desc_datadist_get_iparam( &descAt, 1 );

    CHAMELEON_Ipiv_Create( &descIPIV, side, descAt.mb, K, P, P*Q, NULL );

    CHAMELEON_Perm_Init( dir, descIPIV, PERM );

    /* Call the tile interface */
    ws = CHAMELEON_zlaswp_WS_Alloc( side, &descAt );
    CHAMELEON_zlapmt_Tile_Async( side, dir, &descAt, descIPIV, ws, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInput, ChamUpperLower, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );
    CHAMELEON_zlaswp_WS_Free( ws );

    /* Cleanup the temporary data */
    CHAMELEON_Ipiv_Destroy( &descIPIV );
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );

    chameleon_sequence_destroy( chamctxt, sequence );
    return CHAMELEON_SUCCESS;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  @brief Tile equivalent of CHAMELEON_zlapmt().
 *
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the permutation is done on the rows or the columns.
 *          = ChamLeft:  op(A) = A
 *          = ChamRight: op(A) = A^T
 *
 * @param[in] dir
 *          Specifies the order of the permutation.
 *          = ChamDirForward:  Natural order. P*op(A)
 *          = ChamDirBackward: Reverse order. op(A)*P
 *
 * @param[in,out] A
 *          The M-by-N matrix A.
 *
 * @param[in] IPIV
 *          Vector of pivot indices created with CHAMELEON_Ipiv_Create() and
 *          initialized with CHAMELEON_Perm_Init().
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlapmt
 * @sa CHAMELEON_zlapmt_Tile_Async
 * @sa CHAMELEON_clapmt_Tile
 * @sa CHAMELEON_dlapmt_Tile
 * @sa CHAMELEON_slapmt_Tile
 *
 */
int CHAMELEON_zlapmt_Tile( cham_side_t  side,
                           cham_dir_t   dir,
                           CHAM_desc_t *A,
                           CHAM_ipiv_t *IPIV )
{
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request  = RUNTIME_REQUEST_INITIALIZER;
    int                 status;
    void               *ws;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error("CHAMELEON_zlapmt_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    ws = CHAMELEON_zlaswp_WS_Alloc( side, A );
    CHAMELEON_zlapmt_Tile_Async( side, dir, A, IPIV, ws, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    CHAMELEON_zlaswp_WS_Free( ws );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 *  @brief Non-blocking equivalent of CHAMELEON_zlapmt_Tile().
 *
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether the permutation is done on the rows or the columns.
 *          = ChamLeft:  op(A) = A
 *          = ChamRight: op(A) = A^T
 *
 * @param[in] dir
 *          Specifies the order of the permutation.
 *          = ChamDirForward:  Natural order. P*op(A)
 *          = ChamDirBackward: Reverse order. op(A)*P
 *
 * @param[in,out] A
 *          The M-by-N matrix A.
 *
 * @param[in] IPIV
 *          Vector of pivot indices created with CHAMELEON_Ipiv_Create() and
 *          initialized with CHAMELEON_Perm_Init().
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
 * @sa CHAMELEON_zlapmt
 * @sa CHAMELEON_zlapmt_Tile
 * @sa CHAMELEON_clapmt_Tile_Async
 * @sa CHAMELEON_dlapmt_Tile_Async
 * @sa CHAMELEON_slapmt_Tile_Async
 *
 */
int CHAMELEON_zlapmt_Tile_Async( cham_side_t         side,
                                 cham_dir_t          dir,
                                 CHAM_desc_t        *A,
                                 CHAM_ipiv_t        *IPIV,
                                 void               *user_ws,
                                 RUNTIME_sequence_t *sequence,
                                 RUNTIME_request_t  *request )
{
    CHAM_context_t             *chamctxt;
    struct chameleon_pzlaswp_s *ws;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error("CHAMELEON_zlapmt_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if ( sequence == NULL ) {
        chameleon_fatal_error("CHAMELEON_zlapmt_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if ( request == NULL ) {
        chameleon_fatal_error("CHAMELEON_zlapmt_Tile_Async", "NULL request");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if ( sequence->status == CHAMELEON_SUCCESS ) {
        request->status = CHAMELEON_SUCCESS;
    }
    else {
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED);
    }

    /* Check descriptors for correctness */
    if ( chameleon_desc_check(A) != CHAMELEON_SUCCESS ) {
        chameleon_error("CHAMELEON_zlapmt_Tile_Async", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ( A->mb != A->nb ) {
        chameleon_error("CHAMELEON_zlapmt_Tile_Async", "only matching tile sizes supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if ( chameleon_min( A->m, A->n ) == 0 ) {
        return CHAMELEON_SUCCESS;
    }

    if ( user_ws == NULL ) {
        ws = CHAMELEON_zlaswp_WS_Alloc( side, A );
    }
    else {
        ws = user_ws;
    }

    if ( side == ChamLeft ) {
        chameleon_pzlaswp( ws, ChamDirForward, A, IPIV, sequence, request );
    }
    else {
        chameleon_pzlaswpc( ws, ChamDirForward, A, IPIV, sequence, request );
    }

    if ( user_ws == NULL ) {
        CHAMELEON_Desc_Flush( A, sequence );
        chameleon_sequence_wait( chamctxt, sequence );
        CHAMELEON_zlaswp_WS_Free( ws );
    }

    return CHAMELEON_SUCCESS;
}

