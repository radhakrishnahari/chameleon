/**
 *
 * @file zlaswp.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaswp wrappers
 *
 * @version 1.3.0
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-03-24
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
 *         matrix generated from IPIV.
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
 * @param[in] K1
 *          The first element of IPIV for which an interchange will
 *          be done. Must follow the Fortran numbering standard.
 *
 * @param[in] K2
 *          The last element of ipiv for which an interchange will
 *          be done. Must follow the Fortran numbering standard.
 *
 * @param[in] IPIV
 *          Vector of pivot indices.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlaswp_Tile
 * @sa CHAMELEON_zlaswp_Tile_Async
 * @sa CHAMELEON_claswp
 * @sa CHAMELEON_dlaswp
 * @sa CHAMELEON_slaswp
 *
 */
int CHAMELEON_zlaswp( cham_side_t            side,
                      cham_dir_t             dir,
                      int                    M,
                      int                    N,
                      CHAMELEON_Complex64_t *A,
                      int                    LDA,
                      int                    K1,
                      int                    K2,
                      int                   *IPIV )
{
    int                 status;
    int                 NB;
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request  = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t         descAl, descAt;
    CHAM_ipiv_t        *descIPIV;
    int                 K = ( side == ChamLeft ) ? M : N;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error("CHAMELEON_zlaswp", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( M < 0 ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of M");
        return -2;
    }
    if ( N < 0 ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of N");
        return -3;
    }
    if ( LDA < chameleon_max( 1, M ) ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of LDA");
        return -5;
    }
    if ( ( K1 < 1 ) || ( K1 > K ) ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of K1");
        return -6;
    }
    if ( ( K2 < 1 ) || ( K2 > K ) ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of K2");
        return -7;
    }

    /* Quick return */
    if ( chameleon_min( N, M ) == 0 ) {
        return (double)0.0;
    }

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = chameleon_tune(CHAMELEON_FUNC_ZGEMM, M, N, 0);
    if ( status != CHAMELEON_SUCCESS ) {
        chameleon_error("CHAMELEON_zlaswp", "chameleon_tune() failed");
        return status;
    }

    /* Set NB */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, &descAl, &descAt, ChamDescInput, ChamUpperLower,
                         A, NB, NB, LDA, N, M, N, sequence, &request );
    CHAMELEON_Ipiv_Create( &descIPIV, &descAt, K, IPIV );

    CHAMELEON_Ipiv_Init( &descAt, descIPIV );

    /* Call the tile interface */
    CHAMELEON_zlaswp_Tile_Async( side, dir, &descAt, K1, K2, descIPIV, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInput, ChamUpperLower, sequence, &request );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    CHAMELEON_Ipiv_Destroy( &descIPIV, &descAt );
    chameleon_ztile2lap_cleanup( chamctxt, &descAl, &descAt );

    chameleon_sequence_destroy( chamctxt, sequence );
    return CHAMELEON_SUCCESS;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile
 *
 *  @brief Tile equivalent of CHAMELEON_zlaswp().
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
 * @param[in] K1
 *          The first element of IPIV for which an interchange will
 *          be done. Must follow the Fortran numbering standard.
 *
 * @param[in] K2
 *          The last element of ipiv for which an interchange will
 *          be done. Must follow the Fortran numbering standard.
 *
 * @param[in] IPIV
 *          Vector of pivot indices.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zlaswp
 * @sa CHAMELEON_zlaswp_Tile_Async
 * @sa CHAMELEON_claswp_Tile
 * @sa CHAMELEON_dlaswp_Tile
 * @sa CHAMELEON_slaswp_Tile
 *
 */
int CHAMELEON_zlaswp_Tile( cham_side_t  side,
                           cham_dir_t   dir,
                           CHAM_desc_t *A,
                           int          K1,
                           int          K2,
                           CHAM_ipiv_t *IPIV )
{
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request  = RUNTIME_REQUEST_INITIALIZER;
    int                 status;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error("CHAMELEON_zlaswp_Tile", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if ( ( K1 < 1 ) || ( K1 > A->m ) ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of K1");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    if ( ( K2 < 1 ) || ( K2 > A->m ) ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of K2");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    CHAMELEON_zlaswp_Tile_Async( side, dir, A, K1, K2, IPIV, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Ipiv_Flush( IPIV, sequence );

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
 *  @brief Non-blocking equivalent of CHAMELEON_zlaswp_Tile().
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
 * @param[in] K1
 *          The first element of IPIV for which an interchange will
 *          be done. Must follow the Fortran numbering standard
 *
 * @param[in] K2
 *          The last element of ipiv for which an interchange will
 *          be done. Must follow the Fortran numbering standard.
 *
 * @param[in] IPIV
 *          Vector of pivot indices.
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
 * @sa CHAMELEON_zlaswp
 * @sa CHAMELEON_zlaswp_Tile
 * @sa CHAMELEON_claswp_Tile_Async
 * @sa CHAMELEON_dlaswp_Tile_Async
 * @sa CHAMELEON_slaswp_Tile_Async
 *
 */
int CHAMELEON_zlaswp_Tile_Async( cham_side_t         side,
                                 cham_dir_t          dir,
                                 CHAM_desc_t        *A,
                                 int                 K1,
                                 int                 K2,
                                 CHAM_ipiv_t        *IPIV,
                                 RUNTIME_sequence_t *sequence,
                                 RUNTIME_request_t  *request )
{
    CHAM_context_t             *chamctxt;
    struct chameleon_pzgetrf_s *ws;
    RUNTIME_option_t            options;
    int                         k;
    int                         K = ( side == ChamLeft ) ? A->m : A->n;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error("CHAMELEON_zlaswp_Tile_Async", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if ( ( K1 < 1 ) || ( K1 > K ) ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of K1");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    if ( ( K2 < 1 ) || ( K2 > K ) ) {
        chameleon_error("CHAMELEON_zlaswp", "illegal value of K2");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
    if ( sequence == NULL ) {
        chameleon_fatal_error("CHAMELEON_zlaswp_Tile_Async", "NULL sequence");
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if ( request == NULL ) {
        chameleon_fatal_error("CHAMELEON_zlaswp_Tile_Async", "NULL request");
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
        chameleon_error("CHAMELEON_zlaswp_Tile_Async", "invalid first descriptor");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ( A->mb != A->nb ) {
        chameleon_error("CHAMELEON_zlaswp_Tile_Async", "only matching tile sizes supported");
        return chameleon_request_fail(sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if ( chameleon_min( A->m, A->n ) == 0 ) {
        return CHAMELEON_SUCCESS;
    }

    if ( IPIV->data != NULL ) {
        RUNTIME_options_init( &options, chamctxt, sequence, request );
        if ( side == ChamLeft ) {
            int tempkm, m0;

            for ( k = 0; k < A->mt; k++ ) {
                tempkm = A->get_blkdim( A, k, DIM_m, A->m );
                m0 = k * A->mb;
                INSERT_TASK_ipiv_to_perm( &options, m0, tempkm, tempkm, K1 - 1, K2 - 1,
                                               IPIV, k );
            }
        }
        else {
            int tempkn, n0;

            for ( k = 0; k < A->nt; k++ ) {
                tempkn = A->get_blkdim( A, k, DIM_n, A->n );
                n0 = k * A->nb;
                INSERT_TASK_ipiv_to_perm( &options, n0, tempkn, tempkn, K1 - 1, K2 - 1,
                                               IPIV, k );
            }
        }
        chameleon_sequence_wait( chamctxt, sequence );
    }

    ws = CHAMELEON_zgetrf_WS_Alloc( A );

    if ( side == ChamLeft ) {
        chameleon_pzlaswp( ws, dir, A, IPIV, sequence, request );
    }
    else {
        chameleon_pzlaswpc( ws, dir, A, IPIV, sequence, request );
    }

    CHAMELEON_zgetrf_WS_Free( ws );

    return CHAMELEON_SUCCESS;
}

