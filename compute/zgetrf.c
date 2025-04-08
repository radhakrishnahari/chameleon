/**
 *
 * @file zgetrf.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf wrappers
 *
 * @version 1.3.0
 * @author Omar Zenati
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Matthieu Kuhn
 * @author Lionel Eyraud-Dubois
 * @author Alycia Lisito
 * @author Xavier Lacoste
 * @author Pierre Esterie
 * @date 2024-12-09
 *
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"
#include <limits.h>

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zgetrf_WS_Alloc - Allocate the required workspaces for
 *  asynchronous getrf
 *
 *******************************************************************************
 *
 * @param[in] A
 *          The descriptor of the matrix A.
 *
 *******************************************************************************
 *
 * @retval An allocated opaque pointer to use in CHAMELEON_zgetrf_Tile_Async()
 *         and to free with CHAMELEON_zgetrf_WS_Free().
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrf_Tile_Async
 * @sa CHAMELEON_zgetrf_WS_Free
 *
 */
void *
CHAMELEON_zgetrf_WS_Alloc( const CHAM_desc_t *A )
{
    CHAM_context_t             *chamctxt;
    struct chameleon_pzgetrf_s *ws;
    int                         lookahead, batch_size;
    int                         P = chameleon_desc_datadist_get_iparam( A, 0 );
    int                         Q = chameleon_desc_datadist_get_iparam( A, 1 );

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        return NULL;
    }

    ws = calloc( 1, sizeof( struct chameleon_pzgetrf_s ) );
    ws->alg = ChamGetrfPPiv;
    ws->ib  = CHAMELEON_IB;

    ws->laswp = CHAMELEON_zlaswp_WS_Alloc( ChamLeft, A );

    {
        char *algostr = chameleon_getenv( "CHAMELEON_GETRF_ALGO" );

        if ( algostr != NULL ) {
            if ( strcasecmp( algostr, "nopiv" ) == 0 ) {
                ws->alg = ChamGetrfNoPiv;
            }
            else if ( strcasecmp( algostr, "nopivpercolumn" ) == 0  ) {
                ws->alg = ChamGetrfNoPivPerColumn;
            }
            else if ( strcasecmp( algostr, "ppiv" )  == 0 ) {
                ws->alg = ChamGetrfPPiv;
            }
            else if ( strcasecmp( algostr, "ppivpercolumn" ) == 0  ) {
                ws->alg = ChamGetrfPPivPerColumn;
            }
            else {
                chameleon_error( "CHAMELEON_zgetrf_WS_Alloc", "CHAMELEON_GETRF_ALGO is not one of NoPiv, NoPivPerColumn, PPiv, PPivPerColumn => Switch back to NoPiv\n" );
            }
        }
        chameleon_cleanenv( algostr );
    }

    batch_size = chameleon_getenv_get_value_int( "CHAMELEON_GETRF_BATCH_SIZE", 0 );
    if ( batch_size > CHAMELEON_BATCH_SIZE ) {
        chameleon_warning( "CHAMELEON_BATCH_SIZE", "CHAMELEON_GETRF_BATCH_SIZE must be smaller than CHAMELEON_BATCH_SIZE, please recompile with the right CHAMELEON_BATCH_SIZE, or reduce the CHAMELEON_GETRF_BATCH_SIZE value\n" );
    }
    ws->batch_size_blas2 = chameleon_getenv_get_value_int( "CHAMELEON_GETRF_BATCH_SIZE_BLAS2", batch_size );
    ws->batch_size_blas2 = ( ws->batch_size_blas2 > CHAMELEON_BATCH_SIZE ) ? CHAMELEON_BATCH_SIZE : ws->batch_size_blas2;
    ws->batch_size_blas3 = chameleon_getenv_get_value_int( "CHAMELEON_GETRF_BATCH_SIZE_BLAS3", batch_size );
    ws->batch_size_blas3 = ( ws->batch_size_blas3 > CHAMELEON_BATCH_SIZE ) ? CHAMELEON_BATCH_SIZE : ws->batch_size_blas3;
    ws->batch_size_swap = chameleon_getenv_get_value_int( "CHAMELEON_GETRF_BATCH_SIZE_SWAP", batch_size );
    ws->batch_size_swap = ( ws->batch_size_swap > CHAMELEON_BATCH_SIZE ) ? CHAMELEON_BATCH_SIZE : ws->batch_size_swap;

    ws->ringswitch = chameleon_getenv_get_value_int( "CHAMELEON_GETRF_RINGSWITCH", INT_MAX );

    /* Allocation of U for permutation of the panels */
    if ( ws->alg == ChamGetrfNoPivPerColumn ) {
        chameleon_desc_init( &(ws->U), CHAMELEON_MAT_ALLOC_TILE,
                             ChamComplexDouble, 1, A->nb, A->nb,
                             A->mt, A->nt * A->nb, 0, 0,
                             A->mt, A->nt * A->nb, P, Q,
                             NULL, NULL, A->get_rankof_init, A->get_rankof_init_arg );
    }
    else if ( ( ws->alg == ChamGetrfPPiv )          ||
              ( ws->alg == ChamGetrfPPivPerColumn ) )
    {
        chameleon_desc_init( &(ws->U), CHAMELEON_MAT_ALLOC_TILE,
                             ChamComplexDouble, A->mb, A->nb, A->mb*A->nb,
                             A->m, A->n, 0, 0,
                             A->m, A->n, P, Q,
                             NULL, NULL, A->get_rankof_init, A->get_rankof_init_arg );
        lookahead = chamctxt->lookahead;
        chameleon_desc_init( &(ws->Wl), CHAMELEON_MAT_ALLOC_TILE,
                             ChamComplexDouble, A->mb, A->nb, (A->mb * A->nb),
                             A->mt * A->mb, A->nb * Q * lookahead, 0, 0,
                             A->mt * A->mb, A->nb * Q * lookahead, P, Q,
                             NULL, NULL, A->get_rankof_init, A->get_rankof_init_arg );
    }

    /* Set ib to 1 if per column algorithm */
    if ( ( ws->alg == ChamGetrfNoPivPerColumn ) ||
         ( ws->alg == ChamGetrfPPivPerColumn ) )
    {
        ws->ib = 1;
    }

    /* Allocation of Up for the permutation of the diagonal panel per block */
    if ( ws->alg == ChamGetrfPPiv )
    {
        /* TODO: Should be restricted to diagonal tiles */
        /* Possibly to a single handle with a permutation of the ownership */
        chameleon_desc_init( &(ws->Up), CHAMELEON_MAT_ALLOC_TILE,
                             ChamComplexDouble, ws->ib, A->nb, ws->ib * A->nb,
                             A->mt * ws->ib, A->nt * A->nb, 0, 0,
                             A->mt * ws->ib, A->nt * A->nb, P, Q,
                             NULL, NULL, A->get_rankof_init, A->get_rankof_init_arg );
    }

    return ws;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 * @brief Free the allocated workspaces for asynchronous getrf
 *
 *******************************************************************************
 *
 * @param[in,out] user_ws
 *          On entry, the opaque pointer allocated by
 *          CHAMELEON_zgetrf_WS_Alloc() On exit, all data are freed.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrf_Tile_Async
 * @sa CHAMELEON_zgetrf_WS_Alloc
 *
 */
void
CHAMELEON_zgetrf_WS_Free( void *user_ws )
{
    struct chameleon_pzgetrf_s *ws = (struct chameleon_pzgetrf_s *)user_ws;

    if ( ( ws->alg == ChamGetrfNoPivPerColumn ) ||
         ( ws->alg == ChamGetrfPPiv           ) ||
         ( ws->alg == ChamGetrfPPivPerColumn  ) )
    {
        chameleon_desc_destroy( &(ws->U) );
    }
    if ( ws->alg == ChamGetrfPPiv )
    {
        chameleon_desc_destroy( &(ws->Up) );
    }
    if ( ( ws->alg == ChamGetrfPPiv           ) ||
         ( ws->alg == ChamGetrfPPivPerColumn  ) )
    {
        chameleon_desc_destroy( &(ws->Wl) );
    }
    free( ws );
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 *  CHAMELEON_zgetrf - Computes an LU factorization of a general M-by-N matrix A
 *  using the tile LU algorithm without row pivoting.
 *  WARNING: Don't use this function if you are not sure your matrix is diagonal
 *  dominant.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[out] IPIV
 *          Integer array of dimension min(M,N).
 *          The pivot indices; for 1 <= i <= min(M,N), row i of the
 *          matrix was interchanged with row IPIV(i).
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 * @retval >0 if i, U(i,i) is exactly zero. The factorization has been
 *               completed, but the factor U is exactly singular, and division
 *               by zero will occur if it is used to solve a system of
 *               equations.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrf_Tile
 * @sa CHAMELEON_zgetrf_Tile_Async
 * @sa CHAMELEON_cgetrf
 * @sa CHAMELEON_dgetrf
 * @sa CHAMELEON_sgetrf
 *
 */
int
CHAMELEON_zgetrf( int M, int N, CHAMELEON_Complex64_t *A, int LDA, int *IPIV )
{
    int                 NB;
    int                         status;
    CHAM_desc_t                 descAl, descAt;
    CHAM_ipiv_t                 descIPIV;
    CHAM_context_t             *chamctxt;
    RUNTIME_sequence_t         *sequence = NULL;
    RUNTIME_request_t           request  = RUNTIME_REQUEST_INITIALIZER;
    struct chameleon_pzgetrf_s *ws;
    int                         P, Q;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error( "CHAMELEON_zgetrf", "CHAMELEON not initialized" );
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ( M < 0 ) {
        chameleon_error( "CHAMELEON_zgetrf", "illegal value of M" );
        return -1;
    }
    if ( N < 0 ) {
        chameleon_error( "CHAMELEON_zgetrf", "illegal value of N" );
        return -2;
    }
    if ( LDA < chameleon_max( 1, M ) ) {
        chameleon_error( "CHAMELEON_zgetrf", "illegal value of LDA" );
        return -4;
    }
    /* Quick return */
    if ( chameleon_min( M, N ) == 0 )
        return CHAMELEON_SUCCESS;

    /* Tune NB & IB depending on M, N & NRHS; Set NBNBSIZE */
    status = chameleon_tune( CHAMELEON_FUNC_ZGESV, M, N, 0 );
    if ( status != CHAMELEON_SUCCESS ) {
        chameleon_error( "CHAMELEON_zgetrf", "chameleon_tune() failed" );
        return status;
    }

    /* Set NT & NTRHS */
    NB = CHAMELEON_NB;

    chameleon_sequence_create( chamctxt, &sequence );

    /* Submit the matrix conversion */
    chameleon_zlap2tile( chamctxt, &descAl, &descAt, ChamDescInout, ChamUpperLower,
                         A, NB, NB, LDA, N, M, N, sequence, &request );

    P = chameleon_desc_datadist_get_iparam( &descAt, 0 );
    Q = chameleon_desc_datadist_get_iparam( &descAt, 1 );

    /* Allocate workspace for partial pivoting */
    ws = CHAMELEON_zgetrf_WS_Alloc( &descAt );

    if ( ( ws->alg == ChamGetrfPPivPerColumn ) ||
         ( ws->alg == ChamGetrfPPiv ) )
    {
        chameleon_ipiv_init( &descIPIV, ChamLeft, descAt.mb, chameleon_min( M, N ), P, P*Q, IPIV, chameleon_getrankof_ipiv_2d_diag);
    }

    /* Call the tile interface */
    CHAMELEON_zgetrf_Tile_Async( &descAt, &descIPIV, ws, sequence, &request );

    /* Submit the matrix conversion back */
    chameleon_ztile2lap( chamctxt, &descAl, &descAt,
                         ChamDescInout, ChamUpperLower, sequence, &request );

    if ( ( ws->alg == ChamGetrfPPivPerColumn ) ||
         ( ws->alg == ChamGetrfPPiv ) )
    {
        RUNTIME_ipiv_gather( sequence, &descIPIV, IPIV, 0 );
    }
    chameleon_sequence_wait( chamctxt, sequence );

    /* Cleanup the temporary data */
    if ( ( ws->alg == ChamGetrfPPivPerColumn ) ||
         ( ws->alg == ChamGetrfPPiv ) )
    {
        chameleon_ipiv_destroy( &descIPIV );
    }
    CHAMELEON_zgetrf_WS_Free( ws );
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
 *  CHAMELEON_zgetrf_Tile - Computes the tile LU factorization of a matrix.
 *  Tile equivalent of CHAMELEON_zgetrf().  Operates on matrices stored by
 *  tiles.  All matrices are passed through descriptors.  All dimensions are
 *  taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[in,out] IPIV
 *          On entry, ipiv descriptor associated to A and created with
 *          CHAMELEON_Ipiv_Create().
 *          On exit, it contains the pivot indices associated to the PLU
 *          factorization of A.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval >0 if i, U(i,i) is exactly zero. The factorization has been
 *               completed, but the factor U is exactly singular, and division
 *               by zero will occur if it is used to solve a system of
 *               equations.
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zgetrf
 * @sa CHAMELEON_zgetrf_Tile_Async
 * @sa CHAMELEON_cgetrf_Tile
 * @sa CHAMELEON_dgetrf_Tile
 * @sa CHAMELEON_sgetrf_Tile
 * @sa CHAMELEON_zgetrs_Tile
 *
 */
int
CHAMELEON_zgetrf_Tile( CHAM_desc_t *A, CHAM_ipiv_t *IPIV )
{
    CHAM_context_t     *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t   request  = RUNTIME_REQUEST_INITIALIZER;
    int                 status;
    void               *ws;

    chamctxt = chameleon_context_self();
    if ( chamctxt == NULL ) {
        chameleon_fatal_error( "CHAMELEON_zgetrf_Tile", "CHAMELEON not initialized" );
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    chameleon_sequence_create( chamctxt, &sequence );

    ws = CHAMELEON_zgetrf_WS_Alloc( A );
    CHAMELEON_zgetrf_Tile_Async( A, IPIV, ws, sequence, &request );
    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    CHAMELEON_zgetrf_WS_Free( ws );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t_Tile_Async
 *
 *  CHAMELEON_zgetrf_Tile_Async - Computes the tile LU factorization of a
 *  matrix.  Non-blocking equivalent of CHAMELEON_zgetrf_Tile().  May return
 *  before the computation is finished.  Allows for pipelining of operations ar
 *  runtime.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the tile factors L and U from the factorization.
 *
 * @param[in,out] IPIV
 *          On entry, ipiv descriptor associated to A and created with
 *          CHAMELEON_Ipiv_Create().
 *          On exit, it contains the pivot indices associated to the PLU
 *          factorization of A.
 *
 * @param[in,out] user_ws
 *          The opaque pointer to pre-allocated getrf workspace through
 *          CHAMELEON_zgetrf_WS_Alloc(). If user_ws is NULL, it is automatically
 *          allocated, but BE CAREFULL as it switches the call from asynchronous
 *          to synchronous call.
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
 * @sa CHAMELEON_zgetrf
 * @sa CHAMELEON_zgetrf_Tile
 * @sa CHAMELEON_cgetrf_Tile_Async
 * @sa CHAMELEON_dgetrf_Tile_Async
 * @sa CHAMELEON_sgetrf_Tile_Async
 * @sa CHAMELEON_zgetrs_Tile_Async
 *
 */
int
CHAMELEON_zgetrf_Tile_Async( CHAM_desc_t        *A,
                             CHAM_ipiv_t        *IPIV,
                             void               *user_ws,
                             RUNTIME_sequence_t *sequence,
                             RUNTIME_request_t  *request )
{
    CHAM_context_t             *chamctxt;
    struct chameleon_pzgetrf_s *ws;
    chamctxt = chameleon_context_self();

    if ( chamctxt == NULL ) {
        chameleon_fatal_error( "CHAMELEON_zgetrf_Tile_Async", "CHAMELEON not initialized" );
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if ( chamctxt->scheduler != RUNTIME_SCHED_STARPU ) {
        chameleon_fatal_error( "CHAMELEON_zgetrf_Tile_Async", "CHAMELEON_zgetrf_Tile_Async is only available with StarPU" );
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    if ( sequence == NULL ) {
        chameleon_fatal_error( "CHAMELEON_zgetrf_Tile_Async", "NULL sequence" );
        return CHAMELEON_ERR_UNALLOCATED;
    }
    if ( request == NULL ) {
        chameleon_fatal_error( "CHAMELEON_zgetrf_Tile_Async", "NULL request" );
        return CHAMELEON_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if ( sequence->status == CHAMELEON_SUCCESS ) {
        request->status = CHAMELEON_SUCCESS;
    }
    else {
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_SEQUENCE_FLUSHED );
    }

    /* Check descriptors for correctness */
    if ( chameleon_desc_check( A ) != CHAMELEON_SUCCESS ) {
        chameleon_error( "CHAMELEON_zgetrf_Tile", "invalid first descriptor" );
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE );
    }

    /* Check input arguments */
    if ( A->nb != A->mb ) {
        chameleon_error( "CHAMELEON_zgetrf_Tile", "only square tiles supported" );
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE );
    }
    if ( IPIV->mb != A->mb ) {
        chameleon_error( "CHAMELEON_zgetrf_Tile", "IPIV tiles must have the number of rows as tiles of A" );
        return chameleon_request_fail( sequence, request, CHAMELEON_ERR_ILLEGAL_VALUE );
    }

    if ( user_ws == NULL ) {
        ws = CHAMELEON_zgetrf_WS_Alloc( A );
    }
    else {
        ws = user_ws;
    }

    IPIV->get_rankof = chameleon_getrankof_ipiv_2d_diag;

    chameleon_pzgetrf( ws, A, IPIV, sequence, request );

    if ( user_ws == NULL ) {
        CHAMELEON_Desc_Flush( A, sequence );
        chameleon_sequence_wait( chamctxt, sequence );
        CHAMELEON_zgetrf_WS_Free( ws );
    }
    return CHAMELEON_SUCCESS;
}
