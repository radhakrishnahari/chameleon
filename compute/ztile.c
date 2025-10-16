/**
 *
 * @file ztile.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon auxiliary routines
 *
 * @version 1.3.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Lionel Eyraud-Dubois
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
 * @brief Conversion of CHAMELEON_Complex64_t matrix from LAPACK layout to tile
 *        layout.  Deprecated function, see CHAMELEON_zLap2Desc().
 *
 *******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[in,out] A
 *          Descriptor of the CHAMELEON matrix.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_zLapack_to_Tile( CHAMELEON_Complex64_t *Af77, int LDA, CHAM_desc_t *A )
{
    return CHAMELEON_zLap2Desc( ChamUpperLower, Af77, LDA, A );
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 * @brief Conversion of CHAMELEON_Complex64_t matrix from tile layout to LAPACK
 *        layout.  Deprecated function, see CHAMELEON_zDesc2Lap().
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Descriptor of the CHAMELEON matrix.
 *
 * @param[in,out] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_zTile_to_Lapack( CHAM_desc_t *A, CHAMELEON_Complex64_t *Af77, int LDA )
{
    return CHAMELEON_zDesc2Lap( ChamUpperLower, A, Af77, LDA );
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 * @brief Conversion of a CHAMELEON_Complex64_t matrix from LAPACK layout to
 *        CHAM_desct_t.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of the matrix A:
 *          = ChamUpper: A is upper triangular;
 *          = ChamLower: A is lower triangular;
 *          = ChamUpperLower: A is general.
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[in,out] A
 *          Descriptor of the CHAMELEON matrix.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zDesc2Lap
 * @sa CHAMELEON_cLap2Desc
 * @sa CHAMELEON_dLap2Desc
 * @sa CHAMELEON_sLap2Desc
 *
 */
int CHAMELEON_zLap2Desc( cham_uplo_t uplo, CHAMELEON_Complex64_t *Af77, int LDA, CHAM_desc_t *A )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t B;
    char *lapname;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zLap2Desc", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (chameleon_desc_check( A ) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zLap2Desc", "invalid descriptor");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    /* Create the B descriptor to handle the Lapack format matrix */
    chameleon_asprintf( &lapname, "%slap", A->name );
    status = chameleon_desc_init( &B, lapname, Af77, ChamComplexDouble, A->mb, A->nb,
                                  LDA, A->n, A->m, A->n, 1, 1,
                                  chameleon_getaddr_cm, chameleon_getblkldd_cm, NULL, NULL );
    free( lapname );
    if ( status != CHAMELEON_SUCCESS ) {
        chameleon_error("CHAMELEON_zLap2Desc", "Failed to create the descriptor");
        return status;
    }

    /* Start the computation */
    chameleon_sequence_create( chamctxt, &sequence );

    chameleon_pzlacpy( uplo, &B, A, sequence, &request );

    CHAMELEON_Desc_Flush( &B, sequence );
    CHAMELEON_Desc_Flush( A, sequence );

    chameleon_sequence_wait( chamctxt, sequence );

    /* Destroy temporary B descriptor */
    chameleon_desc_destroy( &B );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup CHAMELEON_Complex64_t
 *
 * @brief Conversion of CHAMELEON_Complex64_t matrix from LAPACK layout to tile
 *        layout. Deprecated function, see CHAMELEON_zDesc2Lap().
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of the matrix A:
 *          = ChamUpper: A is upper triangular;
 *          = ChamLower: A is lower triangular;
 *          = ChamUpperLower: A is general.
 *
 * @param[in] A
 *          Descriptor of the CHAMELEON matrix.
 *
 * @param[in,out] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa CHAMELEON_zLap2Desc
 * @sa CHAMELEON_cDesc2Lap
 * @sa CHAMELEON_dDesc2Lap
 * @sa CHAMELEON_sDesc2Lap
 *
 */
int CHAMELEON_zDesc2Lap( cham_uplo_t uplo, CHAM_desc_t *A, CHAMELEON_Complex64_t *Af77, int LDA )
{
    CHAM_context_t *chamctxt;
    RUNTIME_sequence_t *sequence = NULL;
    RUNTIME_request_t request = RUNTIME_REQUEST_INITIALIZER;
    CHAM_desc_t B;
    char *lapname;
    int status;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_fatal_error("CHAMELEON_zTile_to_Lapack", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (chameleon_desc_check( A ) != CHAMELEON_SUCCESS) {
        chameleon_error("CHAMELEON_zDesc2Lap", "invalid descriptor");
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    /* Create the B descriptor to handle the Lapack format matrix */
    chameleon_asprintf( &lapname, "%slap", A->name );
    status = chameleon_desc_init( &B, lapname, Af77, ChamComplexDouble, A->mb, A->nb,
                                  LDA, A->n, A->m, A->n, 1, 1,
                                  chameleon_getaddr_cm, chameleon_getblkldd_cm, NULL, NULL );
    free( lapname );
    if ( status != CHAMELEON_SUCCESS ) {
        chameleon_error("CHAMELEON_zDesc2Lap", "Failed to create the descriptor");
        return status;
    }

    /* Start the computation */
    chameleon_sequence_create( chamctxt, &sequence );

    chameleon_pzlacpy( uplo, A, &B, sequence, &request );

    CHAMELEON_Desc_Flush( A, sequence );
    CHAMELEON_Desc_Flush( &B, sequence );

    chameleon_sequence_wait( chamctxt, sequence );

    chameleon_desc_destroy( &B );

    status = sequence->status;
    chameleon_sequence_destroy( chamctxt, sequence );
    return status;
}
