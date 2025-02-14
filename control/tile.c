/**
 *
 * @file tile.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon layout conversion wrappers
 *
 * @version 1.3.0
 * @author Jakub Kurzak
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Guillaume Sylvand
 * @author Mathieu Faverge
 * @date 2024-02-18
 *
 ***
 *
 * @defgroup Tile
 * @brief Group routines exposed to users for matrices conversion LAPACK-Tile
 *
 */
#include "control/common.h"
#include "control/auxiliary.h"

/**
 *
 * @ingroup Tile
 *
 * @brief  Conversion from LAPACK layout to tile layout
 *         Deprecated function see CHAMELEON_Lap2Desc().
 *
 ******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[in,out] A
 *          Descriptor of the CHAMELEON matrix in tile layout.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int
CHAMELEON_Lapack_to_Tile(  void *Af77, int LDA, CHAM_desc_t *A )
{
    return CHAMELEON_Lap2Desc( ChamUpperLower, Af77, LDA, A );
}

/**
 *
 * @ingroup Tile
 *
 * @brief Conversion from tile layout to LAPACK layout.
 *        Deprecated function, see CHAMELEON_Desc2Lap().
 *
 ******************************************************************************
 *
 * @param[in] A
 *          Descriptor of the CHAMELEON matrix in tile layout.
 *
 * @param[in,out] Af77
 *          LAPACK matrix (only needed on proc 0).
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int
CHAMELEON_Tile_to_Lapack( CHAM_desc_t *A, void *Af77, int LDA )
{
    return CHAMELEON_Desc2Lap( ChamUpperLower, A, Af77, LDA );
}

/**
 *
 * @ingroup Tile
 *
 * @brief Conversion from LAPACK layout to CHAM_desc_t.
 *
 ******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of the matrix A:
 *          = ChamUpper: A is upper triangular, the lower part is not referenced;
 *          = ChamLower: A is lower triangular, the upper part is not referenced;
 *          = ChamUpperLower: A is general.
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[in,out] A
 *          Descriptor of the CHAMELEON matrix initialized with data from Af77.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Lap2Desc( cham_uplo_t uplo, void *Af77, int LDA, CHAM_desc_t *A )
{
    switch( A->dtyp ) {
#if defined(PRECISION_z)
    case ChamComplexDouble:
        return CHAMELEON_zLap2Desc( uplo, (CHAMELEON_Complex64_t *)Af77, LDA, A );
        break;
#endif
#if defined(PRECISION_c)
    case ChamComplexFloat:
        return CHAMELEON_cLap2Desc( uplo, (CHAMELEON_Complex32_t *)Af77, LDA, A );
        break;
#endif
#if defined(PRECISION_s)
    case ChamRealFloat:
        return CHAMELEON_sLap2Desc( uplo, (float *)Af77, LDA, A );
        break;
#endif
#if defined(PRECISION_d)
    case ChamRealDouble:
        return CHAMELEON_dLap2Desc( uplo, (double *)Af77, LDA, A );
        break;
#endif
    default:
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
}

/**
 *
 * @ingroup Tile
 *
 * @brief Conversion from CHAM_desc_t to LAPACK layout.
 *
 ******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of the matrix A:
 *          = ChamUpper: A is upper triangular, the lower part is not referenced;
 *          = ChamLower: A is lower triangular, the upper part is not referenced;
 *          = ChamUpperLower: A is general.
 *
 * @param[in] A
 *          Descriptor of the CHAMELEON matrix in tile layout.
 *
 * @param[in,out] Af77
 *          LAPACK matrix (only needed on proc 0).
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 ******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Desc2Lap( cham_uplo_t uplo, CHAM_desc_t *A, void *Af77, int LDA )
{
    switch( A->dtyp ) {
#if defined(PRECISION_z)
    case ChamComplexDouble:
        return CHAMELEON_zDesc2Lap( uplo, A, (CHAMELEON_Complex64_t *)Af77, LDA );
        break;
#endif
#if defined(PRECISION_c)
    case ChamComplexFloat:
        return CHAMELEON_cDesc2Lap( uplo, A, (CHAMELEON_Complex32_t *)Af77, LDA );
        break;
#endif
#if defined(PRECISION_s)
    case ChamRealFloat:
        return CHAMELEON_sDesc2Lap( uplo, A, (float *)Af77, LDA );
        break;
#endif
#if defined(PRECISION_d)
    case ChamRealDouble:
        return CHAMELEON_dDesc2Lap( uplo, A, (double *)Af77, LDA );
        break;
#endif
    default:
        return CHAMELEON_ERR_ILLEGAL_VALUE;
    }
}
