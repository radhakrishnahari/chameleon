/**
 *
 * @file chameleon_mf77.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Fortran77 interface for mixed-precision computational routines
 *
 * @version 1.3.0
 * @author Bilel Hadri
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2024-02-18
 *
 */
#include "control/common.h"

#define CHAMELEON_ZCGESV           CHAMELEON_FNAME(zcgesv,   ZCGESV)
#define CHAMELEON_DSGESV           CHAMELEON_FNAME(dsgesv,   DSGESV)
#define CHAMELEON_ZCPOSV           CHAMELEON_FNAME(zcposv,   ZCPOSV)
#define CHAMELEON_DSPOSV           CHAMELEON_FNAME(dsposv,   DSPOSV)
#define CHAMELEON_ZCGELS           CHAMELEON_FNAME(zcgels,   ZCGELS)
#define CHAMELEON_DSGELS           CHAMELEON_FNAME(dsgels,   DSGELS)
#define CHAMELEON_ZCUNGESV         CHAMELEON_FNAME(zcungesv, ZCUNGESV)
#define CHAMELEON_DSUNGESV         CHAMELEON_FNAME(dsungesv, DSUNGESV)

#define CHAMELEON_ZCGESV_TILE       CHAMELEON_TILE_FNAME(zcgesv,   ZCGESV)
#define CHAMELEON_DSGESV_TILE       CHAMELEON_TILE_FNAME(dsgesv,   DSGESV)
#define CHAMELEON_ZCPOSV_TILE       CHAMELEON_TILE_FNAME(zcposv,   ZCPOSV)
#define CHAMELEON_DSPOSV_TILE       CHAMELEON_TILE_FNAME(dsposv,   DSPOSV)
#define CHAMELEON_ZCGELS_TILE       CHAMELEON_TILE_FNAME(zcgels,   ZCGELS)
#define CHAMELEON_DSGELS_TILE       CHAMELEON_TILE_FNAME(dsgels,   DSGELS)
#define CHAMELEON_ZCUNGESV_TILE     CHAMELEON_TILE_FNAME(zcungesv, ZCUNGESV)
#define CHAMELEON_DSUNGESV_TILE     CHAMELEON_TILE_FNAME(dsungesv, DSUNGESV)

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  FORTRAN API - math functions (simple interface)
 */
//void CHAMELEON_ZCGESV(int *N, int *NRHS, CHAMELEON_Complex64_t *A, int *LDA, int *IPIV, CHAMELEON_Complex64_t *B, int *LDB, CHAMELEON_Complex64_t *X, int *LDX, int *ITER, int *INFO)
//{   *INFO = CHAMELEON_zcgesv(*N, *NRHS, A, *LDA, IPIV, B, *LDB, X, *LDX, ITER); }
//
//void CHAMELEON_DSGESV(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, double *X, int *LDX, int *ITER, int *INFO)
//{   *INFO = CHAMELEON_dsgesv(*N, *NRHS, A, *LDA, IPIV, B, *LDB, X, *LDX, ITER); }
//
//void CHAMELEON_ZCPOSV(CHAMELEON_enum *uplo, int *N, int *NRHS, CHAMELEON_Complex64_t *A, int *LDA, CHAMELEON_Complex64_t *B, int *LDB, CHAMELEON_Complex64_t *X, int *LDX, int *ITER, int *INFO)
//{   *INFO = CHAMELEON_zcposv(*uplo, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); }
//
//void CHAMELEON_DSPOSV(CHAMELEON_enum *uplo, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, double *X, int *LDX, int *ITER, int *INFO)
//{   *INFO = CHAMELEON_dsposv(*uplo, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); }
//
//void CHAMELEON_ZCGELS(CHAMELEON_enum *trans, int *M, int *N, int *NRHS, CHAMELEON_Complex64_t *A, int *LDA, CHAMELEON_Complex64_t **T, CHAMELEON_Complex64_t *B, int *LDB, CHAMELEON_Complex64_t *X, int *LDX, int *ITER, int *INFO)
//{   *INFO = CHAMELEON_zcgels(*trans, *M, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); }
//
//void CHAMELEON_DSGELS(CHAMELEON_enum *trans, int *M, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, double *X, int *LDX, int *ITER, int *INFO)
//{   *INFO = CHAMELEON_dsgels(*trans, *M, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); }
//
//void CHAMELEON_ZCUNGESV(CHAMELEON_enum *trans, int *N, int *NRHS, CHAMELEON_Complex64_t *A, int *LDA, CHAMELEON_Complex64_t *B, int *LDB, CHAMELEON_Complex64_t *X, int *LDX, int *ITER, int *INFO)
//{   *INFO = CHAMELEON_zcungesv(*trans, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); }
//
//void CHAMELEON_DSUNGESV(CHAMELEON_enum *trans, int *N, int *NRHS, double *A, int *LDA, double *B, int *LDB, double *X, int *LDX, int *ITER, int *INFO)
//{   *INFO = CHAMELEON_dsungesv(*trans, *N, *NRHS, A, *LDA, B, *LDB, X, *LDX, ITER); }

#ifdef __cplusplus
}
#endif

/**
 *  FORTRAN API - math functions (native interface)
 */
//void CHAMELEON_ZCGESV_TILE(intptr_t *A, int *IPIV, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
//{   *INFO = CHAMELEON_zcgesv_Tile((CHAM_desc_t *)(*A), IPIV, (CHAM_desc_t *)(*B), (CHAM_desc_t *)(*X), ITER); }
//
//void CHAMELEON_DSGESV_TILE(intptr_t *A, int *IPIV, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
//{   *INFO = CHAMELEON_zcgesv_Tile((CHAM_desc_t *)(*A), IPIV, (CHAM_desc_t *)(*B), (CHAM_desc_t *)(*X), ITER); }
//
//void CHAMELEON_ZCPOSV_TILE(CHAMELEON_enum *uplo, intptr_t *A, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
//{   *INFO = CHAMELEON_zcposv_Tile(*uplo, (CHAM_desc_t *)(*A), (CHAM_desc_t *)(*B), (CHAM_desc_t *)(*X), ITER); }
//
//void CHAMELEON_DSPOSV_TILE(CHAMELEON_enum *uplo, intptr_t *A, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
//{   *INFO = CHAMELEON_dsposv_Tile(*uplo, (CHAM_desc_t *)(*A), (CHAM_desc_t *)(*B), (CHAM_desc_t *)(*X), ITER); }
//
//void CHAMELEON_ZCGELS_TILE(CHAMELEON_enum *trans, intptr_t *A, intptr_t *B, intptr_t *T, intptr_t *X, int *ITER, int *INFO)
//{   *INFO = CHAMELEON_zcgels_Tile(*trans, (CHAM_desc_t *)(*A), (CHAM_desc_t *)(*B), (CHAM_desc_t *)(*T), (CHAM_desc_t *)(*X), ITER); }
//
//void CHAMELEON_DSGELS_TILE(CHAMELEON_enum *trans, intptr_t *A, intptr_t *B, intptr_t *T, intptr_t *X, int *ITER, int *INFO)
//{   *INFO = CHAMELEON_dsgels_Tile(*trans, (CHAM_desc_t *)(*A), (CHAM_desc_t *)(*B), (CHAM_desc_t *)(*T), (CHAM_desc_t *)(*X), ITER); }
//
//void CHAMELEON_ZCUNGESV_TILE(CHAMELEON_enum *trans, intptr_t *A, intptr_t *T, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
//{   *INFO = CHAMELEON_zcungesv_Tile(*trans, (CHAM_desc_t *)(*A), (CHAM_desc_t *)(*T), (CHAM_desc_t *)(*B), (CHAM_desc_t *)(*X), ITER); }
//
//void CHAMELEON_DSUNGESV_TILE(CHAMELEON_enum *trans, intptr_t *A, intptr_t *T, intptr_t *B, intptr_t *X, int *ITER, int *INFO)
//{   *INFO = CHAMELEON_dsungesv_Tile(*trans, (CHAM_desc_t *)(*A), (CHAM_desc_t *)(*T), (CHAM_desc_t *)(*B), (CHAM_desc_t *)(*X), ITER); }
