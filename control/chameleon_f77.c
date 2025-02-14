/**
 *
 * @file chameleon_f77.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Fortran77 interface
 *
 * @version 1.3.0
 * @author Bilel Hadri
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Philippe Virouleau
 * @author Lionel Eyraud-Dubois
 * @date 2024-02-18
 *
 */
#include "control/common.h"
#include "chameleon.h"
#include "chameleon_f77.h"

#ifdef __cplusplus
extern "C" {
#endif

    /**
     *  FORTRAN API - auxiliary function prototypes
     */
    void __CHAMELEON_INIT(int *CORES, int *NGPUS, int *INFO)
    {   *INFO = __chameleon_init(*CORES, *NGPUS); }

    void __CHAMELEON_FINALIZE(int *INFO)
    {   *INFO = __chameleon_finalize(); }

    void CHAMELEON_ENABLE(int *lever, int *INFO)
    {   *INFO = CHAMELEON_Enable(*lever); }

    void CHAMELEON_DISABLE(int *lever, int *INFO)
    {   *INFO = CHAMELEON_Disable(*lever); }

    void CHAMELEON_SET(int *param, int *value, int *INFO)
    {   *INFO = CHAMELEON_Set(*param, *value); }

    void CHAMELEON_GET(int *param, int *value, int *INFO)
    {   *INFO = CHAMELEON_Get(*param, value); }

    void CHAMELEON_DEALLOC_HANDLE(size_t *sp, int *INFO)
    {   free((void *)(*sp));
        *INFO = CHAMELEON_SUCCESS; }

    void CHAMELEON_VERSION(int *VER_MAJOR, int *VER_MINOR, int *VER_MICRO, int *INFO)
    {
        *VER_MAJOR = CHAMELEON_VERSION_MAJOR;
        *VER_MINOR = CHAMELEON_VERSION_MINOR;
        *VER_MICRO = CHAMELEON_VERSION_MICRO;
        *INFO = CHAMELEON_SUCCESS;
    }

    /**
     *  FORTRAN API - descriptor allocation and deallocation
     */
    void CHAMELEON_DESC_CREATE(CHAM_desc_t **desc, void *mat, cham_flttype_t *dtyp,
                           int *mb, int *nb, int *bsiz, int *lm, int *ln,
                           int *i, int *j, int *m, int *n, int *p, int *q,
                           int *INFO)
    {   *INFO = CHAMELEON_Desc_Create(desc, mat, *dtyp, *mb, *nb, *bsiz, *lm, *ln, *i, *j, *m, *n, *p, *q); }
    void CHAMELEON_DESC_CREATE_OOC(CHAM_desc_t **desc, cham_flttype_t *dtyp,
                               int *mb, int *nb, int *bsiz, int *lm, int *ln,
                               int *i, int *j, int *m, int *n, int *p, int *q,
                               int *INFO)
    {   *INFO = CHAMELEON_Desc_Create_OOC(desc, *dtyp, *mb, *nb, *bsiz, *lm, *ln, *i, *j, *m, *n, *p, *q); }
    void CHAMELEON_DESC_CREATE_USER(CHAM_desc_t **descptr, void *mat, cham_flttype_t *dtyp,
                                int *mb, int *nb, int *bsiz, int *lm, int *ln,
                                int *i, int *j, int *m, int *n, int *p, int *q,
                                void* (*get_blkaddr)( const CHAM_desc_t*, int, int ),
                                int   (*get_blkldd) ( const CHAM_desc_t*, int      ),
                                int   (*get_rankof) ( const CHAM_desc_t*, int, int ),
                                void* get_rankof_arg,
                                int *INFO)
    {   *INFO = CHAMELEON_Desc_Create_User(descptr, mat, *dtyp, *mb, *nb, *bsiz, *lm, *ln, *i, *j, *m, *n, *p, *q,
                                           get_blkaddr, get_blkldd, get_rankof, get_rankof_arg); }
    void CHAMELEON_DESC_CREATE_OOC_USER(CHAM_desc_t **descptr, cham_flttype_t *dtyp,
                                    int *mb, int *nb, int *bsiz, int *lm, int *ln,
                                    int *i, int *j, int *m, int *n, int *p, int *q,
                                    int (*get_rankof) ( const CHAM_desc_t*, int, int ),
                                    void* get_rankof_arg,
                                    int *INFO)
    {   *INFO = CHAMELEON_Desc_Create_OOC_User(descptr, *dtyp, *mb, *nb, *bsiz, *lm, *ln, *i, *j, *m, *n, *p, *q,
                                               get_rankof, get_rankof_arg); }

    void CHAMELEON_DESC_DESTROY(CHAM_desc_t **desc, int *INFO)
    {   *INFO = CHAMELEON_Desc_Destroy(desc); }

    /**
     *  FORTRAN API - conversion from LAPACK F77 matrix layout to tile layout
     */
    void CHAMELEON_LAPACK_TO_TILE(intptr_t *Af77, int *LDA, intptr_t *A, int *INFO)
    {   *INFO = CHAMELEON_Lap2Desc( ChamUpperLower, (void *)Af77, *LDA, (CHAM_desc_t *)(*A) ); }

    void CHAMELEON_TILE_TO_LAPACK(intptr_t *A, intptr_t *Af77, int *LDA, int *INFO)
    {   *INFO = CHAMELEON_Desc2Lap( ChamUpperLower, (CHAM_desc_t *)(*A), (void *)Af77, *LDA ); }

#ifdef __cplusplus
}
#endif
