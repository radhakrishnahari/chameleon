/**
 *
 * @file context.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon context management routines
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2020-03-03
 *
 ***
 *
 * @defgroup Options
 * @brief Group routines exposed to users to handle options
 *
 */
#include "control/common.h"
#include "control/auxiliary.h"
#include "control/context.h"
#include "chameleon/runtime.h"
#include <stdlib.h>

#if !defined(CHAMELEON_SIMULATION)
#include "coreblas.h"
#endif

/**
 *  Global data
 */
/* master threads context lookup table */
static CHAM_context_t *chameleon_ctxt = NULL;

/**
 *  Create new context
 */
CHAM_context_t *chameleon_context_create()
{
    CHAM_context_t *chamctxt;

    if ( chameleon_ctxt != NULL ) {
        chameleon_error("chameleon_context_create", "a context is already existing\n");
        return NULL;
    }

    chamctxt = (CHAM_context_t*)malloc(sizeof(CHAM_context_t));
    if (chamctxt == NULL) {
        chameleon_error("chameleon_context_create", "malloc() failed");
        return NULL;
    }

    /* These initializations are just in case the user
       disables autotuning and does not set nb and ib */
    chamctxt->nb                 = 320;
    chamctxt->ib                 = 48;
    chamctxt->rhblock            = 4;
    chamctxt->lookahead          = 3;

    chamctxt->nworkers           = 1;
    chamctxt->ncudas             = 0;
    chamctxt->nthreads_per_worker= 1;

    chamctxt->warnings_enabled   = CHAMELEON_TRUE;
    chamctxt->autotuning_enabled = CHAMELEON_FALSE;
    chamctxt->parallel_enabled   = CHAMELEON_FALSE;
    chamctxt->profiling_enabled  = CHAMELEON_FALSE;
    chamctxt->progress_enabled   = CHAMELEON_FALSE;
    chamctxt->generic_enabled    = CHAMELEON_FALSE;

    chamctxt->householder        = ChamFlatHouseholder;
    chamctxt->translation        = ChamOutOfPlace;

    /* Initialize scheduler */
    RUNTIME_context_create(chamctxt);

    chameleon_ctxt = chamctxt;
    return chamctxt;
}


/**
 *  Return context for a thread
 */
CHAM_context_t *chameleon_context_self()
{
    return chameleon_ctxt;
}

/**
 *  Clean the context
 */
int chameleon_context_destroy(){

    RUNTIME_context_destroy(chameleon_ctxt);
    free(chameleon_ctxt);
    chameleon_ctxt = NULL;

    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  CHAMELEON_Enable - Enable CHAMELEON feature.
 *
 *******************************************************************************
 *
 * @param[in] option
 *          Feature to be enabled:
 *          @arg CHAMELEON_WARNINGS   printing of warning messages,
 *          @arg CHAMELEON_AUTOTUNING autotuning for tile size and inner block size.
 *          @arg CHAMELEON_PROFILING_MODE  activate profiling of kernels
 *          @arg CHAMELEON_PROGRESS  activate progress indicator
 *          @arg CHAMELEON_GEMM3M  Use z/cgemm3m for complexe matrix-matrix products
 *          @arg CHAMELEON_GENERIC  enable/disable GEMM3M  Use z/cgemm3m for complexe matrix-matrix products
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Enable(int option)
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Enable", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }

    switch (option)
    {
        case CHAMELEON_WARNINGS:
            chamctxt->warnings_enabled = CHAMELEON_TRUE;
            break;
        case CHAMELEON_AUTOTUNING:
            chamctxt->autotuning_enabled = CHAMELEON_TRUE;
            break;
        case CHAMELEON_PROFILING_MODE:
            RUNTIME_start_profiling();
            break;
        case CHAMELEON_KERNELPROFILE_MODE:
            chamctxt->profiling_enabled = CHAMELEON_TRUE;
            break;
        case CHAMELEON_PROGRESS:
            chamctxt->progress_enabled = CHAMELEON_TRUE;
            break;
        case CHAMELEON_GEMM3M:
#if defined(CBLAS_HAS_ZGEMM3M) && !defined(CHAMELEON_SIMULATION)
            set_coreblas_gemm3m_enabled(1);
#else
            chameleon_error("CHAMELEON_Enable", "cannot enable GEMM3M (not available in cblas)");
#endif
            break;
        /* case CHAMELEON_PARALLEL: */
        /*     chamctxt->parallel_enabled = CHAMELEON_TRUE; */
        /*     break; */
        case CHAMELEON_GENERIC:
            chamctxt->generic_enabled = CHAMELEON_TRUE;
            break;
        default:
            chameleon_error("CHAMELEON_Enable", "illegal parameter value");
            return CHAMELEON_ERR_ILLEGAL_VALUE;
        case CHAMELEON_BOUND:
            break;
    }

    /* Enable at the lower level if required */
    RUNTIME_enable( chamctxt->schedopt, option );

    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  CHAMELEON_Disable - Disable CHAMELEON feature.
 *
 *******************************************************************************
 *
 * @param[in] option
 *          Feature to be disabled:
 *          @arg CHAMELEON_WARNINGS   printing of warning messages,
 *          @arg CHAMELEON_AUTOTUNING autotuning for tile size and inner block size.
 *          @arg CHAMELEON_PROFILING_MODE  deactivate profiling of kernels
 *          @arg CHAMELEON_PROGRESS  deactivate progress indicator
 *          @arg CHAMELEON_GEMM3M  Use z/cgemm3m for complexe matrix-matrix products
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Disable(int option)
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Disable", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    switch ( option )
    {
        case CHAMELEON_WARNINGS:
            chamctxt->warnings_enabled = CHAMELEON_FALSE;
            break;
        case CHAMELEON_AUTOTUNING:
            chamctxt->autotuning_enabled = CHAMELEON_FALSE;
            break;
        case CHAMELEON_PROFILING_MODE:
            RUNTIME_stop_profiling();
            break;
        case CHAMELEON_KERNELPROFILE_MODE:
            chamctxt->profiling_enabled = CHAMELEON_FALSE;
            break;
        case CHAMELEON_PROGRESS:
            chamctxt->progress_enabled = CHAMELEON_FALSE;
            break;
        case CHAMELEON_GEMM3M:
#if defined(CBLAS_HAS_ZGEMM3M) && !defined(CHAMELEON_SIMULATION)
            set_coreblas_gemm3m_enabled(0);
#endif
            break;
        case CHAMELEON_PARALLEL_MODE:
            chamctxt->parallel_enabled = CHAMELEON_FALSE;
            break;
        case CHAMELEON_GENERIC:
            chamctxt->generic_enabled = CHAMELEON_FALSE;
            break;
        default:
            chameleon_error("CHAMELEON_Disable", "illegal parameter value");
            return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    /* Disable at the lower level if required */
    RUNTIME_disable( chamctxt->schedopt, option );

    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  CHAMELEON_Set - Set CHAMELEON parameter.
 *
 *******************************************************************************
 *
 * @param[in] param
 *          Feature to be enabled:
 *          @arg CHAMELEON_TILE_SIZE:        size matrix tile,
 *          @arg CHAMELEON_INNER_BLOCK_SIZE: size of tile inner block,
 *          @arg CHAMELEON_LOOKAHEAD:        depth of the look ahead in algorithms
 *
 * @param[in] value
 *          Value of the parameter.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Set( int param, int value )
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Set", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    switch (param) {
        case CHAMELEON_TILE_SIZE:
            if (value <= 0) {
                chameleon_error("CHAMELEON_Set", "negative tile size");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            chamctxt->nb = value;
            if ( chamctxt->autotuning_enabled ) {
                chamctxt->autotuning_enabled = CHAMELEON_FALSE;
                chameleon_warning("CHAMELEON_Set", "autotuning has been automatically disable\n");
            }
            /* Limit ib to nb */
            chamctxt->ib = chameleon_min( chamctxt->nb, chamctxt->ib );
            break;
        case CHAMELEON_INNER_BLOCK_SIZE:
            if (value <= 0) {
                chameleon_error("CHAMELEON_Set", "negative inner block size");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            if (value > chamctxt->nb) {
                chameleon_error("CHAMELEON_Set", "inner block larger than tile");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            /* if (chamctxt->nb % value != 0) { */
            /*     chameleon_error("CHAMELEON_Set", "inner block does not divide tile"); */
            /*     return CHAMELEON_ERR_ILLEGAL_VALUE; */
            /* } */
            chamctxt->ib = value;

            if ( chamctxt->autotuning_enabled ) {
                chamctxt->autotuning_enabled = CHAMELEON_FALSE;
                chameleon_warning("CHAMELEON_Set", "autotuning has been automatically disable\n");
            }
            break;
        case CHAMELEON_HOUSEHOLDER_MODE:
            if (value != ChamFlatHouseholder && value != ChamTreeHouseholder) {
                chameleon_error("CHAMELEON_Set", "illegal value of CHAMELEON_HOUSEHOLDER_MODE");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            chamctxt->householder = value;
            break;
        case CHAMELEON_HOUSEHOLDER_SIZE:
            if (value <= 0) {
                chameleon_error("CHAMELEON_Set", "negative householder size");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            chamctxt->rhblock = value;
            break;
        case CHAMELEON_TRANSLATION_MODE:
            if (value != ChamInPlace && value != ChamOutOfPlace) {
                chameleon_error("CHAMELEON_Set", "illegal value of CHAMELEON_TRANSLATION_MODE");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            chamctxt->translation = value;
            break;
        case CHAMELEON_LOOKAHEAD:
            if (value < 1) {
                chameleon_error("CHAMELEON_Set", "illegal value of CHAMELEON_LOOKAHEAD");
                return CHAMELEON_ERR_ILLEGAL_VALUE;
            }
            chamctxt->lookahead = value;
            break;
        default:
            chameleon_error("CHAMELEON_Set", "unknown parameter");
            return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    return CHAMELEON_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  CHAMELEON_Get - Get value of CHAMELEON parameter.
 *
 *******************************************************************************
 *
 * @param[in] param
 *          Feature to be enabled:
 *          @arg CHAMELEON_TILE_SIZE:        size matrix tile,
 *          @arg CHAMELEON_INNER_BLOCK_SIZE: size of tile inner block,
 *          @arg CHAMELEON_LOOKAHEAD:        depth of the look ahead in algorithms
 *
 * @param[out] value
 *          Value of the parameter.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int CHAMELEON_Get(int param, int *value)
{
    CHAM_context_t *chamctxt;

    chamctxt = chameleon_context_self();
    if (chamctxt == NULL) {
        chameleon_error("CHAMELEON_Get", "CHAMELEON not initialized");
        return CHAMELEON_ERR_NOT_INITIALIZED;
    }
    switch (param) {
        case CHAMELEON_TILE_SIZE:
            *value = chamctxt->nb;
            return CHAMELEON_SUCCESS;
        case CHAMELEON_INNER_BLOCK_SIZE:
            *value = chamctxt->ib;
            return CHAMELEON_SUCCESS;
        case CHAMELEON_HOUSEHOLDER_MODE:
            *value = chamctxt->householder;
            return CHAMELEON_SUCCESS;
        case CHAMELEON_HOUSEHOLDER_SIZE:
            *value = chamctxt->rhblock;
            return CHAMELEON_SUCCESS;
        case CHAMELEON_TRANSLATION_MODE:
            *value = chamctxt->translation;
            return CHAMELEON_SUCCESS;
        case CHAMELEON_LOOKAHEAD:
            *value = chamctxt->lookahead;
            return CHAMELEON_SUCCESS;
        default:
            chameleon_error("CHAMELEON_Get", "unknown parameter");
            return CHAMELEON_ERR_ILLEGAL_VALUE;
    }

    return CHAMELEON_SUCCESS;
}
