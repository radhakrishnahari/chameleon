/**
 *
 * @file core_zc.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CPU CHAMELEON_Complex64_t kernels header
 *
 * @version 1.3.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Matthieu Kuhn
 * @author Brieuc Nicolas
 * @date 2025-07-09
 * @precisions mixed zc -> ds
 *
 */
/* -*- c-file-style: "GNU" -*- */
#ifndef _REENTRANT
#define _REENTRANT
#endif

#include <eztrace-core/eztrace_config.h>
#include <eztrace-instrumentation/pptrace.h>
#include <eztrace-lib/eztrace.h>
#include <eztrace-lib/eztrace_module.h>

#include "chameleon/struct.h"
#include "chameleon/constants.h"
#include "chameleon/types.h"

/* set to 1 when all the hooks are set.
 * This is usefull in order to avoid recursive calls to mutex_lock for example
 */
static volatile int _coreblas_core_zc_initialized = 0;

#define CURRENT_MODULE coreblas_core_zc
DECLARE_CURRENT_MODULE;

/* pointers to actual coreblas_core_zc functions */


void (*libCORE_clag2z)(int m, int n, const CHAMELEON_Complex32_t * A, int lda, CHAMELEON_Complex64_t * B, int ldb) = NULL;
void (*libCORE_zlag2c)(int m, int n, const CHAMELEON_Complex64_t * A, int lda, CHAMELEON_Complex32_t * B, int ldb, int *info) = NULL;


/* Wrapper functions */



void CORE_clag2z(int m, int n, const CHAMELEON_Complex32_t * A, int lda, CHAMELEON_Complex64_t * B, int ldb) {
  FUNCTION_ENTRY;
   libCORE_clag2z(m, n, A, lda, B, ldb);
  FUNCTION_EXIT;

}



void CORE_zlag2c(int m, int n, const CHAMELEON_Complex64_t * A, int lda, CHAMELEON_Complex32_t * B, int ldb, int *info) {
  FUNCTION_ENTRY;
   libCORE_zlag2c(m, n, A, lda, B, ldb, info);
  FUNCTION_EXIT;

}



PPTRACE_START_INTERCEPT_FUNCTIONS(coreblas_core_zc)
   INTERCEPT3("CORE_clag2z", libCORE_clag2z)
 INTERCEPT3("CORE_zlag2c", libCORE_zlag2c)

PPTRACE_END_INTERCEPT_FUNCTIONS(coreblas_core_zc)
