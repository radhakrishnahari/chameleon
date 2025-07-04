/**
 *
 * @file gpucublas_zc.c
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
#include "gpucublas.h"

/* set to 1 when all the hooks are set.
 * This is usefull in order to avoid recursive calls to mutex_lock for example
 */
static volatile int _gpucublas_zc_initialized = 0;

#define CURRENT_MODULE gpucublas_zc
DECLARE_CURRENT_MODULE;

#if defined(GPUCUBLAS_HAVE_CUDA_TOOLKIT)
/* pointers to actual gpucublas_zc functions */
int (*libCUDA_clag2z)( int m, int n, const cuFloatComplex *A, int lda, cuDoubleComplex *B, int ldb, cublasHandle_t handle ) = NULL;
int (*libCUDA_zlag2c)( int m, int n, const cuDoubleComplex *A, int lda, cuFloatComplex *B, int ldb, cublasHandle_t handle ) = NULL;


/* Wrapper functions */

int CUDA_clag2z(int m, int n, const cuFloatComplex *A, int lda, cuDoubleComplex *B, int ldb, cublasHandle_t handle) {
  FUNCTION_ENTRY;
  int ret = libCUDA_clag2z(m, n, A, lda, B, ldb, handle);
  FUNCTION_EXIT;
  return ret;
}



int CUDA_zlag2c(int m, int n, const cuDoubleComplex *A, int lda, cuFloatComplex *B, int ldb, cublasHandle_t handle) {
  FUNCTION_ENTRY;
  int ret = libCUDA_zlag2c(m, n, A, lda, B, ldb, handle);
  FUNCTION_EXIT;
  return ret;
}
#endif


PPTRACE_START_INTERCEPT_FUNCTIONS(gpucublas_zc)
#if defined(GPUCUBLAS_HAVE_CUDA_TOOLKIT)
 INTERCEPT3("CUDA_clag2z", libCUDA_clag2z)
 INTERCEPT3("CUDA_zlag2c", libCUDA_zlag2c)
 #endif
PPTRACE_END_INTERCEPT_FUNCTIONS(gpucublas_zc)
