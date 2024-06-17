/**
 *
 * @file gpucublas_zc.h
 *
 * @copyright 2023-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon GPU Mixed-precision kernels header
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @date 2023-07-06
 * @precisions mixed zc -> zc ds
 *
 */
#ifndef _gpucublas_zc_h_
#define _gpucublas_zc_h_

#if defined(GPUCUBLAS_HAVE_CUDA_TOOLKIT)
int CUDA_clag2z( int m, int n, const cuFloatComplex *A, int lda, cuDoubleComplex *B, int ldb, cublasHandle_t handle );
int CUDA_zlag2c( int m, int n, const cuDoubleComplex *A, int lda, cuFloatComplex *B, int ldb, cublasHandle_t handle );
#endif

#endif /* _gpucublas_zc_h_ */
