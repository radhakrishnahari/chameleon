/**
 *
 * @file cuda_zlaswp.c
 *
 * @copyright 2012-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon cuda_zgetrf with partial pivoting CPU kernel
 *
 * @version 1.0.0
 * @author Xavier Lacoste
 * @date 2023-09-27
 * @precisions normal z -> c d s
 *
 */
#include "gpucublas.h"


// BLK_X and BLK_Y need to be equal for zlaset_q to deal with diag & offdiag
// when looping over super blocks.
// Formerly, BLK_X and BLK_Y could be different.
#define BLK_X 64
#define BLK_Y BLK_X

/******************************************************************************/
/**
 *   Divides matrix into ceil( m/BLK_X ) x ceil( n/BLK_Y ) blocks.
 *   Each block has BLK_X threads.
 *   Each thread loops across one row, updating BLK_Y entries.
 *
 *   Code imported from magmablas
 *
 *       @author Mark Gates
 *       @author Azzam Haidar
 */
static __device__
void zlacpy_full_device(
                        int m, int n,
                        const CHAMELEON_Complex64_t *dA, int ldda,
                        CHAMELEON_Complex64_t       *dB, int lddb )
{
  int ind = blockIdx.x*BLK_X + threadIdx.x;
  int iby = blockIdx.y*BLK_Y;
  /* check if full block-column */
  bool full = (iby + BLK_Y <= n);
  /* do only rows inside matrix */
  if ( ind < m ) {
    dA += ind + iby*ldda;
    dB += ind + iby*lddb;
    if ( full ) {
      // full block-column
#pragma unroll
      for( int j=0; j < BLK_Y; ++j ) {
        dB[j*lddb] = dA[j*ldda];
      }
    }
    else {
      // partial block-column
      for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
        dB[j*lddb] = dA[j*ldda];
      }
    }
  }
}


/**
 *   Similar to zlacpy_full, but updates only the diagonal and below.
 *   Blocks that are fully above the diagonal exit immediately.
 *
 *   Code imported from magmablas
 *
 *       @author Mark Gates
 *       @author Azzam Haidar
 */
static __device__
void zlacpy_lower_device(
    int m, int n,
    const CHAMELEON_Complex64_t *dA, int ldda,
    CHAMELEON_Complex64_t       *dB, int lddb )
{
  int ind = blockIdx.x*BLK_X + threadIdx.x;
  int iby = blockIdx.y*BLK_Y;
  /* check if full block-column && (below diag) */
  bool full = (iby + BLK_Y <= n && (ind >= iby + BLK_Y));
  /* do only rows inside matrix, and blocks not above diag */
  if ( ind < m && ind + BLK_X > iby ) {
    dA += ind + iby*ldda;
    dB += ind + iby*lddb;
    if ( full ) {
      // full block-column, off-diagonal block
#pragma unroll
      for( int j=0; j < BLK_Y; ++j ) {
        dB[j*lddb] = dA[j*ldda];
      }
    }
    else {
      // either partial block-column or diagonal block
      for( int j=0; j < BLK_Y && iby+j < n && ind >= iby+j; ++j ) {
        dB[j*lddb] = dA[j*ldda];
      }
    }
  }
}


/**
 *   Similar to zlacpy_full, but updates only the diagonal and above.
 *   Blocks that are fully below the diagonal exit immediately.
 *
 *   Code imported from magmablas
 *
 *       @author Mark Gates
 *       @author Azzam Haidar
 */
static __device__
void zlacpy_upper_device(
                         int m, int n,
                         const CHAMELEON_Complex64_t *dA, int ldda,
                         CHAMELEON_Complex64_t       *dB, int lddb )
{
  int ind = blockIdx.x*BLK_X + threadIdx.x;
  int iby = blockIdx.y*BLK_Y;
  /* check if full block-column && (above diag) */
  bool full = (iby + BLK_Y <= n && (ind + BLK_X <= iby));
  /* do only rows inside matrix, and blocks not below diag */
  if ( ind < m && ind < iby + BLK_Y ) {
    dA += ind + iby*ldda;
    dB += ind + iby*lddb;
    if ( full ) {
      // full block-column, off-diagonal block
#pragma unroll
      for( int j=0; j < BLK_Y; ++j ) {
        dB[j*lddb] = dA[j*ldda];
      }
    }
    else {
      // either partial block-column or diagonal block
      for( int j=0; j < BLK_Y && iby+j < n; ++j ) {
        if ( ind <= iby+j ) {
          dB[j*lddb] = dA[j*ldda];
        }
      }
    }
  }
}


/******************************************************************************/
/*
 * kernel wrappers to call the device functions.
 *
 *   Code imported from magmablas
 *
 *       @author Mark Gates
 *       @author Azzam Haidar
 */
__global__
void zlacpy_full_kernel(
                        int m, int n,
                        const CHAMELEON_Complex64_t *dA, int ldda,
                        CHAMELEON_Complex64_t       *dB, int lddb )
{
  zlacpy_full_device(m, n, dA, ldda, dB, lddb);
}

__global__
void zlacpy_lower_kernel(
                         int m, int n,
                         const CHAMELEON_Complex64_t *dA, int ldda,
                         CHAMELEON_Complex64_t       *dB, int lddb )
{
  zlacpy_lower_device(m, n, dA, ldda, dB, lddb);
}

__global__
void zlacpy_upper_kernel(
                         int m, int n,
                         const CHAMELEON_Complex64_t *dA, int ldda,
                         CHAMELEON_Complex64_t       *dB, int lddb )
{
  zlacpy_upper_device(m, n, dA, ldda, dB, lddb);
}



/**
 *******************************************************************************
 *
 * @ingroup CORE_CHAMELEON_Complex64_t
 *
 *  Copies all or part of a two-dimensional matrix A to another matrix B.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - ChamGeneral: entire A,
 *          - ChamUpper:   upper triangle,
 *          - ChamLower:   lower triangle.
 *
 * @param[in] M
 *          The number of rows of the matrices A and B.
 *          m >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrices A and B.
 *          n >= 0.
 *
 * @param[in] A
 *          The m-by-n matrix to copy.
 *
 * @param[in] LDA
 *          The leading dimension of the array A.
 *          lda >= max(1,m).
 *
 * @param[out] B
 *          The m-by-n copy of the matrix A.
 *          On exit, B = A ONLY in the locations specified by uplo.
 *
 * @param[in] LDB
 *          The leading dimension of the array B.
 *          ldb >= max(1,m).
 *
 */
extern "C"
int CUDA_zlacpy(cham_uplo_t uplo, int M, int N,
                const CHAMELEON_Complex64_t *A, int LDA,
                CHAMELEON_Complex64_t *B, int LDB,
                cublasHandle_t handle )
{
#define dA(i_, j_) (dA + (i_) + (j_)*ldda)
#define dB(i_, j_) (dB + (i_) + (j_)*lddb)

  cudaStream_t stream;
  magma_int_t info = 0;
  if ( uplo != ChamLower && uplo != ChamUpper && uplo != ChamUpperLower )
    info = -1;
  else if ( m < 0 )
    info = -2;
  else if ( n < 0 )
    info = -3;
  else if ( ldda < max(1,m))
    info = -5;
  else if ( lddb < max(1,m))
    info = -7;

  if ( info != 0 ) {
    return info;
  }

  if ( m == 0 || n == 0 ) {
    return 0;
  }

      cublasGetStream( handle, &stream );

  assert( BLK_X == BLK_Y );
  const int super_NB = max_blocks*BLK_X;
  dim3 super_grid( magma_ceildiv( m, super_NB ), magma_ceildiv( n, super_NB ) );

  dim3 threads( BLK_X, 1 );
  dim3 grid;

  int mm, nn;
  if ( uplo == ChamLower ) {
    for( unsigned int i=0; i < super_grid.x; ++i ) {
      mm = (i == super_grid.x-1 ? m % super_NB : super_NB);
      grid.x = magma_ceildiv( mm, BLK_X );
      for( unsigned int j=0; j < super_grid.y && j <= i; ++j ) {  // from left to diagonal
        nn = (j == super_grid.y-1 ? n % super_NB : super_NB);
        grid.y = magma_ceildiv( nn, BLK_Y );
        if ( i == j ) {  // diagonal super block
          zlacpy_lower_kernel<<< grid, threads, 0, stream >>>
            ( mm, nn, dA(i*super_NB, j*super_NB), ldda, dB(i*super_NB, j*super_NB), lddb );
        }
        else {           // off diagonal super block
          zlacpy_full_kernel <<< grid, threads, 0, stream >>>
                                 ( mm, nn, dA(i*super_NB, j*super_NB), ldda, dB(i*super_NB, j*super_NB), lddb );
        }
      }
    }
  }
  else if ( uplo == ChamUpper ) {
    for( unsigned int i=0; i < super_grid.x; ++i ) {
      mm = (i == super_grid.x-1 ? m % super_NB : super_NB);
      grid.x = magma_ceildiv( mm, BLK_X );
      for( unsigned int j=i; j < super_grid.y; ++j ) {  // from diagonal to right
        nn = (j == super_grid.y-1 ? n % super_NB : super_NB);
        grid.y = magma_ceildiv( nn, BLK_Y );
        if ( i == j ) {  // diagonal super block
          zlacpy_upper_kernel<<< grid, threads, 0, stream >>>
            ( mm, nn, dA(i*super_NB, j*super_NB), ldda, dB(i*super_NB, j*super_NB), lddb );
        }
        else {           // off diagonal super block
          zlacpy_full_kernel <<< grid, threads, 0, stream >>>
            ( mm, nn, dA(i*super_NB, j*super_NB), ldda, dB(i*super_NB, j*super_NB), lddb );
        }
      }
    }
  }
  else {
    // TODO: use cudaMemcpy or cudaMemcpy2D ?
    for( unsigned int i=0; i < super_grid.x; ++i ) {
      mm = (i == super_grid.x-1 ? m % super_NB : super_NB);
      grid.x = magma_ceildiv( mm, BLK_X );
      for( unsigned int j=0; j < super_grid.y; ++j ) {  // full row
        nn = (j == super_grid.y-1 ? n % super_NB : super_NB);
        grid.y = magma_ceildiv( nn, BLK_Y );
        zlacpy_full_kernel <<< grid, threads, 0, stream >>>
          ( mm, nn, dA(i*super_NB, j*super_NB), ldda, dB(i*super_NB, j*super_NB), lddb );
      }
    }
  }
}
