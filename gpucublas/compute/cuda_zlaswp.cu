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
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Xavier Lacoste
 * @date 2023-09-27
 * @precisions normal z -> c d s
 *
 */
#include "gpucublas.h"

/**
 * Extracts the rows from the tile A into the tile B
 * B[i][:] = A[perm[i]][:] if perm[i] belongs to A.
 *
 * A, and B are stored column wise.
 */
__global__
void cuda_zlaswp_get_kernel(
    int m0, int m, int n, int k,
    const CHAMELEON_Complex64_t *A, int lda,
    CHAMELEON_Complex64_t       *B, int ldb,
    const int *perm )
{
    int i;
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if ( tid < n )
    {
        for( i=0; i<k; i++ )
        {
            int idx = perm[i] - m0;
            if ( ( idx >= 0 ) && (idx < m ) )
            {
                const CHAMELEON_Complex64_t *At = A + idx + tid*lda;
                CHAMELEON_Complex64_t *Bt = B + i + tid*ldb;

                *Bt = *At;
            }
        }
    }
}

/**
 ******************************************************************************
 *
 * @ingroup CUDA_CHAMELEON_Complex64_t
 *
 * CUDA_zlaswp_get extracts the rows from the tile A that have been selected as
 * pivot into the tile B.
 *
 * B[i][:] = A[perm[i]][:] if perm[i] belongs to A.
 *
 *******************************************************************************
 *
 * @param[in] m0
 *         The index of the first row of the tile A into the larger matrix it
 *         belongs to.
 *
 * @param[in] m
 *          The number of rows of the matrix A.
 *
 * @param[in] n
 *         The number of columns of the matrices A and B.
 *
 * @param[in] k
 *         The number of rows of the matrix B. This is the number of potential
 *         pivot that can be extracted from A.
 *
 * @param[in] A
 *          On entry, the matrix A of dimension lda-by-n where to extract the
 *          pivot rows if some are selected in the range m0..m0+m.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,m).
 *
 * @param[inout] B
 *          On entry, a matrix of size ldb-by-n with 0s or already collected
 *          rows.
 *          On exit, B is filled with the selected rows from A, such that for
 *          each row i, B[i] = A[perm[i]-m0-1].
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= max(1,k).
 *
 * @param[in] perm
 *          The permutation array of dimension k.
 *
 * @param[in] handle
 *          The cublas handle
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */
#define NTHREADS 64
int
CUDA_zlaswp_get( int m0, int m, int n, int k,
                 const CHAMELEON_Complex64_t *A, int lda,
                 CHAMELEON_Complex64_t       *B, int ldb,
                 const int *perm,
                 cublasHandle_t handle )
{
    int i;
    cudaStream_t stream;
    dim3 threads( NTHREADS );
    dim3 grid( chameleon_ceil( n, NTHREADS ) );

    /* Check input arguments */
    if (m0 < 0) {
        gpucublas_error(1, "Illegal value of m0");
        return -1;
    }
    if (m < 0) {
        gpucublas_error(2, "Illegal value of m");
        return -2;
    }
    if (n < 0) {
        gpucublas_error(3, "Illegal value of n");
        return -3;
    }
    if (k < 0) {
        gpucublas_error(4, "Illegal value of k");
        return -4;
    }
    if ((lda < chameleon_max(1,m)) && (m > 0)) {
        gpucublas_error(6, "Illegal value of lda");
        return -6;
    }
    if ((ldb < chameleon_max(1,k)) && (k > 0)) {
        gpucublas_error(8, "Illegal value of ldb");
        return -8;
    }

    /* Quick return */
    if ((m == 0) || (n == 0) || (k == 0)) {
        return CHAMELEON_SUCCESS;
    }

    cublasGetStream( handle, &stream );

    for( i=0; i<k; i++ )
    {
        int idx = perm[i] - m0;
        if ( ( idx >= 0 ) && (idx < m ) )
        {
            cublasZcopy(handle, n, A + idx, lda, B + i, ldb);
        }
    }
    /* cuda_zlaswp_get_kernel */
    /*         <<< grid, threads, 0, stream >>> */
    /*     ( m0, m, n, k, A, lda, B, ldb, perm ); */

    return CHAMELEON_SUCCESS;
}

/**
 * Copies the rows from the tile A into the tile B when they are the destination
 * of the pivoted rows.
 *
 * B[invp[i]][:] = A[i][:] if invp[i] belongs to B.
 *
 * A, and B are stored column wise.
 */
__global__
void cuda_zlaswp_set_kernel(
    int m0, int m, int n, int k,
    const CHAMELEON_Complex64_t *A, int lda,
    CHAMELEON_Complex64_t       *B, int ldb,
    const int *invp )
{
    int i;
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if ( tid < n )
    {
        for( i=0; i<k; i++ )
        {
            int idx = invp[i] - m0;

            if ( ( idx >= 0 ) && (idx < m ) )
            {
                const CHAMELEON_Complex64_t *At = A + i + tid*lda;
                CHAMELEON_Complex64_t *Bt = B + idx + tid*ldb;
                *Bt =*At;
            }
        }
    }
}
/**
 ******************************************************************************
 *
 * @ingroup CUDA_CHAMELEON_Complex64_t
 *
 * CUDA_zlaswp_set copies the rows from the tile A into the tile B when they are
 * the destination of the pivoted rows.
 *
 * B[invp[i]][:] = A[i][:] if invp[i] belongs to B.
 *
 *******************************************************************************
 *
 * @param[in] m0
 *         The index of the first row of the tile B into the larger matrix it
 *         belongs to.
 *
 * @param[in] m
 *          The number of rows of the matrix B.
 *
 * @param[in] n
 *         The number of columns of the matrices A and B.
 *
 * @param[in] k
 *         The number of rows of the matrix A. This is the number of potential
 *         pivot that can be inserted into B.
 *
 * @param[in] A
 *          On entry, the matrix A of dimension lda-by-n where to read the
 *          pivoted rows.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,k).
 *
 * @param[inout] B
 *          On entry, a matrix of size ldb-by-n that may require some pivoted rows.
 *          On exit, B is updated with the pivoted rows it needs to receive, such that for
 *          each row i, A[i] = B[invp[i]-m0-1].
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= max(1,m).
 *
 * @param[in] invp
 *          The inverse permutation array of dimension k.
 *
 * @param[in] handle
 *          The cublas handle
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 *
 */
int
CUDA_zlaswp_set( int m0, int m, int n, int k,
                 const CHAMELEON_Complex64_t *A, int lda,
                 CHAMELEON_Complex64_t *B, int ldb,
                 const int *invp,
                 cublasHandle_t handle )
{
    int i;
    cudaStream_t stream;
    dim3 threads( NTHREADS );
    dim3 grid( chameleon_ceil( n, NTHREADS ) );

    /* Check input arguments */
    if (m0 < 0) {
        gpucublas_error(1, "Illegal value of m0");
        return -1;
    }
    if (m < 0) {
        gpucublas_error(2, "Illegal value of m");
        return -2;
    }
    if (n < 0) {
        gpucublas_error(3, "Illegal value of n");
        return -3;
    }
    if (k < 0) {
        gpucublas_error(4, "Illegal value of k");
        return -4;
    }
    if ((lda < chameleon_max(1,k)) && (k > 0)) {
        gpucublas_error(6, "Illegal value of lda");
        return -6;
    }
    if ((ldb < chameleon_max(1,m)) && (m > 0)) {
        gpucublas_error(8, "Illegal value of ldb");
        return -8;
    }

    /* Quick return */
    if ((m == 0) || (n == 0) || (k == 0)) {
        return CHAMELEON_SUCCESS;
    }

    cublasGetStream( handle, &stream );

    for( i=0; i<k; i++ )
    {
        int idx = invp[i] - m0;
        if ( ( idx >= 0 ) && (idx < m ) )
        {
            cublasZcopy(handle, n, A + i, lda, B + idx, ldb);
        }
    }
    /* cuda_zlaswp_set_kernel */
    /*         <<< grid, threads, 0, stream >>> */
    /*     ( m0, m, n, k, A, lda, B, ldb, invp ); */


    return CHAMELEON_SUCCESS;
}
