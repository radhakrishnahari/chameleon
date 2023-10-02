/**
 *
 * @file testing_zlaswp.c
 *
 * @copyright 2019-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaswp testing
 *
 * @version 1.3.0
 * @author Xavier Lacoste
 * @date 2023-28-09
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>
#if !defined(CHAMELEON_SIMULATION)
#include <coreblas.h>
#endif
#ifdef CHAMELEON_USE_CUDA
#include <gpucublas.h>
#include <cuda_runtime.h>
#endif

static cham_fixdbl_t
flops_zlaswp_set( int m0, int m, int n, int k,
                 const int *invp )
{
    int i;
    cham_fixdbl_t flops = 0;
    for( i=0; i<k; i++ )
    {
        int idx = invp[i] - m0;
        if ( ( idx >= 0 ) && (idx < m ) )
        {
            flops += n;
        }
    }

    flops *= sizeof( CHAMELEON_Complex64_t );

    return flops;
}


int
testing_zlaswp_set_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    //int         api   = parameters_getvalue_int( "api" );
    int         nb    = run_arg_get_int( args, "nb", 320 );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         M     = run_arg_get_int( args, "M", N );
    int         LDA   = run_arg_get_int( args, "LDA", M );
    int         LDB   = run_arg_get_int( args, "LDB", M );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    char       *algo  = parameters_getvalue_str( "algo" );
    /* Descriptors */
    CHAMELEON_Complex64_t *A, *B, *Bref, *dA, *dB;
    int i,j;
    int m0 = 0;
    int *invp, *dinvp;
#ifdef CHAMELEON_USE_CUDA
    cublasHandle_t handle;
#endif

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates two different matrices */
    A = malloc( LDA*N*sizeof(CHAMELEON_Complex64_t) );
    B = calloc( LDB*N, sizeof(CHAMELEON_Complex64_t) );
    Bref = calloc( LDB*N, sizeof(CHAMELEON_Complex64_t) );
    invp = malloc( M*sizeof(int));

    /* Fills A with row index and Bref conversly as invp will just ask to copy in reverse order */
    for (i = 0; i < N; i++)
        for (j = 0; j < M; j++)
            A[j+i*LDA] = (double)j;

    for (i = 0; i < M; i++) {
        invp[i] = (i%10)*M + M - 1 - i;
    }

    for (i = 0; i < N; i++)
        for (j = 0; j < M; j++) {
            int idx = invp[j] - m0;
            if ( ( idx >= 0 ) && (idx < M ) ) {
                Bref[idx+i*LDB] = A[j+i*LDA];
            }
        }

    if ( algo != NULL && ((strcmp(algo, "gpu") == 0) || (strcmp(algo, "gpu_zcopy") == 0) )) {
#ifdef CHAMELEON_USE_CUDA
        cudaMalloc( &dA, LDA*N*sizeof(CHAMELEON_Complex64_t));
        cudaMalloc( &dB, LDB*N*sizeof(CHAMELEON_Complex64_t));
        cudaMalloc( &dinvp, M*sizeof(int));
        cudaMemcpy( dA, A, LDA*N*sizeof(CHAMELEON_Complex64_t), cudaMemcpyHostToDevice );
        cudaMemcpy( dB, B, LDB*N*sizeof(CHAMELEON_Complex64_t), cudaMemcpyHostToDevice );
        cudaMemcpy( dinvp, invp, M*sizeof(int), cudaMemcpyHostToDevice );
        cublasCreate(&handle);
#endif
    }
    testing_start( &test_data );

    if ( algo != NULL && (strcmp(algo, "gpu") == 0)) {
#ifdef CHAMELEON_USE_CUDA
        hres = CUDA_zlaswp_set(  m0, M, N, M,
                                 dA, LDA,
                                 dB, LDB,
                                 dinvp,
                                 handle );
#else
        fprintf(stderr, "GPU algo not available\n");
        return 1;
#endif
    } else if ( algo != NULL && (strcmp(algo, "gpu_zcopy") == 0)) {
#ifdef CHAMELEON_USE_CUDA
        hres = CUDA_zlaswp_set_zcopy(  m0, M, N, M,
                                       dA, LDA,
                                       dB, LDB,
                                       invp,
                                       handle );
#else
        fprintf(stderr, "GPU algo not available\n");
        return 1;
#endif
    } else {
        hres = CORE_zlaswp_set( m0, M, N, M,
                                A, LDA,
                                B, LDB,
                                invp );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlaswp_set( m0, M, N, M, invp ) );

    /* Checks their differences */
    if ( check ) {
        if ( algo != NULL && ((strcmp(algo, "gpu") == 0) || (strcmp(algo, "gpu_zcopy") == 0) )) {
#ifdef CHAMELEON_USE_CUDA
            cudaMemcpy( B, dB, LDB*N*sizeof(CHAMELEON_Complex64_t), cudaMemcpyDeviceToHost );
#endif
        }
        hres += check_zmatrices_std( args, uplo, M, N, B, LDB, Bref, LDB );
    }

    free( A );
    free( B );
    free( Bref );
    free( invp );

    return hres;
}

int
testing_zlaswp_get_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    //int         api   = parameters_getvalue_int( "api" );
    int         nb    = run_arg_get_int( args, "nb", 320 );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpperLower );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         M     = run_arg_get_int( args, "M", N );
    int         LDA   = run_arg_get_int( args, "LDA", M );
    int         LDB   = run_arg_get_int( args, "LDB", M );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAMELEON_Complex64_t *A, *B, *Bref;
    int i,j;
    int m0 = 0;
    int *perm;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates two different matrices */
    A = malloc( LDA*N*sizeof(CHAMELEON_Complex64_t) );
    B = calloc( LDB*N,sizeof(CHAMELEON_Complex64_t) );
    Bref = calloc( LDB*N,sizeof(CHAMELEON_Complex64_t) );
    perm = malloc( M*sizeof(int));

    /* Fills A with row index and Bref conversly as perm will just ask to copy in reverse order */
    for (i = 0; i < N; i++)
        for (j = 0; j < M; j++)
            A[j+i*LDA] = j;
    for (i = 0; i < M; i++)
        perm[i] = (i%10*M) + M - 1 - i;

    for (i = 0; i < N; i++)
        for (j = 0; j < M; j++) {
            int idx = perm[j] - m0;
            if ( ( idx >= 0 ) && (idx < M ) )
            Bref[j+i*LDB] = A[idx+i*LDA];
        }

    /* Makes a copy of descA to descB */
    testing_start( &test_data );

    hres = CORE_zlaswp_get( m0, M, N, M,
                            A, LDA,
                            B, LDB,
                            perm );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlaswp_set( m0, M, N, M, perm ) );

    /* Checks their differences */
    if ( check ) {
        hres += check_zmatrices_std( args, uplo, M, N, B, LDB, Bref, LDB );
    }

    free( A );
    free( B );
    free( Bref );
    free( perm );

    return hres;
}

testing_t   test_zlaswp_set;
const char *zlaswp_set_params[] = { "mtxfmt", "nb", "uplo", "m", "n", "lda", "ldb", "seedA", NULL };
const char *zlaswp_set_output[] = { NULL };
const char *zlaswp_set_outchk[] = { "||A||", "||B||", "RETURN", NULL };

testing_t   test_zlaswp_get;
const char *zlaswp_get_params[] = { "mtxfmt", "nb", "uplo", "m", "n", "lda", "ldb", "seedA", NULL };
const char *zlaswp_get_output[] = { NULL };
const char *zlaswp_get_outchk[] = { "||A||", "||B||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlaswp_set_init( void ) __attribute__( ( constructor ) );
void
testing_zlaswp_set_init( void )
{
    test_zlaswp_set.name   = "zlaswp_set";
    test_zlaswp_set.helper = "General matrix row swapping";
    test_zlaswp_set.params = zlaswp_set_params;
    test_zlaswp_set.output = zlaswp_set_output;
    test_zlaswp_set.outchk = zlaswp_set_outchk;
    test_zlaswp_set.fptr_desc = NULL;
    test_zlaswp_set.fptr_std  = testing_zlaswp_set_std;
    test_zlaswp_set.next   = NULL;

    testing_register( &test_zlaswp_set );

    test_zlaswp_get.name   = "zlaswp_get";
    test_zlaswp_get.helper = "General matrix row swapping";
    test_zlaswp_get.params = zlaswp_get_params;
    test_zlaswp_get.output = zlaswp_get_output;
    test_zlaswp_get.outchk = zlaswp_get_outchk;
    test_zlaswp_get.fptr_desc = NULL;
    test_zlaswp_get.fptr_std  = testing_zlaswp_get_std;
    test_zlaswp_get.next   = NULL;

    testing_register( &test_zlaswp_get );
}
