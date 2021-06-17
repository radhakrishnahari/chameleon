/**
 *
 * @file testing_zgemm_batch.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm_batch testing
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>
#include "runtime_codelet_z.h"
#include "runtime_codelets.h"
#include "runtime_energy.h"


static cham_fixdbl_t
flops_zgemm_batch( int nb, int M, int N, int K )
{
    return flops_zgemm( M, N, K ) * nb;
}

int
testing_zgemm_batch( run_arg_list_t *args, int check )
{
    int          Am, An, Bm, Bn;
    int          hres = 0;
    CHAM_desc_t *descA, *descB, *descC;

    /* Read arguments */
    int          nb     = run_arg_get_int( args, "nb", 10 );
    int          ib     = run_arg_get_int( args, "ib", 10 );
    int          P      = parameters_getvalue_int( "P" );
    cham_trans_t transA = run_arg_get_trans( args, "transA", ChamNoTrans );
    cham_trans_t transB = run_arg_get_trans( args, "transB", ChamNoTrans );
    int          N      = run_arg_get_int( args, "N", 320 );
    int          M      = run_arg_get_int( args, "M", N );
    int          K      = run_arg_get_int( args, "K", N );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", random() );
    int                   seedB = run_arg_get_int( args, "seedB", random() );
    int                   seedC = run_arg_get_int( args, "seedC", random() );
    int                   Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zgemm_batch( nb*ib, M, N, K );
    int energy = parameters_getvalue_int( "energy" );

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta",  beta  );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, N );

    /* Calculate the dimensions according to the transposition */
    if ( transA == ChamNoTrans ) {
        Am = M;
        An = K;
    }
    else {
        Am = K;
        An = M;
    }
    if ( transB == ChamNoTrans ) {
        Bm = K;
        Bn = N;
    }
    else {
        Bm = N;
        Bn = K;
    }

    /* Create the matrices */
    CHAMELEON_Desc_Create(
		    &descA, NULL, ChamComplexDouble, Am, An, Am * An,
		    nb * Am, ib * An, 0, 0,
		    nb * Am, ib * An, P, Q );
    CHAMELEON_Desc_Create(
		    &descB, NULL, ChamComplexDouble, Bm, Bn, Bm * Bn,
		    nb * Bm, ib * Bn, 0, 0,
		    nb * Bm, ib * Bn, P, Q );
    CHAMELEON_Desc_Create(
		    &descC, NULL, ChamComplexDouble, M, N, M * N,
		    nb * M,  ib * N,  0, 0,
		    nb * M,  ib * N,  P, Q );

    /* Fill the matrices with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );
    CHAMELEON_zplrnt_Tile( descC, seedC );

    /*Start energy measurement*/
    if ( energy ) {
        RUNTIME_start_energy();
    }

    START_TIMING( t );
    hres = CHAMELEON_zgemm_batch_Tile( transA, transB, alpha, descA, descB, beta, descC );
    STOP_TIMING( t );

    /*Stop energy measurement*/
    if ( energy ) {
        RUNTIME_stop_energy( ChamComplexDouble, TASK_GEMM );
    }

    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );
    CHAMELEON_Desc_Destroy( &descC );

    (void)check;
    return hres;
}

testing_t   test_zgemm_batch;
const char *zgemm_batch_params[] = { "nb", "ib", "transA", "transB", "m", "n", "k",
                                     "alpha", "beta", "seedA", "seedB", "seedC", NULL };
const char *zgemm_batch_output[] = { NULL };
const char *zgemm_batch_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zgemm_batch_init( void ) __attribute__( ( constructor ) );
void
testing_zgemm_batch_init( void )
{
    test_zgemm_batch.name   = "zgemm_batch";
    test_zgemm_batch.helper = "Perform nb*ib general matrix-matrix multiply of size MxNxK";
    test_zgemm_batch.params = zgemm_batch_params;
    test_zgemm_batch.output = zgemm_batch_output;
    test_zgemm_batch.outchk = zgemm_batch_outchk;
    test_zgemm_batch.fptr   = testing_zgemm_batch;
    test_zgemm_batch.next   = NULL;

    testing_register( &test_zgemm_batch );
}
