/**
 *
 * @file testing_zherk_batch.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherk_batch testing
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
#include "runtime_energy.h"
//#include "power_measurement.h"

static cham_fixdbl_t
flops_zherk_batch( int nb, int K, int N )
{
    return flops_zherk( K, N ) * nb;
}

int
testing_zherk_batch( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    test_data.task = TASK_HERK;

    int          Am, An;
    int          hres = 0;
    CHAM_desc_t *descA, *descC;

    /* Read arguments */
    int           nb    = run_arg_get_int( args, "nb", 10 );
    int           ib    = run_arg_get_int( args, "ib", 10 );
    int           P     = parameters_getvalue_int( "P" );
    cham_trans_t  trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_uplo_t   uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int           N     = run_arg_get_int( args, "N", 320 );
    int           K     = run_arg_get_int( args, "K", N );
    double        alpha = testing_dalea();
    double        beta  = testing_dalea();
    double        bump  = testing_dalea();
    int           seedA = run_arg_get_int( args, "seedA", random() );
    int           seedC = run_arg_get_int( args, "seedC", random() );
    int           Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zherk_batch( nb*ib, K, N );
    int energy = parameters_getvalue_int( "energy" );

    alpha = run_arg_get_double( args, "alpha", alpha );
    beta  = run_arg_get_double( args, "beta", beta );
    bump  = run_arg_get_double( args, "bump", bump );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, N );

    /* Calculate the dimensions according to the transposition */
    if ( trans == ChamNoTrans ) {
        Am = N;
        An = K;
    }
    else {
        Am = K;
        An = N;
    }

    /* Create the matrices */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, Am, An, Am * An,
        nb * Am, ib * An, 0, 0,
        nb * Am, ib * An, P, Q );
    CHAMELEON_Desc_Create(
        &descC, NULL, ChamComplexDouble, N, N, N * N,
        nb * N,  ib * N,  0, 0,
        nb * N,  ib * N,  P, Q );

    /* Fill the matrices with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplghe_batch_Tile( bump, descC, seedC );
     
    /* Start measurement */
    testing_start( &test_data );

    hres = CHAMELEON_zherk_batch_Tile( uplo, trans, alpha, descA, beta, descC );
    
    /* Stop measurement */
    testing_stop( &test_data, flops);

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descC );

    (void)check;
    return hres;
}

testing_t   test_zherk_batch;
const char *zherk_batch_params[] = { "nb", "ib", "trans", "uplo", "n", "k",
                                     "alpha", "beta", "seedA", "seedC", "bump", NULL };
const char *zherk_batch_output[] = { NULL };
const char *zherk_batch_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zherk_batch_init( void ) __attribute__( ( constructor ) );
void
testing_zherk_batch_init( void )
{
    test_zherk_batch.name        = "zherk_batch";
    test_zherk_batch.helper      = "Perform nb*ib rank-k updates zherk( uplo, trans, N, K, ... )";
    test_zherk_batch.params      = zherk_batch_params;
    test_zherk_batch.output      = zherk_batch_output;
    test_zherk_batch.outchk      = zherk_batch_outchk;
    test_zherk_batch.fptr_desc   = testing_zherk_batch;
    test_zherk_batch.next        = NULL;

    testing_register( &test_zherk_batch );
}
