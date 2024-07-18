/**
 *
 * @file testing_zherk_batch.c
 *
 * @copyright 2019-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zherk_batch testing
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @date 2024-04-03
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <chameleon_lapack.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>
#if defined(CHAMELEON_TESTINGS_VENDOR) || !defined(CHAMELEON_SIMULATION)
#include <coreblas.h>
#endif

static cham_fixdbl_t
flops_zherk_batch( int nb, int K, int N )
{
    return flops_zherk( K, N ) * nb;
}

int
testing_zherk_batch_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int           async = parameters_getvalue_int( "async" );
    int           nb    = run_arg_get_int( args, "nb", 10 );
    int           ib    = run_arg_get_int( args, "ib", 10 );
    cham_trans_t  trans = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_uplo_t   uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int           N     = run_arg_get_int( args, "N", 320 );
    int           K     = run_arg_get_int( args, "K", N );
    double        alpha = testing_dalea();
    double        beta  = testing_dalea();
    double        bump  = 0.;
    int           seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    int           seedC = run_arg_get_int( args, "seedC", testing_ialea() );

    /* Descriptors */
    int          Am, An;
    CHAM_desc_t *descA, *descC;

    alpha = run_arg_get_double( args, "alpha", alpha );
    beta  = run_arg_get_double( args, "beta",  beta  );
    bump  = run_arg_get_double( args, "bump",  0.    );

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
    parameters_desc_create( "A", &descA, ChamComplexDouble, Am, An, nb * Am, ib * An, nb * Am, ib * An );
    parameters_desc_create( "C", &descC, ChamComplexDouble, N,  N,  nb * N,  ib * N,  nb * N,  ib * N  );

    /* Fill the matrices with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplghe_batch_Tile( bump, descC, seedC );

    /* Start measurement */
    testing_start( &test_data );
    if ( async ) {
        fprintf( stderr, "Async unavailable yet\n" );
    }
    hres = CHAMELEON_zherk_batch_Tile( uplo, trans, alpha, descA, beta, descC );

    /* Stop measurement */
    test_data.hres = hres;
    testing_stop( &test_data, flops_zherk_batch( nb*ib, K, N ) );

    /* Check the solution */
    if ( check ) {
        fprintf( stderr, "Check is not available for zherk_batch\n" );
    }

    parameters_desc_destroy( &descA );
    parameters_desc_destroy( &descC );

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
    test_zherk_batch.name      = "zherk_batch";
    test_zherk_batch.helper    = "Perform nb*ib rank-k updates zherk( uplo, trans, N, K, ... )";
    test_zherk_batch.params    = zherk_batch_params;
    test_zherk_batch.output    = zherk_batch_output;
    test_zherk_batch.outchk    = zherk_batch_outchk;
    test_zherk_batch.fptr_desc = testing_zherk_batch_desc;
    test_zherk_batch.fptr_std  = NULL;
    test_zherk_batch.next      = NULL;

    testing_register( &test_zherk_batch );
}
