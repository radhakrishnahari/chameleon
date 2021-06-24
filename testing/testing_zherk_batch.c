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

    /* PAPI variables*/
    /* int EventSet = PAPI_NULL; */
    /* long long *values; */
    /* int retval; */
    /* values=calloc(N_SOCK * N_EVTS,sizeof(long long)); */
    /* if (values==NULL) { */
    /*     exit(1); */
    /* } */

    /* if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) { */
    /*     perror("unable to initialize PAPI"); */
    /*     exit(1); */
    /* } */

    /* retval = PAPI_create_eventset( &EventSet ); */
    /* if (retval != PAPI_OK) { */
    /*     perror("unable to create eventSet"); */
    /*     exit(1); */
    /* } */

    /* for (int i = 0 ; i < N_SOCK ; i ++ ) */
    /*     add_event(EventSet, i); */

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

    /* Calculate the product */
    //    retval = PAPI_start( EventSet );
     /*Start energy measurement*/
    if ( energy ) {
        RUNTIME_start_energy();
    }

    START_TIMING( t );
    hres = CHAMELEON_zherk_batch_Tile( uplo, trans, alpha, descA, beta, descC );
    STOP_TIMING( t );

    if ( energy ) {
        /* the last parameter is side, which is not used for gemm, let's set it to -1*/
        RUNTIME_stop_energy( ChamComplexDouble, TASK_HERK, -1 );
    }

    /* retval = PAPI_stop( EventSet, values ); */
    /* if (retval != PAPI_OK) { */
    /*     perror("unable to papi stop"); */
    /*     exit(1); */
    /* } */
    gflops = flops * 1.e-9 / t;
    run_arg_add_fixdbl( args, "time", t );
    run_arg_add_fixdbl( args, "gflops", ( hres == CHAMELEON_SUCCESS ) ? gflops : -1. );

    /* for( int s = 0 ; s < N_SOCK ; s ++){ */
    /*     for( int i = 0 ; i < N_EVTS; i++) { */
    /*         printf("%-40s%12.6f J\t(Average Power %.1fW)\n", */
    /*                event_names[i], */
    /*                (double)values[s * N_EVTS + i]/1.0e9, */
    /*                ((double)values[s * N_EVTS + i]/1.0e9)/t); */
    /*     } */
    /* } */

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descC );

    /* retval = PAPI_cleanup_eventset( EventSet ); */
    /* if (retval != PAPI_OK) { */
    /*     perror("unable to cleanup"); */
    /*     exit(1); */
    /* } */

    /* retval = PAPI_destroy_eventset( &EventSet ); */
    /* if (retval != PAPI_OK) { */
    /*     perror("unable to destroy eventset"); */
    /*     exit(1); */
    /* } */

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
    test_zherk_batch.name   = "zherk_batch";
    test_zherk_batch.helper = "Perform nb*ib rank-k updates zherk( uplo, trans, N, K, ... )";
    test_zherk_batch.params = zherk_batch_params;
    test_zherk_batch.output = zherk_batch_output;
    test_zherk_batch.outchk = zherk_batch_outchk;
    test_zherk_batch.fptr   = testing_zherk_batch;
    test_zherk_batch.next   = NULL;

    testing_register( &test_zherk_batch );
}
