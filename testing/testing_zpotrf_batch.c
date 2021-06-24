/**
 *
 * @file testing_zpotrf_batch.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotrf_batch testing
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
flops_zpotrf_batch( int nb, int N )
{
    return flops_zpotrf( N ) * nb;
}

int
testing_zpotrf_batch( run_arg_list_t *args, int check )
{
    int          hres   = 0;
    CHAM_desc_t *descA;

    /* Read arguments */
    int           nb    = run_arg_get_int( args, "nb", 10 );
    int           ib    = run_arg_get_int( args, "ib", 10 );
    int           P     = parameters_getvalue_int( "P" );
    cham_uplo_t   uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int           N     = run_arg_get_int( args, "N", 320 );
    int           seedA = run_arg_get_int( args, "seedA", random() );
    int           Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_zpotrf_batch( nb*ib, N );
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


    CHAMELEON_Set( CHAMELEON_TILE_SIZE, N );

    /* Create the matrices */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, N, N, N * N,
        nb * N,  ib * N,  0, 0,
        nb * N,  ib * N,  P, Q );

    /* Fill the matrices with random values */
    CHAMELEON_zplghe_batch_Tile( (double)N, descA, seedA );

    /* Calculate the product */
    //    retval = PAPI_start( EventSet );
    if ( energy ) {
        RUNTIME_start_energy();
    }

    START_TIMING( t );
    hres = CHAMELEON_zpotrf_batch_Tile( uplo, descA );
    STOP_TIMING( t );

     if ( energy ) {
         RUNTIME_stop_energy( ChamComplexDouble, TASK_POTRF);
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

testing_t   test_zpotrf_batch;
const char *zpotrf_batch_params[] = { "nb", "ib",  "uplo", "n", "seedA", NULL };
const char *zpotrf_batch_output[] = { NULL };
const char *zpotrf_batch_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zpotrf_batch_init( void ) __attribute__( ( constructor ) );
void
testing_zpotrf_batch_init( void )
{
    test_zpotrf_batch.name   = "zpotrf_batch";
    test_zpotrf_batch.helper = "Perform nb*ib Cholesky factorization potrf( uplo, N, ... )";
    test_zpotrf_batch.params = zpotrf_batch_params;
    test_zpotrf_batch.output = zpotrf_batch_output;
    test_zpotrf_batch.outchk = zpotrf_batch_outchk;
    test_zpotrf_batch.fptr   = testing_zpotrf_batch;
    test_zpotrf_batch.next   = NULL;

    testing_register( &test_zpotrf_batch );
}
