/**
 *
 * @file testing_ztrsm_batch.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztrsm_batch testing
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
flops_ztrsm_batch( int nb, cham_side_t side, int M, int N )
{
    return flops_ztrsm( side, M, N ) * nb;
}

int
testing_ztrsm_batch( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    test_data.task = TASK_TRSM;

    int          Am, An;
    int          hres = 0;
    CHAM_desc_t *descA, *descB;

    /* Read arguments */
    int          nb     = run_arg_get_int( args, "nb", 10 );
    int          ib     = run_arg_get_int( args, "ib", 10 );
    int          P      = parameters_getvalue_int( "P" );
    cham_trans_t trans  = run_arg_get_trans( args, "trans", ChamNoTrans );
    cham_side_t  side   = run_arg_get_side( args, "side", ChamLeft );
    cham_uplo_t  uplo   = run_arg_get_uplo( args, "uplo", ChamUpper );
    cham_diag_t  diag   = run_arg_get_diag( args, "diag", ChamNonUnit );
    int          N      = run_arg_get_int( args, "N", 320 );
    int          M      = run_arg_get_int( args, "M", N );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", random() );
    int                   seedB = run_arg_get_int( args, "seedB", random() );
    int                   Q     = parameters_compute_q( P );
    cham_fixdbl_t t, gflops;
    cham_fixdbl_t flops = flops_ztrsm_batch( nb*ib, side, M, N );
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

    alpha = run_arg_get_complex64( args, "alpha", alpha );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, N );

    /* Calculate the dimensions according to the transposition */
    if ( side == ChamLeft ) {
        Am = M;
        An = M;
    }
    else {
        Am = N;
        An = N;
    }

    /* Create the matrices */
    CHAMELEON_Desc_Create(
        &descA, NULL, ChamComplexDouble, Am, An, Am * An,
        nb * Am, ib * An, 0, 0,
        nb * Am, ib * An, P, Q );
    CHAMELEON_Desc_Create(
        &descB, NULL, ChamComplexDouble, M, N, M * N,
        nb * M,  ib * N,  0, 0,
        nb * M,  ib * N,  P, Q );

    /* Fill the matrices with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );

    /* Calculate the product */
    //    retval = PAPI_start( EventSet );
    
    /* Start measurement */
    testing_start( &test_data );

    hres = CHAMELEON_ztrsm_batch_Tile( side, uplo, trans, diag, alpha, descA, descB );
    
    /* Stop measurement */
    testing_stop( &test_data, flops);

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
    CHAMELEON_Desc_Destroy( &descB );

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

testing_t   test_ztrsm_batch;
const char *ztrsm_batch_params[] = { "nb", "ib", "side",  "uplo", "trans", "diag", "m", "n",
                                     "alpha", "seedA", "seedB", NULL };
const char *ztrsm_batch_output[] = { NULL };
const char *ztrsm_batch_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_ztrsm_batch_init( void ) __attribute__( ( constructor ) );
void
testing_ztrsm_batch_init( void )
{
    test_ztrsm_batch.name        = "ztrsm_batch";
    test_ztrsm_batch.helper      = "Perform nb*ib triangular solve trsm( side, uplo, trns, diag, M, N, ... )";
    test_ztrsm_batch.params      = ztrsm_batch_params;
    test_ztrsm_batch.output      = ztrsm_batch_output;
    test_ztrsm_batch.outchk      = ztrsm_batch_outchk;
    test_ztrsm_batch.fptr_desc   = testing_ztrsm_batch;
    test_ztrsm_batch.next        = NULL;

    testing_register( &test_ztrsm_batch );
}
