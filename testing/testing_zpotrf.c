/**
 *
 * @file testing_zpotrf.c
 *
 * @copyright 2019-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zpotrf testing
 *
 * @version 1.3.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Florent Pruvost
 * @date 2023-07-05
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
#include <coreblas/lapacke.h>
#endif
#endif
//#include <papi.h>
//#include "power_measurement.h"

#if !defined(CHAMELEON_TESTINGS_VENDOR)
int
testing_zpotrf_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int         async = parameters_getvalue_int( "async" );
    int         nb    = run_arg_get_nb( args );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAM_desc_t *descA;

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

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, N, N );

    /* Fills the matrix with random values */
    CHAMELEON_zplghe_Tile( (double)N, uplo, descA, seedA );

    /* Calculates the solution */
    //   retval = PAPI_start( EventSet );
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zpotrf_Tile_Async( uplo, descA,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zpotrf_Tile( uplo, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpotrf( N ) );

    /* retval = PAPI_stop( EventSet, values ); */
    /* if (retval != PAPI_OK) { */
    /*     perror("unable to papi stop"); */
    /*     exit(1); */
    /* } */

    /* for( int s = 0 ; s < N_SOCK ; s ++){ */
    /*     for( int i = 0 ; i < N_EVTS; i++) { */
    /*         printf("%-40s%12.6f J\t(Average Power %.1fW)\n", */
    /*                event_names[i], */
    /*                (double)values[s * N_EVTS + i]/1.0e9, */
    /*                ((double)values[s * N_EVTS + i]/1.0e9)/t); */
    /*     } */
    /* } */

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplghe_Tile( (double)N, ChamUpperLower, descA0, seedA );

        hres += check_zxxtrf( args, ChamHermitian, uplo, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    parameters_desc_destroy( &descA );

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

    return hres;
}
#endif

int
testing_zpotrf_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int         api   = parameters_getvalue_int( "api" );
    int         nb    = run_arg_get_nb( args );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         seedA = run_arg_get_int( args, "seedA", testing_ialea() );

    /* Descriptors */
    CHAMELEON_Complex64_t *A;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    A = malloc( LDA*N*sizeof(CHAMELEON_Complex64_t) );

    /* Fills the matrix with random values */
    CHAMELEON_zplghe( (double)N, uplo, N, A, LDA, seedA );

    /* Calculates the solution */
#if defined(CHAMELEON_TESTINGS_VENDOR)
    testing_start( &test_data );
    hres = LAPACKE_zpotrf( LAPACK_COL_MAJOR, chameleon_lapack_const(uplo), N, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpotrf( N ) );
#else
    testing_start( &test_data );
    switch ( api ) {
    case 1:
        hres = CHAMELEON_zpotrf( uplo, N, A, LDA );
        break;
#if !defined(CHAMELEON_SIMULATION)
    case 2:
        CHAMELEON_lapacke_zpotrf( CblasColMajor, chameleon_lapack_const(uplo), N, A, LDA );
        break;
#endif
    default:
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr,
                     "SKIPPED: This function can only be used with the option --api 1 or --api 2.\n" );
        }
        return -1;
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zpotrf( N ) );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAMELEON_Complex64_t *A0 = malloc( LDA*N*sizeof(CHAMELEON_Complex64_t) );
        CHAMELEON_zplghe( (double)N, uplo, N, A0, LDA, seedA );

        hres += check_zxxtrf_std( args, ChamHermitian, uplo, N, N, A0, A, LDA );

        free( A0 );
    }
#endif

    free( A );

    (void)check;
    return hres;
}

testing_t   test_zpotrf;
#if defined(CHAMELEON_TESTINGS_VENDOR)
const char *zpotrf_params[] = { "uplo", "n", "lda", "seedA", NULL };
#else
const char *zpotrf_params[] = { "mtxfmt", "nb", "uplo", "n", "lda", "seedA", NULL };
#endif
const char *zpotrf_output[] = { NULL };
const char *zpotrf_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zpotrf_init( void ) __attribute__( ( constructor ) );
void
testing_zpotrf_init( void )
{
    test_zpotrf.name   = "zpotrf";
    test_zpotrf.helper = "Hermitian positive definite factorization (Cholesky)";
    test_zpotrf.params = zpotrf_params;
    test_zpotrf.output = zpotrf_output;
    test_zpotrf.outchk = zpotrf_outchk;
#if defined(CHAMELEON_TESTINGS_VENDOR)
    test_zpotrf.fptr_desc = NULL;
#else
    test_zpotrf.fptr_desc = testing_zpotrf_desc;
#endif
    test_zpotrf.fptr_std  = testing_zpotrf_std;
    test_zpotrf.next   = NULL;

    testing_register( &test_zpotrf );
}
