/**
 *
 * @file testing_zrecpotrf.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zrecpotrf testing
 *
 * @version 0.9.2
 * @author Lucas Barros de Assis
 * @date 2019-08-12
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zrecpotrf( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Reads arguments */
    int         async  = parameters_getvalue_int( "async" );
    int         nb    = run_arg_get_int( args, "nb", 320 );
    int         l1    = run_arg_get_int( args, "l1", 0 );
    int         l2    = run_arg_get_int( args, "l2", 0 );
    int         l3    = run_arg_get_int( args, "l3", 0 );
    cham_rec_t  rec   = run_arg_get_rec( args, "rec", ChamRecFull );
    int         rarg  = run_arg_get_int( args, "rarg", 1 );
    int         P     = parameters_getvalue_int( "P" );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         LDA   = run_arg_get_int( args, "LDA", N );
    int         seedA = run_arg_get_int( args, "seedA", random() );
    int         Q     = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA;
    int list_nb[] = { nb, l1, l2, l3, 0 };

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    CHAMELEON_Recursive_Desc_Create(
        &descA, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble,
        rec, rarg, list_nb, list_nb, N, N, N, N,
        P, Q, NULL, NULL, NULL, NULL, "A" );

    CHAMELEON_Recursive_Desc_Partition_Submit( descA );

    /* Fills the matrix with random values */
    CHAMELEON_zplghe_Tile( (double)N, uplo, descA, seedA );

    /* Calculates the solution */
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

    CHAMELEON_Recursive_Desc_Unpartition_Submit( descA );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAM_desc_t *descA0;// = CHAMELEON_Desc_Copy( descA, NULL );
        CHAMELEON_Desc_Create(
            &descA0, NULL, ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, N, N, P, Q );
        CHAMELEON_zplghe_Tile( (double)N, uplo, descA0, seedA );
        hres += check_zxxtrf( args, ChamHermitian, uplo, descA0, descA );
        CHAMELEON_Desc_Destroy( &descA0 );
    }

    CHAMELEON_Desc_Destroy( &descA );

    return hres;
}

testing_t   test_zrecpotrf;
const char *zrecpotrf_params[] = { "nb", "l1", "l2", "l3", "rec", "rarg", "uplo", "n", "lda", "seedA", NULL };
const char *zrecpotrf_output[] = { NULL };
const char *zrecpotrf_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zrecpotrf_init( void ) __attribute__( ( constructor ) );
void
testing_zrecpotrf_init( void )
{
    test_zrecpotrf.name        = "zrecpotrf";
    test_zrecpotrf.helper      = "zrecpotrf";
    test_zrecpotrf.params      = zrecpotrf_params;
    test_zrecpotrf.output      = zrecpotrf_output;
    test_zrecpotrf.outchk      = zrecpotrf_outchk;
    test_zrecpotrf.fptr_desc   = testing_zrecpotrf;
    test_zrecpotrf.fptr_std    = NULL;
    test_zrecpotrf.next        = NULL;

    testing_register( &test_zrecpotrf );
}
