/**
 *
 * @file testing_zrecursive.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zrecursive testing
 *
 * @version 0.9.2
 * @author Lucas Barros de Assis
 * @date 2019-07-04
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

static cham_fixdbl_t
flops_zrecursive( cham_uplo_t uplo, int M, int N )
{
    cham_fixdbl_t flops;

    switch ( uplo ) {
        case ChamUpper:
            if ( N > M ) {
                flops = ( M * ( M + 1 ) / 2 ) + M * ( N - M );
            }
            else {
                flops = N * ( N + 1 ) / 2;
            }
            break;
        case ChamLower:
            if ( M > N ) {
                flops = ( N * ( N + 1 ) / 2 ) + N * ( M - N );
            }
            else {
                flops = M * ( M + 1 ) / 2;
            }
            break;
        case ChamUpperLower:
        default:
            flops = M * N;
    }
    flops *= sizeof( CHAMELEON_Complex64_t );

    return flops;
}

int
testing_zrecursive( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Reads arguments */
    int         async = parameters_getvalue_int( "async" );
    int         nb    = run_arg_get_int( args, "nb", 320 );
    int         l1    = run_arg_get_int( args, "l1", 0 );
    int         l2    = run_arg_get_int( args, "l2", 0 );
    int         l3    = run_arg_get_int( args, "l3", 0 );
    cham_rec_t  rec   = run_arg_get_rec( args, "rec", ChamRecFull );
    int         rarg  = run_arg_get_rec( args, "rarg", 1 );
    int         P     = parameters_getvalue_int( "P" );
    cham_uplo_t uplo  = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N     = run_arg_get_int( args, "N", 1000 );
    int         M     = run_arg_get_int( args, "M", N );
    int         LDA   = run_arg_get_int( args, "LDA", M );
    int         LDB   = run_arg_get_int( args, "LDB", M );
    int         seedA = run_arg_get_int( args, "seedA", random() );
    int         Q     = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA, *descB;
    int list_nb[] = { nb, l1, l2, l3, 0 };

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates two different matrices */
    CHAMELEON_Recursive_Desc_Create(
        &descA, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble, rec, rarg, list_nb, list_nb,
        LDA, N, M, N, P, Q, NULL, NULL, NULL, NULL, "A" );
    CHAMELEON_Recursive_Desc_Create(
        &descB, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble, rec, rarg, list_nb, list_nb,
        LDB, N, M, N, P, Q, NULL, NULL, NULL, NULL, "B" );

    /* Fills each matrix with different random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    /* We use seedA + 1, just to create a variation in B */
    CHAMELEON_zplrnt_Tile( descB, seedA + 1 );

    /* Makes a copy of descA to descB */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zlacpy_Tile_Async( uplo, descB, descA,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descB, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zlacpy_Tile( uplo, descB, descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zrecursive( uplo, M, N ) );

    /* Checks their differences */
    if ( check ) {
        hres += check_zmatrices( args, uplo, descA, descB );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );

    return hres;
}

testing_t   test_zrecursive;
const char *zrecursive_params[] = { "nb", "uplo", "m", "n", "lda", "ldb", "seedA", NULL };
const char *zrecursive_output[] = { NULL };
const char *zrecursive_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zrecursive_init( void ) __attribute__( ( constructor ) );
void
testing_zrecursive_init( void )
{
    test_zrecursive.name        = "zrecursive";
    test_zrecursive.helper      = "zrecursive";
    test_zrecursive.params      = zrecursive_params;
    test_zrecursive.output      = zrecursive_output;
    test_zrecursive.outchk      = zrecursive_outchk;
    test_zrecursive.fptr_desc   = testing_zrecursive;
    test_zrecursive.fptr_std    = NULL;
    test_zrecursive.next        = NULL;

    testing_register( &test_zrecursive );
}
