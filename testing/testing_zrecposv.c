/**
 *
 * @file testing_zrecposv.c
 *
 * @copyright 2019-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zrecposv testing
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

static cham_fixdbl_t
flops_zposv( int N, int NRHS )
{
    cham_fixdbl_t flops = flops_zpotrf( N ) + flops_zpotrs( N, NRHS );
    return flops;
}

int
testing_zrecposv( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Reads arguments */
    int         async  = parameters_getvalue_int( "async" );
    /* intptr_t mtxfmt    = parameters_getvalue_int( "mtxfmt" ); */
    int         nb     = run_arg_get_int( args, "nb", 320 );
    int         l1     = run_arg_get_int( args, "l1", 0 );
    int         l2     = run_arg_get_int( args, "l2", 0 );
    int         l3     = run_arg_get_int( args, "l3", 0 );
    cham_rec_t  rec    = run_arg_get_rec( args, "rec", ChamRecFull );
    int         rarg   = run_arg_get_int( args, "rarg", 1 );
    int         P      = parameters_getvalue_int( "P" );
    cham_uplo_t uplo   = run_arg_get_uplo( args, "uplo", ChamUpper );
    int         N      = run_arg_get_int( args, "N", 1000 );
    int         NRHS   = run_arg_get_int( args, "NRHS", 1 );
    int         LDA    = run_arg_get_int( args, "LDA", N );
    int         LDB    = run_arg_get_int( args, "LDB", N );
    int         seedA  = run_arg_get_int( args, "seedA", random() );
    int         seedB  = run_arg_get_int( args, "seedB", random() );
    int         Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA, *descX;
    int list_nb[] = { nb, l1, l2, l3, 0 };

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    CHAMELEON_Recursive_Desc_Create(
        &descA, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble, rec, rarg, list_nb, list_nb,
        LDA, N, N, N, P, Q, NULL, NULL, NULL, NULL, "A" ); /* keeping the macro for now */
    CHAMELEON_Recursive_Desc_Create(
        &descX, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble, rec, rarg, list_nb, list_nb,
        LDB, NRHS, N, NRHS, P, Q, NULL, NULL, NULL, NULL, "A" ); /* keeping the macro for now */

    CHAMELEON_Recursive_Desc_Partition_Submit( descA );
    CHAMELEON_Recursive_Desc_Partition_Submit( descX );

    /* Fills the matrix with random values */
    CHAMELEON_zplghe_Tile( (double)N, uplo, descA, seedA );
    CHAMELEON_zplrnt_Tile( descX, seedB );

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zposv_Tile_Async( uplo, descA, descX,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descX, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zposv_Tile( uplo, descA, descX );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zposv( N, NRHS ) );

    CHAMELEON_Recursive_Desc_Unpartition_Submit( descA );
    CHAMELEON_Recursive_Desc_Unpartition_Submit( descX );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAM_desc_t *descA0, *descB;

        /* Check the factorization */
        descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAMELEON_zplghe_Tile( (double)N, uplo, descA0, seedA );

        hres += check_zxxtrf( args, ChamHermitian, uplo, descA0, descA );

        /* Check the solve */
        descB = CHAMELEON_Desc_Copy( descX, NULL );
        CHAMELEON_zplrnt_Tile( descB, seedB );

        CHAMELEON_zplghe_Tile( (double)N, uplo, descA0, seedA );
        hres += check_zsolve( args, ChamHermitian, ChamNoTrans, uplo, descA0, descX, descB );

        CHAMELEON_Desc_Destroy( &descA0 );
        CHAMELEON_Desc_Destroy( &descB );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descX );

    return hres;
}

testing_t   test_zrecposv;
const char *zrecposv_params[] = { "nb", "l1", "l2", "l3", "rec", "rarg", "uplo", "n", "lda", "seedA", NULL };
const char *zrecposv_output[] = { NULL };
const char *zrecposv_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zrecposv_init( void ) __attribute__( ( constructor ) );
void
testing_zrecposv_init( void )
{
    test_zrecposv.name        = "zrecposv";
    test_zrecposv.helper      = "zrecposv";
    test_zrecposv.params      = zrecposv_params;
    test_zrecposv.output      = zrecposv_output;
    test_zrecposv.outchk      = zrecposv_outchk;
    test_zrecposv.fptr_desc   = testing_zrecposv;
    test_zrecposv.fptr_std    = NULL;
    test_zrecposv.next        = NULL;

    testing_register( &test_zrecposv );
}
