/**
 *
 * @file testing_zrecgetrf_nopiv.c
 *
 * @copyright 2019-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zrecgetrf_nopiv testing
 *
 * @version 1.2.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2022-02-22
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zrecgetrf_nopiv_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int        async  = parameters_getvalue_int( "async" );
    intptr_t   mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int        nb     = run_arg_get_int( args, "nb", 320 );
    int        l1     = run_arg_get_int( args, "l1", 0 );
    int        l2     = run_arg_get_int( args, "l2", 0 );
    int        l3     = run_arg_get_int( args, "l3", 0 );
    cham_rec_t rec    = run_arg_get_rec( args, "rec", ChamRecFull );
    int        rarg   = run_arg_get_int( args, "rarg", 1 );
    int        P      = parameters_getvalue_int( "P" );
    int        N      = run_arg_get_int( args, "N", 1000 );
    int        M      = run_arg_get_int( args, "M", N );
    int        LDA    = run_arg_get_int( args, "LDA", M );
    int        seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    double     bump   = run_arg_get_double( args, "bump", (double)N );
    int        Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t *descA;
    int list_nb[] = { nb, l1, l2, l3, 0 };

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    CHAMELEON_Recursive_Desc_Create(
        &descA, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble,
        rec, rarg, list_nb, list_nb, N, N, N, N,
        P, Q, NULL, NULL, NULL, NULL, "A" );

    /* Fills the matrix with random values */
    CHAMELEON_zplgtr_Tile( 0,    ChamUpper, descA, seedA   );
    CHAMELEON_zplgtr_Tile( bump, ChamLower, descA, seedA+1 );

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgetrf_nopiv_Tile_Async( descA, test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgetrf_nopiv_Tile( descA );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgetrf( M, N ) );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, NULL );
        CHAMELEON_zplgtr_Tile( 0,    ChamUpper, descA0, seedA   );
        CHAMELEON_zplgtr_Tile( bump, ChamLower, descA0, seedA+1 );

        hres += check_zxxtrf( args, ChamGeneral, ChamUpperLower, descA0, descA );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    CHAMELEON_Desc_Destroy( &descA );

    return hres;
}

int
testing_zrecgetrf_nopiv_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int    nb    = run_arg_get_int( args, "nb", 320 );
    int    N     = run_arg_get_int( args, "N", 1000 );
    int    M     = run_arg_get_int( args, "M", N );
    int    LDA   = run_arg_get_int( args, "LDA", M );
    int    seedA = run_arg_get_int( args, "seedA", testing_ialea() );
    double bump  = run_arg_get_double( args, "bump", (double)N );

    /* Descriptors */
    CHAMELEON_Complex64_t *A;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    A = malloc( LDA*N*sizeof(CHAMELEON_Complex64_t) );

    /* Fills the matrix with random values */
    CHAMELEON_zplgtr( 0,    ChamUpper, M, N, A, LDA, seedA   );
    CHAMELEON_zplgtr( bump, ChamLower, M, N, A, LDA, seedA+1 );

    /* Calculates the solution */
    testing_start( &test_data );
    hres = CHAMELEON_zgetrf_nopiv( M, N, A, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgetrf( M, N ) );

    /* Checks the factorisation and residue */
    if ( check ) {
        CHAMELEON_Complex64_t *A0 = malloc( LDA*N*sizeof(CHAMELEON_Complex64_t) );
        CHAMELEON_zplgtr( 0,    ChamUpper, M, N, A0, LDA, seedA   );
        CHAMELEON_zplgtr( bump, ChamLower, M, N, A0, LDA, seedA+1 );

        hres += check_zxxtrf_std( args, ChamGeneral, ChamUpperLower, M, N, A0, A, LDA );

        free( A0 );
    }

    free( A );

    return hres;
}

testing_t   test_zrecgetrf_nopiv;
const char *zrecgetrf_nopiv_params[] = { "mtxfmt", "nb", "l1", "l2", "l3", "rec", "rarg", "m", "n", "lda", "seedA", "bump", NULL };
const char *zrecgetrf_nopiv_output[] = { NULL };
const char *zrecgetrf_nopiv_outchk[] = { "||A||", "||A-fact(A)||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zrecgetrf_nopiv_init( void ) __attribute__( ( constructor ) );
void
testing_zrecgetrf_nopiv_init( void )
{
    test_zrecgetrf_nopiv.name   = "zrecgetrf_nopiv";
    test_zrecgetrf_nopiv.helper = "General factorization (LU without pivoting)";
    test_zrecgetrf_nopiv.params = zrecgetrf_nopiv_params;
    test_zrecgetrf_nopiv.output = zrecgetrf_nopiv_output;
    test_zrecgetrf_nopiv.outchk = zrecgetrf_nopiv_outchk;
    test_zrecgetrf_nopiv.fptr_desc = testing_zrecgetrf_nopiv_desc;
    test_zrecgetrf_nopiv.fptr_std  = testing_zrecgetrf_nopiv_std;
    test_zrecgetrf_nopiv.next   = NULL;

    testing_register( &test_zrecgetrf_nopiv );
}
