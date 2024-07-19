/**
 *
 * @file testing_zungqr_hqr.c
 *
 * @copyright 2019-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zungqr_hqr testing
 *
 * @version 1.3.0
 * @author Lucas Barros de Assis
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2023-07-05
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include <chameleon_lapack.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zungqr_hqr_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int      async  = parameters_getvalue_int( "async" );
    intptr_t mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int      nb     = run_arg_get_nb( args );
    int      ib     = run_arg_get_ib( args );
    int      P      = parameters_getvalue_int( "P" );
    int      N      = run_arg_get_int( args, "N", 1000 );
    int      M      = run_arg_get_int( args, "M", N );
    int      K      = run_arg_get_int( args, "K", chameleon_min( M, N ) );
    int      LDA    = run_arg_get_int( args, "LDA", M );
    int      qr_a   = run_arg_get_int( args, "qra", -1 );
    int      qr_p   = run_arg_get_int( args, "qrp", -1 );
    int      llvl   = run_arg_get_int( args, "llvl", -1 );
    int      hlvl   = run_arg_get_int( args, "hlvl", -1 );
    int      domino = run_arg_get_int( args, "domino", -1 );
    int      seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int      Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAM_desc_t    *descA, *descTS, *descTT, *descQ, *D = NULL;
    libhqr_tree_t   qrtree;
    libhqr_matrix_t matrix;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( N > M ) {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr, "SKIPPED: Incorrect parameters for ungqr_hqr (N > M)\n" );
        }
        return -1;
    }

    if ( K > N ) {
        if ( CHAMELEON_Comm_rank() == 0 ) {
            fprintf( stderr, "SKIPPED: Incorrect parameters for ungqr_hqr (K > N)\n" );
        }
        return -1;
    }

    /* Creates the matrices */
    CHAMELEON_Desc_Create(
        &descA, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, K, 0, 0, M, K, P, Q );
    CHAMELEON_Desc_Create(
        &descQ, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, LDA, N, 0, 0, M, N, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, K, &descTS, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, K, &descTT, P, Q );
#if defined(CHAMELEON_COPY_DIAG)
    {
        D = (CHAM_desc_t*)malloc(sizeof(CHAM_desc_t));
        CHAMELEON_Zdesc_Copy_And_Restrict( descA, D, M, K );
    }
#endif

    /* Initialize matrix tree */
    matrix.mt    = descTS->mt;
    matrix.nt    = descTS->nt;
    matrix.nodes = P * Q;
    matrix.p     = P;

    libhqr_init_hqr( &qrtree, LIBHQR_QR, &matrix, llvl, hlvl, qr_a, qr_p, domino, 0 );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    hres = CHAMELEON_zgeqrf_param_Tile( &qrtree, descA, descTS, descTT );
    if ( hres != CHAMELEON_SUCCESS ) {
        CHAMELEON_Desc_Destroy( &descA  );
        CHAMELEON_Desc_Destroy( &descTS );
        CHAMELEON_Desc_Destroy( &descTT );
        CHAMELEON_Desc_Destroy( &descQ  );
        libhqr_finalize( &qrtree );
        return hres;
    }

    /* Calculates the solution */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zungqr_param_Tile_Async( &qrtree, descA, descTS, descTT, D, descQ,
                                                  test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descTS, test_data.sequence );
        CHAMELEON_Desc_Flush( descTT, test_data.sequence );
        if (D != NULL) {
            CHAMELEON_Desc_Flush( D, test_data.sequence );
        }
        CHAMELEON_Desc_Flush( descQ, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zungqr_param_Tile( &qrtree, descA, descTS, descTT, descQ );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zungqr( M, N, K ) );

    /* Checks the factorisation and orthogonality */
    if ( check ) {
        CHAM_desc_t *descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );
        CHAMELEON_zplrnt_Tile( descA0, seedA );

        hres += check_zortho( args, descQ );
        hres += check_zgeqrf( args, descA0, descA, descQ );

        CHAMELEON_Desc_Destroy( &descA0 );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descTS );
    CHAMELEON_Desc_Destroy( &descTT );
    if (D != NULL) {
        CHAMELEON_Desc_Destroy( &D );
    }
    CHAMELEON_Desc_Destroy( &descQ );
    libhqr_finalize( &qrtree );

    return hres;
}

int
testing_zungqr_hqr_std( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int nb     = run_arg_get_nb( args );
    int ib     = run_arg_get_ib( args );
    int P      = parameters_getvalue_int( "P" );
    int N      = run_arg_get_int( args, "N", 1000 );
    int M      = run_arg_get_int( args, "M", N );
    int K      = run_arg_get_int( args, "K", chameleon_min( M, N ) );
    int LDA    = run_arg_get_int( args, "LDA", M );
    int qr_a   = run_arg_get_int( args, "qra", -1 );
    int qr_p   = run_arg_get_int( args, "qrp", -1 );
    int llvl   = run_arg_get_int( args, "llvl", -1 );
    int hlvl   = run_arg_get_int( args, "hlvl", -1 );
    int domino = run_arg_get_int( args, "domino", -1 );
    int seedA  = run_arg_get_int( args, "seedA", testing_ialea() );
    int Q      = parameters_compute_q( P );

    /* Descriptors */
    CHAMELEON_Complex64_t *A, *Qlap;
    CHAM_desc_t           *descTS, *descTT;
    libhqr_tree_t          qrtree;
    libhqr_matrix_t        matrix;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );
    CHAMELEON_Set( CHAMELEON_INNER_BLOCK_SIZE, ib );

    if ( N > M ) {
        fprintf( stderr, "SKIPPED: Incorrect parameters for ungqr_hqr (N > M)\n" );
        return -1;
    }

    if ( K > N ) {
        fprintf( stderr, "SKIPPED: Incorrect parameters for ungqr_hqr (K > N)\n" );
        return -1;
    }

    /* Creates the matrices */
    A    = malloc( LDA*K*sizeof(CHAMELEON_Complex64_t) );
    Qlap = malloc( LDA*N*sizeof(CHAMELEON_Complex64_t) );
    CHAMELEON_Alloc_Workspace_zgels( M, K, &descTS, P, Q );
    CHAMELEON_Alloc_Workspace_zgels( M, K, &descTT, P, Q );

    /* Initialize matrix tree */
    matrix.mt    = descTS->mt;
    matrix.nt    = descTS->nt;
    matrix.nodes = P * Q;
    matrix.p     = P;

    libhqr_init_hqr( &qrtree, LIBHQR_QR, &matrix, llvl, hlvl, qr_a, qr_p, domino, 0 );

    /* Fills the matrix with random values */
    CHAMELEON_zplrnt( M, K, A, LDA, seedA );
    hres = CHAMELEON_zgeqrf_param( &qrtree, M, K, A, LDA, descTS, descTT );
    if ( hres != CHAMELEON_SUCCESS ) {
        free( A    );
        free( Qlap );
        CHAMELEON_Desc_Destroy( &descTS );
        CHAMELEON_Desc_Destroy( &descTT );
        libhqr_finalize( &qrtree );
        return hres;
    }

    /* Calculates the solution */
    testing_start( &test_data );
    hres = CHAMELEON_zungqr_param( &qrtree, M, N, K, A, LDA, descTS, descTT, Qlap, LDA );
    test_data.hres = hres;
    testing_stop( &test_data, flops_zungqr( M, N, K ) );

    /* Checks the factorisation and orthogonality */
    if ( check ) {
        CHAMELEON_Complex64_t *A0 = malloc( LDA*K*sizeof(CHAMELEON_Complex64_t) );
        CHAMELEON_zplrnt( M, K, A0, LDA, seedA );

        hres += check_zortho_std( args, M, N, Qlap, LDA );
        hres += check_zgeqrf_std( args, M, N, K, A0, A, LDA, Qlap, LDA );

        free( A0 );
    }

    free( A );
    CHAMELEON_Desc_Destroy( &descTS );
    CHAMELEON_Desc_Destroy( &descTT );
    free( Qlap );
    libhqr_finalize( &qrtree );

    return hres;
}

testing_t   test_zungqr_hqr;
const char *zungqr_hqr_params[] = { "mtxfmt", "nb",  "ib",   "m",    "n",      "k",     "lda",
                                    "qra",    "qrp", "llvl", "hlvl", "domino", "seedA", NULL };
const char *zungqr_hqr_output[] = { NULL };
const char *zungqr_hqr_outchk[] = { "||A||", "||I-QQ'||", "||A-fact(A)||", "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zungqr_hqr_init( void ) __attribute__( ( constructor ) );
void
testing_zungqr_hqr_init( void )
{
    test_zungqr_hqr.name   = "zungqr_hqr";
    test_zungqr_hqr.helper = "Q generation with hierarchical reduction trees (QR)";
    test_zungqr_hqr.params = zungqr_hqr_params;
    test_zungqr_hqr.output = zungqr_hqr_output;
    test_zungqr_hqr.outchk = zungqr_hqr_outchk;
    test_zungqr_hqr.fptr_desc = testing_zungqr_hqr_desc;
    test_zungqr_hqr.fptr_std  = testing_zungqr_hqr_std;
    test_zungqr_hqr.next   = NULL;

    testing_register( &test_zungqr_hqr );
}
