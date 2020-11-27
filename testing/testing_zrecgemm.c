/**
 *
 * @file testing_zrecgemm.c
 *
 * @copyright 2019-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zrecgemm testing
 *
 * @version 1.0.0
 * @author Lucas Barros de Assis
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include <chameleon.h>
#include "testings.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>

int
testing_zrecgemm( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int          async  = parameters_getvalue_int( "async" );
    intptr_t     mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int          nb     = run_arg_get_int( args, "nb", 320 );
    int          l1     = run_arg_get_int( args, "l1", 0 );
    int          l2     = run_arg_get_int( args, "l2", 0 );
    int          l3     = run_arg_get_int( args, "l3", 0 );
    cham_rec_t   rec    = run_arg_get_rec( args, "rec", ChamRecFull );
    int          rarg   = run_arg_get_int( args, "rarg", 1 );
    int          P      = parameters_getvalue_int( "P" );
    cham_trans_t transA = run_arg_get_trans( args, "transA", ChamNoTrans );
    cham_trans_t transB = run_arg_get_trans( args, "transB", ChamNoTrans );
    int          N      = run_arg_get_int( args, "N", 1000 );
    int          M      = run_arg_get_int( args, "M", N );
    int          K      = run_arg_get_int( args, "K", N );
    int          LDA    = run_arg_get_int( args, "LDA", ( ( transA == ChamNoTrans ) ? M : K ) );
    int          LDB    = run_arg_get_int( args, "LDB", ( ( transB == ChamNoTrans ) ? K : N ) );
    int          LDC    = run_arg_get_int( args, "LDC", M );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", random() );
    int                   seedB = run_arg_get_int( args, "seedB", random() );
    int                   seedC = run_arg_get_int( args, "seedC", random() );
    int                   Q     = parameters_compute_q( P );

    /* Descriptors */
    int          Am, An, Bm, Bn;
    CHAM_desc_t *descA, *descB, *descC, *descCinit;
    void        *ws = NULL;
    int          list_nb[] = { nb, l1, l2, l3, 0 };

    alpha = run_arg_get_complex64( args, "alpha", alpha );
    beta  = run_arg_get_complex64( args, "beta", beta );

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Calculate the dimensions according to the transposition */
    if ( transA == ChamNoTrans ) {
        Am = M;
        An = K;
    }
    else {
        Am = K;
        An = M;
    }
    if ( transB == ChamNoTrans ) {
        Bm = K;
        Bn = N;
    }
    else {
        Bm = N;
        Bn = K;
    }

    if ( rec == ChamRecSmart ) {
        /* rarg is the percentage of CPU computing power of the node */
        int nt = N/nb;
        float ratio = (float)rarg / 100;
        /* number of RW-tile to partition */
        int ntile = ratio * nt * nt;
        /* Creates the matrices */
        CHAMELEON_Recursive_Desc_Create(
            &descA, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble, ChamRecFull, ntile, list_nb, list_nb,
            LDA, An, Am, An, P, Q, NULL, NULL, NULL, NULL, "A" );
        CHAMELEON_Recursive_Desc_Create(
            &descB, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble, ChamRecFull, ntile, list_nb, list_nb,
            LDB, Bn, Bm, Bn, P, Q, NULL, NULL, NULL, NULL, "B" );
        CHAMELEON_Recursive_Desc_Create(
            &descC, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble, ChamRecSmart, ntile, list_nb, list_nb,
            LDC, N, M, N, P, Q, NULL, NULL, NULL, NULL, "C" );
   }
    else {
        /* Creates the matrices */
        CHAMELEON_Recursive_Desc_Create(
            &descA, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble, rec, rarg, list_nb, list_nb,
            LDA, An, Am, An, P, Q, NULL, NULL, NULL, NULL, "A" );
        CHAMELEON_Recursive_Desc_Create(
            &descB, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble, rec, rarg, list_nb, list_nb,
            LDB, Bn, Bm, Bn, P, Q, NULL, NULL, NULL, NULL, "B" );
        CHAMELEON_Recursive_Desc_Create(
            &descC, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble, rec, rarg, list_nb, list_nb,
            LDC, N, M, N, P, Q, NULL, NULL, NULL, NULL, "C" );
    }

    CHAMELEON_Recursive_Desc_Partition_Submit( descA );
    CHAMELEON_Recursive_Desc_Partition_Submit( descB );
    CHAMELEON_Recursive_Desc_Partition_Submit( descC );

    /* Fill the matrices with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descB, seedB );
    CHAMELEON_zplrnt_Tile( descC, seedC );

    if ( async ) {
        ws = CHAMELEON_zgemm_WS_Alloc( transA, transB, descA, descB, descC );
    }

    /* Calculate the product */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgemm_Tile_Async( transA, transB, alpha, descA, descB, beta, descC, ws,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descB, test_data.sequence );
        CHAMELEON_Desc_Flush( descC, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zgemm_Tile( transA, transB, alpha, descA, descB, beta, descC );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgemm( M, N, K ) );

    if ( ws != NULL ) {
        CHAMELEON_zgemm_WS_Free( ws );
    }

    CHAMELEON_Recursive_Desc_Unpartition_Submit( descA );
    CHAMELEON_Recursive_Desc_Unpartition_Submit( descB );
    CHAMELEON_Recursive_Desc_Unpartition_Submit( descC );

    /* Check the solution */
    if ( check ) {
        CHAMELEON_Desc_Create(
            &descCinit, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, M, N, 0, 0, M, N, P, Q );
        CHAMELEON_zplrnt_Tile( descCinit, seedC );

        hres += check_zgemm( args, transA, transB, alpha, descA, descB, beta, descCinit, descC );

        CHAMELEON_Desc_Destroy( &descCinit );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descB );
    CHAMELEON_Desc_Destroy( &descC );

    return hres;
}

testing_t   test_zrecgemm;
const char *zrecgemm_params[] = { "mtxfmt", "nb", "l1", "l2", "l3", "rec", "rarg", "transA", "transB", "m",     "n",     "k",     "lda", "ldb",
                                  "ldc", "alpha",  "beta",   "seedA", "seedB", "seedC", NULL };
const char *zrecgemm_output[] = { NULL };
const char *zrecgemm_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zrecgemm_init( void ) __attribute__( ( constructor ) );
void
testing_zrecgemm_init( void )
{
    test_zrecgemm.name        = "zrecgemm";
    test_zrecgemm.helper      = "General matrix-matrix multiply";
    test_zrecgemm.params      = zrecgemm_params;
    test_zrecgemm.output      = zrecgemm_output;
    test_zrecgemm.outchk      = zrecgemm_outchk;
    test_zrecgemm.fptr_desc   = testing_zrecgemm;
    test_zrecgemm.fptr_std    = NULL;
    test_zrecgemm.next        = NULL;

    testing_register( &test_zrecgemm );
}
