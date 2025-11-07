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

void
print_descriptor_rec_RW( CHAM_desc_t *desc, char *name, int lvl )
{
    int i;
    if (name) fprintf(stderr,"\n%s (%p) :\n", name, desc);
    for ( i=0; i<desc->lmt*desc->lnt; i++ ) {
        CHAM_tile_t t = desc->tiles[i];
        if ( t.format == CHAMELEON_TILE_DESC ) {
            fprintf(stderr, "[LVL%d] desc-%d (%d, %d)\n", lvl, i, desc->mb, desc->nb);
            print_descriptor_rec_RW( (CHAM_desc_t*)(t.mat), NULL, lvl+1 );
        }
        else if ( t.format == CHAMELEON_TILE_FULLRANK ) {
            /* fprintf(stderr, "[LVL%d] tile-%d\n", lvl, i); */
        }
        else {
            fprintf(stderr, "??\n");
        }
    }
}

int
testing_zrecgemm_RW( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int          async  = parameters_getvalue_int( "async" );
    intptr_t     mtxfmt = parameters_getvalue_int( "mtxfmt" );
    int          N      = run_arg_get_int( args, "N", 1000 );
    int          M      = run_arg_get_int( args, "M", N );
    int          K      = run_arg_get_int( args, "K", N );
    int          nb     = run_arg_get_int( args, "nb", M );
    int          l1     = run_arg_get_int( args, "l1", 0 );
    int          l2     = run_arg_get_int( args, "l2", 0 );
    int          l3     = run_arg_get_int( args, "l3", 0 );
    cham_rec_t   rec    = run_arg_get_rec( args, "rec", ChamRecFull );
    int          rarg   = run_arg_get_int( args, "rarg", 1 );
    int          P      = parameters_getvalue_int( "P" );
    cham_trans_t transA = run_arg_get_trans( args, "transA", ChamNoTrans );
    cham_trans_t transB = run_arg_get_trans( args, "transB", ChamTrans );
    int          LDA    = run_arg_get_int( args, "LDA", ( ( transA == ChamNoTrans ) ? M : K ) );
    int          LDC    = run_arg_get_int( args, "LDC", M );
    CHAMELEON_Complex64_t alpha = testing_zalea();
    CHAMELEON_Complex64_t beta  = testing_zalea();
    int                   seedA = run_arg_get_int( args, "seedA", random() );
    int                   seedC = run_arg_get_int( args, "seedC", random() );
    int                   Q     = parameters_compute_q( P );

    /* Descriptors */
    int          Am, An;
    CHAM_desc_t *descA, *descC, *descCinit;
    int          list_nb[] = { nb, l1, l2, l3, 0 };
    void        *ws = NULL;

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

    /* Creates the matrices */
    CHAMELEON_Recursive_Desc_Create(
        &descA, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble, rec, rarg, list_nb, list_nb,
        LDA, An, Am, An, P, Q, NULL, NULL, NULL, NULL, "A" );
    /* Might need to create descB = tA */
    CHAMELEON_Recursive_Desc_Create_Full(
        &descC, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble, list_nb, list_nb,
        LDC, N, M, N, P, Q, NULL, NULL, NULL, NULL, "C" );

    CHAMELEON_Recursive_Desc_Partition_Submit( descA );
    CHAMELEON_Recursive_Desc_Partition_Submit( descC );

    /* Fill the matrices with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );
    CHAMELEON_zplrnt_Tile( descC, seedC );

    if ( async ) {
        ws = CHAMELEON_zgemm_WS_Alloc( transA, transB, descA, descA, descC );
    }

    /* Calculate the product */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zgemm_Tile_Async( transA, transB, alpha, descA, descA, beta, descC, ws,
                                           test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA, test_data.sequence );
        CHAMELEON_Desc_Flush( descC, test_data.sequence );
    }
    else {
        /* C += A*tA */
        hres = CHAMELEON_zgemm_Tile( transA, transB, alpha, descA, descA, beta, descC );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zgemm( M, N, K ) );

    if ( ws != NULL ) {
        CHAMELEON_zgemm_WS_Free( ws );
    }

    CHAMELEON_Recursive_Desc_Unpartition_Submit( descA );
    CHAMELEON_Recursive_Desc_Unpartition_Submit( descC );

    /* Check the solution */
    if ( check ) {
        CHAMELEON_Desc_Create(
            &descCinit, (void*)(-mtxfmt), ChamComplexDouble, nb, nb, nb * nb, M, N, 0, 0, M, N, P, Q );
        CHAMELEON_zplrnt_Tile( descCinit, seedC );

        hres += check_zgemm( args, transA, transB, alpha, descA, descA, beta, descCinit, descC );

        CHAMELEON_Desc_Destroy( &descCinit );
    }

    CHAMELEON_Desc_Destroy( &descA );
    CHAMELEON_Desc_Destroy( &descC );

    return hres;
}

testing_t   test_zrecgemm_RW;
const char *zrecgemm_RW_params[] = { "mtxfmt", "nb", "l1", "l2", "l3", "rec", "transA", "transB", "m",     "n",     "k",     "lda", "ldb",
                               "ldc", "alpha",  "beta",   "seedA", "seedB", "seedC", NULL };
const char *zrecgemm_RW_output[] = { NULL };
const char *zrecgemm_RW_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zrecgemm_RW_init( void ) __attribute__( ( constructor ) );
void
testing_zrecgemm_RW_init( void )
{
    test_zrecgemm_RW.name        = "zrecgemm_RW";
    test_zrecgemm_RW.helper      = "Inner-product";
    test_zrecgemm_RW.params      = zrecgemm_RW_params;
    test_zrecgemm_RW.output      = zrecgemm_RW_output;
    test_zrecgemm_RW.outchk      = zrecgemm_RW_outchk;
    test_zrecgemm_RW.fptr_desc   = testing_zrecgemm_RW;
    test_zrecgemm_RW.fptr_std    = NULL;
    test_zrecgemm_RW.next        = NULL;

    testing_register( &test_zrecgemm_RW );
}
