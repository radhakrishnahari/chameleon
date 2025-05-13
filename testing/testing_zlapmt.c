/**
 *
 * @file testing_zlapmt.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlapmt testing
 *
 * @version 1.3.0
 * @author Matteo Marcos
 * @date 2025-03-24
 * @precisions normal z -> c d s
 *
 */
#include "chameleon/constants.h"
#include "chameleon/struct.h"
#include "testings.h"
#include "chameleon/chameleon_z.h"
#include "testing_zcheck.h"
#include <chameleon/flops.h>
#include <chameleon/getenv.h>
#include <coreblas/lapacke.h>
#include <chameleon/tasks.h>

static cham_fixdbl_t
flops_zlapmt( cham_fixdbl_t M, cham_fixdbl_t N )
{
    cham_fixdbl_t flops;

    flops = M * N * sizeof( CHAMELEON_Complex64_t );

    return flops;
}

static void
testing_zlapmt_perm_gen( int *PERM,
                         int  k )
{
    int i, j, tmp;

    for ( i = 0; i < k; i++ ) {
        PERM[i] = i + 1;
    }

    for ( i = k - 1; i > 0; i-- ) {
        j = testing_ialea() % ( i + 1 );
        tmp = PERM[i];
        PERM[i] = PERM[j];
        PERM[j] = tmp;
    }

}

int
testing_zlapmt_desc( run_arg_list_t *args, int check )
{
    testdata_t test_data = { .args = args };
    int        hres      = 0;

    /* Read arguments */
    int         async   = parameters_getvalue_int( "async" );
    int         nb      = run_arg_get_nb(  args );
    int         P       = parameters_getvalue_int( "P" );
    cham_side_t side    = run_arg_get_side( args, "side", ChamLeft );
    cham_dir_t  dir     = run_arg_get_dir( args,  "dir", ChamDirForward );
    int         N       = run_arg_get_int( args, "N", 1000 );
    int         M       = run_arg_get_int( args, "M", N );
    int         LDA     = run_arg_get_int( args, "LDA", M );
    int         seedA   = run_arg_get_int( args, "seedA", testing_ialea() );
    int         Q       = parameters_compute_q( P );

    int  K    = ( side == ChamLeft ) ? M : N;
    int  kb   = nb;
    int *PERM = malloc( sizeof(int) * K );

    /* Descriptors */
    CHAM_desc_t *descA;
    CHAM_ipiv_t *descIPIV;
    void        *ws = NULL;

    CHAMELEON_Set( CHAMELEON_TILE_SIZE, nb );

    /* Creates the matrices */
    parameters_desc_create( "A", &descA, ChamComplexDouble, nb, nb, LDA, N, M, N );
    CHAMELEON_Ipiv_Create( &descIPIV, side, kb, K, P, P*Q, NULL );

    /* Fill the matrices with random values */
    CHAMELEON_zplrnt_Tile( descA, seedA );

    /* PERM is initialized with random values that are propagated to the descriptor (should be changed in the future) */
    testing_zlapmt_perm_gen( PERM, K );
    CHAMELEON_Perm_Init( dir, descIPIV, PERM );

    if ( async ) {
        ws = CHAMELEON_zlaswp_WS_Alloc( side, descA );
    }

    /* Calculates the perumtation */
    testing_start( &test_data );
    if ( async ) {
        hres = CHAMELEON_zlapmt_Tile_Async( side, dir, descA, descIPIV, ws,
                                            test_data.sequence, &test_data.request );
        CHAMELEON_Desc_Flush( descA,    test_data.sequence );
        CHAMELEON_Ipiv_Flush( descIPIV, test_data.sequence );
    }
    else {
        hres = CHAMELEON_zlapmt_Tile( side, dir, descA, descIPIV );
    }
    test_data.hres = hres;
    testing_stop( &test_data, flops_zlapmt( M, N ) );

#if !defined(CHAMELEON_SIMULATION)
    if ( check ) {
        CHAM_desc_t *descA0, *descA0c;
        int forwrd = ( dir == ChamDirForward ) ? 1 : 0;

        descA0 = CHAMELEON_Desc_Copy( descA, CHAMELEON_MAT_ALLOC_TILE );

        CHAMELEON_Desc_Create_User(
            &descA0c, (void*)CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble,
            nb, nb, nb*nb, M, N, 0, 0, M, N, 1, 1,
            chameleon_getaddr_cm, chameleon_getblkldd_cm, NULL, NULL );

        CHAMELEON_zplrnt_Tile( descA0c, seedA );

        if ( CHAMELEON_Comm_rank() == 0 ) {
            if ( side == ChamLeft ){
                LAPACKE_zlapmr( LAPACK_COL_MAJOR, forwrd, M, N, descA0c->mat, M, PERM);
            }
            else {
                LAPACKE_zlapmr( LAPACK_ROW_MAJOR, forwrd, N, M, descA0c->mat, M, PERM);
            }
        }

        CHAMELEON_zlacpy_Tile( ChamUpperLower, descA0c, descA0 );
        CHAMELEON_Desc_Destroy( &descA0c );

        hres += check_zmatrices( args, ChamUpperLower, descA, descA0 );

        CHAMELEON_Desc_Destroy( &descA0 );
    }
#endif /* !defined(CHAMELEON_SIMULATION) */

    CHAMELEON_Ipiv_Destroy( &descIPIV );
    parameters_desc_destroy( &descA );
    free( PERM );

    return hres;
}

testing_t   test_zlapmt;
const char *zlapmt_params[] = { "mtxfmt", "nb", "side", "dir", "m", "n", "lda", "seedA", NULL };
const char *zlapmt_output[] = { NULL };
const char *zlapmt_outchk[] = { "RETURN", NULL };

/**
 * @brief Testing registration function
 */
void testing_zlapmt_init( void ) __attribute__( ( constructor ) );
void
testing_zlapmt_init( void )
{
    test_zlapmt.name      = "zlapmt";
    test_zlapmt.helper    = "Row interchange on general matrices";
    test_zlapmt.params    = zlapmt_params;
    test_zlapmt.output    = zlapmt_output;
    test_zlapmt.outchk    = zlapmt_outchk;
    test_zlapmt.fptr_desc = testing_zlapmt_desc;
    test_zlapmt.next      = NULL;

    testing_register( &test_zlapmt );
}

