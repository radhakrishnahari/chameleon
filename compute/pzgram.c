/**
 *
 * @file pzgram.c
 *
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgram parallel algorithm
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2020-03-03
 * @precisions normal z -> s d c z
 *
 */
#include "control/common.h"

#define A( m, n )        A,     (m), (n)
#define W( desc, m, n ) (desc), (m), (n)

static inline void
chameleon_pzgram_internal( cham_uplo_t uplo,
                           CHAM_desc_t *A,
                           CHAM_desc_t *Wcol,
                           CHAM_desc_t *Welt,
                           RUNTIME_option_t *options )
{
    int m, n;
    int MT = A->mt;
    int NT = A->nt;
    int M  = A->m;
    int N  = A->n;
    int P  = Welt->p;
    int Q  = Welt->q;

    /**
     *  1) compute (scl,ssq) over columns in each tile
     */
    for(n = 0; n < NT; n++) {
        int mmin = ( uplo == ChamLower ) ? n                      : 0;
        int mmax = ( uplo == ChamUpper ) ? chameleon_min(n+1, MT) : MT;
        int tempnn = ( n == (NT-1) ) ? N - n * A->nb : A->nb;

        for(m = mmin; m < mmax; m++) {
            int tempmm = ( m == (MT-1) ) ? M - m * A->mb : A->mb;

            if ( n == m ) {
                INSERT_TASK_dsyssq(
                    options, ChamColumnwise, uplo, tempmm,
                    A(m, n), W( Wcol, m, n) );
            }
            else {
                INSERT_TASK_dgessq(
                    options, ChamColumnwise, tempmm, tempnn,
                    A(m, n), W( Wcol, m, n) );
                if ( uplo != ChamUpperLower ) {
                    INSERT_TASK_dgessq(
                        options, ChamRowwise, tempmm, tempnn,
                        A(m, n), W( Wcol, n, m) );
                }
            }
        }
    }

    for(n = 0; n < NT; n++) {
        int tempnn = ( n == (NT-1) ) ? N - n * A->nb : A->nb;

        /**
         *  2) reduce columns (scl,ssq) tiles per processus (between lines)
         */
        for(m = P; m < MT; m++) {
            INSERT_TASK_dplssq(
                options, ChamColumnwise, 1, tempnn,
                W( Wcol, m,   n ),
                W( Wcol, m%P, n ) );
        }

        /**
         *  3) reduce columns (scl,ssq) tiles on the first line of tiles
         */
        for(m = 1; m < P; m++) {
            INSERT_TASK_dplssq(
                options, ChamColumnwise, 1, tempnn,
                W( Wcol, m, n ),
                W( Wcol, 0, n ) );
        }

        /* 4) reduce (scl,ssq) inside each tile of the first line of tiles for the global sum square */
        INSERT_TASK_dplssq(
            options, ChamEltwise, 1, tempnn,
            W( Wcol, 0, n ),
            W( Welt, 0, n ) );

        /* 5) deduce the sum square for each column from the pairs (scl,ssq) -> sqrt(sum) = scl*sqrt(ssq) */
        INSERT_TASK_dplssq2( options, tempnn, W( Wcol, 0, n ) );
    }

    /* 6) reduce global sum squares on each processus (between columns) */
    for(n = Q; n < NT; n++) {
        INSERT_TASK_dplssq( options, ChamEltwise, 1, 1, W( Welt, 0, n), W( Welt, 0, n%Q) );
    }

    /* 7) reduce global sum squares on the first tile (index 0, 0) */
    for(n = 1; n < Q; n++) {
        INSERT_TASK_dplssq(
            options, ChamEltwise, 1, 1, W( Welt, 0, n), W( Welt, 0, 0) );
    }

    /* 8) deduce the global sum square from the pair (scl,ssq) -> sqrt(sum) = scl*sqrt(ssq) */
    INSERT_TASK_dplssq2( options, 1, W( Welt, 0, 0) );

    /* Finally compute Gram matrix coefficients inplace */
    for(n = 0; n < NT; n++) {
        int mmin = ( uplo == ChamLower ) ? n                      : 0;
        int mmax = ( uplo == ChamUpper ) ? chameleon_min(n+1, MT) : MT;
        int tempnn = ( n == (NT-1) ) ? N - n * A->nb : A->nb;

        for(m = mmin; m < mmax; m++) {
            int tempmm = ( m == (MT-1) ) ? M - m * A->mb : A->mb;

            INSERT_TASK_zgram(
                options,
                ( m == n ) ? uplo : ChamUpperLower,
                A->m, A->n, tempmm, tempnn,
                W( Wcol, 0, m ),
                W( Wcol, 0, n ),
                W( Welt, 0, 0 ),
                A( m, n ) );
        }
    }
}

/**
 *
 */
void chameleon_pzgram( cham_uplo_t uplo, CHAM_desc_t *A, RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t Wcol;
    CHAM_desc_t Welt;
    int workmt, worknt;
    int m, n, tempmm, tempnn;

    chamctxt = chameleon_context_self();
    if ( sequence->status != CHAMELEON_SUCCESS ) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    workmt = chameleon_max( A->mt, A->p );
    worknt = chameleon_max( A->nt, A->q );

    RUNTIME_options_ws_alloc( &options, 1, 0 );

    chameleon_desc_init( &Wcol, CHAMELEON_MAT_ALLOC_TILE, ChamRealDouble, 2, A->nb, 2*A->nb,
                         2*workmt, A->n, 0, 0, 2*workmt, A->n, A->p, A->q,
                         NULL, NULL, NULL );

    chameleon_desc_init( &Welt, CHAMELEON_MAT_ALLOC_TILE, ChamRealDouble, 2, 1, 2,
                         2, worknt, 0, 0, 2, worknt, A->p, A->q,
                         NULL, NULL, NULL );

    /* Initialize Wcol */
    for(m = 0; m < Wcol.mt; m++) {
        tempmm = m == Wcol.mt-1 ? Wcol.m-m*Wcol.mb : Wcol.mb;
        for(n = 0; n < Wcol.nt; n++) {
            tempnn = n == Wcol.nt-1 ? Wcol.n-n*Wcol.nb : Wcol.nb;
            INSERT_TASK_dlaset(
                &options,
                ChamUpperLower, tempmm, tempnn,
                -1., -1.,
                W( &Wcol, m, n ) );
        }
    }
    /* Initialize Welt */
    for(m = 0; m < Welt.mt; m++) {
        tempmm = m == Welt.mt-1 ? Welt.m-m*Welt.mb : Welt.mb;
        for(n = 0; n < Welt.nt; n++) {
            tempnn = n == Welt.nt-1 ? Welt.n-n*Welt.nb : Welt.nb;
            INSERT_TASK_dlaset(
                &options,
                ChamUpperLower, tempmm, tempnn,
                -1., -1.,
                W( &Welt, m, n ) );
        }
    }

    chameleon_pzgram_internal( uplo, A, &Wcol, &Welt, &options );

    CHAMELEON_Desc_Flush( &Wcol, sequence );
    CHAMELEON_Desc_Flush( &Welt, sequence );
    CHAMELEON_Desc_Flush( A, sequence );
    RUNTIME_sequence_wait( chamctxt, sequence );

    chameleon_desc_destroy( &Wcol );
    chameleon_desc_destroy( &Welt );

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
