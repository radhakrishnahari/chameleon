/**
 *
 * @file pzgered.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlange parallel algorithm
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Ana Hourcau
 * @date 2024-07-17
 * @precisions normal z -> z d
 *
 */
//ALLOC_WS :  A->mb
//ALLOC_WS :  A->nb
//WS_ADD :  A->mb + A->nb
#include "control/common.h"
#include <coreblas/lapacke.h>

#define A( m, n )        A,     (m), (n)
#define W( desc, m, n ) (desc), (m), (n)

static inline void
chameleon_pzgered_frb( cham_uplo_t uplo,
                       CHAM_desc_t *A, CHAM_desc_t *Wnorm, CHAM_desc_t *Welt,
                       RUNTIME_option_t *options )
{
    double alpha = 1.0;
    double beta  = 0.0;

    int m, n;
    int minMNT = chameleon_min( A->mt, A->nt );
    int minMN  = chameleon_min( A->m,  A->n  );
    int MT = (uplo == ChamUpper) ? minMNT : A->mt;
    int NT = (uplo == ChamLower) ? minMNT : A->nt;
    int M  = (uplo == ChamUpper) ? minMN  : A->m;
    int N  = (uplo == ChamLower) ? minMN  : A->n;
    int P  = Welt->p;
    int Q  = Welt->q;

    /* Initialize workspaces for tile norms */
    for(m = 0; m < Wnorm->mt; m++) {
        int nmin = ( uplo == ChamUpper ) ? m                      : 0;
        int nmax = ( uplo == ChamLower ) ? chameleon_min(m+1, NT) : NT;

        for(n = nmin; n < nmax; n++) {
            INSERT_TASK_dlaset(
                options,
                ChamUpperLower, Wnorm->mb, Wnorm->nb,
                alpha, beta,
                W( Wnorm, m, n ) );
        }
    }

    /* Initialize workspaces */
    for(m = 0; m < Welt->mt; m++) {
        for(n = 0; n < Welt->nt; n++) {
            INSERT_TASK_dlaset(
                options,
                ChamUpperLower, Welt->mb, Welt->nb,
                alpha, beta,
                W( Welt, m, n ) );
        }
    }

    /**
     * Step 1:
     *  For j in [1,Q], Welt(m, j) = reduce( A(m, j+k*Q) )
     */
    for(m = 0; m < MT; m++) {
        int nmin = ( uplo == ChamUpper ) ? m                      : 0;
        int nmax = ( uplo == ChamLower ) ? chameleon_min(m+1, NT) : NT;

        int tempmm = ( m == (MT-1) ) ? M - m * A->mb : A->mb;

        for(n = nmin; n < nmax; n++) {
            int tempnn = ( n == (NT-1) ) ? N - n * A->nb : A->nb;

            if ( (n == m) && (uplo != ChamUpperLower) ) {
                INSERT_TASK_ztrssq(
                    options,
                    uplo, ChamNonUnit, tempmm, tempnn,
                    A(m, n), W( Wnorm, m, n) );
            }
            else {
                INSERT_TASK_zgessq(
                    options,
                    ChamEltwise,
                    tempmm, tempnn,
                    A(m, n), W( Wnorm, m, n) );
            }

            /* Compress the info per line */
            INSERT_TASK_dplssq(
                options, ChamEltwise, 1, 1, W( Wnorm, m, n), W( Welt, m, n%Q) );

            /* Compute the final norm of the tile */
            INSERT_TASK_dplssq2(
                options, 1, W( Wnorm, m, n ) );
        }

        /**
         * Step 2:
         *  For each j, W(m, j) = reduce( Welt(m, 0..Q-1) )
         */
        for(n = 1; n < Q; n++) {
            INSERT_TASK_dplssq(
                options, ChamEltwise, 1, 1, W( Welt, m, n), W( Welt, m, 0) );
        }
    }

    /**
     * Step 3:
     *  For m in 0..P-1, Welt(m, n) = max( Welt(m..mt[P], n ) )
     */
    for(m = P; m < MT; m++) {
        INSERT_TASK_dplssq(
            options, ChamEltwise, 1, 1, W( Welt, m, 0), W( Welt, m%P, 0) );
    }

    /**
     * Step 4:
     *  For each i, Welt(i, n) = max( Welt(0..P-1, n) )
     */
    for(m = 1; m < P; m++) {
        INSERT_TASK_dplssq(
            options, ChamEltwise, 1, 1, W( Welt, m, 0), W( Welt, 0, 0) );
    }

    INSERT_TASK_dplssq2(
        options, 1, W( Welt, 0, 0) );

    /**
     * Broadcast the result
     */
    for(m = 0; m < A->p; m++) {
        for(n = 0; n < A->q; n++) {
            if ( (m != 0) || (n != 0) ) {
                INSERT_TASK_dlacpy(
                    options,
                    ChamUpperLower, 1, 1,
                    W( Welt, 0, 0 ), W( Welt, m, n ) );
            }
        }
    }
}

/**
 *
 */
void chameleon_pzgered( cham_uplo_t uplo, double prec, CHAM_desc_t *A,
                         RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t Wcol;
    CHAM_desc_t Welt;
    double gnorm, lnorm, threshold, eps;

    int workmt, worknt;
    int m, n;

    chamctxt = chameleon_context_self();
    if ( sequence->status != CHAMELEON_SUCCESS ) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    workmt = chameleon_max( A->mt, A->p );
    worknt = chameleon_max( A->nt, A->q );

    RUNTIME_options_ws_alloc( &options, 1, 0 );

    /* Matrix to store the norm of each element */
    chameleon_desc_init( &Wcol, CHAMELEON_MAT_ALLOC_GLOBAL, ChamRealDouble, 2, 1, 2,
                         A->mt * 2, A->nt, 0, 0, A->mt * 2, A->nt, A->p, A->q,
                         NULL, NULL, A->get_rankof_init, A->get_rankof_init_arg );

    /* Matrix to compute the global frobenius norm */
    chameleon_desc_init( &Welt, CHAMELEON_MAT_ALLOC_GLOBAL, ChamRealDouble, 2, 1, 2,
                         workmt*2, worknt, 0, 0, workmt*2, worknt, A->p, A->q,
                         NULL, NULL, NULL, NULL );

    chameleon_pzgered_frb( uplo, A, &Wcol, &Welt, &options );

    CHAMELEON_Desc_Flush( &Wcol, sequence );
    CHAMELEON_Desc_Flush( &Welt, sequence );
    CHAMELEON_Desc_Flush( A,     sequence );

    RUNTIME_sequence_wait( chamctxt, sequence );

    gnorm = *((double *)Welt.get_blkaddr( &Welt, A->myrank / A->q, A->myrank % A->q ));
    chameleon_desc_destroy( &Welt );

    /**
     * Reduce the precision of the tiles if possible
     */
    if ( prec < 0. ) {
#if !defined(CHAMELEON_SIMULATION)
        eps = LAPACKE_dlamch_work('e');
#else
#if defined(PRECISION_z) || defined(PRECISION_d)
        eps = 1.e-15;
#else
        eps = 1.e-7;
#endif
#endif
    }
    else {
        eps = prec;
    }
    threshold = (eps * gnorm) / (double)(chameleon_min(A->mt, A->nt));

#if defined(CHAMELEON_DEBUG_GERED)
    fprintf( stderr,
             "[%2d] The norm of A is:           %e\n"
             "[%2d] The requested precision is: %e\n"
             "[%2d] The computed threshold is:  %e\n",
             A->myrank, gnorm,
             A->myrank, eps,
             A->myrank, threshold );
#endif
    for(m = 0; m < A->mt; m++) {
        int tempmm = ( m == (A->mt-1) ) ? A->m - m * A->mb : A->mb;
        int nmin   = ( uplo == ChamUpper ) ? m                         : 0;
        int nmax   = ( uplo == ChamLower ) ? chameleon_min(m+1, A->nt) : A->nt;

        for(n = nmin; n < nmax; n++) {
            CHAM_tile_t *tile = A->get_blktile( A, m, n );

            int tempnn = ( n == (A->nt-1) ) ? A->n - n * A->nb : A->nb;

            /*
                * u_{high} = 1e-16 (later should be application accuracy)
                * u_{low} = 1e-8
                * ||A_{i,j}||_F  < u_{high} * || A ||_F / (nt * u_{low})
                * ||A_{i,j}||_F  < threshold / u_{low}
                */
            INSERT_TASK_zgered( &options, threshold,
                                tempmm, tempnn, A( m, n ), W( &Wcol, m, n ) );
        }
    }

    CHAMELEON_Desc_Flush( A, sequence );
    RUNTIME_sequence_wait( chamctxt, sequence );

    chameleon_desc_destroy( &Wcol );
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
