/**
 *
 * @file pzgenm2.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2020 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon parallel algorithm to compute 2-norm estimator.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Dalal Sukkari
 * @date 2020-10-06
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m, n)     A,      m, n
#define X(m, n)     &X,     m, n
#define SX(m, n)    &SX,    m, n
#define NRMX( m, n) &NRMX,  m, n
#define NRMSX(m, n) &NRMSX, m, n
#define DROW(m, n)  &DROW,  m, n

void
chameleon_pzgenm2( double tol, const CHAM_desc_t *A, double *result,
                   RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t  *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t X, SX, NRMX, NRMSX, DROW;
    int m, n, k;
    int myp = A->myrank / A->q;
    int myq = A->myrank % A->q;
    int tempmm, tempnn;
    int cnt, maxiter;
    double e0, normx, normsx, beta, scl;

    chamctxt = chameleon_context_self();
    if ( sequence->status != CHAMELEON_SUCCESS ) {
        return;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    /* Initialize the result */
    *result = 0.0;

    /* Init workspace handle for the call to dlange but unused */
    RUNTIME_options_ws_alloc( &options, 1, 0 );

    /**
     * DROW must be allocated with GLOBAL to avoid destruction bewteen the first
     * stage and the second one during Flush.
     * This is the same issue for X and SX to be reused from one iteration to another.
     */
    chameleon_desc_init( &DROW, CHAMELEON_MAT_ALLOC_GLOBAL, ChamRealDouble, 1, A->nb, A->nb,
                         A->p, A->n, 0, 0, A->p, A->n, A->p, A->q,
                         NULL, NULL, NULL );
    /**
     * NRMX must be allocated with GLOBAL to be able to access the norm value
     * after flushing the descriptor.
     * This is the same issue for NRMSX.
     */
    chameleon_desc_init( &NRMX, CHAMELEON_MAT_ALLOC_GLOBAL, ChamRealDouble, 2, 1, 2,
                         A->p * 2, A->q, 0, 0, A->p * 2, A->q, A->p, A->q,
                         NULL, NULL, NULL );

    /**
     * Start by computing the initial vector of the iterative process, and that
     * is used to compute the initial norm.
     *
     * For j=0,n-1, drow[j] = sum( |A_{i,j}|, i=0..m-1 )
     * So drow[j] = sum( S_{p,j}, p=0..P-1 ) with S_{p,j} = sum( |A_{i,j}|, i=0..m-1 \ i%P = p )
     *
     */
    for(n = myq; n < A->nt; n += A->q) {
        tempnn = n == A->nt-1 ? A->n - n * A->nb : A->nb;

        /* Zeroes the local intermediate vector */
        INSERT_TASK_dlaset(
            &options,
            ChamUpperLower, 1, tempnn,
            0., 0.,
            DROW( myp, n ) );

        /* Computes the sums of the local tiles into the local vector */
        for(m = myp; m < A->mt; m += A->p) {
            tempmm = m == A->mt-1 ? A->m - m * A->mb : A->mb;
            INSERT_TASK_dzasum(
                &options,
                ChamColumnwise, ChamUpperLower, tempmm, tempnn,
                A(m, n), DROW( myp, n ) );
        }

        /* Reduce on first row of nodes */
	for(m = 1; m < A->p; m++) {
	    INSERT_TASK_daxpy(
                &options, tempnn, 1.,
                DROW( m, n ), 1,
                DROW( 0, n ), 1 );
	}
    }

    /**
     * Reduce now by columns on processes with (myp == 0)
     */
    if ( myp == 0 )
    {
        INSERT_TASK_dlaset(
            &options,
            ChamUpperLower, NRMX.mb, NRMX.nb,
            1., 0.,
            NRMX( myp, myq ) );

        for( n = myq; n < A->nt; n += A->q ) {
	    tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
	    INSERT_TASK_dgessq(
                &options, ChamEltwise, 1, tempnn,
                DROW( myp, n   ),
                NRMX( myp, myq ) );
	}

        /* Reduce on first row of nodes */
	for(n = 1; n < A->q; n++) {
	    INSERT_TASK_dplssq(
                &options, ChamEltwise, 1, 1,
                NRMX( myp, n ),
                NRMX( myp, 0 ) );
	}

        INSERT_TASK_dplssq2(
            &options, 1, NRMX( myp, 0 ) );
    }

    /* Bcast norm over processes from node (0,0) */
    for(m = 0; m < A->p; m++) {
	for(n = 0; n < A->q; n++) {
            if ( (m != 0) || (n != 0) ) {
                INSERT_TASK_dlacpy(
                    &options,
                    ChamUpperLower, 1, 1, 1,
                    NRMX(0, 0),
                    NRMX(m, n) );
            }
        }
    }

    CHAMELEON_Desc_Flush( &DROW, sequence );
    CHAMELEON_Desc_Flush( &NRMX, sequence );
    CHAMELEON_Desc_Flush( A,     sequence );
    RUNTIME_sequence_wait( chamctxt, sequence );
    *result = *((double *)(NRMX.get_blkaddr( &NRMX, myp, myq )));

    if ( (*result) == 0.0 ) {
	chameleon_desc_destroy( &DROW );
	chameleon_desc_destroy( &NRMX );
	return;
    } else {
	normx = *result;
    }

    chameleon_desc_init( &NRMSX, CHAMELEON_MAT_ALLOC_GLOBAL, ChamRealDouble, 2, 1, 2,
                         A->p * 2, A->q, 0, 0, A->p * 2, A->q, A->p, A->q,
                         NULL, NULL, NULL );
    chameleon_desc_init( &X,  CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble, 1, A->nb, A->nb,
                         A->p, A->n, 0, 0, A->p, A->n, A->p, A->q,
                         NULL, NULL, NULL );
    chameleon_desc_init( &SX, CHAMELEON_MAT_ALLOC_GLOBAL, ChamComplexDouble, A->mb, 1, A->mb,
                         A->m, A->q, 0, 0, A->m, A->q, A->p, A->q,
                         NULL, NULL, NULL );

    cnt = 0;
    e0  = 0.;
    maxiter = chameleon_min( 100, A->n );

    while ( (cnt < maxiter) &&
            (fabs(*result - e0) > (tol * (*result))) )
    {
        e0 = *result;

        /* Initialization of X in the first loop */
        if ( cnt == 0 )
        {
            for (n = myq; n < A->nt; n += A->q) {
                tempnn = n == A->nt-1 ? A->n - n * A->nb : A->nb;

                if ( myp == 0 ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
                    INSERT_TASK_dlag2z(
                        &options,
                        ChamUpperLower, 1, tempnn,
                        DROW( 0, n ),
                        X(    0, n ) );
#else
                    INSERT_TASK_zlacpy(
                        &options,
                        ChamUpperLower, 1, tempnn, tempnn,
                        DROW( 0, n ),
                        X(    0, n ) );
#endif
                }

                /* Broadcast X */
                for (m = 1; m < A->p; m++) {
                    INSERT_TASK_zlacpy(
                        &options,
                        ChamUpperLower, 1, tempnn, tempnn,
                        X( 0, n ),
                        X( m, n ) );
                }
            }
            CHAMELEON_Desc_Flush( &DROW, sequence );
        }

        /**
         * Scale x = x / ||A||
         *
         * We scale the P replica of X, so each row of Q processes possess one
         * copy of the scaled X.
         */
        scl = 1. / e0;
        for (n = myq; n < A->nt; n += A->q) {
            tempnn = n == A->nt-1 ? A->n - n * A->nb : A->nb;

            INSERT_TASK_zlascal(
                &options,
                ChamUpperLower, 1, tempnn, tempnn,
                scl, X( myp, n ) );
        }

        /**
         *  Compute Sx = S * x
         */
        for(m = myp; m < A->mt;  m+=A->p) {
            tempmm = m == A->mt-1 ? A->m - m * A->mb : A->mb;

            for (n = myq; n < A->nt; n += A->q ) {
                tempnn = n == A->nt-1 ? A->n - n * A->nb : A->nb;
                beta   = n == myq ? 0. : 1.;

                INSERT_TASK_zgemv(
                    &options,
                    ChamNoTrans, tempmm, tempnn,
                    1.,   A(  m,   n ),
                          X(  myp, n ), 1,
                    beta, SX( m, myq ), 1 );
            }

            /* Reduce columns */
            for (k = 1; k < chameleon_min( A->q, A->nt ); k++) {
                INSERT_TASK_zaxpy(
                    &options, tempmm, 1.,
                    SX( m, k ), 1,
                    SX( m, 0 ), 1 );
            }
            /* Broadcast SX to ease the following gemv */
            for (k = 1; k < A->q; k++) {
                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpperLower, tempmm, 1, tempmm,
                    SX( m, 0 ),
                    SX( m, k ) );
            }
        }

        /**
         *  Compute x = S' * S * x = S' * Sx
         */
        for ( n = myq; n < A->nt; n += A->q ) {
            tempnn = n == A->nt-1 ? A->n - n * A->nb : A->nb;

            for( m = myp; m < A->mt;  m += A->p ) {
                tempmm = m == A->mt-1 ? A->m - m * A->mb : A->mb;
                beta   = m == myp ? 0. : 1.;

                INSERT_TASK_zgemv(
                    &options,
                    ChamConjTrans, tempmm, tempnn,
                    1.,   A(  m,   n ),
                          SX( m, myq ), 1,
                    beta, X(  myp, n ), 1 );
            }

            /* Reduce rows */
            for (k = 1; k < chameleon_min( A->p, A->mt ); k++) {
                INSERT_TASK_zaxpy(
                    &options, tempnn, 1.,
                    X( k, n ), 1,
                    X( 0, n ), 1 );
            }
            /* Broadcast */
            for (k = 1; k < A->p; k++) {
                INSERT_TASK_zlacpy(
                    &options,
                    ChamUpperLower, 1, tempnn, tempnn,
                    X( 0, n ),
                    X( k, n ) );
            }
        }

        /**
         * Compute ||x||
         *
         * All rows of Q nodes compute the same thing in parallel due to replication.
         */
        {
            INSERT_TASK_dlaset(
                &options,
                ChamUpperLower, NRMX.mb, NRMX.nb,
                1., 0.,
                NRMX( myp, myq ) );

            for( n = myq; n < A->nt; n += A->q ) {
                tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;

                INSERT_TASK_zgessq(
                    &options, ChamEltwise, 1, tempnn,
                    X(    myp, n   ),
                    NRMX( myp, myq ) );
            }

            /* Reduce columns  */
            for(n = 1; n < chameleon_min( A->q, A->nt ); n++) {
                INSERT_TASK_dplssq(
                    &options, ChamEltwise, 1, 1,
                    NRMX( myp, n ),
                    NRMX( myp, 0 ) );
            }

            INSERT_TASK_dplssq2(
                &options, 1, NRMX( myp, 0 ) );

            /* Broadcast the results to processes in the same row */
            for(n = 1; n < A->q; n++) {
                INSERT_TASK_dlacpy(
                    &options,
                    ChamUpperLower, 1, 1, 1,
                    NRMX( myp, 0 ),
                    NRMX( myp, n ) );
            }
        }

        /**
         * Compute ||Sx||
         *
         * All columns of P nodes compute the same thing in parallel due to replication.
         */
        {
            INSERT_TASK_dlaset(
                &options,
                ChamUpperLower, NRMSX.mb, NRMSX.nb,
                1., 0.,
                NRMSX( myp, myq ) );

            for( m = myp; m < A->mt; m += A->p ) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                INSERT_TASK_zgessq(
                    &options, ChamEltwise, tempmm, 1,
                    SX(    m,   myq ),
                    NRMSX( myp, myq ) );
            }

            /* Reduce rows */
            for( m = 1; m < chameleon_min( A->p, A->mt ); m++ ) {
                INSERT_TASK_dplssq(
                    &options, ChamEltwise, 1, 1,
                    NRMSX( m, myq ),
                    NRMSX( 0, myq ) );
            }

            INSERT_TASK_dplssq2(
                &options, 1, NRMSX( 0, myq ) );

            /* Broadcast the results to processes in the same column */
            for(m = 1; m < A->p; m++) {
                INSERT_TASK_dlacpy(
                    &options,
                    ChamUpperLower, 1, 1, 1,
                    NRMSX( 0, myq ),
                    NRMSX( m, myq ) );
            }
        }

        CHAMELEON_Desc_Flush( &NRMX,  sequence );
        CHAMELEON_Desc_Flush( &NRMSX, sequence );
        CHAMELEON_Desc_Flush( &X,     sequence );
        CHAMELEON_Desc_Flush( &SX,    sequence );
        CHAMELEON_Desc_Flush( A,      sequence );
        RUNTIME_sequence_wait(chamctxt, sequence);
        normx  = *((double *)(NRMX.get_blkaddr(  &NRMX,  myp, myq )));
        normsx = *((double *)(NRMSX.get_blkaddr( &NRMSX, myp, myq )));

        *result = normx / normsx;
        cnt++;
    }

    if ( (cnt >= maxiter) &&
         (fabs((*result) - e0) > (tol * (*result))) ) {
        *result = INFINITY;
    }

    RUNTIME_options_ws_free( &options );
    RUNTIME_options_finalize( &options, chamctxt );

    chameleon_desc_destroy( &DROW );
    chameleon_desc_destroy( &NRMX );
    chameleon_desc_destroy( &NRMSX );
    chameleon_desc_destroy( &X );
    chameleon_desc_destroy( &SX );

    return;
}
