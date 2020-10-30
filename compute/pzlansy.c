/**
 *
 * @file pzlansy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlansy parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for CHAMELEON 0.9.2
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
//ALLOC_WS :  A->mb
//WS_ADD :  A->mb

#include <stdlib.h>
#include <math.h>
#include "control/common.h"

#define A( m, n )        A,     (m), (n)
#define W( desc, m, n ) (desc), (m), (n)

static inline void
chameleon_pzlansy_inf( cham_uplo_t uplo, CHAM_desc_t *A,
                       CHAM_desc_t *Wcol, CHAM_desc_t *Welt,
                       RUNTIME_option_t *options)
{
    int m, n;
    int MT = A->mt;
    int NT = A->nt;
    int M  = A->m;
    int N  = A->n;
    int P  = Welt->p;
    int Q  = Welt->q;

    /**
     * Step 1:
     *  For j in [1,Q], Wcol(m, j) = reduce( A(m, j+k*Q) )
     */
    for(m = 0; m < MT; m++) {
        int nmin = ( uplo == ChamUpper ) ? m                      : 0;
        int nmax = ( uplo == ChamLower ) ? chameleon_min(m+1, NT) : NT;

        int tempmm = ( m == (MT-1) ) ? M - m * A->mb : A->mb;

        for(n = nmin; n < nmax; n++) {
            int tempnn = ( n == (NT-1) ) ? N - n * A->nb : A->nb;

            if ( n == m ) {
                INSERT_TASK_dzasum(
                    options,
                    ChamRowwise, uplo, tempmm, tempnn,
                    A(m, n), W( Wcol, m, n) );
            }
            else {
                INSERT_TASK_dzasum(
                    options,
                    ChamRowwise, ChamUpperLower, tempmm, tempnn,
                    A(m, n), W( Wcol, m, n) );

                INSERT_TASK_dzasum(
                    options,
                    ChamColumnwise, ChamUpperLower, tempmm, tempnn,
                    A(m, n), W( Wcol, n, m) );
            }
        }
    }

    for(m = 0; m < MT; m++) {
        int tempmm = ( m == (MT-1) ) ? M - m * A->mb : A->mb;

        for(n = Q; n < NT; n++) {
            INSERT_TASK_daxpy(
                options, tempmm, 1.,
                W( Wcol, m, n   ), 1,
                W( Wcol, m, n%Q ), 1 );
        }

        /**
         * Step 2:
         *  For each j, W(m, j) = reduce( Wcol(m, 0..Q-1) )
         */
        for(n = 1; n < Q; n++) {
            INSERT_TASK_daxpy(
                options, tempmm, 1.,
                W( Wcol, m, n ), 1,
                W( Wcol, m, 0 ), 1 );
        }

        INSERT_TASK_dlange(
            options,
            ChamMaxNorm, tempmm, 1, A->nb,
            W( Wcol, m, 0), W( Welt, m, 0));
    }

    /**
     * Step 3:
     *  For m in 0..P-1, W( Welt, m, n) = max( Wcol(m..mt[P], n ) )
     */
    for(m = P; m < MT; m++) {
        INSERT_TASK_dlange_max(
            options,
            W( Welt, m, 0), W( Welt, m%P, 0) );
    }

    /**
     * Step 4:
     *  For each i, W( Welt, i, n) = max( W( Welt, 0..P-1, n) )
     */
    for(m = 1; m < P; m++) {
        INSERT_TASK_dlange_max(
            options,
            W( Welt, m, 0), W( Welt, 0, 0) );
    }
}

static inline void
chameleon_pzlansy_max( cham_trans_t trans, cham_uplo_t uplo, CHAM_desc_t *A,
                       CHAM_desc_t *Welt, RUNTIME_option_t *options)
{
    int m, n;
    int MT = A->mt;
    int NT = A->nt;
    int M  = A->m;
    int N  = A->n;
    int P  = Welt->p;
    int Q  = Welt->q;

    /**
     * Step 1:
     *  For j in [1,Q], Welt(m, j) = reduce( A(m, j+k*Q) )
     */
    for(m = 0; m < MT; m++) {
        int nmin = (uplo == ChamUpper ) ? m                      : 0;
        int nmax = (uplo == ChamLower ) ? chameleon_min(m+1, NT) : NT;

        int tempmm = ( m == (MT-1) ) ? M - m * A->mb : A->mb;

        for(n = nmin; n < nmax; n++) {
            int tempnn = ( n == (NT-1) ) ? N - n * A->nb : A->nb;

            if ( n == m ) {
                if ( trans == ChamConjTrans) {
                    INSERT_TASK_zlanhe(
                        options,
                        ChamMaxNorm, uplo, tempmm, A->nb,
                        A(m, n), W( Welt, m, n));
                }
                else {
                    INSERT_TASK_zlansy(
                        options,
                        ChamMaxNorm, uplo, tempmm, A->nb,
                        A(m, n), W( Welt, m, n));
                }
            }
            else {
                INSERT_TASK_zlange(
                    options,
                    ChamMaxNorm, tempmm, tempnn, A->nb,
                    A(m, n), W( Welt, m, n));
            }

            if ( n >= Q ) {
                INSERT_TASK_dlange_max(
                    options,
                    W( Welt, m, n), W( Welt, m, n%Q) );
            }
        }

        /**
         * Step 2:
         *  For each j, W(m, j) = reduce( Welt(m, 0..Q-1) )
         */
        for(n = 1; n < Q; n++) {
            INSERT_TASK_dlange_max(
                options,
                W( Welt, m, n), W( Welt, m, 0) );
        }
    }

    /**
     * Step 3:
     *  For m in 0..P-1, Welt(m, n) = max( Welt(m..mt[P], n ) )
     */
    for(m = P; m < MT; m++) {
        INSERT_TASK_dlange_max(
            options,
            W( Welt, m, 0), W( Welt, m%P, 0) );
    }

    /**
     * Step 4:
     *  For each i, Welt(i, n) = max( Welt(0..P-1, n) )
     */
    for(m = 1; m < P; m++) {
        INSERT_TASK_dlange_max(
            options,
            W( Welt, m, 0), W( Welt, 0, 0) );
    }
}

static inline void
chameleon_pzlansy_frb( cham_trans_t trans, cham_uplo_t uplo,
                       CHAM_desc_t *A, CHAM_desc_t *Welt,
                       RUNTIME_option_t *options)
{
    int m, n;
    int MT = A->mt;
    int NT = A->nt;
    int M  = A->m;
    int N  = A->n;
    int P  = Welt->p;
    int Q  = Welt->q;

    /**
     * Step 1:
     *  For j in [1,Q], Welt(m, j) = reduce( A(m, j+k*Q) )
     */
    for(m = 0; m < MT; m++) {
        int nmin = (uplo == ChamUpper ) ? m                      : 0;
        int nmax = (uplo == ChamLower ) ? chameleon_min(m+1, NT) : NT;

        int tempmm = ( m == (MT-1) ) ? M - m * A->mb : A->mb;

        for(n = nmin; n < nmax; n++) {
            int tempnn = ( n == (NT-1) ) ? N - n * A->nb : A->nb;

            if ( n == m ) {
                if ( trans == ChamConjTrans) {
                    INSERT_TASK_zhessq(
                        options, ChamEltwise, uplo, tempmm,
                        A(m, n), W( Welt, m, n) );
                }
                else {
                    INSERT_TASK_zsyssq(
                        options, ChamEltwise, uplo, tempmm,
                        A(m, n), W( Welt, m, n) );
                }
            }
            else {
                INSERT_TASK_zgessq(
                    options, ChamEltwise, tempmm, tempnn,
                    A(m, n), W( Welt, m, n) );
                INSERT_TASK_zgessq(
                    options, ChamEltwise, tempmm, tempnn,
                    A(m, n), W( Welt, n, m) );
            }
        }
    }

    for(m = 0; m < MT; m++) {
        for(n = Q; n < NT; n++) {
            INSERT_TASK_dplssq(
                options, ChamEltwise, 1, 1, W( Welt, m, n), W( Welt, m, n%Q) );
        }

        /**
         * Step 2:
         *  For each j, W(m, j) = reduce( W( Welt, m, 0..Q-1) )
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
}

/**
 *
 */
void chameleon_pzlansy_generic( cham_normtype_t norm, cham_uplo_t uplo, cham_trans_t trans,
                                CHAM_desc_t *A, double *result,
                                RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t Wcol;
    CHAM_desc_t Welt;
    double alpha = 0.0;
    double beta  = 0.0;

    int workmt, worknt;
    int m, n, wcol_init = 0;

    chamctxt = chameleon_context_self();
    if ( sequence->status != CHAMELEON_SUCCESS ) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    *result = 0.0;

    workmt = chameleon_max( A->mt, A->p );
    worknt = chameleon_max( A->nt, A->q );

    switch ( norm ) {
    case ChamOneNorm:
    case ChamInfNorm:
        RUNTIME_options_ws_alloc( &options, 1, 0 );

        chameleon_desc_init( &Wcol, CHAMELEON_MAT_ALLOC_TILE, ChamRealDouble, A->mb, 1, A->mb,
                             workmt * A->mb, worknt, 0, 0, workmt * A->mb, worknt, A->p, A->q,
                             NULL, NULL, NULL );
        wcol_init = 1;

        /*
         * Use the global allocator for Welt, otherwise flush may free the data before the result is read.
         */
        chameleon_desc_init( &Welt, CHAMELEON_MAT_ALLOC_GLOBAL, ChamRealDouble, 1, 1, 1,
                             workmt, A->q, 0, 0, workmt, A->q, A->p, A->q,
                             NULL, NULL, NULL );
        break;

        /*
         *  ChamFrobeniusNorm
         */
    case ChamFrobeniusNorm:
        RUNTIME_options_ws_alloc( &options, 1, 0 );

        alpha = 1.;
        chameleon_desc_init( &Welt, CHAMELEON_MAT_ALLOC_GLOBAL, ChamRealDouble, 2, 1, 2,
                             workmt*2, worknt, 0, 0, workmt*2, worknt, A->p, A->q,
                             NULL, NULL, NULL );
        break;

        /*
         *  ChamMaxNorm
         */
    case ChamMaxNorm:
    default:
        RUNTIME_options_ws_alloc( &options, 1, 0 );

        chameleon_desc_init( &Welt, CHAMELEON_MAT_ALLOC_GLOBAL, ChamRealDouble, 1, 1, 1,
                             workmt, worknt, 0, 0, workmt, worknt, A->p, A->q,
                             NULL, NULL, NULL );
    }

    /* Initialize workspaces */
    if ( (norm == ChamInfNorm) ||
         (norm == ChamOneNorm) )
    {
        /* Initialize Wcol tile */
        for(m = 0; m < Wcol.mt; m++) {
            for(n = 0; n < Wcol.nt; n++) {
                INSERT_TASK_dlaset(
                    &options,
                    ChamUpperLower, Wcol.mb, Wcol.nb,
                    alpha, beta,
                    W( &Wcol, m, n ) );
            }
        }
    }
    for(m = 0; m < Welt.mt; m++) {
        for(n = 0; n < Welt.nt; n++) {
            INSERT_TASK_dlaset(
                &options,
                ChamUpperLower, Welt.mb, Welt.nb,
                alpha, beta,
                W( &Welt, m, n ) );
        }
    }

    switch ( norm ) {
    case ChamOneNorm:
    case ChamInfNorm:
        chameleon_pzlansy_inf( uplo, A, &Wcol, &Welt, &options );
        CHAMELEON_Desc_Flush( &Wcol, sequence );
        break;

    case ChamFrobeniusNorm:
        chameleon_pzlansy_frb( trans, uplo, A, &Welt, &options );
        break;

    case ChamMaxNorm:
    default:
        chameleon_pzlansy_max( trans, uplo, A, &Welt, &options );
    }

    /**
     * Broadcast the result
     */
    for(m = 0; m < A->p; m++) {
        for(n = 0; n < A->q; n++) {
            if ( (m != 0) || (n != 0) ) {
                INSERT_TASK_dlacpy(
                    &options,
                    ChamUpperLower, 1, 1, 1,
                    W( &Welt, 0, 0 ), W( &Welt, m, n ));
            }
        }
    }

    if ( wcol_init ) {
        CHAMELEON_Desc_Flush( &Wcol, sequence );
    }
    CHAMELEON_Desc_Flush( &Welt, sequence );
    CHAMELEON_Desc_Flush( A, sequence );
    RUNTIME_sequence_wait(chamctxt, sequence);

    *result = *(double *)Welt.get_blkaddr( &Welt, A->myrank / A->q, A->myrank % A->q );

    if ( wcol_init ) {
        chameleon_desc_destroy( &Wcol );
    }
    chameleon_desc_destroy( &Welt );

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}
