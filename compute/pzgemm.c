/**
 *
 * @file pzgemm.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgemm parallel algorithm
 *
 * @version 1.2.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Alycia Lisito
 * @date 2022-02-22
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(  _m_, _n_ ) A,  (_m_), (_n_)
#define B(  _m_, _n_ ) B,  (_m_), (_n_)
#define C(  _m_, _n_ ) C,  (_m_), (_n_)
#define WA( _m_, _n_ ) WA, (_m_), (_n_)
#define WB( _m_, _n_ ) WB, (_m_), (_n_)

/**
 *  Parallel tile matrix-matrix multiplication.
 *  Generic algorithm for any data distribution with a stationnary A.
 */
static inline void
chameleon_pzgemm_Astat( CHAM_context_t *chamctxt, cham_trans_t transA, cham_trans_t transB,
                        CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                        CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                        RUNTIME_option_t *options )
{
    const CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;
    RUNTIME_sequence_t *sequence = options->sequence;
    int                 m, n, k;
    int                 tempmm, tempnn, tempkn, tempkm;
    int                 myrank  = RUNTIME_comm_rank( chamctxt );
    int8_t             *reduceC = calloc( C->mt * C->nt, sizeof(int8_t) );

    /* Set C tiles to redux mode. */
    for (n = 0; n < C->nt; n++) {
        for (m = 0; m < C->mt; m++) {
            /* The node owns the C tile. */
            if ( C->get_rankof( C(m, n) ) == myrank ) {
                reduceC[ n * C->mt + m ] = 1;
                RUNTIME_zgersum_set_methods( C(m, n) );
                continue;
            }

            /*
             * The node owns the A tile that will define the locality of the
             * computations.
             */
            if ( transA == ChamNoTrans ) {
                for (k = 0; k < A->nt; k++) {
                    if ( A->get_rankof( A(m, k) ) == myrank ) {
                        reduceC[ n * C->mt + m ] = 1;
                        RUNTIME_zgersum_set_methods( C(m, n) );
                        break;
                    }
                }
            }
            else {
                for (k = 0; k < A->mt; k++) {
                    if ( A->get_rankof( A(k, m) ) == myrank ) {
                        reduceC[ n * C->mt + m ] = 1;
                        RUNTIME_zgersum_set_methods( C(m, n) );
                        break;
                    }
                }
            }
        }
    }

    for (n = 0; n < C->nt; n++) {
        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
        for (m = 0; m < C->mt; m++) {
            tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

            /* Scale C */
            options->forcesub = 0;
            INSERT_TASK_zlascal( options, ChamUpperLower, tempmm, tempnn, C->mb,
                                 beta, C, m, n );
            options->forcesub = reduceC[ n * C->mt + m ];

            /*
             *  A: ChamNoTrans / B: ChamNoTrans
             */
            if (transA == ChamNoTrans) {
                if (transB == ChamNoTrans) {
                    for (k = 0; k < A->nt; k++) {
                        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

                        INSERT_TASK_zgemm_Astat(
                            options,
                            transA, transB,
                            tempmm, tempnn, tempkn, A->mb,
                            alpha, A(m, k),  /* lda * Z */
                                   B(k, n),  /* ldb * Y */
                            zone,  C(m, n)); /* ldc * Y */
                    }
                }
                /*
                 *  A: ChamNoTrans / B: Cham[Conj]Trans
                 */
                else {
                    for (k = 0; k < A->nt; k++) {
                        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;

                        INSERT_TASK_zgemm_Astat(
                            options,
                            transA, transB,
                            tempmm, tempnn, tempkn, A->mb,
                            alpha, A(m, k),  /* lda * Z */
                                   B(n, k),  /* ldb * Z */
                            zone,  C(m, n)); /* ldc * Y */
                    }
                }
            }
            /*
             *  A: Cham[Conj]Trans / B: ChamNoTrans
             */
            else {
                if (transB == ChamNoTrans) {
                    for (k = 0; k < A->mt; k++) {
                        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;

                        INSERT_TASK_zgemm_Astat(
                            options,
                            transA, transB,
                            tempmm, tempnn, tempkm, A->mb,
                            alpha, A(k, m),  /* lda * X */
                                   B(k, n),  /* ldb * Y */
                            zone,  C(m, n)); /* ldc * Y */
                    }
                }
                /*
                 *  A: Cham[Conj]Trans / B: Cham[Conj]Trans
                 */
                else {
                    for (k = 0; k < A->mt; k++) {
                        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;

                        INSERT_TASK_zgemm_Astat(
                            options,
                            transA, transB,
                            tempmm, tempnn, tempkm, A->mb,
                            alpha, A(k, m),  /* lda * X */
                                   B(n, k),  /* ldb * Z */
                            zone,  C(m, n)); /* ldc * Y */
                    }
                }
            }
            RUNTIME_zgersum_submit_tree( options, C(m, n) );
            RUNTIME_data_flush( sequence, C(m, n) );
        }
    }
    options->forcesub = 0;
    free( reduceC );

    (void)chamctxt;
}

/**
 *  Parallel tile matrix-matrix multiplication
 *  SUMMA algorithm for 2D block-cyclic data distribution.
 */
static inline void
chameleon_pzgemm_summa( CHAM_context_t *chamctxt, cham_trans_t transA, cham_trans_t transB,
                        CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                        CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                        CHAM_desc_t *WA, CHAM_desc_t *WB,
                        RUNTIME_option_t *options )
{
    RUNTIME_sequence_t *sequence = options->sequence;
    int m, n, k, p, q, KT, K, lp, lq;
    int tempmm, tempnn, tempkk;
    int lookahead, myp, myq;

    CHAMELEON_Complex64_t zbeta;
    CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;

    lookahead = chamctxt->lookahead;

    KT  = transA == ChamNoTrans ? A->nt : A->mt;
    K   = transA == ChamNoTrans ? A->n  : A->m;
    myp = C->myrank / C->q;
    myq = C->myrank % C->q;

    /*
     *  A: ChamNoTrans / B: ChamNoTrans
     */
    for (k = 0; k < KT; k++ ) {
        lp = (k % lookahead) * C->p;
        lq = (k % lookahead) * C->q;
        tempkk = k == KT - 1 ? K - k * A->nb : A->nb;
        zbeta = k == 0 ? beta : zone;

        /* Transfert ownership of the k column of A */
        for (m = 0; m < C->mt; m ++ ) {
            tempmm = m == C->mt-1 ? C->m - m * C->mb : C->mb;

            if ( transA == ChamNoTrans ) {
                INSERT_TASK_zlacpy(
                    options,
                    ChamUpperLower, tempmm, tempkk,
                    A(  m,  k ),
                    WA( m, (k % C->q) + lq ) );

                RUNTIME_data_flush( sequence, A( m, k ) );

                for ( q=1; q < C->q; q++ ) {
                    INSERT_TASK_zlacpy(
                        options,
                        ChamUpperLower, tempmm, tempkk,
                        WA( m, ((k+q-1) % C->q) + lq ),
                        WA( m, ((k+q)   % C->q) + lq ) );
                }
            }
            else {
                INSERT_TASK_zlacpy(
                    options,
                    ChamUpperLower, tempkk, tempmm,
                    A(  k,  m ),
                    WA( m, (m % C->q) + lq ) );

                RUNTIME_data_flush( sequence, A( k, m ) );

                for ( q=1; q < C->q; q++ ) {
                    INSERT_TASK_zlacpy(
                        options,
                        ChamUpperLower, tempkk, tempmm,
                        WA( m, ((m+q-1) % C->q) + lq ),
                        WA( m, ((m+q)   % C->q) + lq ) );
                }
            }
        }

        /* Transfert ownership of the k row of B */
        for (n = 0; n < C->nt; n++) {
            tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

            if ( transB == ChamNoTrans ) {
                INSERT_TASK_zlacpy(
                    options,
                    ChamUpperLower, tempkk, tempnn,
                    B(   k,              n ),
                    WB( (k % C->p) + lp, n ) );

                RUNTIME_data_flush( sequence, B( k, n ) );

                for ( p=1; p < C->p; p++ ) {
                    INSERT_TASK_zlacpy(
                        options,
                        ChamUpperLower, tempkk, tempnn,
                        WB( ((k+p-1) % C->p) + lp, n ),
                        WB( ((k+p)   % C->p) + lp, n ) );
                }
            }
            else {
                INSERT_TASK_zlacpy(
                    options,
                    ChamUpperLower, tempnn, tempkk,
                    B(   n,              k ),
                    WB( (n % C->p) + lp, n ) );

                RUNTIME_data_flush( sequence, B( n, k ) );

                for ( p=1; p < C->p; p++ ) {
                    INSERT_TASK_zlacpy(
                        options,
                        ChamUpperLower, tempnn, tempkk,
                        WB( ((n+p-1) % C->p) + lp, n ),
                        WB( ((n+p)   % C->p) + lp, n ) );
                }
            }
        }

        for (m = myp; m < C->mt; m+=C->p) {
            tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

            for (n = myq; n < C->nt; n+=C->q) {
                tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                INSERT_TASK_zgemm(
                    options,
                    transA, transB,
                    tempmm, tempnn, tempkk, A->mb,
                    alpha, WA( m,        myq + lq ),  /* lda * Z */
                           WB( myp + lp, n        ),  /* ldb * Y */
                    zbeta, C(  m,        n        ) ); /* ldc * Y */
            }
        }
    }

    CHAMELEON_Desc_Flush( WA, sequence );
    CHAMELEON_Desc_Flush( WB, sequence );
}

/**
 *  Parallel tile matrix-matrix multiplication.
 *  Generic algorithm for any data distribution.
 */
static inline void
chameleon_pzgemm_generic( CHAM_context_t *chamctxt, cham_trans_t transA, cham_trans_t transB,
                          CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                          CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                          RUNTIME_option_t *options )
{
    RUNTIME_sequence_t *sequence = options->sequence;

    int m, n, k;
    int tempmm, tempnn, tempkn, tempkm;

    CHAMELEON_Complex64_t zbeta;
    CHAMELEON_Complex64_t zone = (CHAMELEON_Complex64_t)1.0;

    for (m = 0; m < C->mt; m++) {
        tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
        for (n = 0; n < C->nt; n++) {
            tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
            /*
             *  A: ChamNoTrans / B: ChamNoTrans
             */
            if (transA == ChamNoTrans) {
                if (transB == ChamNoTrans) {
                    for (k = 0; k < A->nt; k++) {
                        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                        zbeta = k == 0 ? beta : zone;
                        INSERT_TASK_zgemm(
                            options,
                            transA, transB,
                            tempmm, tempnn, tempkn, A->mb,
                            alpha, A(m, k),  /* lda * Z */
                            B(k, n),  /* ldb * Y */
                            zbeta, C(m, n)); /* ldc * Y */
                    }
                }
                /*
                 *  A: ChamNoTrans / B: Cham[Conj]Trans
                 */
                else {
                    for (k = 0; k < A->nt; k++) {
                        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                        zbeta = k == 0 ? beta : zone;
                        INSERT_TASK_zgemm(
                            options,
                            transA, transB,
                            tempmm, tempnn, tempkn, A->mb,
                            alpha, A(m, k),  /* lda * Z */
                            B(n, k),  /* ldb * Z */
                            zbeta, C(m, n)); /* ldc * Y */
                    }
                }
            }
            /*
             *  A: Cham[Conj]Trans / B: ChamNoTrans
             */
            else {
                if (transB == ChamNoTrans) {
                    for (k = 0; k < A->mt; k++) {
                        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                        zbeta = k == 0 ? beta : zone;
                        INSERT_TASK_zgemm(
                            options,
                            transA, transB,
                            tempmm, tempnn, tempkm, A->mb,
                            alpha, A(k, m),  /* lda * X */
                            B(k, n),  /* ldb * Y */
                            zbeta, C(m, n)); /* ldc * Y */
                    }
                }
                /*
                 *  A: Cham[Conj]Trans / B: Cham[Conj]Trans
                 */
                else {
                    for (k = 0; k < A->mt; k++) {
                        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                        zbeta = k == 0 ? beta : zone;
                        INSERT_TASK_zgemm(
                            options,
                            transA, transB,
                            tempmm, tempnn, tempkm, A->mb,
                            alpha, A(k, m),  /* lda * X */
                            B(n, k),  /* ldb * Z */
                            zbeta, C(m, n)); /* ldc * Y */
                    }
                }
            }
            RUNTIME_data_flush( sequence, C(m, n) );
        }
        if (transA == ChamNoTrans) {
            for (k = 0; k < A->nt; k++) {
                RUNTIME_data_flush( sequence, A(m, k) );
            }
        } else {
            for (k = 0; k < A->mt; k++) {
                RUNTIME_data_flush( sequence, A(k, m) );
            }
        }
    }

    (void)chamctxt;
}

/**
 *  Parallel tile matrix-matrix multiplication wrapper.
 */
void
chameleon_pzgemm( struct chameleon_pzgemm_s *ws,
                  cham_trans_t transA, cham_trans_t transB,
                  CHAMELEON_Complex64_t alpha, CHAM_desc_t *A, CHAM_desc_t *B,
                  CHAMELEON_Complex64_t beta,  CHAM_desc_t *C,
                  RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    cham_gemm_t alg = (ws != NULL) ? ws->alg : ChamGemmAlgGeneric;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    switch( alg ) {
    case ChamGemmAlgAuto:
    case ChamGemmAlgSummaB: /* Switch back to generic since it does not exist yet. */
    case ChamGemmAlgGeneric:
        chameleon_pzgemm_generic( chamctxt, transA, transB, alpha, A, B, beta, C, &options );
        break;

    case ChamGemmAlgSummaC:
        chameleon_pzgemm_summa( chamctxt, transA, transB, alpha, A, B, beta, C,
                                &(ws->WA), &(ws->WB), &options );
        break;

    case ChamGemmAlgSummaA:
        chameleon_pzgemm_Astat( chamctxt, transA, transB, alpha, A, B, beta, C, &options );
        break;
    }

    RUNTIME_options_finalize( &options, chamctxt );
}
