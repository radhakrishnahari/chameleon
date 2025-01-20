/**
 *
 * @file pzplrnk.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplrnk parallel algorithm
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @author Lionel Eyraud-Dubois
 * @date 2023-07-05
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define WA(m, n) WA, m,  n
#define WB(m, n) WB, m,  n
#define C(m, n)  C,  m,  n

/**
 *  chameleon_pzplrnk - Generate a random rank-k matrix by tiles.
 */
static inline void
chameleon_pzplrnk_generic( CHAM_context_t         *chamctxt,
                           int                     K,
                           CHAM_desc_t            *WA,
                           CHAM_desc_t            *WB,
                           CHAM_desc_t            *C,
                           unsigned long long int  seedA,
                           unsigned long long int  seedB,
                           RUNTIME_option_t       *options )
{
    RUNTIME_sequence_t   *sequence = options->sequence;
    CHAMELEON_Complex64_t zbeta;
    int  m, n, k, KT;
    int  tempmm, tempnn, tempkk;
    int  myrank  = RUNTIME_comm_rank( chamctxt );
    int  initA;
    int *initB = malloc( C->nt * sizeof(int) );

    KT = (K + C->mb - 1) / C->mb;

    for (k = 0; k < KT; k++) {
        tempkk = k == KT-1 ? K - k * WA->nb : WA->nb;
        zbeta  = k == 0 ? 0. : 1.;

        memset( initB, 0, C->nt * sizeof(int) );

        for (m = 0; m < C->mt; m++) {
            tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

            initA = 0;

            for (n = 0; n < C->nt; n++) {
                tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                if ( C->get_rankof( C(m, n) ) == myrank ) {
                    if ( !initA ) {
                        INSERT_TASK_zplrnt(
                            options,
                            tempmm, tempkk, WA(m, myrank),
                            WA->m, m * WA->mb, k * WA->nb, seedA );
                        initA = 1;
                    }
                    if ( !initB[n] ) {
                        INSERT_TASK_zplrnt(
                            options,
                            tempkk, tempnn, WB(myrank, n),
                            WB->m, k * WB->mb, n * WB->nb, seedB );
                        initB[n] = 1;
                    }

                    INSERT_TASK_zgemm(
                        options,
                        ChamNoTrans, ChamNoTrans,
                        tempmm, tempnn, tempkk, C->mb,
                        1.,    WA(m, myrank),
                               WB(myrank, n),
                        zbeta, C(m, n));
                }
            }
            if ( initA ) {
                RUNTIME_data_flush( sequence, WA(m, myrank) );
            }
        }
        for (n = 0; n < C->nt; n++) {
            if ( initB[n] ) {
                RUNTIME_data_flush( sequence, WB(myrank, n) );
            }
        }
    }

    free( initB );
}

/**
 *  chameleon_pzplrnk - Generate a random rank-k matrix by tiles on a 2dbc grid.
 */
static inline void
chameleon_pzplrnk_2dbc( CHAM_context_t         *chamctxt,
                        int                     K,
                        CHAM_desc_t            *WA,
                        CHAM_desc_t            *WB,
                        CHAM_desc_t            *C,
                        unsigned long long int  seedA,
                        unsigned long long int  seedB,
                        RUNTIME_option_t       *options )
{
    RUNTIME_sequence_t   *sequence = options->sequence;
    CHAMELEON_Complex64_t zbeta;
    int m, n, k, KT;
    int tempmm, tempnn, tempkk;
    int p, q, myp, myq;

    KT  = (K + C->mb - 1) / C->mb;
    p   = chameleon_desc_datadist_get_iparam( C, 0 );
    q   = chameleon_desc_datadist_get_iparam( C, 1 );
    myp = C->myrank / q;
    myq = C->myrank % q;

    for (k = 0; k < KT; k++) {
        tempkk = k == KT-1 ? K - k * WA->nb : WA->nb;
        zbeta  = k == 0 ? 0. : 1.;

        for (n = myq; n < C->nt; n += q) {
            tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

            INSERT_TASK_zplrnt(
                options,
                tempkk, tempnn, WB(myp, n),
                WB->m, k * WB->mb, n * WB->nb, seedB );
        }

        for (m = myp; m < C->mt; m += p) {
            tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

            INSERT_TASK_zplrnt(
                options,
                tempmm, tempkk, WA(m, myq),
                WA->m, m * WA->mb, k * WA->nb, seedA );

            for (n = myq; n < C->nt; n+=q) {
                tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                INSERT_TASK_zgemm(
                    options,
                    ChamNoTrans, ChamNoTrans,
                    tempmm, tempnn, tempkk, C->mb,
                    1.,    WA(m, myq),
                           WB(myp, n),
                    zbeta,  C(m, n));
            }
            RUNTIME_data_flush( sequence, WA(m, myq) );
        }
        for (n = myq; n < C->nt; n+=q) {
            RUNTIME_data_flush( sequence, WB(myp, n) );
        }
    }
}

/**
 * Rank-k matrix generator.
 */
void
chameleon_pzplrnk( int                         K,
                   CHAM_desc_t                *C,
                   unsigned long long int      seedA,
                   unsigned long long int      seedB,
                   RUNTIME_sequence_t         *sequence,
                   RUNTIME_request_t          *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    CHAM_desc_t WA, WB;
    int p, q;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    p = chameleon_desc_datadist_get_iparam( C, 0 );
    q = chameleon_desc_datadist_get_iparam( C, 1 );
    if ( ( chamctxt->generic_enabled != CHAMELEON_TRUE )  &&
         ( C->get_rankof_init == chameleon_getrankof_2d ) &&
         ( (chameleon_desc_datadist_get_iparam(C, 0) != 1) ||
           (chameleon_desc_datadist_get_iparam(C, 1) != 1) ) )
    {
        chameleon_desc_init( &WA, CHAMELEON_MAT_ALLOC_TILE,
                             ChamComplexDouble, C->mb, C->nb, (C->mb * C->nb),
                             C->m, C->nb * q, 0, 0,
                             C->m, C->nb * q, p, q,
                             NULL, NULL, NULL, NULL );
        chameleon_desc_init( &WB, CHAMELEON_MAT_ALLOC_TILE,
                             ChamComplexDouble, C->mb, C->nb, (C->mb * C->nb),
                             C->mb * p, C->n, 0, 0,
                             C->mb * p, C->n, p, q,
                             NULL, NULL, NULL, NULL );

        chameleon_pzplrnk_2dbc( chamctxt, K, &WA, &WB, C, seedA, seedB, &options );
    }
    else {
        int np = p * q;
        chameleon_desc_init( &WA, CHAMELEON_MAT_ALLOC_TILE,
                             ChamComplexDouble, C->mb, C->nb, (C->mb * C->nb),
                             C->m, C->nb * np, 0, 0,
                             C->m, C->nb * np, 1, np,
                             NULL, NULL, NULL, NULL );
        chameleon_desc_init( &WB, CHAMELEON_MAT_ALLOC_TILE,
                             ChamComplexDouble, C->mb, C->nb, (C->mb * C->nb),
                             C->mb * np, C->n, 0,  0,
                             C->mb * np, C->n, np, 1,
                             NULL, NULL, NULL, NULL );

        chameleon_pzplrnk_generic( chamctxt, K, &WA, &WB, C, seedA, seedB, &options );
    }

    RUNTIME_desc_flush( &WA, sequence );
    RUNTIME_desc_flush( &WB, sequence );
    RUNTIME_desc_flush(  C,  sequence );

    chameleon_sequence_wait( chamctxt, sequence );
    chameleon_desc_destroy( &WA );
    chameleon_desc_destroy( &WB );

    RUNTIME_options_finalize( &options, chamctxt );
}
