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

#define C(m, n) C,  m,  n
#define WA(m, n) &WA,  m,  n
#define WB(m, n) &WB,  m,  n

/**
 *  chameleon_pzplrnk - Generate a random rank-k matrix by tiles.
 */
void chameleon_pzplrnk( int K, CHAM_desc_t *C,
                        unsigned long long int seedA,
                        unsigned long long int seedB,
                        RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    int m, n, k, KT;
    int tempmm, tempnn, tempkk;
    int myp, myq;
    CHAMELEON_Complex64_t zbeta;
    CHAM_desc_t WA, WB;

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init( &options, chamctxt, sequence, request );

    chameleon_desc_init( &WA, CHAMELEON_MAT_ALLOC_TILE,
                         ChamComplexDouble, C->mb, C->nb, (C->mb * C->nb),
                         C->mt * C->mb, C->nb * chameleon_desc_datadist_get_iparam(C, 1), 0, 0,
                         C->mt * C->mb, C->nb * chameleon_desc_datadist_get_iparam(C, 1),
                         chameleon_desc_datadist_get_iparam(C, 0),
                         chameleon_desc_datadist_get_iparam(C, 1),
                         NULL, NULL, NULL, NULL );
    chameleon_desc_init( &WB, CHAMELEON_MAT_ALLOC_TILE,
                         ChamComplexDouble, C->mb, C->nb, (C->mb * C->nb),
                         C->mb * chameleon_desc_datadist_get_iparam(C, 0), C->nt * C->nb, 0, 0,
                         C->mb * chameleon_desc_datadist_get_iparam(C, 0), C->nt * C->nb,
                         chameleon_desc_datadist_get_iparam(C, 0),
                         chameleon_desc_datadist_get_iparam(C, 1),
                         NULL, NULL, NULL, NULL );

    KT = (K + C->mb - 1) / C->mb;
    myp = C->myrank / chameleon_desc_datadist_get_iparam(C, 1);
    myq = C->myrank % chameleon_desc_datadist_get_iparam(C, 1);

    for (k = 0; k < KT; k++) {
        tempkk = k == KT-1 ? K - k * WA.nb : WA.nb;
        zbeta  = k == 0 ? 0. : 1.;

        for (n = myq; n < C->nt; n+=chameleon_desc_datadist_get_iparam(C, 1)) {
            tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

            INSERT_TASK_zplrnt(
                &options,
                tempkk, tempnn, WB(myp, n),
                WB.m, k * WB.mb, n * WB.nb, seedB );
        }

        for (m = myp; m < C->mt; m+=chameleon_desc_datadist_get_iparam(C, 0)) {
            tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;

            INSERT_TASK_zplrnt(
                &options,
                tempmm, tempkk, WA(m, myq),
                WA.m, m * WA.mb, k * WA.nb, seedA );

            for (n = myq; n < C->nt; n+=chameleon_desc_datadist_get_iparam(C, 1)) {
                tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;

                INSERT_TASK_zgemm(
                    &options,
                    ChamNoTrans, ChamNoTrans,
                    tempmm, tempnn, tempkk, C->mb,
                    1.,    WA(m, myq),
                           WB(myp, n),
                    zbeta,  C(m, n));
            }
            RUNTIME_data_flush( sequence, WA(m, 0) );
        }
        for (n = myq; n < C->nt; n+=chameleon_desc_datadist_get_iparam(C, 1)) {
            RUNTIME_data_flush( sequence, WB(0, n) );
        }
    }

    RUNTIME_desc_flush( &WA, sequence );
    RUNTIME_desc_flush( &WB, sequence );
    RUNTIME_desc_flush(  C,  sequence );
    chameleon_sequence_wait( chamctxt, sequence );
    chameleon_desc_destroy( &WA );
    chameleon_desc_destroy( &WB );

    RUNTIME_options_finalize( &options, chamctxt );
}
