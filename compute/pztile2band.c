/**
 *
 * @file pztile2band.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon ztile2band parallel algorithm
 *
 * @version 1.0.0
 * @author Azzam Haidar
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2020-03-03
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n)   A,  m, n
#define B(m, n)  B, m, n

/**
 *  Parallel copy of a band matrix from full NxN tile storage to band storage (LDABxN).
 */
void chameleon_pztile2band(cham_uplo_t uplo, CHAM_desc_t *A, CHAM_desc_t *B,
                       RUNTIME_sequence_t *sequence, RUNTIME_request_t *request)
{
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    int j;
    int tempjm, tempjn;
    int minmnt = chameleon_min(A->mt, A->nt);

    chamctxt = chameleon_context_self();
    if (sequence->status != CHAMELEON_SUCCESS) {
        return;
    }
    RUNTIME_options_init(&options, chamctxt, sequence, request);

    /* The code is actually incorrect due to the removal of the ld (Need new insert_task dedicated) */
    assert( 0 );

    /*
     *  ChamLower => Lower Band
     */
    if ( uplo == ChamLower ) {
       for (j = 0; j < minmnt; j++){
           /* Compute dimension on N with B since it is dimensioned with chameleon_min(A->m, A->n) */
           assert( A->m == B->n );
           assert( A->n >= B->n );

           tempjm = j == A->mt-1 ? A->m - j * A->mb : A->mb;
           tempjn = j == B->nt-1 ? B->n - j * B->nb : B->nb;

           INSERT_TASK_zlaset(
               &options,
               ChamUpperLower, B->mb, tempjn,
               0., 0.,
               B(0, j) );

           INSERT_TASK_zlacpy(
               &options,
               ChamLower, tempjm, tempjn, A->nb,
               A(j, j),
               B(0, j) );

           if( j<minmnt-1 ){
               tempjm = (j+1) == A->mt-1 ? A->m-(j+1)*A->mb : A->mb;
               INSERT_TASK_zlacpyx(
                   &options,
                   ChamUpper, tempjm, tempjn, A->nb,
                   0,     A(j+1, j),
                   A->nb, B(0,   j));
           }
       }
    }
    else if ( uplo == ChamUpper ) {
       for (j = 0; j < minmnt; j++){
           /* Compute dimension on M with B since it is dimensioned with chameleon_min(A->m, A->n) */
           assert( A->n == B->n );
           assert( A->m >= B->n );
           tempjn = j == A->nt-1 ? A->n - j * A->nb : A->nb;

           INSERT_TASK_zlaset(
               &options,
               ChamUpperLower, B->mb, tempjn,
               0., 0.,
               B(0, j) );

           if(j > 0){
               INSERT_TASK_zlacpy(
                   &options,
                   ChamLower, A->mb, tempjn, A->nb,
                   A(j-1, j),
                   B(0,   j));
           }

           tempjm = j == B->nt-1 ? B->n - j * B->nb : B->nb;
           INSERT_TASK_zlacpyx(
               &options,
               ChamUpper, tempjm, tempjn, A->nb,
               0,     A(j, j),
               A->nb, B(0, j));
       }
    }
    RUNTIME_options_finalize(&options, chamctxt);
}
#undef B
#undef A
