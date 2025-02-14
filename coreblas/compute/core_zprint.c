/**
 *
 * @file core_zprint.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zprint CPU kernel
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "coreblas/lapacke.h"
#include "coreblas.h"

void
CORE_zprint( FILE *file, const char *header,
             cham_uplo_t uplo, int M, int N, int Am, int An,
             const CHAMELEON_Complex64_t *A, int lda )
{
    FILE *output = (file == NULL) ? stdout : file;
    int   i, j;

    fflush( output );
    fprintf( output, "--- %10s (%2d, %2d) / %p, %d:\n", header, Am, An, (void*)A, lda );
    for(i=0; i<M; i++) {
        fprintf( output, "  " );
        for(j=0; j<N; j++) {
#if defined(PRECISION_z) || defined(PRECISION_c)
            fprintf( output, " (% e, % e)",
                     creal( A[j*lda + i] ),
                     cimag( A[j*lda + i] ));
#else
            fprintf( output, " % e", A[j*lda + i] );
#endif
        }
        fprintf( output, "\n" );
    }
    fprintf( output, "-------------------------\n" );
    fflush( output );
}
