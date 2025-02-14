/**
 *
 * @file parsec/codelet_zipiv_allreduce.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon parsec codelets to do the reduction
 *
 * @version 1.3.0
 * @author Alycia Lisito
 * @date 2024-06-11
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"

void INSERT_TASK_zipiv_allreduce( const RUNTIME_option_t *options,
                                  CHAM_desc_t            *A,
                                  CHAM_ipiv_t            *ipiv,
                                  int                     k,
                                  int                     h,
                                  int                     n,
                                  void                   *ws )
{
    (void)options;
    (void)A;
    (void)ipiv;
    (void)k;
    (void)h;
    (void)n;
    (void)ws;
}
