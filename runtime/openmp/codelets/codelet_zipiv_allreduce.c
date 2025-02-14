/**
 *
 * @file openmp/codelet_zipiv_allreduce.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon openmp codelets to do the reduction
 *
 * @version 1.3.0
 * @author Alycia Lisito
 * @date 2024-06-11
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_openmp.h"

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
