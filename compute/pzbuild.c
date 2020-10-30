/**
 *
 * @file pzbuild.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zbuild parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for CHAMELEON 0.9.2
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Guillaume Sylvand
 * @date 2020-03-03
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m, n) A,  m,  n
/**
 *  Parallel tile matrix generation
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *            = ChamUpperLower: All the matrix A
 *            = ChamUpper: Upper triangular part
 *            = ChamLower: Lower triangular part
 *
 * @param[in] A
 *          On exit, The matrix A generated.
 *
 * @param[in] user_data
 *          The data used in the matrix generation.
 *
 * @param[in] user_build_callback
 *          The function called by the codelet to fill the tiles
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 */
void chameleon_pzbuild( cham_uplo_t uplo, CHAM_desc_t *A, void *user_data, void* user_build_callback,
                    RUNTIME_sequence_t *sequence, RUNTIME_request_t *request )
{
  CHAM_context_t *chamctxt;
  RUNTIME_option_t options;

  int m, n;

  chamctxt = chameleon_context_self();
  if (sequence->status != CHAMELEON_SUCCESS)
    return;
  RUNTIME_options_init(&options, chamctxt, sequence, request);

  for (m = 0; m < A->mt; m++) {
    for (n = 0; n < A->nt; n++) {

      if ( ( uplo == ChamUpper && m <= n ) ||
           ( uplo == ChamLower && m >= n ) ||
           ( uplo == ChamUpperLower ) )
        INSERT_TASK_zbuild(
              &options,
              A(m, n),
              user_data, user_build_callback );
    }
  }

  RUNTIME_options_finalize( &options, chamctxt);
}
