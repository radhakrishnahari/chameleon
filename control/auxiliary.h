/**
 *
 * @file auxiliary.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon auxiliary header
 *
 * @version 1.3.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Loris Lucido
 * @date 2024-02-18
 *
 */
#ifndef _chameleon_auxiliary_h_
#define _chameleon_auxiliary_h_

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "chameleon/struct.h"
#include "chameleon/tasks.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  Internal routines
 */
void chameleon_warning      (const char *func_name, const char* msg_text);
void chameleon_error        (const char *func_name, const char* msg_text);
void chameleon_fatal_error  (const char *func_name, const char* msg_text);
int  chameleon_rank         (CHAM_context_t *chamctxt);
int  chameleon_tune         (cham_tasktype_t func, int M, int N, int NRHS);

#ifdef __cplusplus
}
#endif

#endif /* _chameleon_auxiliary_h_ */
