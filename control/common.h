/**
 *
 * @file common.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon common header file
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2020-03-03
 *
 */
/**
 *  CHAMELEON facilities of interest to both CHAMELEON core developer
 *  and also of interest to CHAMELEON community contributor.
 */
#ifndef _chameleon_common_h_
#define _chameleon_common_h_


#if defined( _WIN32 ) || defined( _WIN64 )
#include <io.h>
#else
#include <unistd.h>
#endif

/**
 * Implementation headers
 */
#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#if defined(CHAMELEON_USE_CUBLAS_V2)
#include <cublas.h>
#include <cublas_v2.h>
#else
#include <cublas.h>
#endif
#endif

#if defined(CHAMELEON_USE_MPI)
#include <mpi.h>
#endif

/**
 *  Line to avoid conflict with other linear algebra libraries, because, we
 *  don't know why but lapacke provide a wrong interface of lapack in fortran
 */
#ifndef LAPACK_NAME
#define LAPACK_NAME(a, b) lapackef77_##a
#endif

/**
 *  Chameleon header files
 */
#include "chameleon.h"

#include "control/global.h"
#include "control/auxiliary.h"
#include "control/context.h"
#include "control/descriptor.h"
#include "control/async.h"

/**
 *  Global shortcuts
 */
#define CHAMELEON_RANK        chameleon_rank(chamctxt)
#define CHAMELEON_NB          chamctxt->nb
#define CHAMELEON_IB          chamctxt->ib
#define CHAMELEON_RHBLK       chamctxt->rhblock
#define CHAMELEON_TRANSLATION chamctxt->translation
#define CHAMELEON_PARALLEL    chamctxt->parallel_enabled
#define CHAMELEON_PROFILING   chamctxt->profiling_enabled

/**
 *  IPT internal define
 */
#define ChamIPT_NoDep   0
#define ChamIPT_Panel   1
#define ChamIPT_All     2

/**
 *  Global array of LAPACK constants
 *
 * TODO: this is needed to be able to call LAPACKE in Eigen values algorithms.
 * As soons as LAPACKE is not needed anymore in these algorithms this
 * declaration can be removed. This will also avoid a direct dependency between
 * chameleon lib and coreblas one which should not exist.
 */
extern char *chameleon_lapack_constants[];
#define chameleon_lapack_const(chameleon_const) chameleon_lapack_constants[chameleon_const][0]

#ifdef __cplusplus
extern "C" {
#endif

void chameleon_pmap( cham_uplo_t uplo, CHAM_desc_t *A,
                     cham_unary_operator_t operator, void *op_args,
                     RUNTIME_sequence_t *sequence, RUNTIME_request_t *request );

#include "control/compute_s.h"
#include "control/compute_d.h"
#include "control/compute_c.h"
#include "control/compute_z.h"

/*
void chameleon_pdlag2s(CHAM_context_t *chamctxt);
void chameleon_pzlag2c(CHAM_context_t *chamctxt);
void chameleon_pslag2d(CHAM_context_t *chamctxt);
void chameleon_pclag2z(CHAM_context_t *chamctxt);
*/

#ifdef __cplusplus
}
#endif

#endif /* _chameleon_common_h_ */
