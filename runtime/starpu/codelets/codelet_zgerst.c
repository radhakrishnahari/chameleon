/**
 *
 * @file starpu/codelet_zgerst.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgerst StarPU codelet
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @date 2024-07-17
 * @precisions normal z -> d
 *
 */
#include "chameleon_starpu.h"
#include <coreblas/lapacke.h>
#include "runtime_codelet_zc.h"
#include "runtime_codelet_z.h"

void INSERT_TASK_zgerst( const RUNTIME_option_t *options,
                         int m, int n,
                         const CHAM_desc_t *A, int Am, int An )
{
    CHAM_tile_t          *tileA;
    int64_t               mm, nn;
#if defined(CHAMELEON_USE_MPI)
    int                   tag;
#endif
    starpu_data_handle_t *handleAin;
    starpu_data_handle_t  handleAout;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    tileA = A->get_blktile( A, Am, An );
    if ( tileA->flttype == ChamComplexDouble ) {
        return;
    }

    /* Get the Input handle */
    mm = Am + (A->i / A->mb);
    nn = An + (A->j / A->nb);
    handleAin = A->schedopt;
    handleAin += ((int64_t)A->lmt) * nn + mm;

    assert( *handleAin != NULL );

#if defined(CHAMELEON_USE_MPI)
    tag = starpu_mpi_data_get_tag( *handleAin );
#endif /* defined(CHAMELEON_USE_MPI) */

    starpu_cham_tile_register( &handleAout, -1, tileA, ChamComplexDouble );

    switch( tileA->flttype ) {
#if defined(CHAMELEON_USE_CUDA) && (CUDA_VERSION >= 7500)
#if defined(PRECISION_d)
    /*
     * Restore from half precision
     */
    case ChamComplexHalf:
#if defined(CHAMELEON_DEBUG_GERED)
        fprintf( stderr,
                 "[%2d] Convert back the tile ( %d, %d ) from half precision\n",
                 A->myrank, Am, An );
#endif
        rt_shm_starpu_insert_task(
            &cl_hlag2d,
            STARPU_VALUE,    &m,                 sizeof(int),
            STARPU_VALUE,    &n,                 sizeof(int),
            STARPU_R,        *handleAin,
            STARPU_W,         handleAout,
            STARPU_PRIORITY,  options->priority,
            STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "hlag2d",
#endif
            0);
        break;
#endif
#endif

    case ChamComplexFloat:
#if defined(CHAMELEON_DEBUG_GERED)
        fprintf( stderr,
                 "[%2d] Convert back the tile ( %d, %d ) from half precision\n",
                 A->myrank, Am, An );
#endif
        rt_shm_starpu_insert_task(
            &cl_clag2z,
            STARPU_VALUE,    &m,                 sizeof(int),
            STARPU_VALUE,    &n,                 sizeof(int),
            STARPU_R,        *handleAin,
            STARPU_W,         handleAout,
            STARPU_PRIORITY,  options->priority,
            STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "clag2z",
#endif
            0);
        break;

    default:
        fprintf( stderr, "ERROR: Unknonw input datatype" );
    }

    starpu_data_unregister_submit( *handleAin );
    *handleAin = handleAout;
    tileA->flttype = ChamComplexDouble;
#if defined(CHAMELEON_USE_MPI)
    starpu_mpi_data_register( handleAout, tag, tileA->rank );
#endif
}
