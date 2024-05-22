/**
 *
 * @file starpu/codelet_zgered.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgered StarPU codelet
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

void INSERT_TASK_zgered( const RUNTIME_option_t *options,
                         double threshold, double Anorm, int m, int n,
                         const CHAM_desc_t *A, int Am, int An )
{
    CHAM_tile_t          *tileA;
    double                u_low;
    int64_t               mm, nn;
#if defined(CHAMELEON_USE_MPI)
    int                   tag;
#endif
    starpu_data_handle_t *handleAin;
    starpu_data_handle_t  handleAout;

    CHAMELEON_BEGIN_ACCESS_DECLARATION;
    CHAMELEON_ACCESS_RW(A, Am, An);
    CHAMELEON_END_ACCESS_DECLARATION;

    /* Get the Input handle */
    mm = Am + (A->i / A->mb);
    nn = An + (A->j / A->nb);
    handleAin = A->schedopt;
    handleAin += ((int64_t)A->lmt) * nn + mm;

    assert( *handleAin != NULL );

    /*
     * Lets convert the tile precision based on the following criteria:
     *
     * ||A_{i,j}||_F  < u_{high} * || A ||_F / (nt * u_{low})
     * ||A_{i,j}||_F  < u_{high} * || A ||_F / nt *  1/ u_{low}
     * ||A_{i,j}||_F  < threshold / u_{low}
     */

    tileA = A->get_blktile( A, Am, An );
#if defined(CHAMELEON_USE_MPI)
    tag = starpu_mpi_data_get_tag( *handleAin );
#endif /* defined(CHAMELEON_USE_MPI) */

#if defined(CHAMELEON_USE_CUDA) && (CUDA_VERSION >= 7500)
#if defined(PRECISION_d)
    if ( options->withcuda ) {
        /*
         * Check for half precision
         */
        u_low = 1.e-4;
        if ( Anorm < (threshold / u_low) ) {
#if defined(CHAMELEON_DEBUG_GERED)
            fprintf( stderr,
                     "[%2d] Convert the tile ( %d, %d ) to half precision\n",
                     A->myrank, Am, An );
#endif
            starpu_cham_tile_register( &handleAout, -1, tileA, ChamComplexHalf );

            rt_starpu_insert_task(
                &cl_dlag2h,
                STARPU_VALUE,    &m,                 sizeof(int),
                STARPU_VALUE,    &n,                 sizeof(int),
                STARPU_R,        *handleAin,
                STARPU_W,         handleAout,
                STARPU_PRIORITY,  options->priority,
                STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "dlag2h",
#endif
                0);

            starpu_data_unregister_submit( *handleAin );
            *handleAin = handleAout;
            tileA->flttype = ChamComplexHalf;
#if defined(CHAMELEON_USE_MPI)
            starpu_mpi_data_register( handleAout, tag, tileA->rank );
#endif
            return;
        }
    }
#endif
#endif

    /*
     * Check for single precision
     */
#if !defined(CHAMELEON_SIMULATION)
    u_low = LAPACKE_slamch_work('e');
#else
    u_low = 1e-8;
#endif
    if ( Anorm < (threshold / u_low) ) {
#if defined(CHAMELEON_DEBUG_GERED)
        fprintf( stderr,
                 "[%2d] Convert the tile ( %d, %d ) to single precision\n",
                 A->myrank, Am, An );
#endif
        starpu_cham_tile_register( &handleAout, -1, tileA, ChamComplexFloat );

        rt_starpu_insert_task(
            &cl_zlag2c,
            STARPU_VALUE,    &m,                 sizeof(int),
            STARPU_VALUE,    &n,                 sizeof(int),
            STARPU_R,        *handleAin,
            STARPU_W,         handleAout,
            STARPU_PRIORITY,  options->priority,
            STARPU_EXECUTE_ON_WORKER, options->workerid,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "zlag2c",
#endif
            0);

        starpu_data_unregister_submit( *handleAin );
        *handleAin = handleAout;
        tileA->flttype = ChamComplexFloat;
#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_data_register( *handleAin, tag, tileA->rank );
#endif
        return;
    }
}
