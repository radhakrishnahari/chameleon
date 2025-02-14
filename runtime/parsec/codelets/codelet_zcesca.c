/**
 *
 * @file parsec/codelet_zcesca.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zcesca PaRSEC codelet
 *
 * @version 1.3.0
 * @author Florent Pruvost
 * @date 2024-02-18
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zcesca_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int center, scale;
    cham_store_t axis;
    int m, n, mt, nt;
    CHAMELEON_Complex64_t *Gi;
    int ldgi;
    CHAMELEON_Complex64_t *Gj;
    int ldgj;
    CHAMELEON_Complex64_t *G;
    double *Di;
    int lddi;
    double *Dj;
    int lddj;
    CHAMELEON_Complex64_t *A;
    int lda;

    parsec_dtd_unpack_args(
        this_task, &center, &scale, &axis, &m, &n, &mt, &nt, &Gi, &ldgi, &Gj, &ldgj, &G, &Di, &lddi, &Dj, &lddj, &A, &lda );

    CORE_zcesca( center, scale, axis,
                 m, n, mt, nt,
                 Gi, ldgi,
                 Gj, ldgj,
                 G,
                 Di, lddi,
                 Dj, lddj,
                 A, lda);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}


void INSERT_TASK_zcesca( const RUNTIME_option_t *options,
                         int center, int scale, cham_store_t axis,
                         int m, int n, int mt, int nt,
                         const CHAM_desc_t *Gi, int Gim, int Gin,
                         const CHAM_desc_t *Gj, int Gjm, int Gjn,
                         const CHAM_desc_t *G, int Gm, int Gn,
                         const CHAM_desc_t *Di, int Dim, int Din,
                         const CHAM_desc_t *Dj, int Djm, int Djn,
                         CHAM_desc_t *A, int Am, int An)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    CHAM_tile_t *tileGi = Gi->get_blktile( Gi, Gim, Gin );
    CHAM_tile_t *tileGj = Gj->get_blktile( Gj, Gjm, Gjn );
    CHAM_tile_t *tileDi = Di->get_blktile( Di, Dim, Din );
    CHAM_tile_t *tileDj = Dj->get_blktile( Dj, Djm, Djn );
    CHAM_tile_t *tileA  = A->get_blktile( A, Am, An );
    parsec_dtd_insert_task(
        PARSEC_dtd_taskpool, CORE_zcesca_parsec, options->priority, PARSEC_DEV_CPU, "cesca",
        sizeof(int),   &center, PARSEC_VALUE,
        sizeof(int),   &scale, PARSEC_VALUE,
        sizeof(int),   &axis, PARSEC_VALUE,
        sizeof(int),   &m,    PARSEC_VALUE,
        sizeof(int),   &n,    PARSEC_VALUE,
        sizeof(int),   &mt,   PARSEC_VALUE,
        sizeof(int),   &nt,   PARSEC_VALUE,
        PASSED_BY_REF, RTBLKADDR( Gi, CHAMELEON_Complex64_t, Gim, Gin ), chameleon_parsec_get_arena_index( Gi ) | PARSEC_INPUT,
        sizeof(int), &(tileGi->ld), PARSEC_VALUE,
        PASSED_BY_REF, RTBLKADDR( Gj, CHAMELEON_Complex64_t, Gjm, Gjn ), chameleon_parsec_get_arena_index( Gj ) | PARSEC_INPUT,
        sizeof(int), &(tileGj->ld), PARSEC_VALUE,
        PASSED_BY_REF, RTBLKADDR( G, CHAMELEON_Complex64_t, Gm, Gn ), chameleon_parsec_get_arena_index( G ) | PARSEC_INPUT,
        PASSED_BY_REF, RTBLKADDR( Di, double, Dim, Din ), chameleon_parsec_get_arena_index( Di ) | PARSEC_INPUT,
        sizeof(int), &(tileDi->ld), PARSEC_VALUE,
        PASSED_BY_REF, RTBLKADDR( Dj, double, Djm, Djn ), chameleon_parsec_get_arena_index( Dj ) | PARSEC_INPUT,
        sizeof(int), &(tileDj->ld), PARSEC_VALUE,
        PASSED_BY_REF, RTBLKADDR( A, CHAMELEON_Complex64_t, Am, An ), chameleon_parsec_get_arena_index( A ) | PARSEC_INOUT | PARSEC_AFFINITY,
        sizeof(int), &(tileA->ld), PARSEC_VALUE,
        PARSEC_DTD_ARG_END );
}
