/**
 *
 * @file quark/codelet_zlaswp.c
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark codelets to apply zlaswp on a panel
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @date 2025-03-24
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks_z.h"
#include "coreblas/coreblas_ztile.h"

static void CORE_zlaswp_get_quark( Quark *quark )
{
    int          m0, k, *perm;
    CHAM_tile_t *A, *B;

    quark_unpack_args_5( quark, m0, k, perm, A, B );

    TCORE_zlaswp_get( m0, A->m, A->n, k, A, B, perm );
}

void INSERT_TASK_zlaswp_get( const RUNTIME_option_t *options,
                             cham_dir_t dir, int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *U, int Um, int Un )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LASWP;

    QUARK_Insert_Task(
        opt->quark, CORE_zlaswp_get_quark, (Quark_Task_Flags*)opt,
        sizeof(int),          &m0, VALUE,
        sizeof(int),          &k,  VALUE,
        sizeof(int*),         RUNTIME_perm_getaddr( ipiv, ipivk ),     INPUT,
        sizeof(CHAM_tile_t*), RTBLKADDR(A, ChamComplexDouble, Am, An), INPUT,
        sizeof(CHAM_tile_t*), RTBLKADDR(U, ChamComplexDouble, Um, Un), INOUT,
        0 );

    (void)dir;
    (void)m;
    (void)n;
}

static void CORE_zlaswp_set_quark( Quark *quark )
{
    int          m0, k, *invp;
    CHAM_tile_t *A, *B;

    quark_unpack_args_5( quark, m0, k, invp, A, B );

    TCORE_zlaswp_set( m0, A->m, A->n, k, A, B, invp );
}

void INSERT_TASK_zlaswp_set( const RUNTIME_option_t *options,
                             cham_dir_t dir, int m0, int m, int n, int k,
                             const CHAM_ipiv_t *ipiv, int ipivk,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LASWP;

    QUARK_Insert_Task(
        opt->quark, CORE_zlaswp_set_quark, (Quark_Task_Flags*)opt,
        sizeof(int),          &m0, VALUE,
        sizeof(int),          &k,  VALUE,
        sizeof(int*),         RUNTIME_invp_getaddr( ipiv, ipivk ),     INPUT,
        sizeof(CHAM_tile_t*), RTBLKADDR(A, ChamComplexDouble, Am, An), INPUT,
        sizeof(CHAM_tile_t*), RTBLKADDR(B, ChamComplexDouble, Bm, Bn), INOUT,
        0 );

    (void)dir;
    (void)m;
    (void)n;
}
