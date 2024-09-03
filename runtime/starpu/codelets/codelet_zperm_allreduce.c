/**
 *
 * @file starpu/codelet_zperm_allreduce.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelets to do the reduction
 *
 * @version 1.3.0
 * @author Alycia Lisito
 * @date 2024-06-11
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu_internal.h"
#include "runtime_codelet_z.h"
#include <coreblas/cblas_wrapper.h>

#if defined(CHAMELEON_USE_MPI)
void
INSERT_TASK_zperm_allreduce( const RUNTIME_option_t *options,
                             const CHAM_desc_t      *A,
                             CHAM_ipiv_t            *ipiv,
                             int                     ipivk,
                             int                     k,
                             int                     n,
                             CHAM_desc_t            *U,
                             int                     Um,
                             int                     Un,
                             void                   *ws )
{
    struct chameleon_pzgetrf_s *tmp = (struct chameleon_pzgetrf_s *)ws;
}

void
INSERT_TASK_zperm_allreduce_send_A( const RUNTIME_option_t *options,
                                    CHAM_desc_t            *A,
                                    int                     Am,
                                    int                     An,
                                    int                     myrank,
                                    int                     np,
                                    int                    *proc_involved )
{
    int p, rank;

    for ( p = 0; p < np; p ++ ) {
        if ( proc_involved[ p ] == myrank ) {
            continue;
        }
        starpu_mpi_get_data_on_node_detached( options->sequence->comm,
                                              RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
                                              proc_involved[ p ], NULL, NULL );
    }
}

void
INSERT_TASK_zperm_allreduce_send_perm( const RUNTIME_option_t *options,
                                       CHAM_ipiv_t            *ipiv,
                                       int                     ipivk,
                                       int                     myrank,
                                       int                     np,
                                       int                    *proc_involved )
{
    int p;

    for ( p = 0; p < np; p++ ) {
        if ( proc_involved[ p ] == myrank ) {
            continue;
        }
        starpu_mpi_get_data_on_node_detached( options->sequence->comm,
                                              RUNTIME_perm_getaddr( ipiv, ipivk ),
                                              proc_involved[ p ], NULL, NULL );
    }
}

void
INSERT_TASK_zperm_allreduce_send_invp( const RUNTIME_option_t *options,
                                       CHAM_ipiv_t            *ipiv,
                                       int                     ipivk,
                                       const CHAM_desc_t      *A,
                                       int                     k,
                                       int                     n )
{
    int b, rank;

    for ( b = k+1; (b < A->mt) && ((b-(k+1)) < A->p); b ++ ) {
        rank = A->get_rankof( A, b, n );
        if ( rank == A->myrank ) {
            continue;
        }
        starpu_mpi_get_data_on_node_detached( options->sequence->comm,
                                              RUNTIME_invp_getaddr( ipiv, ipivk ),
                                              rank, NULL, NULL );
    }
}
#else
void
INSERT_TASK_zperm_allreduce_send_A( const RUNTIME_option_t *options,
                                    CHAM_desc_t            *A,
                                    int                     Am,
                                    int                     An,
                                    int                     myrank,
                                    int                     np,
                                    int                    *proc_involved  )
{
    (void)options;
    (void)A;
    (void)Am;
    (void)An;
    (void)myrank;
    (void)np;
    (void)proc_involved;
}

void
INSERT_TASK_zperm_allreduce_send_perm( const RUNTIME_option_t *options,
                                       CHAM_ipiv_t            *ipiv,
                                       int                     ipivk,
                                       int                     myrank,
                                       int                     np,
                                       int                    *proc_involved  )
{
    (void)options;
    (void)ipiv;
    (void)ipivk;
    (void)myrank;
    (void)np;
    (void)proc_involved;
}

void
INSERT_TASK_zperm_allreduce_send_invp( const RUNTIME_option_t *options,
                                       CHAM_ipiv_t            *ipiv,
                                       int                     ipivk,
                                       const CHAM_desc_t      *A,
                                       int                     k,
                                       int                     n )
{
    (void)options;
    (void)ipiv;
    (void)ipivk;
    (void)A;
    (void)k;
    (void)n;
}

void
INSERT_TASK_zperm_allreduce( const RUNTIME_option_t *options,
                             const CHAM_desc_t      *A,
                             CHAM_ipiv_t            *ipiv,
                             int                     ipivk,
                             int                     k,
                             int                     n,
                             CHAM_desc_t            *U,
                             int                     Um,
                             int                     Un,
                             void                   *ws )
{
    (void)options;
    (void)A;
    (void)ipiv;
    (void)ipivk;
    (void)k;
    (void)n;
    (void)U;
    (void)Um;
    (void)Un;
    (void)ws;
}
#endif
