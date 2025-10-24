/**
 *
 * @file openmp/codelet_ipiv.c
 *
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon OpenMP codelets to convert pivot to permutations
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-03-24
 *
 */
#include "chameleon_openmp.h"
#include "chameleon/tasks.h"
#include "coreblas.h"

void INSERT_TASK_ipiv_init( const RUNTIME_option_t *options,
                            CHAM_ipiv_t *ipiv )
{
    assert( 0 );
    (void)options;
    (void)ipiv;
}

void INSERT_TASK_ipiv_init_data( const RUNTIME_option_t *options,
                                 CHAM_ipiv_t            *ipiv )
{
    assert( 0 );
    (void)options;
    (void)ipiv;
}

void INSERT_TASK_ipiv_reducek( const RUNTIME_option_t *options,
                               CHAM_desc_pivot_t *pivot, int k, int h, int rank )
{
    assert( 0 );
    (void)options;
    (void)pivot;
    (void)k;
    (void)h;
    (void)rank;
}

void INSERT_TASK_ipiv_to_perm( const RUNTIME_option_t *options,
                               int m0, int m, int k, int K1, int K2,
                               const CHAM_ipiv_t *ipivdesc, int ipivk )
{
    int *ipiv = NULL; // get pointer from ipivdesc
    int *perm = NULL; // get pointer from ipivdesc
    int *invp = NULL; // get pointer from ipivdesc

#pragma omp task firstprivate( m0, m, k ) depend( in:ipiv[0] ) depend( inout:perm[0] ) depend( inout:invp[0] )
    {
        CORE_ipiv_to_perm( m0, m, k, 1, m, ipiv, perm, invp );
    }

    (void)options;
    (void)K1;
    (void)K2;
    (void)ipivk;
}

void INSERT_TASK_perm_init( const RUNTIME_option_t *options,
                            cham_dir_t              dir,
                            CHAM_ipiv_t            *ipiv,
                            int                    *PERM )
{
    assert( 0 );
    (void)options;
    (void)dir;
    (void)ipiv;
    (void)PERM;
}
