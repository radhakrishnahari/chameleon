/**
 *
 * @file parsec/codelet_ipiv.c
 *
 * @copyright 2023-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Parsec codelets to convert pivot to permutations
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Matteo Marcos
 * @date 2025-03-24
 *
 */
#include "chameleon_parsec.h"
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
                                 CHAM_ipiv_t *ipiv )
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

static inline int
CORE_ipiv_to_perm_parsec( parsec_execution_stream_t *context,
                          parsec_task_t             *this_task )
{
    int m0, m, k;
    int *ipiv, *perm, *invp;

    parsec_dtd_unpack_args(
        this_task, &m0, &m, &k, &ipiv, &perm, &invp );

    CORE_ipiv_to_perm( m0, m, k, 1, m, ipiv, perm, invp );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void INSERT_TASK_ipiv_to_perm( const RUNTIME_option_t *options,
                               int m0, int m, int k, int K1, int K2,
                               const CHAM_ipiv_t *ipivdesc, int ipivk )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_ipiv_to_perm_parsec, options->priority, "ipiv_to_perm",
        sizeof(int),         &m0,           VALUE,
        sizeof(int),         &m,            VALUE,
        sizeof(int),         &k,            VALUE,
        PASSED_BY_REF, RUNTIME_ipiv_getaddr( ipivdesc, ipivk ), chameleon_parsec_get_arena_index_ipiv( ipivdesc ) | INPUT,
        PASSED_BY_REF, RUNTIME_perm_getaddr( ipivdesc, ipivk ), chameleon_parsec_get_arena_index_perm( ipivdesc ) | OUTPUT,
        PASSED_BY_REF, RUNTIME_invp_getaddr( ipivdesc, ipivk ), chameleon_parsec_get_arena_index_invp( ipivdesc ) | OUTPUT,
        PARSEC_DTD_ARG_END );

    (void)K1;
    (void)K2;
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
