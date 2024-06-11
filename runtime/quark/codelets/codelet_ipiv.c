/**
 *
 * @file quark/codelet_ipiv.c
 *
 * @copyright 2023-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark codelets to convert pivot to permutations
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @date 2024-03-16
 *
 */
#include "chameleon_quark.h"
#include "chameleon/tasks.h"
#include "coreblas.h"

void INSERT_TASK_ipiv_init( const RUNTIME_option_t *options,
                            CHAM_ipiv_t *ipiv )
{
    assert( 0 );
    (void)options;
    (void)ipiv;
}

void INSERT_TASK_ipiv_reducek( const RUNTIME_option_t *options,
                               CHAM_ipiv_t *ipiv, int k, int h, int rank )
{
    assert( 0 );
    (void)options;
    (void)ipiv;
    (void)k;
    (void)h;
    (void)rank;
}

static inline void
CORE_ipiv_to_perm_quark( Quark *quark )
{
    int m0, m, k;
    int *ipiv, *perm, *invp;

    quark_unpack_args_6( quark, m0, m, k, ipiv, perm, invp );

    CORE_ipiv_to_perm( m0, m, k, ipiv, perm, invp );
}

void INSERT_TASK_ipiv_to_perm( const RUNTIME_option_t *options,
                               int m0, int m, int k,
                               const CHAM_ipiv_t *ipivdesc, int ipivk )
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);

    QUARK_Insert_Task(
        opt->quark, CORE_ipiv_to_perm_quark, (Quark_Task_Flags*)opt,
        sizeof(int),  &m0,  VALUE,
        sizeof(int),  &m,   VALUE,
        sizeof(int),  &k,   VALUE,
        sizeof(int*), RUNTIME_ipiv_getaddr( ipivdesc, ipivk ), INPUT,
        sizeof(int*), RUNTIME_perm_getaddr( ipivdesc, ipivk ), OUTPUT,
        sizeof(int*), RUNTIME_invp_getaddr( ipivdesc, ipivk ), OUTPUT,
        0 );
}
