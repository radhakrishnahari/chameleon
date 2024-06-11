/**
 *
 * @file quark/runtime_descriptor_ipiv.c
 *
 * @copyright 2022-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Quark descriptor routines
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Florent Pruvost
 * @date 2024-03-16
 *
 */
#include "chameleon_quark.h"

void RUNTIME_ipiv_create( CHAM_ipiv_t *ipiv,
                          const CHAM_desc_t *desc )
{
    assert( 0 );
    (void)ipiv;
    (void)desc;
}

void RUNTIME_ipiv_destroy( CHAM_ipiv_t *ipiv,
                           const CHAM_desc_t *desc )
{
    assert( 0 );
    (void)ipiv;
    (void)desc;
}

void *RUNTIME_ipiv_getaddr( const CHAM_ipiv_t *ipiv, int m )
{
    assert( 0 );
    (void)ipiv;
    (void)m;
    return NULL;
}

void *RUNTIME_nextpiv_getaddr( const CHAM_ipiv_t *ipiv, int rank, int m, int h )
{
    assert( 0 );
    (void)ipiv;
    (void)rank;
    (void)m;
    (void)h;
    return NULL;
}

void *RUNTIME_prevpiv_getaddr( const CHAM_ipiv_t *ipiv, int rank, int m, int h )
{
    assert( 0 );
    (void)ipiv;
    (void)rank;
    (void)m;
    (void)h;
    return NULL;
}

void *RUNTIME_perm_getaddr( const CHAM_ipiv_t *ipiv, int k )
{
    assert( 0 );
    (void)ipiv;
    (void)k;
    return NULL;
}

void *RUNTIME_invp_getaddr( const CHAM_ipiv_t *ipiv, int k )
{
    assert( 0 );
    (void)ipiv;
    (void)k;
    return NULL;
}

void RUNTIME_ipiv_flushk( const RUNTIME_sequence_t *sequence,
                          const CHAM_ipiv_t *ipiv, int m )
{
    assert( 0 );
    (void)sequence;
    (void)ipiv;
    (void)m;
}

void RUNTIME_ipiv_flush( const RUNTIME_sequence_t *sequence,
                         const CHAM_ipiv_t        *ipiv )
{
    assert( 0 );
    (void)ipiv;
    (void)sequence;
}

void RUNTIME_perm_flushk( const RUNTIME_sequence_t *sequence,
                          const CHAM_ipiv_t *ipiv, int m )
{
    assert( 0 );
    (void)sequence;
    (void)ipiv;
    (void)m;
}

void RUNTIME_ipiv_gather( const RUNTIME_sequence_t *sequence,
                          CHAM_ipiv_t *desc, int *ipiv, int node )
{
    assert( 0 );
    (void)sequence;
    (void)desc;
    (void)ipiv;
    (void)node;
}
