/**
 *
 * @file openmp/runtime_descriptor_ipiv.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon OpenMP descriptor routines
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Florent Pruvost
 * @author Matteo Marcos
 * @date 2025-10-15
 *
 */
#include "chameleon_openmp.h"

void RUNTIME_ipiv_create( CHAM_ipiv_t *ipiv )
{
    assert( 0 );
    (void)ipiv;
}

void RUNTIME_pivot_create( CHAM_desc_pivot_t *pivot )
{
    assert( 0 );
    (void)pivot;
}

void RUNTIME_ipiv_destroy( CHAM_ipiv_t *ipiv )
{
    assert( 0 );
    (void)ipiv;
}

void RUNTIME_pivot_destroy_submit( const RUNTIME_sequence_t *sequence,
                                   CHAM_desc_pivot_t        *pivot )
{
    assert( 0 );
    (void)sequence;
    (void)pivot;
}

void RUNTIME_pivot_destroy( CHAM_desc_pivot_t *pivot )
{
    assert( 0 );
    (void)pivot;
}

void *RUNTIME_ipiv_getaddr( const CHAM_ipiv_t *ipiv, int m )
{
    assert( 0 );
    (void)ipiv;
    (void)m;
    return NULL;
}

void *RUNTIME_nextpiv_getaddr( const CHAM_desc_pivot_t *pivot, int rank, int k, int h )
{
    assert( 0 );
    (void)pivot;
    (void)rank;
    (void)k;
    (void)h;
    return NULL;
}

void *RUNTIME_prevpiv_getaddr( const CHAM_desc_pivot_t *pivot, int rank, int k, int h )
{
    assert( 0 );
    (void)pivot;
    (void)rank;
    (void)k;
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

void RUNTIME_pivot_flushk( const RUNTIME_sequence_t *sequence,
                           const CHAM_desc_pivot_t *pivot, int rank )
{
    assert( 0 );
    (void)sequence;
    (void)pivot;
    (void)rank;
}

void RUNTIME_pivot_flush( const RUNTIME_sequence_t *sequence,
                          const CHAM_desc_pivot_t  *pivot )
{
    assert( 0 );
    (void)pivot;
    (void)sequence;
}

void RUNTIME_ipiv_flush( const RUNTIME_sequence_t *sequence,
                         const CHAM_ipiv_t        *ipiv )
{
    assert( 0 );
    (void)sequence;
    (void)ipiv;
}

void RUNTIME_ipiv_flushk( const RUNTIME_sequence_t *sequence,
                          const CHAM_ipiv_t *ipiv, int m )
{
    assert( 0 );
    (void)sequence;
    (void)ipiv;
    (void)m;
}

void RUNTIME_ipiv_invalidate( CHAM_desc_pivot_t *pivot,
                              int                k,
                              int                h,
                              int                myrank )
{
    (void)pivot;
    (void)k;
    (void)h;
    (void)myrank;
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

void RUNTIME_cpui_create ( CHAM_perm_t *ws )
{
    assert( 0 );
    (void)ws;
}

void *RUNTIME_cpui_getaddr( const CHAM_perm_t *ws,
                            int                m,
                            int                n )
{
    assert( 0 );
    (void)ws;
    (void)m;
    (void)n;
    return NULL;
}

void RUNTIME_cpui_destroy( CHAM_perm_t *ws )
{
    assert( 0 );
    (void)ws;
}

void RUNTIME_cpui_flushk( const RUNTIME_sequence_t *sequence,
                          int                       rank,
                          const CHAM_perm_t        *ws,
                          int                       m,
                          int                       n )
{
    assert( 0 );
    (void)sequence;
    (void)ws;
    (void)rank;
    (void)m;
    (void)n;
}

