/**
 *
 * @file starpu/runtime_descriptor_ipiv.c
 *
 * @copyright 2022-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU descriptor routines
 *
 * @version 1.3.0
 * @author Mathieu Faverge
 * @author Matthieu Kuhn
 * @author Alycia Lisito
 * @author Florent Pruvost
 * @author Pierre Esterie
 * @author Matteo Marcos
 * @author Samuel Thibault
 * @date 2025-07-15
 *
 */
#include "chameleon_starpu_internal.h"

/**
 *  Create ws_pivot runtime structures
 */
void RUNTIME_ipiv_create( CHAM_ipiv_t *ipiv )
{
    assert( ipiv );
    size_t                nbhandles = 3 * ipiv->mt;
    starpu_data_handle_t *handles   = calloc( nbhandles, sizeof(starpu_data_handle_t) );
    ipiv->ipiv    = handles;
    handles += ipiv->mt;
    ipiv->perm    = handles;
    handles += ipiv->mt;
    ipiv->invp    = handles;
#if defined(CHAMELEON_USE_MPI)
    /*
     * Book the number of tags required to describe pivot structure
     * One per handle type
     */
    {
        chameleon_starpu_tag_init();
        ipiv->mpitag_ipiv = chameleon_starpu_tag_book( nbhandles );
        if ( ipiv->mpitag_ipiv == -1 ) {
            chameleon_fatal_error("RUNTIME_ipiv_create", "Can't pursue computation since no more tags are available for ipiv structure");
            return;
        }
        ipiv->mpitag_perm    = ipiv->mpitag_ipiv + ipiv->mt;
        ipiv->mpitag_invp    = ipiv->mpitag_perm + ipiv->mt;
    }
#endif
}

/**
 *  Create ws_pivot runtime structures
 */
void RUNTIME_pivot_create( CHAM_desc_pivot_t *pivot )
{
    assert( pivot );
    size_t                nbhandles = 2 * pivot->P;
    starpu_data_handle_t *handles   = calloc( nbhandles, sizeof(starpu_data_handle_t) );
    pivot->nextpiv = handles;
    handles += pivot->P;
    pivot->prevpiv = handles;
#if defined(CHAMELEON_USE_MPI)
    /*
     * Book the number of tags required to describe pivot structure
     * One per handle type
     */
    {
        chameleon_starpu_tag_init();
        pivot->mpitag_nextpiv = chameleon_starpu_tag_book( nbhandles );
        if ( pivot->mpitag_nextpiv == -1 ) {
            chameleon_fatal_error("RUNTIME_pivot_create", "Can't pursue computation since no more tags are available for pivot structure");
            return;
        }
        pivot->mpitag_prevpiv = pivot->mpitag_nextpiv + pivot->P;
    }
#endif
}

/**
 *  Destroy ws_pivot runtime structures
 */
void RUNTIME_ipiv_destroy( CHAM_ipiv_t *ipiv )
{
    int                   i;
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(ipiv->ipiv);
    size_t                nbhandles = 3 * ipiv->mt;

    if ( handle ) {
        for(i=0; i<nbhandles; i++) {
            if ( *handle != NULL ) {
                starpu_data_unregister( *handle );
                *handle = NULL;
            }
            handle++;
        }

        free( ipiv->ipiv    );
        ipiv->ipiv    = NULL;
        ipiv->perm    = NULL;
        ipiv->invp    = NULL;
        chameleon_starpu_tag_release( ipiv->mpitag_ipiv );
    }
}

/**
 *  Asynchronously destroy ws_pivot runtime structures
 */
void RUNTIME_pivot_destroy_submit( const RUNTIME_sequence_t *sequence,
                                   CHAM_desc_pivot_t        *pivot )
{
    int                   i;
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(pivot->nextpiv);
    size_t                nbhandles = 2 * pivot->P;

    if ( !handle ) {
        return;
    }

    for ( i = 0; i < nbhandles; i++ ) {
        if ( *handle != NULL ) {
            starpu_data_unregister_submit( *handle );
            *handle = NULL;
        }
        handle++;
    }

    free( pivot->nextpiv );
    pivot->nextpiv = NULL;
    pivot->prevpiv = NULL;
    (void)sequence;
}

/**
 *  Destroy ws_pivot runtime structures
 */
void RUNTIME_pivot_destroy( CHAM_desc_pivot_t *pivot )
{
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(pivot->nextpiv);

    if ( handle ) {
        int    i;
        size_t nbhandles = 2 * pivot->P;

        for ( i = 0; i < nbhandles; i++ ) {
            if ( *handle != NULL ) {
                starpu_data_unregister( *handle );
                *handle = NULL;
            }
            handle++;
        }

        free( pivot->nextpiv );
        pivot->nextpiv = NULL;
        pivot->prevpiv = NULL;
    }
    chameleon_starpu_tag_release( pivot->mpitag_nextpiv );
}

void *RUNTIME_ipiv_getaddr( const CHAM_ipiv_t *ipiv, int m )
{
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(ipiv->ipiv);
    int64_t mm = m + (ipiv->i / ipiv->mb);

    handle += mm;
    assert( handle );

    if ( *handle != NULL ) {
        return (void*)(*handle);
    }

    int ncols = (mm == (ipiv->mt-1)) ? ipiv->m - mm * ipiv->mb : ipiv->mb;

    starpu_vector_data_register( handle, -1, (uintptr_t)NULL, ncols, sizeof(int) );

#if defined(CHAMELEON_USE_MPI)
    {
        int                owner = ipiv->get_rankof( ipiv, m, m );
        int64_t            tag   = ipiv->mpitag_ipiv + mm;
        starpu_mpi_data_register( *handle, tag, owner );
    }
#endif /* defined(CHAMELEON_USE_MPI) */

    assert( *handle );
    return (void*)(*handle);
}

void *RUNTIME_nextpiv_getaddr( const CHAM_desc_pivot_t *pivot, int rank, int k, int h )
{
    starpu_data_handle_t *nextpiv = (starpu_data_handle_t*)(pivot->nextpiv);
    int                   Q       = pivot->Q;

    nextpiv += rank/Q;
    assert( nextpiv );

    if ( *nextpiv != NULL ) {
        return (void*)(*nextpiv);
    }
    int     owner = rank;
    int     ncols = pivot->nb;
    int64_t tag   = pivot->mpitag_nextpiv + owner/Q;

    cppi_register( nextpiv, pivot->dtyp, ncols, tag, owner );

    assert( *nextpiv );
    (void)h;
    return (void*)(*nextpiv);
}

void *RUNTIME_prevpiv_getaddr( const CHAM_desc_pivot_t *pivot, int rank, int k, int h )
{
    starpu_data_handle_t *prevpiv = (starpu_data_handle_t*)(pivot->prevpiv);
    int                   Q       = pivot->Q;

    prevpiv += rank/Q;
    assert( prevpiv );

    if ( *prevpiv != NULL ) {
        return (void*)(*prevpiv);
    }

    int     owner = rank;
    int     ncols = pivot->nb;
    int64_t tag   = pivot->mpitag_prevpiv + owner/Q;

    cppi_register( prevpiv, pivot->dtyp, ncols, tag, owner );

    assert( *prevpiv );
    (void)h;
    return (void*)(*prevpiv);
}

void *RUNTIME_perm_getaddr( const CHAM_ipiv_t *ipiv, int m )
{
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(ipiv->perm);
    int64_t mm = m + (ipiv->i / ipiv->mb);

    handle += mm;
    assert( handle );

    if ( *handle != NULL ) {
        return (void*)(*handle);
    }

    int ncols = ipiv->mb;

    starpu_vector_data_register( handle, -1, (uintptr_t)NULL, ncols, sizeof(int) );

#if defined(CHAMELEON_USE_MPI)
    {
        int                owner = ipiv->get_rankof( ipiv, m, m );
        int64_t            tag   = ipiv->mpitag_perm + mm;
        starpu_mpi_data_register( *handle, tag, owner );
    }
#endif /* defined(CHAMELEON_USE_MPI) */

    assert( *handle );
    return (void*)(*handle);
}

void *RUNTIME_invp_getaddr( const CHAM_ipiv_t *ipiv, int m )
{
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(ipiv->invp);
    int64_t mm = m + (ipiv->i / ipiv->mb);

    handle += mm;
    assert( handle );

    if ( *handle != NULL ) {
        return (void*)(*handle);
    }

    int ncols = ipiv->mb;

    starpu_vector_data_register( handle, -1, (uintptr_t)NULL, ncols, sizeof(int) );

#if defined(CHAMELEON_USE_MPI)
    {
        int                owner = ipiv->get_rankof( ipiv, m, m );
        int64_t            tag   = ipiv->mpitag_invp + mm;
        starpu_mpi_data_register( *handle, tag, owner );
    }
#endif /* defined(CHAMELEON_USE_MPI) */

    assert( *handle );
    return (void*)(*handle);
}

void RUNTIME_pivot_flushk( const RUNTIME_sequence_t *sequence,
                           const CHAM_desc_pivot_t *pivot, int rank )
{
    starpu_data_handle_t *handle;
    int                   Q = pivot->Q;

    handle = (starpu_data_handle_t*)(pivot->nextpiv);
    handle += rank/Q;

    if ( *handle != NULL ) {
#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_cache_flush( sequence->comm, *handle );
        if ( starpu_mpi_data_get_rank( *handle ) == rank )
#endif
        {
            chameleon_starpu_data_wont_use( *handle );
        }
    }

    handle = (starpu_data_handle_t*)(pivot->prevpiv);
    handle += rank/Q;

    if ( *handle != NULL ) {
#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_cache_flush( sequence->comm, *handle );
        if ( starpu_mpi_data_get_rank( *handle ) == rank )
#endif
        {
            chameleon_starpu_data_wont_use( *handle );
        }
    }

    (void)sequence;
    (void)pivot;
    (void)rank;
}

void RUNTIME_pivot_flush( const RUNTIME_sequence_t *sequence,
                          const CHAM_desc_pivot_t  *pivot )
{
    int m;

    for (m = 0; m < pivot->Q; m++)
    {
        RUNTIME_pivot_flushk( sequence, pivot, m );
    }
}

void RUNTIME_ipiv_flush( const RUNTIME_sequence_t *sequence,
                         const CHAM_ipiv_t        *ipiv )
{
    int m;

    for (m = 0; m < ipiv->mt; m++)
    {
        RUNTIME_ipiv_flushk( sequence, ipiv, m );
    }
}

void RUNTIME_ipiv_flushk( const RUNTIME_sequence_t *sequence,
                          const CHAM_ipiv_t *ipiv, int m )
{
    starpu_data_handle_t *handle;
    int64_t mm = m + ( ipiv->i / ipiv->mb );

    handle = (starpu_data_handle_t*)(ipiv->ipiv);
    handle += mm;

    if ( *handle != NULL ) {
#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_cache_flush( sequence->comm, *handle );
        if ( starpu_mpi_data_get_rank( *handle ) == ipiv->myrank )
#endif
        {
            chameleon_starpu_data_wont_use( *handle );
        }
    }

    (void)sequence;
    (void)ipiv;
    (void)m;
}

void RUNTIME_perm_flushk( const RUNTIME_sequence_t *sequence,
                          const CHAM_ipiv_t *ipiv, int m )
{
    starpu_data_handle_t *handle;
    int64_t mm = m + ( ipiv->i / ipiv->mb );

    handle = (starpu_data_handle_t*)(ipiv->perm);
    handle += mm;

    if ( *handle != NULL ) {
#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_cache_flush( sequence->comm, *handle );
        if ( starpu_mpi_data_get_rank( *handle ) == ipiv->myrank )
#endif
        {
            chameleon_starpu_data_wont_use( *handle );
        }
    }

    handle = (starpu_data_handle_t*)(ipiv->invp);
    handle += mm;

    if ( *handle != NULL ) {
#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_cache_flush( sequence->comm, *handle );
        if ( starpu_mpi_data_get_rank( *handle ) == ipiv->myrank )
#endif
        {
            chameleon_starpu_data_wont_use( *handle );
        }
    }

    (void)sequence;
    (void)ipiv;
    (void)m;
}

void RUNTIME_ipiv_gather( const RUNTIME_sequence_t *sequence,
                          CHAM_ipiv_t *desc, int *ipiv, int node )
{
    int64_t mt   = desc->mt;
    int64_t mb   = desc->mb;
    int64_t tag  = chameleon_starpu_tag_book( (int64_t)(desc->mt) );
    int     rank = CHAMELEON_Comm_rank();
    int     m;

    for (m = 0; m < mt; m++, ipiv += mb) {
        starpu_data_handle_t ipiv_src = RUNTIME_ipiv_getaddr( desc, m );

#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_get_data_on_node( sequence->comm, ipiv_src, node );
        if ( rank == node )
#endif
        {
            starpu_data_handle_t ipiv_dst;
            int       ncols     = (m == (mt-1)) ? desc->m - m * mb : mb;
            uintptr_t ipivptr   = (rank == node) ? (uintptr_t)ipiv : 0;
            int       home_node = (rank == node) ? STARPU_MAIN_RAM : -1;

            starpu_vector_data_register( &ipiv_dst, home_node, ipivptr, ncols, sizeof(int) );

#if defined(CHAMELEON_USE_MPI)
            starpu_mpi_data_register( ipiv_dst, tag + m, 0 );
#endif /* defined(CHAMELEON_USE_MPI) */

            assert( ipiv_dst );

            starpu_data_cpy( ipiv_dst, ipiv_src, 0, NULL, NULL );
            starpu_data_unregister( ipiv_dst );
        }
    }

    chameleon_starpu_tag_release( tag );
}

void RUNTIME_cpui_create( CHAM_perm_t *ws )
{
    size_t                nbhandles = ( ws->side == ChamLeft ) ? ws->nt * ws->NP :
                                                                 ws->mt * ws->NP;
    starpu_data_handle_t *handles = calloc( nbhandles, sizeof(starpu_data_handle_t) );

    ws->ws = handles;

#if defined(CHAMELEON_USE_MPI)
    /*
     * Book the number of tags required to describe workspace structure
     * One per handle type
     */
    {
        chameleon_starpu_tag_init();
        ws->mpitag_ws = chameleon_starpu_tag_book( nbhandles );
        if ( ws->mpitag_ws == -1 ) {
            chameleon_fatal_error("RUNTIME_cpui_create", "Can't pursue computation since no more tags are available for workspace structure");
            return;
        }
    }
#endif
}

void *
RUNTIME_cpui_getaddr( const CHAM_perm_t *ws,
                      int                m,
                      int                n )
{
    starpu_data_handle_t *ptr_ws = (starpu_data_handle_t*)(ws->ws);
    int                   ws_idx = ( ws->side == ChamLeft) ? m + n * ws->NP :
                                                             n + m * ws->NP;
    cham_side_t           side   = ws->side;
    int                   ncols, mrows;

    ptr_ws += ws_idx;
    assert( ptr_ws );

    if ( *ptr_ws != NULL ) {
        return (void*)(*ptr_ws);
    }

    int owner = ( side == ChamLeft ) ? m : n;

    if ( side == ChamLeft ) {
        ncols = ( n == ws->nt - 1 ) ? ws->n - n * ws->nb : ws->nb ;
        mrows = ws->mb;
    }
    else {
        ncols = ( m == ws->mt - 1 ) ? ws->m - m * ws->mb : ws->mb ;
        mrows = ws->nb;
    }

    int64_t tag   = ws->mpitag_ws + ws_idx;

    cpui_register( ptr_ws, side, ws->dtyp, mrows, ncols, tag, owner );

    assert( *ptr_ws );
    return (void*)(*ptr_ws);
}

/**
 *  Destroy workspace runtime structures
 */
void RUNTIME_cpui_destroy( CHAM_perm_t *ws )
{
    size_t                i;
    starpu_data_handle_t *handle    = (starpu_data_handle_t*)(ws->ws);
    size_t                nbhandles = ( ws->side == ChamLeft ) ? ws->nt * ws->NP :
                                                                 ws->mt * ws->NP;

    for( i = 0; i < nbhandles; i++ ) {
        if ( *handle != NULL ) {
            starpu_data_unregister_submit( *handle );
            *handle = NULL;
        }
        handle++;
    }

    free( ws->ws );
    ws->ws = NULL;
    chameleon_starpu_tag_release( ws->mpitag_ws );
}

void RUNTIME_cpui_flushk( const RUNTIME_sequence_t *sequence,
                          int                       rank,
                          const CHAM_perm_t        *ws,
                          int                       m,
                          int                       n )
{
    starpu_data_handle_t *handle;
    int                   ws_idx = ( ws->side == ChamLeft) ? m + n * ws->NP :
                                                             n + m * ws->NP;

    handle = (starpu_data_handle_t*)(ws->ws);
    handle += ws_idx;

    if ( *handle != NULL ) {
#if defined(CHAMELEON_USE_MPI)
        starpu_mpi_cache_flush( sequence->comm, *handle );
        if ( starpu_mpi_data_get_rank( *handle ) == rank )
#endif
        {
            chameleon_starpu_data_wont_use( *handle );
        }
    }
}

