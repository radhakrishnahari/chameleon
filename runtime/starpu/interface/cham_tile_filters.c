#include "chameleon_starpu.h"

void starpu_CHAM_tile_filter_square_block( void *father_interface, void *child_interface,
                                           struct starpu_data_filter *f,
                                           unsigned id, unsigned nchunks)
{
    starpu_cham_tile_interface_t *father = (starpu_cham_tile_interface_t *)father_interface;
    starpu_cham_tile_interface_t *child = (starpu_cham_tile_interface_t *)child_interface;
    CHAM_tile_t *father_tile = &(father->tile);
    CHAM_desc_t *desc = (CHAM_desc_t*)(f->filter_arg_ptr);

    size_t elemsize = CHAMELEON_Element_Size(father->flttype);

    STARPU_ASSERT_MSG( nchunks <= (unsigned)(desc->lmt * desc->lnt),
                       "cannot split %u elements in %u parts",
                       (unsigned)(desc->lmt * desc->lnt), nchunks );

    int    child_m, child_n;

    child->id         = father->id;
    child->flttype    = father->flttype;
    child->dev_handle = 0;
    child->allocsize  = -1;
    child->tilesize   = 0;

    /*
     * Copy the tile information from the desc that is always defined on the CPU
     * side.  Note that the father_tile may not be a CHAMELEON_TILE_DESC after
     * copy on a remote device.
     */
    child->tile = desc->tiles[id];

    child_m = child->tile.m;
    child_n = child->tile.n;
    assert( child_m <= desc->mb );
    assert( child_n <= desc->nb );
    child->tilesize  = child_m * child_n * elemsize;
    child->allocsize = child_m * child_n * elemsize;

    if (father->dev_handle)
    {
        child->dev_handle = father->dev_handle;
    }

    if ( father_tile->format & CHAMELEON_TILE_DESC ) {
        assert( CHAM_tile_get_ptr( &(child->tile) ) == desc->get_blkaddr( desc, id % desc->lmt, id / desc->lmt ) );
    }
    else {
        assert( father_tile->format & CHAMELEON_TILE_FULLRANK );
        size_t offset = chameleon_getaddr_cm_offset( desc, id % desc->lmt, id / desc->lmt, father_tile->ld );
        child->tile.mat = ((char *)(father_tile->mat)) + offset * elemsize;
        child->tile.ld  = father_tile->ld;
    }
}
