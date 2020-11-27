#include "chameleon_starpu.h"

void starpu_CHAM_tile_filter_square_block(void *father_interface, void *child_interface,
                                          STARPU_ATTRIBUTE_UNUSED struct starpu_data_filter *f,
                                          unsigned id, unsigned nchunks)
{
    starpu_cham_tile_interface_t *father = (starpu_cham_tile_interface_t *)father_interface;
    starpu_cham_tile_interface_t *child = (starpu_cham_tile_interface_t *)child_interface;
    CHAM_tile_t *father_tile = &(father->tile);

    assert( father_tile->format & CHAMELEON_TILE_DESC );

    //CHAM_desc_t *desc = (CHAM_desc_t*)(father_tile->mat);
    CHAM_desc_t *desc = (CHAM_desc_t*)(f->filter_arg_ptr);

    size_t elemsize = CHAMELEON_Element_Size(father->flttype);

    STARPU_ASSERT_MSG( nchunks <= (unsigned)(desc->lmt * desc->lnt),
                       "cannot split %u elements in %u parts",
                       (unsigned)(desc->lmt * desc->lnt), nchunks );

    int    child_m, child_n;

    child->id         = father->id;
    child->dev_handle = 0;
    child->flttype    = father->flttype;
    child->allocsize  = -1;
    child->tilesize   = 0;

    child->tile = desc->tiles[id];
    child_m = child->tile.m;
    child_n = child->tile.n;
    assert( child_m <= desc->mb );
    assert( child_n <= desc->nb );

    /* } */
    /* else */
    /* { */
    /*     /\* Should be col major for chameleon (row major rn) *\/ */
    /*     /\* Utiliser le champ de starpu_data_filter, filter_arg pour avoir le côté du bloc (nchunks doit être ajusté en initialisant le filtre) *\/ */
    /*     unsigned nchunks_col_line = (unsigned)(sqrt((double)nchunks)); */
    /*     starpu_filter_nparts_compute_chunk_size_and_offset(n, nchunks_col_line, elemsize, id%nchunks_col_line, father_tile.ld, &child_n, &offset); */
    /*     size_t mem_offset = offset; */
    /*     starpu_filter_nparts_compute_chunk_size_and_offset(m, nchunks_col_line, elemsize, id/nchunks_col_line, 1, &child_m, &offset); */
    /*     offset += mem_offset; */

    /*     child->tile.m = child_m; */
    /*     child->tile.n = child_n; */
    /* } */
    child->tilesize  = child_m * child_n * elemsize;
    child->allocsize = child_m * child_n * elemsize;

    /* STARPU_ASSERT_MSG(matrix_father->allocsize == matrix_father->m * matrix_father->n * matrix_father->elemsize, */
    /*                   "partitioning matrix with non-trivial allocsize not supported yet, patch welcome"); */

    if (father->dev_handle)
    {
        //child->tile.mat   = father_tile.mat + offset;
        //child->tile.ld    = father_tile.ld;
        child->dev_handle = (intptr_t)(child->tile.mat);
    }
}
