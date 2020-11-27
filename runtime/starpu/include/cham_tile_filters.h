#ifndef CHAM_TILE_FILTER_H
#define CHAM_TILE_FILTER_H

#include <starpu.h>
#include "cham_tile_interface.h"

void starpu_CHAM_tile_filter_square_block(void *father_interface, void *child_interface,
                                          STARPU_ATTRIBUTE_UNUSED struct starpu_data_filter *f,
                                          unsigned id, unsigned nchunks);

#endif
