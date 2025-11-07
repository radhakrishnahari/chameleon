#ifndef DEBUG_H
#define DEBUG_H

#include <starpu.h>

static inline
void print_descriptor( CHAM_desc_t *desc, char *name, int lvl ) {
    int i;
    starpu_data_handle_t *handles = (starpu_data_handle_t*)(desc->schedopt);
    if (name) fprintf(stderr,"%s (%d,%d) (%d,%d) (%p) :\n",
                      name, desc->lm, desc->ln, desc->mb, desc->nb, desc);
    for ( i=0; i<desc->lmt*desc->lnt; i++ ) {
        CHAM_tile_t t = desc->tiles[i];
        if ( t.format == CHAMELEON_TILE_DESC ) {
            CHAM_desc_t *subdesc = (CHAM_desc_t*)(t.mat);
#if defined(CHAMELEON_KERNELS_TRACE)
            (void)handles;
            fprintf(stderr, "[LVL%d] desc-%d (%d,%d) (%d,%d) %s %p\n",
                    lvl, i, subdesc->lm, subdesc->ln, subdesc->mb, subdesc->nb, t.name, handles[i]);
#else
            fprintf(stderr, "[LVL%d] desc-%d (%d,%d) (%d,%d) | H=%p | mat=%p\n",
                    lvl, i, subdesc->lm, subdesc->ln, subdesc->mb, subdesc->nb, subdesc->mat, handles[i]);
#endif
            print_descriptor( subdesc, NULL, lvl+1 );
        }
        else if ( t.format == CHAMELEON_TILE_FULLRANK ) {
#if defined(CHAMELEON_KERNELS_TRACE)
            fprintf(stderr, "[LVL%d] tile-%d %s %p\n", lvl, i, t.name, handles[i]);
#else
            fprintf(stderr, "[LVL%d] tile-%d %p %p\n", lvl, i, t.mat, handles[i]);
#endif
        }
        else {
            fprintf(stderr, "??\n");
        }
    }
    if (!lvl) fprintf(stderr, "\n");
}

#endif  /* DEBUG_H */
