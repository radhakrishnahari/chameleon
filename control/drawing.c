#include "chameleon.h"
#include "rapack.h"

static void path_cell(FILE *fp, float step_x, float step_y, float x, float y) {
    fprintf(fp, "newpath\n");
    fprintf(fp, "%f %f moveto\n", x, y);
    fprintf(fp, "%f %f lineto\n", x+step_x, y);
    fprintf(fp, "%f %f lineto\n", x+step_x, y+step_y);
    fprintf(fp, "%f %f lineto\n", x, y+step_y);
    fprintf(fp, "%f %f lineto\n", x, y);
    fprintf(fp, "closepath\n");
}

static void moveto_cell(FILE *fp, float step_x, float step_y, float x, float y) {
    fprintf(fp, "%f %f moveto\n", x+step_x/2.0, y+step_y/2.0);
}

void chameleon_tile2ps( cham_uplo_t uplo, CHAM_desc_t *descAt, FILE *fp ) {
    int i, j;
    CHAM_tile_t *tile;

    int im_height = 1080;
    int im_width = 1080;
    float margin_x = 40;
    float margin_y = 40;
    float step_x = ((float)im_height - 2 * margin_x) / (float)descAt->lnt;
    float step_y = ((float)im_width - 2 * margin_y) / (float)descAt->lmt;

    fprintf(fp, "%%!PS\n");
    fprintf(fp, "%%%%BoundingBox: 0 0 %d %d\n", im_width, im_height);
    fprintf(fp, "%% [chameleon_descriptor]\n");
    fprintf(fp, "%% lm = %d\n", descAt->lm);
    fprintf(fp, "%% ln = %d\n", descAt->ln);
    fprintf(fp, "%% lmt = %d\n", descAt->lmt);
    fprintf(fp, "%% lnt = %d\n", descAt->lnt);
    fprintf(fp, "%% mb = %d\n", descAt->mb);
    fprintf(fp, "%% nb = %d\n", descAt->nb);

    for( j=0; j < descAt->nt; j++ ) {
        for( i=0; i < descAt->mt; i++ ) {
            float x = margin_x + (float)j * step_x;
            float y = margin_y + (descAt->lmt - 1 - (float)i) * step_y;
            tile = descAt->get_blktile( descAt, i, j );
            fprintf(fp, "%% drawing cell %d %d\n", i, j);

            path_cell(fp, step_x, step_y, x, y);
            switch(tile->format) {
                case CHAMELEON_TILE_FULLRANK:
                    fprintf(fp, "1.0 0.0 0.0 setrgbcolor\n");
                    break;
                case CHAMELEON_TILE_LOWRANK:
                    {
	                    rpk_matrix_t *Tra = (rpk_matrix_t*) tile->mat;
                        if (!Tra) {
                            continue;
                            fprintf(fp, "0.1 0.1 0.1 setrgbcolor\n");
                        }
                        else if (Tra->rk == -1) {
                            fprintf(fp, "1.0 0.0 0.0 setrgbcolor\n");
                        }
                        else if (Tra->rk == 0) {
                            fprintf(fp, "0.8 setgray\n");
                        }
                        else {
                            fprintf(fp, "0.0 1.0 0.0 setrgbcolor\n");
                        }
                    }
                    break;
            }
            fprintf(fp, "fill\n");

            if (tile->format == CHAMELEON_TILE_LOWRANK) {
                rpk_matrix_t *Tra = (rpk_matrix_t*) tile->mat;
                if ( Tra && (Tra->rk > 0) ) {
                    fprintf(fp, "/Times-Bold findfont 10 scalefont setfont\n");
                    moveto_cell(fp, step_x, step_y, x, y);
                    fprintf(fp, "0.0 setgray\n");
                    fprintf(fp, "(%d) show\n", Tra->rk);
                }
            }
        }
    }

	tile = descAt->tiles;
    for(j=0; j < descAt->lnt; j++) {
        for(i=0; i < descAt->lmt; i++, tile++) {
            float x = margin_x + (float)j * step_x;
            float y = margin_y + (descAt->lmt - 1 - (float)i) * step_y;
            path_cell(fp, step_x, step_y, x, y);
            fprintf(fp, "0.0 0.0 0.0 setrgbcolor\n");
            fprintf(fp, "stroke\n");
        }
    }

    fprintf(fp, "showpage\n");
}
