/**
 *
 * @file zge2ps.c
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 ***
 *
 * @brief Chameleon auxiliary routines for testing structures
 *
 * @version 1.3.0
 * @author Abel Calluaud
 * @date 2023-03-29
 * @precisions normal z -> c d s
 *
 */
#include "chameleon.h"
#include "coreblas.h"
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include "control/common.h"

typedef struct zge2ps_args_s {
    double font_size;
    double offset_m;
    double offset_n;
    int bbox_width;
    int bbox_height;
    int rkmax;
    double m_step;
    double n_step;
    char **fill_tile_ps;
    char **stroke_tile_ps;
} zge2ps_args_t;

static int ge2ps_id = 0;

/* Build the rank matrix */
int zge2ps_build_rkmat_op_cpu( void *args, cham_uplo_t uplo, int m, int n, int ndata,
                   const CHAM_desc_t *desc, CHAM_tile_t *tile, ... ) {
    CHAM_tile_t *tileR;
    CHAM_desc_t *descRk;
    CHAM_tile_t *tileSz;
    CHAM_desc_t *descSz;
    va_list vargs;
    int rank;
    double *rank_mat;
    double *size_mat;

    va_start(vargs, tile);
    descRk = va_arg(vargs, CHAM_desc_t*);
    (void)descRk;
    tileR = va_arg(vargs, CHAM_tile_t*);
    
    descSz = va_arg(vargs, CHAM_desc_t*);
    (void)descSz;
    tileSz = va_arg(vargs, CHAM_tile_t*);
    va_end(vargs);

    rank = TCORE_zlrcrk( tile );
    rank_mat = tileR->mat;
    size_mat = tileSz->mat;
    *rank_mat = (double)rank;

    if (rank == -1) {
        *size_mat = tile->m * tile->n;
    }
    else {
        *size_mat = ((double)rank) * (tile->m + tile->n);
    }

    return CHAMELEON_SUCCESS;
}

/* Format the rank matrix to postscript */
int zge2ps_draw_op_cpu( void *args, cham_uplo_t uplo, int m, int n, int ndata,
                   const CHAM_desc_t *descRk, CHAM_tile_t *tileR, ... ) {
    CHAM_tile_t *tilePs;
    va_list vargs;
    char *stroke_tile_script = NULL;
    char *fill_tile_script = NULL;
    size_t buffer_size = 4096;
    int rc;
    double *rank_mat = (double*)CHAM_tile_get_ptr( tileR );
    zge2ps_args_t *zge2ps_args = (zge2ps_args_t*)args;
    double rank = *rank_mat;
    double norm_rank = rank / zge2ps_args->rkmax;
    double moff, noff;
    char color[256];

    if (norm_rank > 1.0) {
        norm_rank = 1.0;
    }

    if ( rank == -1 ) {
        snprintf(color, 256, "1.0 0.2 0.2");
    }
    else {
        double r_start = 0.2;
        double g_start = 0.2;
        double b_start = 1.0;
        double r_end = 0.2;
        double g_end = 1.0;
        double b_end = 0.2;
        double r_val = r_start + (r_end - r_start) * norm_rank;
        double g_val = g_start + (g_end - g_start) * norm_rank;
        double b_val = b_start + (b_end - b_start) * norm_rank;
        snprintf(color, 256, "%.2f %.2f %.2f", r_val, g_val, b_val);
    }

    moff = m * zge2ps_args->m_step,
    noff = zge2ps_args->offset_m + (descRk->nt - n - 1) * zge2ps_args->n_step;

    stroke_tile_script = malloc(buffer_size);
    
    if ( rank >= 0.0 ) {
        rc = snprintf(
            stroke_tile_script,
            buffer_size,
            "0.0 0.0 0.0 setrgbcolor\n%d %d moveto (%d) show\n",
            (int)moff + 4, (int)noff + 4,
            (int)rank
        );
        if (rc < 0) {
            fprintf(stderr, "Error while formatting the rank matrix\n");
            return CHAMELEON_ERR_UNEXPECTED;
        }
    }
    else {
        rc = 0;
    }
    rc = snprintf(
        stroke_tile_script + rc,
        buffer_size,
        "newpath\n%d %d moveto\n%d %d lineto\n%d %d lineto\n%d %d lineto\nclosepath\n0.0 setgray\nstroke\n",
        (int)moff, (int)noff,
        (int)(moff + zge2ps_args->m_step), (int)noff,
        (int)(moff + zge2ps_args->m_step), (int)(noff + zge2ps_args->n_step),
        (int)moff, (int)(noff + zge2ps_args->n_step)
    );
    if (rc < 0) {
        fprintf(stderr, "Error while formatting the rank matrix\n");
        return CHAMELEON_ERR_UNEXPECTED;
    }
    zge2ps_args->stroke_tile_ps[m * descRk->nt + n] = stroke_tile_script;
    
    fill_tile_script = malloc(buffer_size); 
    rc = snprintf(
        fill_tile_script,
        buffer_size,
        "newpath\n%d %d moveto\n%d %d lineto\n%d %d lineto\n%d %d lineto\nclosepath\n%s setrgbcolor\nfill\n",
        (int)moff, (int)noff,
        (int)(moff + zge2ps_args->m_step), (int)noff,
        (int)(moff + zge2ps_args->m_step), (int)(noff + zge2ps_args->n_step),
        (int)moff, (int)(noff + zge2ps_args->n_step),
        color
    );
    if (rc < 0) {
        fprintf(stderr, "Error while formatting the rank matrix\n");
        return CHAMELEON_ERR_UNEXPECTED;
    }
    zge2ps_args->fill_tile_ps[m * descRk->nt + n] = fill_tile_script;

    return CHAMELEON_SUCCESS;
}

int
CHAMELEON_zge2ps_Tile( FILE *psfile, cham_uplo_t uplo, CHAM_desc_t *descAt ) {
    CHAM_context_t *chamctxt;
    int status = CHAMELEON_SUCCESS;
    CHAM_desc_t *descRk;
    CHAM_desc_t *descSz;
    cham_map_operator_t zge2ps_op = {0};
    cham_map_data_t zge2ps_data[3];
    zge2ps_args_t zge2ps_args;
    int mprank = CHAMELEON_Comm_rank();
    int mpsize = CHAMELEON_Comm_size();
    FILE* psfile_mp;
    char *psfile_mp_path;
    int psfile_mp_path_len;
    int i;
    int j;
    double rkmax;
    chamctxt = chameleon_context_self();
    ge2ps_id++;

    zge2ps_args.offset_n = 0;
    zge2ps_args.bbox_width = 720;
    zge2ps_args.bbox_height = ((double)descAt->lnt / (double)descAt->lmt) * zge2ps_args.bbox_width;
    zge2ps_args.stroke_tile_ps = malloc( descAt->mt * descAt->nt * sizeof(char*) );
    zge2ps_args.fill_tile_ps = malloc( descAt->mt * descAt->nt * sizeof(char*) );
    zge2ps_args.m_step = (double)zge2ps_args.bbox_width / (double)descAt->mt;
    zge2ps_args.n_step = (double)zge2ps_args.bbox_height / (double)descAt->nt;
    zge2ps_args.rkmax = (double)descAt->nb / 4.0;
    zge2ps_args.font_size = (double)((zge2ps_args.n_step >= zge2ps_args.m_step)? zge2ps_args.n_step : zge2ps_args.m_step) / 3.0;
    zge2ps_args.offset_m = 1.5 * 8.0 * zge2ps_args.font_size;

    for(i = 0; i < descAt->mt; i++) {
        for(j = 0; j < descAt->nt; j++) {
            zge2ps_args.fill_tile_ps[i * descAt->nt + j] = NULL;
            zge2ps_args.stroke_tile_ps[i * descAt->nt + j] = NULL;
        }
    }

    CHAMELEON_Desc_Create(
        &descRk,
        CHAMELEON_MAT_ALLOC_GLOBAL,
        ChamRealDouble,
        1, 1, /* tile size */
        1, /* number of elements in each tile */
        descAt->lmt, descAt->lnt, /* size of the full matrix */
        0, 0, /* indexes of the submatrix */
        descAt->lmt, descAt->lnt, /* sizes of the submatrix */
        1, 1 /* only master process own the rank matrix */
    );

    CHAMELEON_Desc_Create(
        &descSz,
        CHAMELEON_MAT_ALLOC_GLOBAL,
        ChamRealDouble,
        1, 1, /* tile size */
        1, /* number of elements in each tile */
        descAt->lmt, descAt->lnt, /* size of the full matrix */
        0, 0, /* indexes of the submatrix */
        descAt->lmt, descAt->lnt, /* sizes of the submatrix */
        1, 1 /* only master process own the rank matrix */
    );

    zge2ps_op.name = "zge2ps_build_rk_mat";
    zge2ps_op.cpufunc = zge2ps_build_rkmat_op_cpu;
    zge2ps_data[0].access = ChamR;
    zge2ps_data[0].desc = descAt;
    zge2ps_data[1].access = ChamW;
    zge2ps_data[1].desc = descRk;
    zge2ps_data[2].access = ChamW;
    zge2ps_data[2].desc = descSz;
    CHAMELEON_mapv_Tile( uplo, 3, zge2ps_data, &zge2ps_op, NULL);

    rkmax = CHAMELEON_zlange_Tile( ChamMaxNorm, descRk );
    zge2ps_args.rkmax = 2 * rkmax;

    zge2ps_op.name = "zge2ps_stroke";
    zge2ps_op.cpufunc = zge2ps_draw_op_cpu;
    zge2ps_data[0].access = ChamR;
    zge2ps_data[0].desc = descRk;
    CHAMELEON_mapv_Tile( uplo, 1, zge2ps_data, &zge2ps_op, (void*)&zge2ps_args );

    psfile_mp_path_len = 1024;
    psfile_mp_path = malloc( psfile_mp_path_len * sizeof(char) );
    snprintf(psfile_mp_path, psfile_mp_path_len, "ge2ps%d.%d", ge2ps_id, mprank);
    psfile_mp = fopen(psfile_mp_path, "w");

    {
        for (i = 0; i < descRk->mt; i++) {
            for (int j = 0; j < descRk->nt; j++) {
                char* fill_tile_script = zge2ps_args.fill_tile_ps[i * descRk->nt + j];
                if (fill_tile_script) {
                    fprintf(psfile_mp, "\n%%tile %d,%d\n%s", i, j, fill_tile_script);
                }
            }
        }
        for (i = 0; i < descRk->mt; i++) {
            for (int j = 0; j < descRk->nt; j++) {
                char* stroke_tile_script = zge2ps_args.stroke_tile_ps[i * descRk->nt + j];
                if (stroke_tile_script) {
                    fprintf(psfile_mp, "\n%%tile %d,%d\n%s", i, j, stroke_tile_script);
                }
            }
        }
    }
    free(psfile_mp_path);
    fclose(psfile_mp);

    RUNTIME_barrier(chamctxt);

    if ( mprank == 0 ) {
        double comp_size = 0.0;
        double full_size = (double)descAt->lm * descAt->ln;
        double comp_rate;

        fprintf(psfile, "%%!PS\n");
        fprintf(psfile, "%%%%BoundingBox 0 0 %d %d\n", zge2ps_args.bbox_width, (int)(zge2ps_args.bbox_height + zge2ps_args.offset_m));
        fprintf(psfile, "/Courier\n");
        fprintf(psfile, "%.2f selectfont\n", zge2ps_args.font_size);
        fprintf(psfile, "2.0 setlinewidth\n");

        for (i = 0; i < mpsize; i++ ) {
            psfile_mp_path = malloc( psfile_mp_path_len * sizeof(char) );
            snprintf(psfile_mp_path, psfile_mp_path_len, "ge2ps%d.%d", ge2ps_id, i);
            psfile_mp = fopen(psfile_mp_path, "r");
            if (psfile_mp) {
                char c;
                while ((c = fgetc(psfile_mp)) != EOF) {
                    fputc(c, psfile);
                }
                fclose(psfile_mp);
            }
            free(psfile_mp_path);
        }
        
        /*
         * Compute compression rate
         */
        for( i = 0; i < descAt->mt; i++ ) {
            for( j = 0; j < descAt->nt; j++ ) {
                comp_size += ((double*)descSz->mat)[j * descSz->mt + i];
            }
        }
        comp_rate = ( full_size - comp_size ) / ( full_size );

        /* Print matrix informations */
        fprintf(psfile, "0 0 0 setrgbcolor\n");
        
        fprintf(psfile, "%.2f %.2f moveto\n", 12.0, 12.0);
        fprintf(psfile, "(Dimensions: %d x %d) show\n", descAt->lm, descAt->ln);

        fprintf(psfile, "%.2f %.2f moveto\n", 12.0, 12.0 + 1.5 * zge2ps_args.font_size);
        fprintf(psfile, "(Tile size: %d x %d) show\n", descAt->mb, descAt->nb);
       
        fprintf(psfile, "%.2f %.2f moveto\n", 12.0, 12.0 + 2 * 1.5 * zge2ps_args.font_size);
        fprintf(psfile, "(Compression rate: %.2f%%) show\n", comp_rate * 100.0);

        fprintf(psfile, "%.2f %.2f moveto\n", 12.0, 12.0 + 3 * 1.5 * zge2ps_args.font_size);
        fprintf(psfile, "(Rank max: %.2f) show\n", rkmax);

        fprintf(psfile, "showpage\n");
    }

    free(zge2ps_args.fill_tile_ps);
    free(zge2ps_args.stroke_tile_ps);
    
    RUNTIME_barrier(chamctxt);

    chameleon_desc_destroy( descRk );
    chameleon_desc_destroy( descSz );

    return status;
}
