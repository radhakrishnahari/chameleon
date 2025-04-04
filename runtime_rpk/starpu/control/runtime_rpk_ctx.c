/**
 *
 * @file runtime_rpk_ctx.c
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon rapack context configuration
 *
 * @version 1.2.0
 * @author Abel Calluaud
 * @date 2024-03-27
 *
 */

#include "chameleon/struct.h"
#include "chameleon/struct_context.h"
#include "../control/context.h"
#include "runtime_rpk.h"
#include "runtime_lrmem.h"
#include <starpu.h>
#include "cham_tile_interface.h"
#include <rapack.h>

/* Array of RAPACK contextes for each worker/precision combination */
static rpk_ctx_t *cham_rpk_ctx;

/* Array of CHAMALEON-StarPU-RAPACK contextes for each worker/precision combination */
static cham_starpu_rpk_ctx_t *cham_starpu_rpk_ctx;

/*
 * Rapack locking function for synchronization on C when doing a gemm
 */
static void
cham_rpk_lock(const rpk_ctx_t *ctx, void* arg) {
    (void)ctx;
    (void)arg;
}

/*
 * Rapack function to set rank limit for using low rank representation
 */
static rpk_int_t
cham_rpkx_getrklimit( const rpk_ctx_t *ctx, rpk_int_t m, rpk_int_t n) {
    return rpk_imax(1, ( ctx->minratio * chameleon_min(m,n) ) / 4);
}

void
runtime_rpk_ctx_init( void ) {
	int nworkers = starpu_cpu_worker_get_count();
	int wid;
	int flttype;
	cham_rpk_ctx = malloc( sizeof( rpk_ctx_t ) * nworkers * ChamComplexDoubleMixed );
	cham_starpu_rpk_ctx = malloc( sizeof( cham_starpu_rpk_ctx_t ) * nworkers * ChamComplexDoubleMixed );

	for (flttype = 0; flttype < (int)ChamComplexDoubleMixed; flttype++) {
		for (wid = 0; wid < nworkers; wid++) {
			cham_rpk_ctx[flttype * nworkers + wid].minratio = 1.0;
			cham_rpk_ctx[flttype * nworkers + wid].orthmeth = 0;
			cham_rpk_ctx[flttype * nworkers + wid].use_reltol = 0;
			cham_rpk_ctx[flttype * nworkers + wid].tolerance = 0;
			cham_rpk_ctx[flttype * nworkers + wid].get_rklimit = cham_rpkx_getrklimit;
			cham_rpk_ctx[flttype * nworkers + wid].rpk_rradd = NULL;
			cham_rpk_ctx[flttype * nworkers + wid].rpk_ge2lr = NULL;
			cham_rpk_ctx[flttype * nworkers + wid].lock = cham_rpk_lock;
			cham_rpk_ctx[flttype * nworkers + wid].unlock = cham_rpk_lock;
           
            /* 
            cham_rpk_ctx[flttype * nworkers + wid].malloc = runtime_rpk_malloc;
			cham_rpk_ctx[flttype * nworkers + wid].free = runtime_rpk_free;
            cham_rpk_ctx[flttype * nworkers + wid].lrmem_alloc = runtime_rpk_lrmem_alloc;
			cham_rpk_ctx[flttype * nworkers + wid].lrmem_free = runtime_rpk_lrmem_free;
			cham_rpk_ctx[flttype * nworkers + wid].lrmem_set = runtime_rpk_lrmem_set;
			cham_rpk_ctx[flttype * nworkers + wid].hook_data = &cham_starpu_rpk_ctx[flttype * nworkers + wid];
			cham_starpu_rpk_ctx[flttype * nworkers + wid].node = -1;
			cham_starpu_rpk_ctx[flttype * nworkers + wid].handle = 0;
            */

			/* TODO: remove this when the runtime functions are fixed */
			cham_rpk_ctx[flttype * nworkers + wid].malloc = rpkx_malloc_default;
			cham_rpk_ctx[flttype * nworkers + wid].free = rpkx_free_default;
			cham_rpk_ctx[flttype * nworkers + wid].lrmem_alloc = rpkx_lrmem_alloc_default;
			cham_rpk_ctx[flttype * nworkers + wid].lrmem_free = rpkx_lrmem_free_default;
			cham_rpk_ctx[flttype * nworkers + wid].lrmem_set = rpkx_lrmem_set_default;
		}
	}

	for (wid = 0; wid < nworkers; wid++) {
		cham_rpk_ctx[ChamRealFloat * nworkers + wid].rpk_rradd = NULL;
		cham_rpk_ctx[ChamRealFloat * nworkers + wid].rpk_ge2lr = NULL;
	
		cham_rpk_ctx[ChamRealDouble * nworkers + wid].rpk_rradd = NULL;
		cham_rpk_ctx[ChamRealDouble * nworkers + wid].rpk_ge2lr = NULL;

		cham_rpk_ctx[ChamComplexFloat * nworkers + wid].rpk_rradd = NULL;
		cham_rpk_ctx[ChamComplexFloat * nworkers + wid].rpk_ge2lr = NULL;
		
		cham_rpk_ctx[ChamComplexDouble * nworkers + wid].rpk_rradd = NULL;
		cham_rpk_ctx[ChamComplexDouble * nworkers + wid].rpk_ge2lr = NULL;
	}
}

void
runtime_rpk_ctx_exit( void ) {
	free(cham_starpu_rpk_ctx);
	free(cham_rpk_ctx);
}

const rpk_ctx_t *
runtime_rpk_ctx_get( cham_flttype_t type ) {
    int wid = 0;
	int nworkers = starpu_cpu_worker_get_count();
	rpk_ctx_t *ctx = &cham_rpk_ctx[type * nworkers + wid];
    
    rpk_coeftype_t rpk_coeftype;
    switch (type) {
        case ChamRealFloat:
            rpk_coeftype = RapackFloat;
            break;
        case ChamRealDouble:
            rpk_coeftype = RapackDouble;
            break;
        case ChamComplexFloat:
            rpk_coeftype = RapackComplex32;
            break;
        case ChamComplexDouble:
            rpk_coeftype = RapackComplex64;
            break; 
        default:
            assert(0);
            exit(1);
    }

    CHAM_context_t *chamctxt = chameleon_context_self();
    rpk_compmeth_t rpk_compmeth;
    switch ( chamctxt->lrmeth ) {
        case ChamLRMethodSVD:
	        rpk_compmeth = RapackCompressMethodSVD; 
            break;
        case ChamLRMethodPQRCP:
	        rpk_compmeth = RapackCompressMethodPQRCP; 
	        break;
        case ChamLRMethodRQRCP:
	        rpk_compmeth = RapackCompressMethodRQRCP;
            break;
        case ChamLRMethodTQRCP:
	        rpk_compmeth = RapackCompressMethodTQRCP;
            break;
        case ChamLRMethodRQRRT:
	        rpk_compmeth = RapackCompressMethodRQRRT;
            break;
        default:
            fprintf(stderr, "Unknown low rank compression method %d\n", (int)chamctxt->lrmeth);
	        exit(2);
            break;
    }

    ctx->tolerance = chamctxt->accuracy;
    ctx->rpk_ge2lr = rpk_ge2lr_functions[rpk_compmeth][rpk_coeftype - 2];
    ctx->rpk_rradd = rpk_rradd_functions[rpk_compmeth][rpk_coeftype - 2];

	return ctx;
}


#if 0
const rpk_ctx_t *
runtime_rpk_ctx_get( cham_flttype_t type, int wid ) {
	struct starpu_task* task = starpu_task_get_current();
	int hid;
	starpu_data_handle_t rpk_handle = 0;
	
	if (task) {
		for (hid = 0; hid < task->nbuffers; hid++) {
			starpu_data_handle_t handle = STARPU_TASK_GET_HANDLE( task, hid );
			struct starpu_data_interface_ops *ops = starpu_data_get_interface_ops( handle );
			
			if ( ops->interfaceid != STARPU_CHAM_TILE_INTERFACE_ID ) {
				continue;
			}

			if ( task->modes[hid] & STARPU_W ) {
				if ( rpk_handle ) {
					fprintf(stderr, "error: Task accessing multiple low rank tiles in write mode is not supported yet\n");
					exit(1);
				}
				rpk_handle = handle;
				break;
			} 
		}
	}

	runtime_rpk_handle_set( type, wid, rpk_handle );

	int nworkers = starpu_cpu_worker_get_count();
	rpk_ctx_t *ctx = &cham_rpk_ctx[type * nworkers + wid];
    
    /*
     * TODO: this could be done only once per chameleon_init config
     */

    rpk_coeftype_t rpk_coeftype;
    switch (type) {
        case ChamRealFloat:
            rpk_coeftype = RapackFloat;
            break;
        case ChamRealDouble:
            rpk_coeftype = RapackDouble;
            break;
        case ChamComplexFloat:
            rpk_coeftype = RapackComplex32;
            break;
        case ChamComplexDouble:
            rpk_coeftype = RapackComplex64;
            break; 
        default:
            assert(0);
            exit(1);
    }

    CHAM_context_t *chamctxt = chameleon_context_self();
    rpk_compmeth_t rpk_compmeth;
    switch ( chamctxt->lrmeth ) {
        case ChamLRMethodSVD:
	        rpk_compmeth = RapackCompressMethodSVD; 
            break;
        case ChamLRMethodPQRCP:
	        rpk_compmeth = RapackCompressMethodPQRCP; 
	        break;
        case ChamLRMethodRQRCP:
	        rpk_compmeth = RapackCompressMethodRQRCP;
            break;
        case ChamLRMethodTQRCP:
	        rpk_compmeth = RapackCompressMethodTQRCP;
            break;
        case ChamLRMethodRQRRT:
	        rpk_compmeth = RapackCompressMethodRQRRT;
            break;
        default:
            fprintf(stderr, "Unknown low rank compression method %d\n", (int)chamctxt->lrmeth);
	        exit(2);
            break;
    }

    ctx->tolerance = chamctxt->accuracy;
    ctx->rpk_ge2lr = rpk_ge2lr_functions[rpk_compmeth][rpk_coeftype - 2];
    ctx->rpk_rradd = rpk_rradd_functions[rpk_compmeth][rpk_coeftype - 2];

    return ctx;
}

void
runtime_rpk_handle_set( cham_flttype_t flttype, int wid, starpu_data_handle_t handle ) {
	int nworkers = starpu_cpu_worker_get_count();
	cham_starpu_rpk_ctx_t *starpu_rpk = &cham_starpu_rpk_ctx[flttype * nworkers + wid];
	starpu_rpk->handle = handle;
}

void
runtime_rpk_node_set( cham_flttype_t flttype, int wid, unsigned node ) {
	int nworkers = starpu_cpu_worker_get_count();
	cham_starpu_rpk_ctx_t *starpu_rpk = &cham_starpu_rpk_ctx[flttype * nworkers + wid];
	starpu_rpk->node = node;
}
#endif
