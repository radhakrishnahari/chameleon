/**
 *
 * @file starpu/coreblas_runtime/rpk_lrmem.c
 *
 * @copyright 2024-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU RAPACK interoperability runtime functions
 *
 * @version 1.3.0
 * @author Abel Calluaud
 * @date 2024-04-18
 *
 */
#include <starpu.h>
#include "runtime_lrmem.h"
#include "cham_tile_interface.h"
#include "runtime_rpk.h"
#include "rapack.h"

void *
runtime_rpk_malloc( const rpk_ctx_t *ctx, size_t size )
{
#if defined(CHAMELEON_SIMULATION) && !defined(STARPU_MALLOC_SIMULATION_FOLDED) && !defined(CHAMELEON_USE_MPI)
    return (void*) 1;
#else
    void *ptr;

	assert(size > 0);

	cham_starpu_rpk_ctx_t *starpu_rpk = (cham_starpu_rpk_ctx_t*)ctx->hook_data;

	ptr = (void*)starpu_malloc_on_node_flags(
		starpu_rpk->node, size, STARPU_MALLOC_PINNED | STARPU_MALLOC_COUNT
	);
    return ptr;
#endif
}

void
runtime_rpk_free( const rpk_ctx_t *ctx, void *ptr, size_t size )
{
#if defined(CHAMELEON_SIMULATION) && !defined(STARPU_MALLOC_SIMULATION_FOLDED) && !defined(CHAMELEON_USE_MPI)
    (void)ptr; (void)size;
    return;
#else
	assert(size > 0);
	cham_starpu_rpk_ctx_t *starpu_rpk = (cham_starpu_rpk_ctx_t*)ctx->hook_data;
    
	starpu_free_on_node_flags(
		starpu_rpk->node, (uintptr_t)ptr, size, STARPU_MALLOC_PINNED | STARPU_MALLOC_COUNT
	);
#endif
}


#if 0
/*
 * Allocate the memory buffers for low rank matrix representation in factorized
 * form A = U * V
 */
rpk_int_t
runtime_rpk_lrmem_alloc(
    const rpk_ctx_t *ctx,
    void           **new_u,
    void           **new_v,
    size_t          u_size,
    size_t          v_size,
    size_t          elem_size
) {
    struct starpu_task* task = starpu_task_get_current();
    
    /* Get the data handle id associated to the rapack matrix */
	cham_starpu_rpk_ctx_t *starpu_rpk = (cham_starpu_rpk_ctx_t*)ctx->hook_data;

    if ( !starpu_rpk->handle ) {
        /*
         * No data handle associated to this matrix.
         * Let's use the runtime allocator.
         */
        if ( u_size > 0 ) {
            *new_u = runtime_rpk_malloc( ctx, u_size * elem_size );
        }
        else {
            *new_u = NULL;
        }

        if ( v_size > 0 ) {
            *new_v = runtime_rpk_malloc( ctx, v_size * elem_size );
        }
        else {
            *new_v = NULL;
        }
        return 0;
    }
    
	assert( task );

    /* Get the chameleon tile interface associated to the data handle */
    starpu_cham_tile_interface_t *tile_interface = cti_handle_get_interface( starpu_rpk->handle );/

    /* Get the current worker id */
    unsigned workerid = starpu_worker_get_id_check();

    /* Get the memory node of the current worker */
	unsigned dst_node = starpu_worker_get_memory_node( workerid );

    /*
     * TODO: allocate U and V in a single buffer
     */

    /* Allocate the new memory buffers */
    if ( u_size > 0 ) {
        *new_u = (void*)starpu_malloc_on_node_flags(
            dst_node,
            u_size * elem_size,
            STARPU_MALLOC_PINNED | STARPU_MALLOC_COUNT | STARPU_MEMORY_OVERFLOW
        );
    }
    else {
        *new_u = NULL;
    }

    if ( v_size > 0 ) {
        *new_v = (void*)starpu_malloc_on_node_flags(
            dst_node,
            v_size * elem_size,
            STARPU_MALLOC_PINNED | STARPU_MALLOC_COUNT | STARPU_MEMORY_OVERFLOW
        );
    }
    else {
        *new_v = 0;
    }

    return 0;
}

/*
 * Free the memory buffers for low rank matrix representation
 */
rpk_int_t
runtime_rpk_lrmem_free(
    const rpk_ctx_t *ctx,
    void            *original_u,
    void            *original_v,
    size_t           u_size,
    size_t           v_size,
    size_t           elem_size
) {
    struct starpu_task* task = starpu_task_get_current();

    /* Get the data handle id associated to the rapack matrix */
	cham_starpu_rpk_ctx_t *starpu_rpk = (cham_starpu_rpk_ctx_t*)ctx->hook_data;

    /*
     * Case of no data handle associated
     */
    if ( !starpu_rpk->handle ) {
        if ( original_u ) {
            runtime_rpk_free( ctx, original_u, u_size * elem_size );
        }
        if ( original_v ) {
            runtime_rpk_free( ctx, original_v, v_size * elem_size );
        }
        return 0;
    }

	/* Get the chameleon tile interface associated to the data handle */
    starpu_cham_tile_interface_t *tile_interface = cti_handle_get_interface( starpu_rpk->handle );

    /* Get the current worker id */
    unsigned workerid = starpu_worker_get_id_check();

    /* Get the memory node of the current worker */
    unsigned dst_node = starpu_worker_get_memory_node( workerid );

    /* Free the u memory buffer */
    if ( u_size > 0 ) {
        starpu_free_on_node_flags(
            dst_node,
            (uintptr_t)original_u,
            u_size * elem_size,
            STARPU_MALLOC_PINNED | STARPU_MALLOC_COUNT | STARPU_MEMORY_OVERFLOW
        );
    }

    /* Free the v memory buffer */
    if ( v_size > 0 ) {
        starpu_free_on_node_flags(
            dst_node,
            (uintptr_t)original_v,
            v_size * elem_size,
            STARPU_MALLOC_PINNED | STARPU_MALLOC_COUNT | STARPU_MEMORY_OVERFLOW
        );
    }
    return 0;
}

/*
 * Set a low rank representation of a rapack matrix.
 */
rpk_int_t
runtime_rpk_lrmem_set(
    const rpk_ctx_t *ctx,
    rpk_matrix_t    *Ara,
    void            *new_u,
    void            *new_v
) {
    struct starpu_task* task = starpu_task_get_current();

	cham_starpu_rpk_ctx_t *starpu_rpk = (cham_starpu_rpk_ctx_t*)ctx->hook_data;

    int rc = 0;

    /*
     * If we are updating a buffer associated with a runtime-managed 
     * RAPACK matrix, we need to update the tile interface
     */
    if ( starpu_rpk->handle ) {
    	assert( task );

        /* Get the chameleon tile interface associated to the data handle */
        starpu_cham_tile_interface_t *tile_interface = cti_handle_get_interface( starpu_rpk->handle );

		rpk_matrix_t *tile_Ara = (rpk_matrix_t*)(tile_interface->tile.mat);
		tile_Ara->u = new_u;
		tile_Ara->v = new_v;
    }

    Ara->u = new_u;
    Ara->v = new_v;

    return rc;
}
#endif
