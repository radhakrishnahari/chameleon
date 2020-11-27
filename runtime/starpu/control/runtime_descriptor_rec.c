/**
 *
 * @file starpu/runtime_descriptor.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU recursive descriptor routines
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Gwenole Lucas
 * @date 2020-03-03
 *
 */
#include "chameleon_starpu.h"

#if defined(CHAMELEON_USE_BUBBLE)
static void
runtime_recdesc_partition_submit_rec( CHAM_desc_t *desc )
{
    int m, n;
    for ( m=0; m<desc->mt; m++ ) {
        for ( n=0; n<desc->nt; n++ ) {
            CHAM_tile_t *tile = desc->get_blktile( desc, m, n );
            starpu_data_handle_t  handle = (starpu_data_handle_t)RUNTIME_data_getaddr(desc, m, n);

            if ( !(tile->format & CHAMELEON_TILE_DESC) ) {
                continue;
            }

            CHAM_desc_t *tiledesc = (CHAM_desc_t*)(tile->mat);
            int nparts = tiledesc->nt * tiledesc->mt;
            struct starpu_data_filter f =
                {
                    .filter_func = starpu_CHAM_tile_filter_square_block,
                    .nchildren = nparts,
                    .filter_arg_ptr = tiledesc
                };

            starpu_data_partition_plan( handle, &f, tiledesc->schedopt );

            /* starpu_data_partition_submit_sequential_consistency( handle, nparts, tiledesc->schedopt, 0 ); */
            RUNTIME_recdesc_partition_submit( tiledesc );
        }
    }
}

void
RUNTIME_recdesc_partition_submit( CHAM_desc_t *desc )
{
    runtime_recdesc_partition_submit_rec(desc);
    starpu_task_wait_for_all();
}

static void
runtime_recdesc_unpartition_submit_rec( CHAM_desc_t *desc )
{
    int m, n;
    for ( m=0; m<desc->mt; m++ ) {
        for ( n=0; n<desc->nt; n++ ) {
            CHAM_tile_t *tile = desc->get_blktile( desc, m, n );

            if ( !(tile->format & CHAMELEON_TILE_DESC) ) {
                continue;
            }

            CHAM_desc_t *tiledesc = (CHAM_desc_t*)(tile->mat);
            RUNTIME_recdesc_unpartition_submit( tiledesc );

            starpu_data_handle_t handle = (starpu_data_handle_t)RUNTIME_data_getaddr(desc, m, n);
            int nparts = tiledesc->nt * tiledesc->mt;

            /* Because of automatic partitionning we need to first check */
            /* if the data needs to be unpartitionned <- This fix is not enough */
            /* starpu_data_unpartition_submit( handle, nparts, tiledesc->schedopt, STARPU_MAIN_RAM ); */
            starpu_data_partition_clean( handle, nparts, tiledesc->schedopt );

            /*
             * Make sure the handles are set to NULL to avoid double
             * starpu_data_unregister_submit() by RUNTIME_Destroy() and
             * starpu_data_partition_clean()
             */
            memset( tiledesc->schedopt, 0, nparts * sizeof(starpu_data_handle_t) );
        }
    }
}

void
RUNTIME_recdesc_unpartition_submit( CHAM_desc_t *desc )
{
    runtime_recdesc_unpartition_submit_rec(desc);
    starpu_task_wait_for_all();
}
#else

void
RUNTIME_recdesc_partition_submit( CHAM_desc_t *desc )
{
    (void)desc;
}

void
RUNTIME_recdesc_unpartition_submit( CHAM_desc_t *desc )
{
    (void)desc;
}

#endif /* defined(CHAMELEON_USE_BUBBLE) */
