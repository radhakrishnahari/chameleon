/**
 *
 * @file starpu/runtime_profiling.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU profiling routines
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Samuel Thibault
 * @date 2020-03-03
 *
 */
#include <math.h>
#include "chameleon_starpu.h"
#include "starpu.h"
#include "runtime_codelet_z.h"
#include "starpu_perfmodel.h"
#include <starpu_fxt.h>
#include "chameleon/config.h"
#include "runtime_codelet_profile.h"
#include "runtime_codelets.h"

struct starpu_task *__chameleon_starpu_energy_task  = NULL;
static int          __chameleon_starpu_energy_ntask = 0;

/* static struct starpu_perfmodel cl_##cl_name##_model = */
/* { */
/*     .type = STARPU_HISTORY_BASED, */
/*     .symbol = #cl_name */
/* }; */

/* static struct starpu_perfmodel cl_##cl_name##_energy_model = */
/* { */
/*     .type = STARPU_HISTORY_BASED, */
/*     .symbol = #cl_name"_energy" */
/* }; */

void RUNTIME_start_energy() {
    starpu_energy_start(-1, STARPU_CPU_WORKER);

    assert( __chameleon_starpu_energy_task == NULL );
    assert( __chameleon_starpu_energy_ntask == 0 );
}

int chameleon_starpu_register_energy_task()
{
    if ( __chameleon_starpu_energy_task == NULL ) {
        __chameleon_starpu_energy_task = starpu_task_create();
        __chameleon_starpu_energy_task->cl          = NULL;
        __chameleon_starpu_energy_task->synchronous = 1;
        __chameleon_starpu_energy_task->destroy     = 0;

        __chameleon_starpu_energy_ntask = 1;
        return 1;
    }
    else {
        __chameleon_starpu_energy_ntask++;
        return 0;
    }
}

/* void create_fake_task_and_count_total_tasks( const char* name, struct starpu_codelet *cl ) */
/* { */
/*     struct starpu_task *task = starpu_task_create(); */

/*     task->cl         = cl_zgemm; */
/*     task->synchronous = 1; */
/*     task->destroy     = 0; */
/*     task->flops      = flops_zgemm( m, n, k ); */

/*     static int times = 0; */
/*     if(times == 0) */
/*     { */
/*         task = starpu_task_create(); */
/*         task->cl = cl; */
/*         if( strcmp(name, "zgemm") ){ */
/*             energy = &cl_zgemm_energy_model; */
/*         } */
/*         else if ( strcmp(name, "ztrsm") ){ */
/*             //task->cl = &cl_ztrsm; */
/*             energy = &cl_ztrsm_energy_model; */
/*         } */
/*         else if ( strcmp(name, "zherk") ){ */
/*             //task->cl = &cl_zherk; */
/*             energy = &cl_zherk_energy_model; */
/*         } */
/*         else if ( strcmp(name, "zpotrf") ){ */
/*             //task->cl = &cl_zpotrf; */
/*             energy = &cl_zpotrf_energy_model; */
/*         } */
/*         else fprintf(stderr, "unhandled kernel"); */
/*         times = 1; */
/*     } */
/*     ntasks ++; */
/* } */

void RUNTIME_stop_energy( cham_flttype_t flttype, cham_tasktype_t kernel )
{
    struct starpu_task *task = __chameleon_starpu_energy_task;

    if ( task == NULL ) {
        fprintf( stderr, "RUNTIME_stop_energy: energy task not initialized\n");
        return;
    }

    switch( flttype ) {
    case ChamComplexDouble:
        RUNTIME_zget_codelet( kernel );
        break;

    case ChamComplexFloat:
        RUNTIME_cget_codelet( kernel );
        break;

    case ChamRealDouble:
        RUNTIME_dget_codelet( kernel );
        break;

    case ChamRealFloat:
        RUNTIME_sget_codelet( kernel );
        break;

    default:
        fprintf( stderr, "RUNTIME_stop_energy: Unknown floating point arithmetic\n");
        return;
    }

    if ( task->cl == NULL ) {
        fprintf( stderr, "RUNTIME_stop_energy: Could not find the associated codelet\n");
        goto cleanup;
    }

    starpu_energy_stop( task->cl->energy_model, task, 0,
                        __chameleon_starpu_energy_ntask, -1, STARPU_CPU_WORKER );

  cleanup:
    starpu_task_destroy(task);
    __chameleon_starpu_energy_task = NULL;
    __chameleon_starpu_energy_ntask = 0;
}
