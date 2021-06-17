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
#include "runtime_energy.h"

struct starpu_task * task = NULL;
struct starpu_perfmodel *energy;
static int ntasks = 0;

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

void RUNTIME_start_energy(){

    unsigned worker;
    starpu_energy_start(-1, STARPU_CPU_WORKER);
}

void create_fake_task_and_count_total_tasks(const char* name, struct starpu_codelet *cl)
{
    static int times = 0;
    if(times == 0)
    {
        task = starpu_task_create();
        task->cl = cl;
        if( strcmp(name, "zgemm") ){
            energy = &cl_zgemm_energy_model;
        }
        else if ( strcmp(name, "ztrsm") ){
            //task->cl = &cl_ztrsm;
            energy = &cl_ztrsm_energy_model;
        }
        else if ( strcmp(name, "zherk") ){
            //task->cl = &cl_zherk;
            energy = &cl_zherk_energy_model;
        }
        else if ( strcmp(name, "zpotrf") ){
            //task->cl = &cl_zpotrf;
            energy = &cl_zpotrf_energy_model;
        }
        else fprintf(stderr, "unhandled kernel");
        times = 1;
    }
    ntasks ++;
}

void RUNTIME_stop_energy(){

    if(task == NULL)
    {
        fprintf(stderr, "task for energy measurement not created\n");
        exit(0);
    }
    starpu_energy_stop(energy, task , 0, ntasks, -1, STARPU_CPU_WORKER);
    //starpu_task_destroy (task);
}


