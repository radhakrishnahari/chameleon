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
#include "runtime_profiling.h"

static struct starpu_perfmodel cl_zgemm_model =
{
    .type = STARPU_HISTORY_BASED,
    .symbol = "zgemm"
};

static struct starpu_perfmodel cl_zgemm_energy_model =
{
    .type = STARPU_HISTORY_BASED,
    .symbol = "zgemm_energy"
};

void RUNTIME_start_energy(){

    unsigned worker;
    
    starpu_energy_start(-1, STARPU_CPU_WORKER);
}

void RUNTIME_stop_energy_gemm(){

  
    struct starpu_task *task = starpu_task_create(); 
   
    int ntasks=starpu_worker_get_count_by_type(STARPU_CPU_WORKER);
    //callback = options->profiling ? cl_zgemm_callback : NULL;

    task->cl = &cl_zgemm;
    printf("nbre tasks %d",ntasks);    
    starpu_energy_stop(&cl_zgemm_energy_model, task , 0, ntasks, -1, STARPU_CPU_WORKER);
    //starpu_task_destroy (task);
}


