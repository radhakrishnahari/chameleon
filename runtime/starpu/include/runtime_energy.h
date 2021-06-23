/**
 *
 * @file starpu/runtime_codelet_profile.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU codelet profiling header
 *
 * @version 1.1.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2020-12-01
 *
 */
#ifndef _runtime_energy_h_
#define _runtime_energy_h_

/* #define CHAMELEON_CL_ENERGY(name)                               \ */
/*     extern struct starpu_perfmodel cl_##name##_energy_model; */

extern struct starpu_task *__chameleon_starpu_energy_task;

int starpu_energy_start(int workerid, enum starpu_worker_archtype archi);
int starpu_energy_stop(struct starpu_perfmodel *model, struct starpu_task *task, unsigned nimpl, unsigned ntasks, int workerid, enum starpu_worker_archtype archi);

void RUNTIME_start_energy();
void RUNTIME_stop_energy( cham_flttype_t flttype, cham_tasktype_t kernel );

int chameleon_starpu_register_energy_task();

#endif /* _runtime_energy_h */
