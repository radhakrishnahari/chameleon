/**
 *
 * @file starpu/runtime_zprofiling.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU CHAMELEON_Complex64_t kernel progiling
 *
 * @version 1.2.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2022-02-22
 * @precisions normal z -> s d c
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

void RUNTIME_zdisplay_allprofile()
{

    /* BLAS 3 */
    profiling_display_zgemm_info();
#if defined(PRECISION_z) || defined(PRECISION_c)
    profiling_display_zhemm_info();
    profiling_display_zher2k_info();
    profiling_display_zherk_info();
    profiling_display_zsytrf_nopiv_info();
#endif
    profiling_display_zsymm_info();
    profiling_display_zsyr2k_info();
    profiling_display_zsyrk_info();
    profiling_display_ztrmm_info();
    profiling_display_ztrsm_info();

    /* Lapack */
    profiling_display_zgelqt_info();
    profiling_display_zgeqrt_info();
    profiling_display_zgessm_info();
    profiling_display_zgetrf_incpiv_info();
    profiling_display_zgetrf_info();
    profiling_display_zgetrf_nopiv_info();
    profiling_display_zlauum_info();
    profiling_display_zpotrf_info();
    profiling_display_zssssm_info();
    profiling_display_ztplqt_info();
    profiling_display_ztpmlqt_info();
    profiling_display_ztpmqrt_info();
    profiling_display_ztpqrt_info();
    profiling_display_ztrtri_info();
    profiling_display_ztstrf_info();
    profiling_display_zunmlq_info();
    profiling_display_zunmqr_info();

    profiling_display_zlange_info();
}

void RUNTIME_zdisplay_oneprofile( cham_tasktype_t kernel )
{
    switch( kernel ) {
        /* Blas 3 */
    case TASK_GEMM:         profiling_display_zgemm_info();         break;
#if defined(PRECISION_z) || defined(PRECISION_c)
    case TASK_HEMM:         profiling_display_zhemm_info();         break;
    case TASK_HER2K:        profiling_display_zher2k_info();        break;
    case TASK_HERK:         profiling_display_zherk_info();         break;
    case TASK_SYTRF_NOPIV:  profiling_display_zsytrf_nopiv_info();  break;
#endif
    case TASK_SYMM:         profiling_display_zsymm_info();         break;
    case TASK_SYR2K:        profiling_display_zsyr2k_info();        break;
    case TASK_SYRK:         profiling_display_zsyrk_info();         break;
    case TASK_TRMM:         profiling_display_ztrmm_info();         break;
    case TASK_TRSM:         profiling_display_ztrsm_info();         break;

        /* Lapack */
    case TASK_GELQT:        profiling_display_zgelqt_info();        break;
    case TASK_GEQRT:        profiling_display_zgeqrt_info();        break;
    case TASK_GESSM:        profiling_display_zgessm_info();        break;
    case TASK_GETRF:        profiling_display_zgetrf_info();        break;
    case TASK_GETRF_INCPIV: profiling_display_zgetrf_incpiv_info(); break;
    case TASK_GETRF_NOPIV:  profiling_display_zgetrf_nopiv_info();  break;
    case TASK_LAUUM:        profiling_display_zlauum_info();        break;
    case TASK_POTRF:        profiling_display_zpotrf_info();        break;
    case TASK_SSSSM:        profiling_display_zssssm_info();        break;
    case TASK_TRTRI:        profiling_display_ztrtri_info();        break;
    case TASK_TSTRF:        profiling_display_ztstrf_info();        break;

    case TASK_TPLQT:        profiling_display_ztplqt_info();        break;
    case TASK_TPMLQT:       profiling_display_ztpmlqt_info();        break;
    case TASK_TPMQRT:       profiling_display_ztpmqrt_info();        break;
    case TASK_TPQRT:        profiling_display_ztpqrt_info();        break;

    case TASK_UNMLQ:        profiling_display_zunmlq_info();        break;
    case TASK_UNMQR:        profiling_display_zunmqr_info();        break;

    case TASK_LANGE:        profiling_display_zlange_info();        break;

    default:
        return;
    }
}

void
RUNTIME_zget_codelet( cham_tasktype_t kernel )
{
    struct starpu_task *task = __chameleon_starpu_energy_task;

    switch( kernel ) {
        /* Blas 3 */
    case TASK_GEMM:
        task->cl    = &cl_zgemm;
        task->flops = flops_zgemm( cti_handle_get_m(task->handles[2]),
                                   cti_handle_get_n(task->handles[2]),
                                   cti_handle_get_n(task->handles[0]) );
        break;

/* #if defined(PRECISION_z) || defined(PRECISION_c) */
/*     case TASK_HEMM:         return &cl_zhemm; */
/*     case TASK_HER2K:        return &cl_zher2k; */
/*     case TASK_HERK:         return &cl_zherk; */
/*     case TASK_SYTRF_NOPIV:  return &zsytrf_nopiv; */
/* #endif */
/*     case TASK_SYMM:         return &cl_zsymm; */
/*     case TASK_SYR2K:        return &cl_zsyr2k; */
/*     case TASK_SYRK:         return &cl_zsyrk; */
/*     case TASK_TRMM:         return &cl_ztrmm; */
/*     case TASK_TRSM:         return &cl_ztrsm; */

/*         /\* Lapack *\/ */
/*     case TASK_GELQT:        return &cl_zgelqt; */
/*     case TASK_GEQRT:        return &cl_zgeqrt; */
/*     case TASK_GESSM:        return &cl_zgessm; */
/*     case TASK_GETRF:        return &cl_zgetrf; */
/*     case TASK_GETRF_INCPIV: return &cl_zgetrf_incpiv; */
/*     case TASK_GETRF_NOPIV:  return &cl_zgetrf_nopiv; */
/*     case TASK_LAUUM:        return &cl_zlauum; */
/*     case TASK_POTRF:        return &cl_zpotrf; */
/*     case TASK_SSSSM:        return &cl_zssssm; */
/*     case TASK_TRTRI:        return &cl_ztrtri; */
/*     case TASK_TSTRF:        return &cl_ztstrf; */

/*     case TASK_TPLQT:        return &cl_ztplqt; */
/*     case TASK_TPMLQT:       return &cl_ztpmlqt; */
/*     case TASK_TPMQRT:       return &cl_ztpmqrt; */
/*     case TASK_TPQRT:        return &cl_ztpqrt; */

/*     case TASK_UNMLQ:        return &cl_zunmlq; */
/*     case TASK_UNMQR:        return &cl_zunmqr; */

/*     case TASK_LANGE:        return &cl_zlange; */

    default:
        return NULL;
    }
}
