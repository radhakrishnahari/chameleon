/**
 *
 * @file runtime_fake.c
 *
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon runtime routines
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2020-10-28
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>

// #include "runtime_fake.h"
#include "control/auxiliary.h"
#include "control/common.h"
#include "chameleon/runtime.h"

/* declaration of fake functions */
void (*_RUNTIME_context_create)( CHAM_context_t *ctxt ) = NULL;
void (*_RUNTIME_context_destroy)( CHAM_context_t *ctxt ) = NULL;
void (*_RUNTIME_enable)( void *runtime_ctxt, int option ) = NULL;
void (*_RUNTIME_disable)( void *runtime_ctxt, int option ) = NULL;
int (*_RUNTIME_init)( CHAM_context_t *ctxt, int ncpus, int ncudas, int nthreads_per_worker ) = NULL;
void (*_RUNTIME_finalize)( CHAM_context_t *ctxt ) = NULL;
void (*_RUNTIME_pause)( CHAM_context_t *ctxt ) = NULL;
void (*_RUNTIME_resume)( CHAM_context_t *ctxt ) = NULL;
void (*_RUNTIME_barrier)( CHAM_context_t *ctxt ) = NULL;
void (*_RUNTIME_progress)( CHAM_context_t *ctxt ) = NULL;
int (*_RUNTIME_thread_rank)( CHAM_context_t *ctxt ) = NULL;
int (*_RUNTIME_thread_size)( CHAM_context_t *ctxt ) = NULL;
int (*_RUNTIME_comm_rank)( CHAM_context_t *ctxt ) = NULL;
int (*_RUNTIME_comm_size)( CHAM_context_t *ctxt ) = NULL;
void (*_RUNTIME_comm_set_tag_sizes)( int user_tag_width, int user_tag_sep ) = NULL;
int (*_RUNTIME_sequence_create)( CHAM_context_t *ctxt, RUNTIME_sequence_t *sequence ) = NULL;
int (*_RUNTIME_sequence_destroy)( CHAM_context_t *ctxt, RUNTIME_sequence_t *sequence) = NULL;
int (*_RUNTIME_sequence_wait)( CHAM_context_t *ctxt, RUNTIME_sequence_t *sequence ) = NULL;
void (*_RUNTIME_sequence_flush)( CHAM_context_t *ctxt,
                                 RUNTIME_sequence_t *sequence,
                                 RUNTIME_request_t  *request,
                                 int                 status ) = NULL;
int (*_RUNTIME_request_create)( CHAM_context_t *ctxt, RUNTIME_request_t *request ) = NULL;
int (*_RUNTIME_request_destroy)( CHAM_context_t *ctxt, RUNTIME_request_t *request) = NULL;
int (*_RUNTIME_request_set)( CHAM_context_t  *chamctxt,
                             RUNTIME_request_t *request,
                             int param, int value ) = NULL;
void *(*_RUNTIME_malloc)( size_t size ) = NULL;
void (*_RUNTIME_free)( void *ptr, size_t size ) = NULL;
void (*_RUNTIME_desc_create)( CHAM_desc_t *desc ) = NULL;
void (*_RUNTIME_desc_destroy)( CHAM_desc_t *desc ) = NULL;
int (*_RUNTIME_desc_acquire)( const CHAM_desc_t *desc ) = NULL;
int (*_RUNTIME_desc_release)( const CHAM_desc_t *desc ) = NULL;
void (*_RUNTIME_desc_flush)( const CHAM_desc_t *desc, const RUNTIME_sequence_t *sequence ) = NULL;
void (*_RUNTIME_flush)( ) = NULL;
void (*_RUNTIME_data_flush)( const RUNTIME_sequence_t *sequence,
                             const CHAM_desc_t *A, int Am, int An ) = NULL;
#if defined(CHAMELEON_USE_MIGRATE)
void (*_RUNTIME_data_migrate)( const RUNTIME_sequence_t *sequence,
                               const CHAM_desc_t *A, int Am, int An, int new_rank ) = NULL;
#endif
void (*_RUNTIME_options_init)( RUNTIME_option_t   *options,
                               CHAM_context_t  *context,
                               RUNTIME_sequence_t *sequence,
                               RUNTIME_request_t  *request ) = NULL;
void (*_RUNTIME_options_finalize)( RUNTIME_option_t  *options,
                                   CHAM_context_t *context ) = NULL;
int (*_RUNTIME_options_ws_alloc)( RUNTIME_option_t *options,
                                  size_t wsize, size_t hsize ) = NULL;
int (*_RUNTIME_options_ws_free)( RUNTIME_option_t *options ) = NULL;
void (*_RUNTIME_zlocality_allrestore) () = NULL;
void (*_RUNTIME_clocality_allrestore) () = NULL;
void (*_RUNTIME_dlocality_allrestore) () = NULL;
void (*_RUNTIME_slocality_allrestore) () = NULL;
void (*_RUNTIME_zlocality_allrestrict)(uint32_t device) = NULL;
void (*_RUNTIME_zlocality_onerestrict)(cham_tasktype_t task, uint32_t device) = NULL;
void (*_RUNTIME_zlocality_onerestore) (cham_tasktype_t task) = NULL;
void (*_RUNTIME_clocality_allrestrict)(uint32_t device) = NULL;
void (*_RUNTIME_clocality_onerestrict)(cham_tasktype_t task, uint32_t device) = NULL;
void (*_RUNTIME_clocality_onerestore) (cham_tasktype_t task) = NULL;
void (*_RUNTIME_dlocality_allrestrict)(uint32_t device) = NULL;
void (*_RUNTIME_dlocality_onerestrict)(cham_tasktype_t task, uint32_t device) = NULL;
void (*_RUNTIME_dlocality_onerestore) (cham_tasktype_t task) = NULL;
void (*_RUNTIME_slocality_allrestrict)(uint32_t device) = NULL;
void (*_RUNTIME_slocality_onerestrict)(cham_tasktype_t task, uint32_t device) = NULL;
void (*_RUNTIME_slocality_onerestore) (cham_tasktype_t task) = NULL;
void   (*_RUNTIME_schedprofile_display) () = NULL;
void   (*_RUNTIME_kernelprofile_display)() = NULL;
double (*_RUNTIME_get_time)() = NULL;
void  (*_RUNTIME_iteration_push)     (CHAM_context_t*, unsigned long iteration) = NULL;
void  (*_RUNTIME_iteration_pop)      (CHAM_context_t*) = NULL;
void (*_RUNTIME_start_profiling)() = NULL;
void (*_RUNTIME_stop_profiling)() = NULL;
void (*_RUNTIME_start_stats)() = NULL;
void (*_RUNTIME_stop_stats)() = NULL;
void (*_RUNTIME_zdisplay_allprofile) () = NULL;
void (*_RUNTIME_zdisplay_oneprofile) (cham_tasktype_t task) = NULL;
void (*_RUNTIME_cdisplay_allprofile) () = NULL;
void (*_RUNTIME_cdisplay_oneprofile) (cham_tasktype_t task) = NULL;
void (*_RUNTIME_ddisplay_allprofile) () = NULL;
void (*_RUNTIME_ddisplay_oneprofile) (cham_tasktype_t task) = NULL;
void (*_RUNTIME_sdisplay_allprofile) () = NULL;
void (*_RUNTIME_sdisplay_oneprofile) (cham_tasktype_t task) = NULL;

void (*_INSERT_TASK_map)( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, const CHAM_desc_t *A, int Am, int An,
                          cham_unary_operator_t op_fct, void *op_args ) = NULL;

void (*_INSERT_TASK_slag2c)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int m, int n,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_scasum)( const RUNTIME_option_t *options,
                             cham_store_t storev, cham_uplo_t uplo, int M, int N,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_caxpy)( const RUNTIME_option_t *options,
                            int M, CHAMELEON_Complex32_t alpha,
                            const CHAM_desc_t *A, int Am, int An, int incA,
                            const CHAM_desc_t *B, int Bm, int Bn, int incB ) = NULL;

void (*_INSERT_TASK_cbuild)( const RUNTIME_option_t *options,
                             const CHAM_desc_t *A, int Am, int An,
                             void *user_data, void* user_build_callback ) = NULL;

void (*_INSERT_TASK_cgeadd)( const RUNTIME_option_t *options,
                             cham_trans_t trans, int m, int n, int nb,
                             CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                             CHAMELEON_Complex32_t beta,  const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_cgelqt)( const RUNTIME_option_t *options,
                             int m, int n, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_cgemv) ( const RUNTIME_option_t *options,
                             cham_trans_t trans, int m, int n,
                             CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *X, int Xm, int Xn, int incX,
                             CHAMELEON_Complex32_t beta,  const CHAM_desc_t *Y, int Ym, int Yn, int incY ) = NULL;

void (*_INSERT_TASK_cgemm)( const RUNTIME_option_t *options,
                            cham_trans_t transA, cham_trans_t transB,
                            int m, int n, int k, int nb,
                            CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn,
                            CHAMELEON_Complex32_t beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_cgeqrt)( const RUNTIME_option_t *options,
                             int m, int n, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_cgessm)( const RUNTIME_option_t *options,
                             int m, int n, int k, int ib, int nb,
                             int *IPIV,
                             const CHAM_desc_t *L, int Lm, int Ln,
                             const CHAM_desc_t *D, int Dm, int Dn,
                             const CHAM_desc_t *A, int Am, int An ) = NULL;

void (*_INSERT_TASK_cgessq)( const RUNTIME_option_t *options,
                             cham_store_t storev, int m, int n,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) = NULL;

void (*_INSERT_TASK_cgetrf)( const RUNTIME_option_t *options,
                             int m, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             int *IPIV,
                             cham_bool_t check_info, int iinfo ) = NULL;

void (*_INSERT_TASK_cgetrf_incpiv)( const RUNTIME_option_t *options,
                                    int m, int n, int ib, int nb,
                                    const CHAM_desc_t *A, int Am, int An,
                                    const CHAM_desc_t *L, int Lm, int Ln,
                                    int *IPIV,
                                    cham_bool_t check_info, int iinfo ) = NULL;

void (*_INSERT_TASK_cgetrf_nopiv)( const RUNTIME_option_t *options,
                                   int m, int n, int ib, int nb,
                                   const CHAM_desc_t *A, int Am, int An, int iinfo ) = NULL;

void (*_INSERT_TASK_che2ge)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo,
                             int m, int n, int mb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_chemm)( const RUNTIME_option_t *options,
                            cham_side_t side, cham_uplo_t uplo,
                            int m, int n, int nb,
                            CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn,
                            CHAMELEON_Complex32_t beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_cher2k)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_trans_t trans,
                             int n, int k, int nb,
                             CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn,
                             float beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_cherfb)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo,
                             int n, int k, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn,
                             const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_cherk)( const RUNTIME_option_t *options,
                            cham_uplo_t uplo, cham_trans_t trans,
                            int n, int k, int nb,
                            float alpha, const CHAM_desc_t *A, int Am, int An,
                            float beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_chessq)( const RUNTIME_option_t *options,
                             cham_store_t storev, cham_uplo_t uplo, int n,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) = NULL;

void (*_INSERT_TASK_clacpy)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int m, int n, int mb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_clacpyx)( const RUNTIME_option_t *options,
                              cham_uplo_t uplo, int m, int n, int mb,
                              int displA, const CHAM_desc_t *A, int Am, int An,
                              int displB, const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_clange)( const RUNTIME_option_t *options,
                             cham_normtype_t norm, int M, int N, int NB,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_clange_max)( const RUNTIME_option_t *options,
                                 const CHAM_desc_t *A, int Am, int An,
                                 const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_clanhe)( const RUNTIME_option_t *options,
                             cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_clansy)( const RUNTIME_option_t *options,
                             cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_clantr)( const RUNTIME_option_t *options,
                             cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag,
                             int M, int N, int NB,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_clascal)( const RUNTIME_option_t *options,
                              cham_uplo_t uplo,
                              int m, int n, int nb,
                              CHAMELEON_Complex32_t alpha,
                              const CHAM_desc_t *A, int Am, int An ) = NULL;

void (*_INSERT_TASK_claset)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int n1, int n2,
                             CHAMELEON_Complex32_t alpha, CHAMELEON_Complex32_t beta,
                             const CHAM_desc_t *tileA, int tileAm, int tileAn ) = NULL;

void (*_INSERT_TASK_claset2)( const RUNTIME_option_t *options,
                              cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex32_t alpha,
                              const CHAM_desc_t *tileA, int tileAm, int tileAn ) = NULL;

void (*_INSERT_TASK_clatro)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_trans_t trans, int m, int n, int mb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_clauum)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An ) = NULL;

void (*_INSERT_TASK_cplghe)( const RUNTIME_option_t *options,
                             float bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                             int bigM, int m0, int n0, unsigned long long int seed ) = NULL;

void (*_INSERT_TASK_cplgsy)( const RUNTIME_option_t *options,
                             CHAMELEON_Complex32_t bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                             int bigM, int m0, int n0, unsigned long long int seed ) = NULL;

void (*_INSERT_TASK_cplrnt)( const RUNTIME_option_t *options,
                             int m, int n, const CHAM_desc_t *A, int Am, int An,
                             int bigM, int m0, int n0, unsigned long long int seed ) = NULL;

void (*_INSERT_TASK_cplssq)( const RUNTIME_option_t *options,
                             cham_store_t storev, int M, int N,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                             const CHAM_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn ) = NULL;

void (*_INSERT_TASK_cplssq2)( const RUNTIME_option_t *options, int N,
                              const CHAM_desc_t *RESULT, int RESULTm, int RESULTn ) = NULL;

void (*_INSERT_TASK_cpotrf)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             int iinfo ) = NULL;

void (*_INSERT_TASK_cssssm)( const RUNTIME_option_t *options,
                             int m1, int n1, int m2, int n2, int k, int ib, int nb,
                             const CHAM_desc_t *A1, int A1m, int A1n,
                             const CHAM_desc_t *A2, int A2m, int A2n,
                             const CHAM_desc_t *L1, int L1m, int L1n,
                             const CHAM_desc_t *L2, int L2m, int L2n,
                             const int *IPIV ) = NULL;

void (*_INSERT_TASK_csymm)( const RUNTIME_option_t *options,
                            cham_side_t side, cham_uplo_t uplo,
                            int m, int n, int nb,
                            CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn,
                            CHAMELEON_Complex32_t beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_csyr2k)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_trans_t trans,
                             int n, int k, int nb,
                             CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn,
                             CHAMELEON_Complex32_t beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_csyrk)( const RUNTIME_option_t *options,
                            cham_uplo_t uplo, cham_trans_t trans,
                            int n, int k, int nb,
                            CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                            CHAMELEON_Complex32_t beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_csyssq)( const RUNTIME_option_t *options,
                             cham_store_t storev, cham_uplo_t uplo, int n,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) = NULL;

void (*_INSERT_TASK_csytrf_nopiv)( const RUNTIME_option_t *options,
                                   cham_uplo_t uplo, int n, int nb,
                                   const CHAM_desc_t *A, int Am, int An,
                                   int iinfo ) = NULL;

void (*_INSERT_TASK_ctplqt)( const RUNTIME_option_t *options,
                             int m, int n, int l, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn,
                             const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_ctpmlqt)( const RUNTIME_option_t *options,
                              cham_side_t side, cham_trans_t trans,
                              int M, int N, int K, int L, int ib, int nb,
                              const CHAM_desc_t *V, int Vm, int Vn,
                              const CHAM_desc_t *T, int Tm, int Tn,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_ctpmqrt)( const RUNTIME_option_t *options,
                              cham_side_t side, cham_trans_t trans,
                              int m, int n, int k, int l, int ib, int nb,
                              const CHAM_desc_t *V, int Vm, int Vn,
                              const CHAM_desc_t *T, int Tm, int Tn,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_ctpqrt)( const RUNTIME_option_t *options,
                             int m, int n, int l, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn,
                             const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_ctradd)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                             CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                             CHAMELEON_Complex32_t beta,  const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_ctrasm)( const RUNTIME_option_t *options,
                             cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_ctrmm)( const RUNTIME_option_t *options,
                            cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                            int m, int n, int nb,
                            CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_ctrsm)( const RUNTIME_option_t *options,
                            cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                            int m, int n, int nb,
                            CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_ctrssq)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_diag_t diag,
                             int m, int n,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) = NULL;

void (*_INSERT_TASK_ctrtri)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_diag_t diag, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             int iinfo ) = NULL;

void (*_INSERT_TASK_ctsmlq_hetra1)( const RUNTIME_option_t *options,
                                    cham_side_t side, cham_trans_t trans,
                                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                    const CHAM_desc_t *A1, int A1m, int A1n,
                                    const CHAM_desc_t *A2, int A2m, int A2n,
                                    const CHAM_desc_t *V, int Vm, int Vn,
                                    const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_ctsmqr_hetra1)( const RUNTIME_option_t *options,
                                    cham_side_t side, cham_trans_t trans,
                                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                    const CHAM_desc_t *A1, int A1m, int A1n,
                                    const CHAM_desc_t *A2, int A2m, int A2n,
                                    const CHAM_desc_t *V, int Vm, int Vn,
                                    const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_ctstrf)( const RUNTIME_option_t *options,
                             int m, int n, int ib, int nb,
                             const CHAM_desc_t *U, int Um, int Un,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *L, int Lm, int Ln,
                             int *IPIV,
                             cham_bool_t check_info, int iinfo ) = NULL;

void (*_INSERT_TASK_cunmlq)( const RUNTIME_option_t *options,
                             cham_side_t side, cham_trans_t trans,
                             int m, int n, int ib,  int nb, int k,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn,
                             const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_cunmqr)( const RUNTIME_option_t *options,
                             cham_side_t side, cham_trans_t trans,
                             int m, int n, int k, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn,
                             const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_cgram)( const RUNTIME_option_t *options,
                            cham_uplo_t uplo,
                            int m, int n, int mt, int nt,
                            const CHAM_desc_t *Di, int Dim, int Din,
                            const CHAM_desc_t *Dj, int Djm, int Djn,
                            const CHAM_desc_t *D, int Dm, int Dn,
                            CHAM_desc_t *A, int Am, int An) = NULL;

void (*_INSERT_TASK_dlag2z)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int m, int n,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_dasum)( const RUNTIME_option_t *options,
                            cham_store_t storev, cham_uplo_t uplo, int M, int N,
                            const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_daxpy)( const RUNTIME_option_t *options,
                            int M, double alpha,
                            const CHAM_desc_t *A, int Am, int An, int incA,
                            const CHAM_desc_t *B, int Bm, int Bn, int incB ) = NULL;

void (*_INSERT_TASK_dbuild)( const RUNTIME_option_t *options,
                             const CHAM_desc_t *A, int Am, int An,
                             void *user_data, void* user_build_callback ) = NULL;

void (*_INSERT_TASK_dgeadd)( const RUNTIME_option_t *options,
                             cham_trans_t trans, int m, int n, int nb,
                             double alpha, const CHAM_desc_t *A, int Am, int An,
                             double beta,  const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_dgelqt)( const RUNTIME_option_t *options,
                             int m, int n, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_dgemv) ( const RUNTIME_option_t *options,
                             cham_trans_t trans, int m, int n,
                             double alpha, const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *X, int Xm, int Xn, int incX,
                             double beta,  const CHAM_desc_t *Y, int Ym, int Yn, int incY ) = NULL;

void (*_INSERT_TASK_dgemm)( const RUNTIME_option_t *options,
                            cham_trans_t transA, cham_trans_t transB,
                            int m, int n, int k, int nb,
                            double alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn,

                            double beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;
void (*_INSERT_TASK_dgeqrt)( const RUNTIME_option_t *options,
                             int m, int n, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_dgessm)( const RUNTIME_option_t *options,
                             int m, int n, int k, int ib, int nb,
                             int *IPIV,
                             const CHAM_desc_t *L, int Lm, int Ln,
                             const CHAM_desc_t *D, int Dm, int Dn,
                             const CHAM_desc_t *A, int Am, int An ) = NULL;

void (*_INSERT_TASK_dgessq)( const RUNTIME_option_t *options,
                             cham_store_t storev, int m, int n,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) = NULL;

void (*_INSERT_TASK_dgetrf)( const RUNTIME_option_t *options,
                             int m, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             int *IPIV,
                             cham_bool_t check_info, int iinfo ) = NULL;

void (*_INSERT_TASK_dgetrf_incpiv)( const RUNTIME_option_t *options,
                                    int m, int n, int ib, int nb,
                                    const CHAM_desc_t *A, int Am, int An,
                                    const CHAM_desc_t *L, int Lm, int Ln,
                                    int *IPIV,
                                    cham_bool_t check_info, int iinfo ) = NULL;

void (*_INSERT_TASK_dgetrf_nopiv)( const RUNTIME_option_t *options,
                                   int m, int n, int ib, int nb,
                                   const CHAM_desc_t *A, int Am, int An, int iinfo ) = NULL;

void (*_INSERT_TASK_dhe2ge)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo,
                             int m, int n, int mb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_dlacpy)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int m, int n, int mb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_dlacpyx)( const RUNTIME_option_t *options,
                              cham_uplo_t uplo, int m, int n, int mb,
                              int displA, const CHAM_desc_t *A, int Am, int An,
                              int displB, const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_dlange)( const RUNTIME_option_t *options,
                             cham_normtype_t norm, int M, int N, int NB,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_dlange_max)( const RUNTIME_option_t *options,
                                 const CHAM_desc_t *A, int Am, int An,
                                 const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_dlansy)( const RUNTIME_option_t *options,
                             cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_dlantr)( const RUNTIME_option_t *options,
                             cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag,
                             int M, int N, int NB,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_dlascal)( const RUNTIME_option_t *options,
                              cham_uplo_t uplo,
                              int m, int n, int nb,
                              double alpha,
                              const CHAM_desc_t *A, int Am, int An ) = NULL;

void (*_INSERT_TASK_dlaset)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int n1, int n2,
                             double alpha, double beta,
                             const CHAM_desc_t *tileA, int tileAm, int tileAn ) = NULL;

void (*_INSERT_TASK_dlaset2)( const RUNTIME_option_t *options,
                              cham_uplo_t uplo, int n1, int n2, double alpha,
                              const CHAM_desc_t *tileA, int tileAm, int tileAn ) = NULL;

void (*_INSERT_TASK_dlatro)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_trans_t trans, int m, int n, int mb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_dlauum)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An ) = NULL;

void (*_INSERT_TASK_dplgsy)( const RUNTIME_option_t *options,
                             double bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                             int bigM, int m0, int n0, unsigned long long int seed ) = NULL;

void (*_INSERT_TASK_dplrnt)( const RUNTIME_option_t *options,
                             int m, int n, const CHAM_desc_t *A, int Am, int An,
                             int bigM, int m0, int n0, unsigned long long int seed ) = NULL;

void (*_INSERT_TASK_dplssq)( const RUNTIME_option_t *options,
                             cham_store_t storev, int M, int N,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                             const CHAM_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn ) = NULL;

void (*_INSERT_TASK_dplssq2)( const RUNTIME_option_t *options, int N,
                              const CHAM_desc_t *RESULT, int RESULTm, int RESULTn ) = NULL;

void (*_INSERT_TASK_dpotrf)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             int iinfo ) = NULL;

void (*_INSERT_TASK_dssssm)( const RUNTIME_option_t *options,
                             int m1, int n1, int m2, int n2, int k, int ib, int nb,
                             const CHAM_desc_t *A1, int A1m, int A1n,
                             const CHAM_desc_t *A2, int A2m, int A2n,
                             const CHAM_desc_t *L1, int L1m, int L1n,
                             const CHAM_desc_t *L2, int L2m, int L2n,
                             const int *IPIV ) = NULL;

void (*_INSERT_TASK_dsymm)( const RUNTIME_option_t *options,
                            cham_side_t side, cham_uplo_t uplo,
                            int m, int n, int nb,
                            double alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn,
                            double beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_dsyr2k)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_trans_t trans,
                             int n, int k, int nb,
                             double alpha, const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn,
                             double beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_dsyrfb)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo,
                             int n, int k, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn,
                             const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_dsyrk)( const RUNTIME_option_t *options,
                            cham_uplo_t uplo, cham_trans_t trans,
                            int n, int k, int nb,
                            double alpha, const CHAM_desc_t *A, int Am, int An,
                            double beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_dsyssq)( const RUNTIME_option_t *options,
                             cham_store_t storev, cham_uplo_t uplo, int n,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) = NULL;

void (*_INSERT_TASK_dsytrf_nopiv)( const RUNTIME_option_t *options,
                                   cham_uplo_t uplo, int n, int nb,
                                   const CHAM_desc_t *A, int Am, int An,
                                   int iinfo ) = NULL;

void (*_INSERT_TASK_dtplqt)( const RUNTIME_option_t *options,
                             int m, int n, int l, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn,
                             const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_dtpmlqt)( const RUNTIME_option_t *options,
                              cham_side_t side, cham_trans_t trans,
                              int M, int N, int K, int L, int ib, int nb,
                              const CHAM_desc_t *V, int Vm, int Vn,
                              const CHAM_desc_t *T, int Tm, int Tn,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_dtpmqrt)( const RUNTIME_option_t *options,
                              cham_side_t side, cham_trans_t trans,
                              int m, int n, int k, int l, int ib, int nb,
                              const CHAM_desc_t *V, int Vm, int Vn,
                              const CHAM_desc_t *T, int Tm, int Tn,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_dtpqrt)( const RUNTIME_option_t *options,
                             int m, int n, int l, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn,
                             const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_dtradd)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                             double alpha, const CHAM_desc_t *A, int Am, int An,
                             double beta,  const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_dtrasm)( const RUNTIME_option_t *options,
                             cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_dtrmm)( const RUNTIME_option_t *options,
                            cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                            int m, int n, int nb,
                            double alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_dtrsm)( const RUNTIME_option_t *options,
                            cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                            int m, int n, int nb,
                            double alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_dtrssq)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_diag_t diag,
                             int m, int n,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) = NULL;

void (*_INSERT_TASK_dtrtri)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_diag_t diag, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             int iinfo ) = NULL;

void (*_INSERT_TASK_dtsmlq_hetra1)( const RUNTIME_option_t *options,
                                    cham_side_t side, cham_trans_t trans,
                                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                    const CHAM_desc_t *A1, int A1m, int A1n,
                                    const CHAM_desc_t *A2, int A2m, int A2n,
                                    const CHAM_desc_t *V, int Vm, int Vn,
                                    const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_dtsmqr_hetra1)( const RUNTIME_option_t *options,
                                    cham_side_t side, cham_trans_t trans,
                                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                    const CHAM_desc_t *A1, int A1m, int A1n,
                                    const CHAM_desc_t *A2, int A2m, int A2n,
                                    const CHAM_desc_t *V, int Vm, int Vn,
                                    const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_dtstrf)( const RUNTIME_option_t *options,
                             int m, int n, int ib, int nb,
                             const CHAM_desc_t *U, int Um, int Un,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *L, int Lm, int Ln,
                             int *IPIV,
                             cham_bool_t check_info, int iinfo ) = NULL;

void (*_INSERT_TASK_dormlq)( const RUNTIME_option_t *options,
                             cham_side_t side, cham_trans_t trans,
                             int m, int n, int ib,  int nb, int k,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn,
                             const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_dormqr)( const RUNTIME_option_t *options,
                             cham_side_t side, cham_trans_t trans,
                             int m, int n, int k, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn,
                             const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_dgram)( const RUNTIME_option_t *options,
                            cham_uplo_t uplo,
                            int m, int n, int mt, int nt,
                            const CHAM_desc_t *Di, int Dim, int Din,
                            const CHAM_desc_t *Dj, int Djm, int Djn,
                            const CHAM_desc_t *D, int Dm, int Dn,
                            CHAM_desc_t *A, int Am, int An) = NULL;

void (*_INSERT_TASK_sasum)( const RUNTIME_option_t *options,
                            cham_store_t storev, cham_uplo_t uplo, int M, int N,
                            const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_saxpy)( const RUNTIME_option_t *options,
                            int M, float alpha,
                            const CHAM_desc_t *A, int Am, int An, int incA,
                            const CHAM_desc_t *B, int Bm, int Bn, int incB ) = NULL;

void (*_INSERT_TASK_sbuild)( const RUNTIME_option_t *options,
                             const CHAM_desc_t *A, int Am, int An,
                             void *user_data, void* user_build_callback ) = NULL;

void (*_INSERT_TASK_sgeadd)( const RUNTIME_option_t *options,
                             cham_trans_t trans, int m, int n, int nb,
                             float alpha, const CHAM_desc_t *A, int Am, int An,
                             float beta,  const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_sgelqt)( const RUNTIME_option_t *options,
                             int m, int n, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_sgemv) ( const RUNTIME_option_t *options,
                             cham_trans_t trans, int m, int n,
                             float alpha, const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *X, int Xm, int Xn, int incX,
                             float beta,  const CHAM_desc_t *Y, int Ym, int Yn, int incY ) = NULL;

void (*_INSERT_TASK_sgemm)( const RUNTIME_option_t *options,
                            cham_trans_t transA, cham_trans_t transB,
                            int m, int n, int k, int nb,
                            float alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn,
                            float beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_sgeqrt)( const RUNTIME_option_t *options,
                             int m, int n, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_sgessm)( const RUNTIME_option_t *options,
                             int m, int n, int k, int ib, int nb,
                             int *IPIV,
                             const CHAM_desc_t *L, int Lm, int Ln,
                             const CHAM_desc_t *D, int Dm, int Dn,
                             const CHAM_desc_t *A, int Am, int An ) = NULL;

void (*_INSERT_TASK_sgessq)( const RUNTIME_option_t *options,
                             cham_store_t storev, int m, int n,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) = NULL;

void (*_INSERT_TASK_sgetrf)( const RUNTIME_option_t *options,
                             int m, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             int *IPIV,
                             cham_bool_t check_info, int iinfo ) = NULL;

void (*_INSERT_TASK_sgetrf_incpiv)( const RUNTIME_option_t *options,
                                    int m, int n, int ib, int nb,
                                    const CHAM_desc_t *A, int Am, int An,
                                    const CHAM_desc_t *L, int Lm, int Ln,
                                    int *IPIV,
                                    cham_bool_t check_info, int iinfo ) = NULL;

void (*_INSERT_TASK_sgetrf_nopiv)( const RUNTIME_option_t *options,
                                   int m, int n, int ib, int nb,
                                   const CHAM_desc_t *A, int Am, int An, int iinfo ) = NULL;

void (*_INSERT_TASK_she2ge)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo,
                             int m, int n, int mb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_slacpy)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int m, int n, int mb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_slacpyx)( const RUNTIME_option_t *options,
                              cham_uplo_t uplo, int m, int n, int mb,
                              int displA, const CHAM_desc_t *A, int Am, int An,
                              int displB, const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_slange)( const RUNTIME_option_t *options,
                             cham_normtype_t norm, int M, int N, int NB,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_slange_max)( const RUNTIME_option_t *options,
                                 const CHAM_desc_t *A, int Am, int An,
                                 const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_slansy)( const RUNTIME_option_t *options,
                             cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_slantr)( const RUNTIME_option_t *options,
                             cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag,
                             int M, int N, int NB,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_slascal)( const RUNTIME_option_t *options,
                              cham_uplo_t uplo,
                              int m, int n, int nb,
                              float alpha,
                              const CHAM_desc_t *A, int Am, int An ) = NULL;

void (*_INSERT_TASK_slaset)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int n1, int n2,
                             float alpha, float beta,
                             const CHAM_desc_t *tileA, int tileAm, int tileAn ) = NULL;

void (*_INSERT_TASK_slaset2)( const RUNTIME_option_t *options,
                              cham_uplo_t uplo, int n1, int n2, float alpha,
                              const CHAM_desc_t *tileA, int tileAm, int tileAn ) = NULL;

void (*_INSERT_TASK_slatro)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_trans_t trans, int m, int n, int mb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_slauum)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An ) = NULL;

void (*_INSERT_TASK_splgsy)( const RUNTIME_option_t *options,
                             float bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                             int bigM, int m0, int n0, unsigned long long int seed ) = NULL;

void (*_INSERT_TASK_splrnt)( const RUNTIME_option_t *options,
                             int m, int n, const CHAM_desc_t *A, int Am, int An,
                             int bigM, int m0, int n0, unsigned long long int seed ) = NULL;

void (*_INSERT_TASK_splssq)( const RUNTIME_option_t *options,
                             cham_store_t storev, int M, int N,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                             const CHAM_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn ) = NULL;

void (*_INSERT_TASK_splssq2)( const RUNTIME_option_t *options, int N,
                              const CHAM_desc_t *RESULT, int RESULTm, int RESULTn ) = NULL;

void (*_INSERT_TASK_spotrf)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             int iinfo ) = NULL;

void (*_INSERT_TASK_sssssm)( const RUNTIME_option_t *options,
                             int m1, int n1, int m2, int n2, int k, int ib, int nb,
                             const CHAM_desc_t *A1, int A1m, int A1n,
                             const CHAM_desc_t *A2, int A2m, int A2n,
                             const CHAM_desc_t *L1, int L1m, int L1n,
                             const CHAM_desc_t *L2, int L2m, int L2n,
                             const int *IPIV ) = NULL;

void (*_INSERT_TASK_ssymm)( const RUNTIME_option_t *options,
                            cham_side_t side, cham_uplo_t uplo,
                            int m, int n, int nb,
                            float alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn,
                            float beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_ssyr2k)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_trans_t trans,
                             int n, int k, int nb,
                             float alpha, const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn,
                             float beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_ssyrfb)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo,
                             int n, int k, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn,
                             const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_ssyrk)( const RUNTIME_option_t *options,
                            cham_uplo_t uplo, cham_trans_t trans,
                            int n, int k, int nb,
                            float alpha, const CHAM_desc_t *A, int Am, int An,
                            float beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_ssyssq)( const RUNTIME_option_t *options,
                             cham_store_t storev, cham_uplo_t uplo, int n,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) = NULL;

void (*_INSERT_TASK_ssytrf_nopiv)( const RUNTIME_option_t *options,
                                   cham_uplo_t uplo, int n, int nb,
                                   const CHAM_desc_t *A, int Am, int An,
                                   int iinfo ) = NULL;

void (*_INSERT_TASK_stplqt)( const RUNTIME_option_t *options,
                             int m, int n, int l, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn,
                             const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_stpmlqt)( const RUNTIME_option_t *options,
                              cham_side_t side, cham_trans_t trans,
                              int M, int N, int K, int L, int ib, int nb,
                              const CHAM_desc_t *V, int Vm, int Vn,
                              const CHAM_desc_t *T, int Tm, int Tn,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_stpmqrt)( const RUNTIME_option_t *options,
                              cham_side_t side, cham_trans_t trans,
                              int m, int n, int k, int l, int ib, int nb,
                              const CHAM_desc_t *V, int Vm, int Vn,
                              const CHAM_desc_t *T, int Tm, int Tn,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_stpqrt)( const RUNTIME_option_t *options,
                             int m, int n, int l, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn,
                             const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_stradd)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                             float alpha, const CHAM_desc_t *A, int Am, int An,
                             float beta,  const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_strasm)( const RUNTIME_option_t *options,
                             cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_strmm)( const RUNTIME_option_t *options,
                            cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                            int m, int n, int nb,
                            float alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_strsm)( const RUNTIME_option_t *options,
                            cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                            int m, int n, int nb,
                            float alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_strssq)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_diag_t diag,
                             int m, int n,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) = NULL;

void (*_INSERT_TASK_strtri)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_diag_t diag, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             int iinfo ) = NULL;

void (*_INSERT_TASK_stsmlq_hetra1)( const RUNTIME_option_t *options,
                                    cham_side_t side, cham_trans_t trans,
                                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                    const CHAM_desc_t *A1, int A1m, int A1n,
                                    const CHAM_desc_t *A2, int A2m, int A2n,
                                    const CHAM_desc_t *V, int Vm, int Vn,
                                    const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_stsmqr_hetra1)( const RUNTIME_option_t *options,
                                    cham_side_t side, cham_trans_t trans,
                                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                    const CHAM_desc_t *A1, int A1m, int A1n,
                                    const CHAM_desc_t *A2, int A2m, int A2n,
                                    const CHAM_desc_t *V, int Vm, int Vn,
                                    const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_ststrf)( const RUNTIME_option_t *options,
                             int m, int n, int ib, int nb,
                             const CHAM_desc_t *U, int Um, int Un,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *L, int Lm, int Ln,
                             int *IPIV,
                             cham_bool_t check_info, int iinfo ) = NULL;

void (*_INSERT_TASK_sormlq)( const RUNTIME_option_t *options,
                             cham_side_t side, cham_trans_t trans,
                             int m, int n, int ib,  int nb, int k,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn,
                             const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_sormqr)( const RUNTIME_option_t *options,
                             cham_side_t side, cham_trans_t trans,
                             int m, int n, int k, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn,
                             const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_sgram)( const RUNTIME_option_t *options,
                            cham_uplo_t uplo,
                            int m, int n, int mt, int nt,
                            const CHAM_desc_t *Di, int Dim, int Din,
                            const CHAM_desc_t *Dj, int Djm, int Djn,
                            const CHAM_desc_t *D, int Dm, int Dn,
                            CHAM_desc_t *A, int Am, int An) = NULL;

void (*_INSERT_TASK_dzasum)( const RUNTIME_option_t *options,
                             cham_store_t storev, cham_uplo_t uplo, int M, int N,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_zaxpy)( const RUNTIME_option_t *options,
                            int M, CHAMELEON_Complex64_t alpha,
                            const CHAM_desc_t *A, int Am, int An, int incA,
                            const CHAM_desc_t *B, int Bm, int Bn, int incB ) = NULL;

void (*_INSERT_TASK_zbuild)( const RUNTIME_option_t *options,
                             const CHAM_desc_t *A, int Am, int An,
                             void *user_data, void* user_build_callback ) = NULL;

void (*_INSERT_TASK_zgeadd)( const RUNTIME_option_t *options,
                             cham_trans_t trans, int m, int n, int nb,
                             CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                             CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_zgelqt)( const RUNTIME_option_t *options,
                             int m, int n, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_zgemv) ( const RUNTIME_option_t *options,
                             cham_trans_t trans, int m, int n,
                             CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *X, int Xm, int Xn, int incX,
                             CHAMELEON_Complex64_t beta,  const CHAM_desc_t *Y, int Ym, int Yn, int incY ) = NULL;

void (*_INSERT_TASK_zgemm)( const RUNTIME_option_t *options,
                            cham_trans_t transA, cham_trans_t transB,
                            int m, int n, int k, int nb,
                            CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn,
                            CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_zgeqrt)( const RUNTIME_option_t *options,
                             int m, int n, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_zgessm)( const RUNTIME_option_t *options,
                             int m, int n, int k, int ib, int nb,
                             int *IPIV,
                             const CHAM_desc_t *L, int Lm, int Ln,
                             const CHAM_desc_t *D, int Dm, int Dn,
                             const CHAM_desc_t *A, int Am, int An ) = NULL;

void (*_INSERT_TASK_zgessq)( const RUNTIME_option_t *options,
                             cham_store_t storev, int m, int n,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) = NULL;

void (*_INSERT_TASK_zgetrf)( const RUNTIME_option_t *options,
                             int m, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             int *IPIV,
                             cham_bool_t check_info, int iinfo ) = NULL;

void (*_INSERT_TASK_zgetrf_incpiv)( const RUNTIME_option_t *options,
                                    int m, int n, int ib, int nb,
                                    const CHAM_desc_t *A, int Am, int An,
                                    const CHAM_desc_t *L, int Lm, int Ln,
                                    int *IPIV,
                                    cham_bool_t check_info, int iinfo ) = NULL;

void (*_INSERT_TASK_zgetrf_nopiv)( const RUNTIME_option_t *options,
                                   int m, int n, int ib, int nb,
                                   const CHAM_desc_t *A, int Am, int An, int iinfo ) = NULL;

void (*_INSERT_TASK_zhe2ge)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo,
                             int m, int n, int mb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_zhemm)( const RUNTIME_option_t *options,
                            cham_side_t side, cham_uplo_t uplo,
                            int m, int n, int nb,
                            CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn,
                            CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_zher2k)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_trans_t trans,
                             int n, int k, int nb,
                             CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn,
                             double beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_zherfb)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo,
                             int n, int k, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn,
                             const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_zherk)( const RUNTIME_option_t *options,
                            cham_uplo_t uplo, cham_trans_t trans,
                            int n, int k, int nb,
                            double alpha, const CHAM_desc_t *A, int Am, int An,
                            double beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_zhessq)( const RUNTIME_option_t *options,
                             cham_store_t storev, cham_uplo_t uplo, int n,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) = NULL;

void (*_INSERT_TASK_zlacpy)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int m, int n, int mb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_zlacpyx)( const RUNTIME_option_t *options,
                              cham_uplo_t uplo, int m, int n, int mb,
                              int displA, const CHAM_desc_t *A, int Am, int An,
                              int displB, const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_zlange)( const RUNTIME_option_t *options,
                             cham_normtype_t norm, int M, int N, int NB,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_zlange_max)( const RUNTIME_option_t *options,
                                 const CHAM_desc_t *A, int Am, int An,
                                 const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_zlanhe)( const RUNTIME_option_t *options,
                             cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_zlansy)( const RUNTIME_option_t *options,
                             cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_zlantr)( const RUNTIME_option_t *options,
                             cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag,
                             int M, int N, int NB,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_zlascal)( const RUNTIME_option_t *options,
                              cham_uplo_t uplo,
                              int m, int n, int nb,
                              CHAMELEON_Complex64_t alpha,
                              const CHAM_desc_t *A, int Am, int An ) = NULL;

void (*_INSERT_TASK_zlaset)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int n1, int n2,
                             CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta,
                             const CHAM_desc_t *tileA, int tileAm, int tileAn ) = NULL;

void (*_INSERT_TASK_zlaset2)( const RUNTIME_option_t *options,
                              cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha,
                              const CHAM_desc_t *tileA, int tileAm, int tileAn ) = NULL;

void (*_INSERT_TASK_zlatro)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_trans_t trans, int m, int n, int mb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_zlauum)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An ) = NULL;
void (*_INSERT_TASK_zplghe)( const RUNTIME_option_t *options,
                             double bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                             int bigM, int m0, int n0, unsigned long long int seed ) = NULL;

void (*_INSERT_TASK_zplgsy)( const RUNTIME_option_t *options,
                             CHAMELEON_Complex64_t bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                             int bigM, int m0, int n0, unsigned long long int seed ) = NULL;

void (*_INSERT_TASK_zplrnt)( const RUNTIME_option_t *options,
                             int m, int n, const CHAM_desc_t *A, int Am, int An,
                             int bigM, int m0, int n0, unsigned long long int seed ) = NULL;

void (*_INSERT_TASK_zplssq)( const RUNTIME_option_t *options,
                             cham_store_t storev, int M, int N,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                             const CHAM_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn ) = NULL;

void (*_INSERT_TASK_zplssq2)( const RUNTIME_option_t *options, int N,
                              const CHAM_desc_t *RESULT, int RESULTm, int RESULTn ) = NULL;

void (*_INSERT_TASK_zpotrf)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             int iinfo ) = NULL;

void (*_INSERT_TASK_zssssm)( const RUNTIME_option_t *options,
                             int m1, int n1, int m2, int n2, int k, int ib, int nb,
                             const CHAM_desc_t *A1, int A1m, int A1n,
                             const CHAM_desc_t *A2, int A2m, int A2n,
                             const CHAM_desc_t *L1, int L1m, int L1n,
                             const CHAM_desc_t *L2, int L2m, int L2n,
                             const int *IPIV ) = NULL;

void (*_INSERT_TASK_zsymm)( const RUNTIME_option_t *options,
                            cham_side_t side, cham_uplo_t uplo,
                            int m, int n, int nb,
                            CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn,
                            CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_zsyr2k)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_trans_t trans,
                             int n, int k, int nb,
                             CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn,
                             CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_zsyrk)( const RUNTIME_option_t *options,
                            cham_uplo_t uplo, cham_trans_t trans,
                            int n, int k, int nb,
                            CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                            CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_zsyssq)( const RUNTIME_option_t *options,
                             cham_store_t storev, cham_uplo_t uplo, int n,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) = NULL;

void (*_INSERT_TASK_zsytrf_nopiv)( const RUNTIME_option_t *options,
                                   cham_uplo_t uplo, int n, int nb,
                                   const CHAM_desc_t *A, int Am, int An,
                                   int iinfo ) = NULL;

void (*_INSERT_TASK_ztplqt)( const RUNTIME_option_t *options,
                             int m, int n, int l, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn,
                             const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_ztpmlqt)( const RUNTIME_option_t *options,
                              cham_side_t side, cham_trans_t trans,
                              int M, int N, int K, int L, int ib, int nb,
                              const CHAM_desc_t *V, int Vm, int Vn,
                              const CHAM_desc_t *T, int Tm, int Tn,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_ztpmqrt)( const RUNTIME_option_t *options,
                              cham_side_t side, cham_trans_t trans,
                              int m, int n, int k, int l, int ib, int nb,
                              const CHAM_desc_t *V, int Vm, int Vn,
                              const CHAM_desc_t *T, int Tm, int Tn,
                              const CHAM_desc_t *A, int Am, int An,
                              const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_ztpqrt)( const RUNTIME_option_t *options,
                             int m, int n, int l, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn,
                             const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_ztradd)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                             CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                             CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_ztrasm)( const RUNTIME_option_t *options,
                             cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_ztrmm)( const RUNTIME_option_t *options,
                            cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                            int m, int n, int nb,
                            CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_ztrsm)( const RUNTIME_option_t *options,
                            cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                            int m, int n, int nb,
                            CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                            const CHAM_desc_t *B, int Bm, int Bn ) = NULL;

void (*_INSERT_TASK_ztrssq)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_diag_t diag,
                             int m, int n,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) = NULL;

void (*_INSERT_TASK_ztrtri)( const RUNTIME_option_t *options,
                             cham_uplo_t uplo, cham_diag_t diag, int n, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             int iinfo ) = NULL;

void (*_INSERT_TASK_ztsmlq_hetra1)( const RUNTIME_option_t *options,
                                    cham_side_t side, cham_trans_t trans,
                                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                    const CHAM_desc_t *A1, int A1m, int A1n,
                                    const CHAM_desc_t *A2, int A2m, int A2n,
                                    const CHAM_desc_t *V, int Vm, int Vn,
                                    const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_ztsmqr_hetra1)( const RUNTIME_option_t *options,
                                    cham_side_t side, cham_trans_t trans,
                                    int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                    const CHAM_desc_t *A1, int A1m, int A1n,
                                    const CHAM_desc_t *A2, int A2m, int A2n,
                                    const CHAM_desc_t *V, int Vm, int Vn,
                                    const CHAM_desc_t *T, int Tm, int Tn ) = NULL;

void (*_INSERT_TASK_ztstrf)( const RUNTIME_option_t *options,
                             int m, int n, int ib, int nb,
                             const CHAM_desc_t *U, int Um, int Un,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *L, int Lm, int Ln,
                             int *IPIV,
                             cham_bool_t check_info, int iinfo ) = NULL;

void (*_INSERT_TASK_zunmlq)( const RUNTIME_option_t *options,
                             cham_side_t side, cham_trans_t trans,
                             int m, int n, int ib,  int nb, int k,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn,
                             const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_zunmqr)( const RUNTIME_option_t *options,
                             cham_side_t side, cham_trans_t trans,
                             int m, int n, int k, int ib, int nb,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *T, int Tm, int Tn,
                             const CHAM_desc_t *C, int Cm, int Cn ) = NULL;

void (*_INSERT_TASK_zgram)( const RUNTIME_option_t *options,
                            cham_uplo_t uplo,
                            int m, int n, int mt, int nt,
                            const CHAM_desc_t *Di, int Dim, int Din,
                            const CHAM_desc_t *Dj, int Djm, int Djn,
                            const CHAM_desc_t *D, int Dm, int Dn,
                            CHAM_desc_t *A, int Am, int An) = NULL;

void (*_INSERT_TASK_slag2d)( const RUNTIME_option_t *options,
                             int m, int n, int nb,
                             CHAM_desc_t *A, int Am, int An, int lda,
                             CHAM_desc_t *B, int Bm, int Bn, int ldb ) = NULL;

void (*_INSERT_TASK_dlag2s)( const RUNTIME_option_t *options,
                             int m, int n, int nb,
                             CHAM_desc_t *A, int Am, int An, int lda,
                             CHAM_desc_t *B, int Bm, int Bn, int ldb ) = NULL;

void (*_INSERT_TASK_clag2z)( const RUNTIME_option_t *options,
                             int m, int n, int nb,
                             CHAM_desc_t *A, int Am, int An, int lda,
                             CHAM_desc_t *B, int Bm, int Bn, int ldb ) = NULL;

void (*_INSERT_TASK_zlag2c)( const RUNTIME_option_t *options,
                             int m, int n, int nb,
                             CHAM_desc_t *A, int Am, int An, int lda,
                             CHAM_desc_t *B, int Bm, int Bn, int ldb ) = NULL;
/**
 * @brief Macro to call dlsym for every runtime symbols
 */
#define RUNTIME_DLSYM( name )                                           \
    *(void**)(&_##name) = dlsym( _runtime_handle, #name );              \
    if ( (error = dlerror()) != NULL ){                                 \
        chameleon_fatal_error("runtime_fake", "dlsym has failed");      \
    }


/* dlopen handle */
void *_runtime_handle = NULL;

int RUNTIME_init(CHAM_context_t *ctxt,
                 int ncpus,
                 int ncudas,
                 int nthreads_per_worker) {
    return _RUNTIME_init( ctxt, ncpus, ncudas, nthreads_per_worker );
}

void RUNTIME_finalize(CHAM_context_t *ctxt)
{
    return _RUNTIME_finalize( ctxt );
}

void RUNTIME_context_create( CHAM_context_t *ctxt ) {
    char *error;
    char *cham_runtime_lib_s;

    /* chameleon runtime library must be provided */
    cham_runtime_lib_s = getenv("CHAMELEON_RUNTIME_LIBRARY");
    if ( cham_runtime_lib_s == NULL ){
        chameleon_fatal_error("runtime_fake", "CHAMELEON_RUNTIME_LIBRARY environment variable is not set.\n"
                              "Please choose a runtime system to use by giving the absolute path to a runtime chameleon library,"
                              " e.g. export CHAMELEON_RUNTIME_LIBRARY=$CHAMELEON_DIR/lib/chameleon_starpu.so.");
    }

    /* dlopen the given library */
    _runtime_handle = dlopen (cham_runtime_lib_s, RTLD_NOW);
    if ( !_runtime_handle ){
        fprintf(stderr, "%s\n", dlerror());
        fprintf(stderr, "CHAMELEON ERROR: %s dlopen fails\n", cham_runtime_lib_s);
        chameleon_fatal_error("runtime_fake", "dlopen has failed");
    }
    dlerror();    /* Clear any existing error */

    /* load symbols with dlsym */
    RUNTIME_DLSYM( INSERT_TASK_slag2c            )
    RUNTIME_DLSYM( INSERT_TASK_scasum            )
    RUNTIME_DLSYM( INSERT_TASK_caxpy             )
    RUNTIME_DLSYM( INSERT_TASK_cbuild            )
    RUNTIME_DLSYM( INSERT_TASK_cgeadd            )
    RUNTIME_DLSYM( INSERT_TASK_cgelqt            )
    RUNTIME_DLSYM( INSERT_TASK_cgemm             )
    RUNTIME_DLSYM( INSERT_TASK_cgemv             )
    RUNTIME_DLSYM( INSERT_TASK_cgeqrt            )
    RUNTIME_DLSYM( INSERT_TASK_cgessm            )
    RUNTIME_DLSYM( INSERT_TASK_cgessq            )
    RUNTIME_DLSYM( INSERT_TASK_cgetrf            )
    RUNTIME_DLSYM( INSERT_TASK_cgetrf_incpiv     )
    RUNTIME_DLSYM( INSERT_TASK_cgetrf_nopiv      )
    RUNTIME_DLSYM( INSERT_TASK_cgram             )
    RUNTIME_DLSYM( INSERT_TASK_che2ge            )
    RUNTIME_DLSYM( INSERT_TASK_chemm             )
    RUNTIME_DLSYM( INSERT_TASK_cher2k            )
    RUNTIME_DLSYM( INSERT_TASK_cherfb            )
    RUNTIME_DLSYM( INSERT_TASK_cherk             )
    RUNTIME_DLSYM( INSERT_TASK_chessq            )
    RUNTIME_DLSYM( INSERT_TASK_clacpy            )
    RUNTIME_DLSYM( INSERT_TASK_clacpyx           )
    RUNTIME_DLSYM( INSERT_TASK_clange            )
    RUNTIME_DLSYM( INSERT_TASK_clange_max        )
    RUNTIME_DLSYM( INSERT_TASK_clanhe            )
    RUNTIME_DLSYM( INSERT_TASK_clansy            )
    RUNTIME_DLSYM( INSERT_TASK_clantr            )
    RUNTIME_DLSYM( INSERT_TASK_clascal           )
    RUNTIME_DLSYM( INSERT_TASK_claset            )
    RUNTIME_DLSYM( INSERT_TASK_claset2           )
    RUNTIME_DLSYM( INSERT_TASK_clatro            )
    RUNTIME_DLSYM( INSERT_TASK_clauum            )
    RUNTIME_DLSYM( INSERT_TASK_cplghe            )
    RUNTIME_DLSYM( INSERT_TASK_cplgsy            )
    RUNTIME_DLSYM( INSERT_TASK_cplrnt            )
    RUNTIME_DLSYM( INSERT_TASK_cplssq            )
    RUNTIME_DLSYM( INSERT_TASK_cplssq2           )
    RUNTIME_DLSYM( INSERT_TASK_cpotrf            )
    RUNTIME_DLSYM( INSERT_TASK_cssssm            )
    RUNTIME_DLSYM( INSERT_TASK_csymm             )
    RUNTIME_DLSYM( INSERT_TASK_csyr2k            )
    RUNTIME_DLSYM( INSERT_TASK_csyrk             )
    RUNTIME_DLSYM( INSERT_TASK_csyssq            )
    RUNTIME_DLSYM( INSERT_TASK_csytrf_nopiv      )
    RUNTIME_DLSYM( INSERT_TASK_ctplqt            )
    RUNTIME_DLSYM( INSERT_TASK_ctpmlqt           )
    RUNTIME_DLSYM( INSERT_TASK_ctpmqrt           )
    RUNTIME_DLSYM( INSERT_TASK_ctpqrt            )
    RUNTIME_DLSYM( INSERT_TASK_ctradd            )
    RUNTIME_DLSYM( INSERT_TASK_ctrasm            )
    RUNTIME_DLSYM( INSERT_TASK_ctrmm             )
    RUNTIME_DLSYM( INSERT_TASK_ctrsm             )
    RUNTIME_DLSYM( INSERT_TASK_ctrssq            )
    RUNTIME_DLSYM( INSERT_TASK_ctrtri            )
    RUNTIME_DLSYM( INSERT_TASK_ctsmlq_hetra1     )
    RUNTIME_DLSYM( INSERT_TASK_ctsmqr_hetra1     )
    RUNTIME_DLSYM( INSERT_TASK_ctstrf            )
    RUNTIME_DLSYM( INSERT_TASK_cunmlq            )
    RUNTIME_DLSYM( INSERT_TASK_cunmqr            )
    RUNTIME_DLSYM( INSERT_TASK_dlag2z            )
    RUNTIME_DLSYM( INSERT_TASK_dasum             )
    RUNTIME_DLSYM( INSERT_TASK_daxpy             )
    RUNTIME_DLSYM( INSERT_TASK_dbuild            )
    RUNTIME_DLSYM( INSERT_TASK_dgeadd            )
    RUNTIME_DLSYM( INSERT_TASK_dgelqt            )
    RUNTIME_DLSYM( INSERT_TASK_dgemm             )
    RUNTIME_DLSYM( INSERT_TASK_dgemv             )
    RUNTIME_DLSYM( INSERT_TASK_dgeqrt            )
    RUNTIME_DLSYM( INSERT_TASK_dgessm            )
    RUNTIME_DLSYM( INSERT_TASK_dgessq            )
    RUNTIME_DLSYM( INSERT_TASK_dgetrf            )
    RUNTIME_DLSYM( INSERT_TASK_dgetrf_incpiv     )
    RUNTIME_DLSYM( INSERT_TASK_dgetrf_nopiv      )
    RUNTIME_DLSYM( INSERT_TASK_dgram             )
    RUNTIME_DLSYM( INSERT_TASK_dhe2ge            )
    RUNTIME_DLSYM( INSERT_TASK_dlacpy            )
    RUNTIME_DLSYM( INSERT_TASK_dlacpyx           )
    RUNTIME_DLSYM( INSERT_TASK_dlange            )
    RUNTIME_DLSYM( INSERT_TASK_dlange_max        )
    RUNTIME_DLSYM( INSERT_TASK_dlansy            )
    RUNTIME_DLSYM( INSERT_TASK_dlantr            )
    RUNTIME_DLSYM( INSERT_TASK_dlascal           )
    RUNTIME_DLSYM( INSERT_TASK_dlaset            )
    RUNTIME_DLSYM( INSERT_TASK_dlaset2           )
    RUNTIME_DLSYM( INSERT_TASK_dlatro            )
    RUNTIME_DLSYM( INSERT_TASK_dlauum            )
    RUNTIME_DLSYM( INSERT_TASK_dormlq            )
    RUNTIME_DLSYM( INSERT_TASK_dormqr            )
    RUNTIME_DLSYM( INSERT_TASK_dplgsy            )
    RUNTIME_DLSYM( INSERT_TASK_dplrnt            )
    RUNTIME_DLSYM( INSERT_TASK_dplssq            )
    RUNTIME_DLSYM( INSERT_TASK_dplssq2           )
    RUNTIME_DLSYM( INSERT_TASK_dpotrf            )
    RUNTIME_DLSYM( INSERT_TASK_dssssm            )
    RUNTIME_DLSYM( INSERT_TASK_dsymm             )
    RUNTIME_DLSYM( INSERT_TASK_dsyr2k            )
    RUNTIME_DLSYM( INSERT_TASK_dsyrfb            )
    RUNTIME_DLSYM( INSERT_TASK_dsyrk             )
    RUNTIME_DLSYM( INSERT_TASK_dsyssq            )
    RUNTIME_DLSYM( INSERT_TASK_dtplqt            )
    RUNTIME_DLSYM( INSERT_TASK_dtpmlqt           )
    RUNTIME_DLSYM( INSERT_TASK_dtpmqrt           )
    RUNTIME_DLSYM( INSERT_TASK_dtpqrt            )
    RUNTIME_DLSYM( INSERT_TASK_dtradd            )
    RUNTIME_DLSYM( INSERT_TASK_dtrasm            )
    RUNTIME_DLSYM( INSERT_TASK_dtrmm             )
    RUNTIME_DLSYM( INSERT_TASK_dtrsm             )
    RUNTIME_DLSYM( INSERT_TASK_dtrssq            )
    RUNTIME_DLSYM( INSERT_TASK_dtrtri            )
    RUNTIME_DLSYM( INSERT_TASK_dtsmlq_hetra1     )
    RUNTIME_DLSYM( INSERT_TASK_dtsmqr_hetra1     )
    RUNTIME_DLSYM( INSERT_TASK_dtstrf            )
    RUNTIME_DLSYM( INSERT_TASK_dzasum            )
    RUNTIME_DLSYM( INSERT_TASK_slag2c            )
    RUNTIME_DLSYM( INSERT_TASK_sasum             )
    RUNTIME_DLSYM( INSERT_TASK_saxpy             )
    RUNTIME_DLSYM( INSERT_TASK_sbuild            )
    RUNTIME_DLSYM( INSERT_TASK_sgeadd            )
    RUNTIME_DLSYM( INSERT_TASK_sgelqt            )
    RUNTIME_DLSYM( INSERT_TASK_sgemm             )
    RUNTIME_DLSYM( INSERT_TASK_sgemv             )
    RUNTIME_DLSYM( INSERT_TASK_sgeqrt            )
    RUNTIME_DLSYM( INSERT_TASK_sgessm            )
    RUNTIME_DLSYM( INSERT_TASK_sgessq            )
    RUNTIME_DLSYM( INSERT_TASK_sgetrf            )
    RUNTIME_DLSYM( INSERT_TASK_sgetrf_incpiv     )
    RUNTIME_DLSYM( INSERT_TASK_sgetrf_nopiv      )
    RUNTIME_DLSYM( INSERT_TASK_sgram             )
    RUNTIME_DLSYM( INSERT_TASK_she2ge            )
    RUNTIME_DLSYM( INSERT_TASK_slacpy            )
    RUNTIME_DLSYM( INSERT_TASK_slacpyx           )
    RUNTIME_DLSYM( INSERT_TASK_slange            )
    RUNTIME_DLSYM( INSERT_TASK_slange_max        )
    RUNTIME_DLSYM( INSERT_TASK_slansy            )
    RUNTIME_DLSYM( INSERT_TASK_slantr            )
    RUNTIME_DLSYM( INSERT_TASK_slascal           )
    RUNTIME_DLSYM( INSERT_TASK_slaset            )
    RUNTIME_DLSYM( INSERT_TASK_slaset2           )
    RUNTIME_DLSYM( INSERT_TASK_slatro            )
    RUNTIME_DLSYM( INSERT_TASK_slauum            )
    RUNTIME_DLSYM( INSERT_TASK_sormlq            )
    RUNTIME_DLSYM( INSERT_TASK_sormqr            )
    RUNTIME_DLSYM( INSERT_TASK_splgsy            )
    RUNTIME_DLSYM( INSERT_TASK_splrnt            )
    RUNTIME_DLSYM( INSERT_TASK_splssq            )
    RUNTIME_DLSYM( INSERT_TASK_splssq2           )
    RUNTIME_DLSYM( INSERT_TASK_spotrf            )
    RUNTIME_DLSYM( INSERT_TASK_sssssm            )
    RUNTIME_DLSYM( INSERT_TASK_ssymm             )
    RUNTIME_DLSYM( INSERT_TASK_ssyr2k            )
    RUNTIME_DLSYM( INSERT_TASK_ssyrfb            )
    RUNTIME_DLSYM( INSERT_TASK_ssyrk             )
    RUNTIME_DLSYM( INSERT_TASK_ssyssq            )
    RUNTIME_DLSYM( INSERT_TASK_stplqt            )
    RUNTIME_DLSYM( INSERT_TASK_stpmlqt           )
    RUNTIME_DLSYM( INSERT_TASK_stpmqrt           )
    RUNTIME_DLSYM( INSERT_TASK_stpqrt            )
    RUNTIME_DLSYM( INSERT_TASK_stradd            )
    RUNTIME_DLSYM( INSERT_TASK_strasm            )
    RUNTIME_DLSYM( INSERT_TASK_strmm             )
    RUNTIME_DLSYM( INSERT_TASK_strsm             )
    RUNTIME_DLSYM( INSERT_TASK_strssq            )
    RUNTIME_DLSYM( INSERT_TASK_strtri            )
    RUNTIME_DLSYM( INSERT_TASK_stsmlq_hetra1     )
    RUNTIME_DLSYM( INSERT_TASK_stsmqr_hetra1     )
    RUNTIME_DLSYM( INSERT_TASK_ststrf            )
    RUNTIME_DLSYM( INSERT_TASK_dzasum            )
    RUNTIME_DLSYM( INSERT_TASK_zaxpy             )
    RUNTIME_DLSYM( INSERT_TASK_zbuild            )
    RUNTIME_DLSYM( INSERT_TASK_zgeadd            )
    RUNTIME_DLSYM( INSERT_TASK_zgelqt            )
    RUNTIME_DLSYM( INSERT_TASK_zgemm             )
    RUNTIME_DLSYM( INSERT_TASK_zgemv             )
    RUNTIME_DLSYM( INSERT_TASK_zgeqrt            )
    RUNTIME_DLSYM( INSERT_TASK_zgessm            )
    RUNTIME_DLSYM( INSERT_TASK_zgessq            )
    RUNTIME_DLSYM( INSERT_TASK_zgetrf            )
    RUNTIME_DLSYM( INSERT_TASK_zgetrf_incpiv     )
    RUNTIME_DLSYM( INSERT_TASK_zgetrf_nopiv      )
    RUNTIME_DLSYM( INSERT_TASK_zgram             )
    RUNTIME_DLSYM( INSERT_TASK_zhe2ge            )
    RUNTIME_DLSYM( INSERT_TASK_zhemm             )
    RUNTIME_DLSYM( INSERT_TASK_zher2k            )
    RUNTIME_DLSYM( INSERT_TASK_zherfb            )
    RUNTIME_DLSYM( INSERT_TASK_zherk             )
    RUNTIME_DLSYM( INSERT_TASK_zhessq            )
    RUNTIME_DLSYM( INSERT_TASK_zlacpy            )
    RUNTIME_DLSYM( INSERT_TASK_zlacpyx           )
    RUNTIME_DLSYM( INSERT_TASK_zlange            )
    RUNTIME_DLSYM( INSERT_TASK_zlange_max        )
    RUNTIME_DLSYM( INSERT_TASK_zlanhe            )
    RUNTIME_DLSYM( INSERT_TASK_zlansy            )
    RUNTIME_DLSYM( INSERT_TASK_zlantr            )
    RUNTIME_DLSYM( INSERT_TASK_zlascal           )
    RUNTIME_DLSYM( INSERT_TASK_zlaset            )
    RUNTIME_DLSYM( INSERT_TASK_zlaset2           )
    RUNTIME_DLSYM( INSERT_TASK_zlatro            )
    RUNTIME_DLSYM( INSERT_TASK_zlauum            )
    RUNTIME_DLSYM( INSERT_TASK_zplghe            )
    RUNTIME_DLSYM( INSERT_TASK_zplgsy            )
    RUNTIME_DLSYM( INSERT_TASK_zplrnt            )
    RUNTIME_DLSYM( INSERT_TASK_zplssq            )
    RUNTIME_DLSYM( INSERT_TASK_zplssq2           )
    RUNTIME_DLSYM( INSERT_TASK_zpotrf            )
    RUNTIME_DLSYM( INSERT_TASK_zssssm            )
    RUNTIME_DLSYM( INSERT_TASK_zsymm             )
    RUNTIME_DLSYM( INSERT_TASK_zsyr2k            )
    RUNTIME_DLSYM( INSERT_TASK_zsyrk             )
    RUNTIME_DLSYM( INSERT_TASK_zsyssq            )
    RUNTIME_DLSYM( INSERT_TASK_zsytrf_nopiv      )
    RUNTIME_DLSYM( INSERT_TASK_ztplqt            )
    RUNTIME_DLSYM( INSERT_TASK_ztpmlqt           )
    RUNTIME_DLSYM( INSERT_TASK_ztpmqrt           )
    RUNTIME_DLSYM( INSERT_TASK_ztpqrt            )
    RUNTIME_DLSYM( INSERT_TASK_ztradd            )
    RUNTIME_DLSYM( INSERT_TASK_ztrasm            )
    RUNTIME_DLSYM( INSERT_TASK_ztrmm             )
    RUNTIME_DLSYM( INSERT_TASK_ztrsm             )
    RUNTIME_DLSYM( INSERT_TASK_ztrssq            )
    RUNTIME_DLSYM( INSERT_TASK_ztrtri            )
    RUNTIME_DLSYM( INSERT_TASK_ztsmlq_hetra1     )
    RUNTIME_DLSYM( INSERT_TASK_ztsmqr_hetra1     )
    RUNTIME_DLSYM( INSERT_TASK_ztstrf            )
    RUNTIME_DLSYM( INSERT_TASK_zunmlq            )
    RUNTIME_DLSYM( INSERT_TASK_zunmqr            )
    // RUNTIME_DLSYM( INSERT_TASK_slag2d            )
    // RUNTIME_DLSYM( INSERT_TASK_dlag2s            )
    // RUNTIME_DLSYM( INSERT_TASK_clag2z            )
    // RUNTIME_DLSYM( INSERT_TASK_zlag2c            )
    RUNTIME_DLSYM( INSERT_TASK_map               )
    RUNTIME_DLSYM( RUNTIME_barrier               )
    RUNTIME_DLSYM( RUNTIME_comm_rank             )
    RUNTIME_DLSYM( RUNTIME_comm_set_tag_sizes    )
    RUNTIME_DLSYM( RUNTIME_comm_size             )
    RUNTIME_DLSYM( RUNTIME_context_create        )
    RUNTIME_DLSYM( RUNTIME_context_destroy       )
    RUNTIME_DLSYM( RUNTIME_data_flush            )
#if defined(CHAMELEON_USE_MIGRATE)
    RUNTIME_DLSYM( RUNTIME_data_migrate          )
#endif
    RUNTIME_DLSYM( RUNTIME_desc_acquire          )
    RUNTIME_DLSYM( RUNTIME_desc_create           )
    RUNTIME_DLSYM( RUNTIME_desc_destroy          )
    RUNTIME_DLSYM( RUNTIME_desc_flush            )
    RUNTIME_DLSYM( RUNTIME_desc_release          )
    RUNTIME_DLSYM( RUNTIME_disable               )
    RUNTIME_DLSYM( RUNTIME_enable                )
    RUNTIME_DLSYM( RUNTIME_finalize              )
    RUNTIME_DLSYM( RUNTIME_flush                 )
    RUNTIME_DLSYM( RUNTIME_free                  )
    RUNTIME_DLSYM( RUNTIME_get_time              )
    RUNTIME_DLSYM( RUNTIME_init                  )
    RUNTIME_DLSYM( RUNTIME_iteration_pop         )
    RUNTIME_DLSYM( RUNTIME_iteration_push        )
    RUNTIME_DLSYM( RUNTIME_kernelprofile_display )
    RUNTIME_DLSYM( RUNTIME_malloc                )
    RUNTIME_DLSYM( RUNTIME_options_finalize      )
    RUNTIME_DLSYM( RUNTIME_options_init          )
    RUNTIME_DLSYM( RUNTIME_options_ws_alloc      )
    RUNTIME_DLSYM( RUNTIME_options_ws_free       )
    RUNTIME_DLSYM( RUNTIME_pause                 )
    RUNTIME_DLSYM( RUNTIME_progress              )
    RUNTIME_DLSYM( RUNTIME_request_create        )
    RUNTIME_DLSYM( RUNTIME_request_destroy       )
    RUNTIME_DLSYM( RUNTIME_request_set           )
    RUNTIME_DLSYM( RUNTIME_resume                )
    RUNTIME_DLSYM( RUNTIME_schedprofile_display  )
    RUNTIME_DLSYM( RUNTIME_sequence_create       )
    RUNTIME_DLSYM( RUNTIME_sequence_destroy      )
    RUNTIME_DLSYM( RUNTIME_sequence_flush        )
    RUNTIME_DLSYM( RUNTIME_sequence_wait         )
    RUNTIME_DLSYM( RUNTIME_start_profiling       )
    RUNTIME_DLSYM( RUNTIME_start_stats           )
    RUNTIME_DLSYM( RUNTIME_stop_profiling        )
    RUNTIME_DLSYM( RUNTIME_stop_stats            )
    RUNTIME_DLSYM( RUNTIME_thread_rank           )
    RUNTIME_DLSYM( RUNTIME_thread_size           )
    RUNTIME_DLSYM( RUNTIME_zlocality_allrestrict )

    return _RUNTIME_context_create( ctxt );
}

void RUNTIME_context_destroy( CHAM_context_t *ctxt ) {
    _RUNTIME_context_destroy( ctxt );
    if ( _runtime_handle ) {
        dlclose( _runtime_handle );
        _runtime_handle = NULL;
    }
}

void RUNTIME_enable( void *runtime_ctxt, int option ) {
    return _RUNTIME_enable( runtime_ctxt, option );
}

void RUNTIME_disable( void *runtime_ctxt, int option ) {
    return _RUNTIME_disable( runtime_ctxt, option );
}

void RUNTIME_pause( CHAM_context_t *ctxt ) {
    return _RUNTIME_pause( ctxt );
}

void RUNTIME_resume( CHAM_context_t *ctxt ) {
    return _RUNTIME_resume( ctxt );
}

void RUNTIME_barrier( CHAM_context_t *ctxt ) {
    return _RUNTIME_barrier( ctxt );
}

void RUNTIME_progress( CHAM_context_t *ctxt ) {
    return _RUNTIME_progress( ctxt );
}

int RUNTIME_thread_rank( CHAM_context_t *ctxt ) {
    return _RUNTIME_thread_rank( ctxt );
}

int RUNTIME_thread_size( CHAM_context_t *ctxt ) {
    return _RUNTIME_thread_size( ctxt );
}

int RUNTIME_comm_rank( CHAM_context_t *ctxt ) {
    return _RUNTIME_comm_rank( ctxt );
}

int RUNTIME_comm_size( CHAM_context_t *ctxt ) {
    return _RUNTIME_comm_size( ctxt );
}

void RUNTIME_comm_set_tag_sizes( int user_tag_width, int user_tag_sep ) {
    return _RUNTIME_comm_set_tag_sizes( user_tag_width, user_tag_sep );
}

int RUNTIME_sequence_create( CHAM_context_t *ctxt, RUNTIME_sequence_t *sequence ) {
    return _RUNTIME_sequence_create( ctxt, sequence );
}

int RUNTIME_sequence_destroy( CHAM_context_t *ctxt, RUNTIME_sequence_t *sequence ) {
    return _RUNTIME_sequence_destroy( ctxt, sequence );
}

int RUNTIME_sequence_wait( CHAM_context_t *ctxt, RUNTIME_sequence_t *sequence ) {
    return _RUNTIME_sequence_wait( ctxt, sequence );
}

void RUNTIME_sequence_flush( CHAM_context_t *ctxt,
                             RUNTIME_sequence_t *sequence,
                             RUNTIME_request_t  *request,
                             int                 status ) {
    return _RUNTIME_sequence_flush( ctxt, sequence, request, status );
}

int RUNTIME_request_create( CHAM_context_t *ctxt, RUNTIME_request_t *request ) {
    return _RUNTIME_request_create( ctxt, request );
}

int RUNTIME_request_destroy( CHAM_context_t *ctxt, RUNTIME_request_t *request ) {
    return _RUNTIME_request_destroy( ctxt, request );
}

int RUNTIME_request_set( CHAM_context_t  *ctxt,
                         RUNTIME_request_t *request,
                         int param, int value ) {
    return _RUNTIME_request_set( ctxt, request, param, value );
}

void *RUNTIME_malloc( size_t size ) {
    return _RUNTIME_malloc( size );
}

void RUNTIME_free( void *ptr, size_t size ) {
    return _RUNTIME_free( ptr, size );
}

void RUNTIME_desc_create( CHAM_desc_t *desc ) {
    return _RUNTIME_desc_create( desc );
}

void RUNTIME_desc_destroy( CHAM_desc_t *desc ) {
    return _RUNTIME_desc_destroy( desc );
}

int RUNTIME_desc_acquire( const CHAM_desc_t *desc ) {
    return _RUNTIME_desc_acquire( desc );
}

int RUNTIME_desc_release( const CHAM_desc_t *desc ) {
    return _RUNTIME_desc_release( desc );
}

void RUNTIME_desc_flush( const CHAM_desc_t *desc, const RUNTIME_sequence_t *sequence ) {
    return _RUNTIME_desc_flush( desc, sequence );
}

void RUNTIME_flush( ) {
    return _RUNTIME_flush( );
}

void RUNTIME_data_flush( const RUNTIME_sequence_t *sequence,
                         const CHAM_desc_t *A, int Am, int An ) {
    return _RUNTIME_data_flush( sequence, A, Am, An );
}

#if defined(CHAMELEON_USE_MIGRATE)
void RUNTIME_data_migrate( const RUNTIME_sequence_t *sequence,
                           const CHAM_desc_t *A, int Am, int An, int new_rank ) {
    return _RUNTIME_data_migrate( sequence, A, Am, An, new_rank );
}
#endif

void RUNTIME_options_init( RUNTIME_option_t *options,
                           CHAM_context_t *context,
                           RUNTIME_sequence_t *sequence,
                           RUNTIME_request_t *request ) {
    return _RUNTIME_options_init( options, context, sequence, request );
}

void RUNTIME_options_finalize( RUNTIME_option_t  *options,
                               CHAM_context_t *context ) {
    return _RUNTIME_options_finalize( options, context );
}

int RUNTIME_options_ws_alloc( RUNTIME_option_t *options,
                              size_t wsize, size_t hsize ) {
    return _RUNTIME_options_ws_alloc( options, wsize, hsize );
}

int RUNTIME_options_ws_free( RUNTIME_option_t *options ) {
    return _RUNTIME_options_ws_free( options );
}

void RUNTIME_zlocality_allrestore () {
    return _RUNTIME_zlocality_allrestore();
}
void RUNTIME_clocality_allrestore () {
    return _RUNTIME_clocality_allrestore();
}
void RUNTIME_dlocality_allrestore () {
    return _RUNTIME_dlocality_allrestore();
}
void RUNTIME_slocality_allrestore () {
    return _RUNTIME_slocality_allrestore();
}

void RUNTIME_zlocality_allrestrict(uint32_t device) {
    return _RUNTIME_zlocality_allrestrict( device );
}
void RUNTIME_zlocality_onerestrict(cham_tasktype_t task, uint32_t device) {
    return _RUNTIME_zlocality_onerestrict( task, device );
}
void RUNTIME_zlocality_onerestore(cham_tasktype_t task) {
    return _RUNTIME_zlocality_onerestore( task );
}

void RUNTIME_clocality_allrestrict(uint32_t device) {
    return _RUNTIME_clocality_allrestrict( device );
}
void RUNTIME_clocality_onerestrict(cham_tasktype_t task, uint32_t device) {
    return _RUNTIME_clocality_onerestrict( task, device );
}
void RUNTIME_clocality_onerestore(cham_tasktype_t task) {
    return _RUNTIME_clocality_onerestore( task );
}

void RUNTIME_dlocality_allrestrict(uint32_t device) {
    return _RUNTIME_dlocality_allrestrict( device );
}
void RUNTIME_dlocality_onerestrict(cham_tasktype_t task, uint32_t device) {
    return _RUNTIME_dlocality_onerestrict( task, device );
}
void RUNTIME_dlocality_onerestore (cham_tasktype_t task) {
    return _RUNTIME_dlocality_onerestore( task );
}

void RUNTIME_slocality_allrestrict(uint32_t device) {
    return _RUNTIME_slocality_allrestrict( device );
}
void RUNTIME_slocality_onerestrict(cham_tasktype_t task, uint32_t device) {
    return _RUNTIME_slocality_onerestrict( task, device );
}
void RUNTIME_slocality_onerestore(cham_tasktype_t task) {
    return _RUNTIME_slocality_onerestore( task );
}

void   RUNTIME_schedprofile_display () {
    return _RUNTIME_schedprofile_display();
}
void   RUNTIME_kernelprofile_display() {
    return _RUNTIME_kernelprofile_display();
}
double RUNTIME_get_time() {
    return _RUNTIME_get_time();
}

void  RUNTIME_iteration_push     (CHAM_context_t* ctxt, unsigned long iteration) {
    return _RUNTIME_iteration_push( ctxt, iteration);
}
void  RUNTIME_iteration_pop      (CHAM_context_t* ctxt) {
    return _RUNTIME_iteration_pop( ctxt );
}

void RUNTIME_start_profiling() {
    return _RUNTIME_start_profiling();
}
void RUNTIME_stop_profiling() {
    return _RUNTIME_stop_profiling();
}

void RUNTIME_start_stats() {
    return _RUNTIME_start_stats();
}
void RUNTIME_stop_stats() {
    return _RUNTIME_stop_stats();
}

void RUNTIME_zdisplay_allprofile () {
    RUNTIME_zdisplay_allprofile ();
}
void RUNTIME_zdisplay_oneprofile (cham_tasktype_t task) {
    return _RUNTIME_zdisplay_oneprofile( task );
}

void RUNTIME_cdisplay_allprofile () {
    return _RUNTIME_cdisplay_allprofile ();
}
void RUNTIME_cdisplay_oneprofile (cham_tasktype_t task){
    return _RUNTIME_cdisplay_oneprofile( task );
}

void RUNTIME_ddisplay_allprofile () {
    return _RUNTIME_ddisplay_allprofile ();
}
void RUNTIME_ddisplay_oneprofile (cham_tasktype_t task){
    return _RUNTIME_ddisplay_oneprofile ( task );
}

void RUNTIME_sdisplay_allprofile () {
    return _RUNTIME_sdisplay_allprofile ();
}
void RUNTIME_sdisplay_oneprofile (cham_tasktype_t task) {
    return _RUNTIME_sdisplay_oneprofile( task );
}

void INSERT_TASK_map( const RUNTIME_option_t *options,
                      cham_uplo_t uplo, const CHAM_desc_t *A, int Am, int An,
                      cham_unary_operator_t op_fct, void *op_args ) {
    return _INSERT_TASK_map( options, uplo, A, Am, An, op_fct, op_args );
}

void INSERT_TASK_scasum( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int M, int N,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_scasum( options, storev, uplo, M, N, A, Am, An, B, Bm, Bn );
}
void INSERT_TASK_caxpy( const RUNTIME_option_t *options,
                        int M, CHAMELEON_Complex32_t alpha,
                        const CHAM_desc_t *A, int Am, int An, int incA,
                        const CHAM_desc_t *B, int Bm, int Bn, int incB ) {
    return _INSERT_TASK_caxpy( options, M, alpha, A, Am, An, incA, B, Bm, Bn, incB );
}
void INSERT_TASK_cbuild( const RUNTIME_option_t *options,
                         const CHAM_desc_t *A, int Am, int An,
                         void *user_data, void* user_build_callback ) {
    return _INSERT_TASK_cbuild( options,
                                A, Am, An,
                                user_data, user_build_callback );
}
void INSERT_TASK_cgeadd( const RUNTIME_option_t *options,
                         cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex32_t beta,  const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_cgeadd( options,
                                trans, m, n, nb,
                                alpha, A, Am, An,
                                beta,  B, Bm, Bn );
}
void INSERT_TASK_cgelqt( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_cgelqt( options,
                                m, n, ib, nb,
                                A, Am, An,
                                T, Tm, Tn );
}
void INSERT_TASK_cgemv( const RUNTIME_option_t *options,
                        cham_trans_t trans, int m, int n,
                        CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *X, int Xm, int Xn, int incX,
                        CHAMELEON_Complex32_t beta,  const CHAM_desc_t *Y, int Ym, int Yn, int incY ) {
    return _INSERT_TASK_cgemv( options,
                               trans, m, n,
                               alpha, A, Am, An,
                               X, Xm, Xn, incX,
                               beta,  Y, Ym, Yn, incY );
}
void INSERT_TASK_cgemm( const RUNTIME_option_t *options,
                        cham_trans_t transA, cham_trans_t transB,
                        int m, int n, int k, int nb,
                        CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Complex32_t beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_cgemm( options,
                               transA, transB,
                               m, n, k, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn,
                               beta, C, Cm, Cn );
}
void INSERT_TASK_cgeqrt( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_cgeqrt( options,
                                m, n, ib, nb,
                                A, Am, An,
                                T, Tm, Tn );
}
void INSERT_TASK_cgessm( const RUNTIME_option_t *options,
                         int m, int n, int k, int ib, int nb,
                         int *IPIV,
                         const CHAM_desc_t *L, int Lm, int Ln,
                         const CHAM_desc_t *D, int Dm, int Dn,
                         const CHAM_desc_t *A, int Am, int An ) {
    return _INSERT_TASK_cgessm( options,
                                m, n, k, ib, nb,
                                IPIV,
                                L, Lm, Ln,
                                D, Dm, Dn,
                                A, Am, An );
}
void INSERT_TASK_cgessq( const RUNTIME_option_t *options,
                         cham_store_t storev, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) {
    return _INSERT_TASK_cgessq( options,
                                storev, m, n,
                                A, Am, An,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
}
void INSERT_TASK_cgetrf( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo ) {
    return _INSERT_TASK_cgetrf( options,
                                m, n, nb,
                                A, Am, An,
                                IPIV,
                                check_info, iinfo );
}
void INSERT_TASK_cgetrf_incpiv( const RUNTIME_option_t *options,
                                int m, int n, int ib, int nb,
                                const CHAM_desc_t *A, int Am, int An,
                                const CHAM_desc_t *L, int Lm, int Ln,
                                int *IPIV,
                                cham_bool_t check_info, int iinfo ) {
    return _INSERT_TASK_cgetrf_incpiv( options,
                                       m, n, ib, nb,
                                       A, Am, An,
                                       L, Lm, Ln,
                                       IPIV,
                                       check_info, iinfo );
}
void INSERT_TASK_cgetrf_nopiv( const RUNTIME_option_t *options,
                               int m, int n, int ib, int nb,
                               const CHAM_desc_t *A, int Am, int An, int iinfo ) {
    return _INSERT_TASK_cgetrf_nopiv( options,
                                      m, n, ib, nb,
                                      A, Am, An, iinfo );
}
void INSERT_TASK_che2ge( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_che2ge( options,
                                uplo,
                                m, n, mb,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_chemm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo,
                        int m, int n, int nb,
                        CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Complex32_t beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_chemm( options,
                               side, uplo,
                               m, n, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn,
                               beta, C, Cm, Cn );
}
void INSERT_TASK_cher2k( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int n, int k, int nb,
                         CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         float beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_cher2k( options,
                                uplo, trans,
                                n, k, nb,
                                alpha, A, Am, An,
                                B, Bm, Bn,
                                beta, C, Cm, Cn );
}
void INSERT_TASK_cherfb( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_cherfb( options,
                                uplo,
                                n, k, ib, nb,
                                A, Am, An,
                                T, Tm, Tn,
                                C, Cm, Cn );
}
void INSERT_TASK_cherk( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_trans_t trans,
                        int n, int k, int nb,
                        float alpha, const CHAM_desc_t *A, int Am, int An,
                        float beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_cherk( options,
                               uplo, trans,
                               n, k, nb,
                               alpha, A, Am, An,
                               beta, C, Cm, Cn );
}
void INSERT_TASK_chessq( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) {
    return _INSERT_TASK_chessq( options,
                                storev, uplo, n,
                                A, Am, An,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
}
void INSERT_TASK_clacpy( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_clacpy( options,
                                uplo, m, n, mb,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_clacpyx( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n, int mb,
                          int displA, const CHAM_desc_t *A, int Am, int An,
                          int displB, const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_clacpyx( options,
                                 uplo, m, n, mb,
                                 displA, A, Am, An,
                                 displB, B, Bm, Bn );
}
void INSERT_TASK_clange( const RUNTIME_option_t *options,
                         cham_normtype_t norm, int M, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_clange( options,
                                norm, M, N, NB,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_clange_max( const RUNTIME_option_t *options,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_clange_max( options,
                                    A, Am, An,
                                    B, Bm, Bn );
}
void INSERT_TASK_clanhe( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_clanhe( options,
                                norm, uplo, N, NB,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_clansy( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_clansy( options,
                                norm, uplo, N, NB,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_clantr( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag,
                         int M, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_clantr( options,
                                norm, uplo, diag,
                                M, N, NB,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_clascal( const RUNTIME_option_t *options,
                          cham_uplo_t uplo,
                          int m, int n, int nb,
                          CHAMELEON_Complex32_t alpha,
                          const CHAM_desc_t *A, int Am, int An )  {
    return _INSERT_TASK_clascal( options,
                                 uplo,
                                 m, n, nb,
                                 alpha,
                                 A, Am, An );
}
void INSERT_TASK_claset( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n1, int n2,
                         CHAMELEON_Complex32_t alpha, CHAMELEON_Complex32_t beta,
                         const CHAM_desc_t *tileA, int tileAm, int tileAn ) {
    return _INSERT_TASK_claset( options,
                                uplo, n1, n2,
                                alpha, beta,
                                tileA, tileAm, tileAn );
}
void INSERT_TASK_claset2( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex32_t alpha,
                          const CHAM_desc_t *tileA, int tileAm, int tileAn ) {
    return _INSERT_TASK_claset2( options,
                                 uplo, n1, n2, alpha,
                                 tileA, tileAm, tileAn );
}
void INSERT_TASK_clatro( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_clatro( options,
                                uplo, trans, m, n, mb,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_clauum( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An ) {
    return _INSERT_TASK_clauum( options,
                                uplo, n, nb,
                                A, Am, An );
}
void INSERT_TASK_cplghe( const RUNTIME_option_t *options,
                         float bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed ) {
    return _INSERT_TASK_cplghe( options,
                                bump, m, n, A, Am, An,
                                bigM, m0, n0, seed );
}
void INSERT_TASK_cplgsy( const RUNTIME_option_t *options,
                         CHAMELEON_Complex32_t bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed ) {
    return _INSERT_TASK_cplgsy( options,
                                bump, m, n, A, Am, An,
                                bigM, m0, n0, seed );
}
void INSERT_TASK_cplrnt( const RUNTIME_option_t *options,
                         int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed ) {
    return _INSERT_TASK_cplrnt( options,
                                m, n, A, Am, An,
                                bigM, m0, n0, seed );
}
void INSERT_TASK_cplssq( const RUNTIME_option_t *options,
                         cham_store_t storev, int M, int N,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                         const CHAM_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn ) {
    return _INSERT_TASK_cplssq( options,
                                storev, M, N,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn,
                                SCLSSQ,     SCLSSQm,     SCLSSQn );
}
void INSERT_TASK_cplssq2( const RUNTIME_option_t *options, int N,
                          const CHAM_desc_t *RESULT, int RESULTm, int RESULTn ) {
    return _INSERT_TASK_cplssq2( options, N,
                                 RESULT, RESULTm, RESULTn );
}
void INSERT_TASK_cpotrf( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int iinfo ) {
    return _INSERT_TASK_cpotrf( options,
                                uplo, n, nb,
                                A, Am, An,
                                iinfo );
}
void INSERT_TASK_cssssm( const RUNTIME_option_t *options,
                         int m1, int n1, int m2, int n2, int k, int ib, int nb,
                         const CHAM_desc_t *A1, int A1m, int A1n,
                         const CHAM_desc_t *A2, int A2m, int A2n,
                         const CHAM_desc_t *L1, int L1m, int L1n,
                         const CHAM_desc_t *L2, int L2m, int L2n,
                         const int *IPIV ) {
    return _INSERT_TASK_cssssm( options,
                                m1, n1, m2, n2, k, ib, nb,
                                A1, A1m, A1n,
                                A2, A2m, A2n,
                                L1, L1m, L1n,
                                L2, L2m, L2n,
                                IPIV );
}
void INSERT_TASK_csymm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo,
                        int m, int n, int nb,
                        CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Complex32_t beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_csymm( options,
                               side, uplo,
                               m, n, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn,
                               beta, C, Cm, Cn );
}
void INSERT_TASK_csyr2k( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int n, int k, int nb,
                         CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         CHAMELEON_Complex32_t beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_csyr2k( options,
                                uplo, trans,
                                n, k, nb,
                                alpha, A, Am, An,
                                B, Bm, Bn,
                                beta, C, Cm, Cn );
}
void INSERT_TASK_csyrk( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_trans_t trans,
                        int n, int k, int nb,
                        CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                        CHAMELEON_Complex32_t beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_csyrk( options,
                               uplo, trans,
                               n, k, nb,
                               alpha, A, Am, An,
                               beta, C, Cm, Cn );
}
void INSERT_TASK_csyssq( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) {
    return _INSERT_TASK_csyssq( options,
                                storev, uplo, n,
                                A, Am, An,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
}
void INSERT_TASK_csytrf_nopiv( const RUNTIME_option_t *options,
                               cham_uplo_t uplo, int n, int nb,
                               const CHAM_desc_t *A, int Am, int An,
                               int iinfo ) {
    return _INSERT_TASK_csytrf_nopiv( options,
                                      uplo, n, nb,
                                      A, Am, An,
                                      iinfo );
}
void INSERT_TASK_ctplqt( const RUNTIME_option_t *options,
                         int m, int n, int l, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_ctplqt( options,
                                m, n, l, ib, nb,
                                A, Am, An,
                                B, Bm, Bn,
                                T, Tm, Tn );
}
void INSERT_TASK_ctpmlqt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int M, int N, int K, int L, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_ctpmlqt( options,
                                 side, trans,
                                 M, N, K, L, ib, nb,
                                 V, Vm, Vn,
                                 T, Tm, Tn,
                                 A, Am, An,
                                 B, Bm, Bn );
}
void INSERT_TASK_ctpmqrt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int m, int n, int k, int l, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_ctpmqrt( options,
                                 side, trans,
                                 m, n, k, l, ib, nb,
                                 V, Vm, Vn,
                                 T, Tm, Tn,
                                 A, Am, An,
                                 B, Bm, Bn );
}
void INSERT_TASK_ctpqrt( const RUNTIME_option_t *options,
                         int m, int n, int l, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_ctpqrt( options,
                                m, n, l, ib, nb,
                                A, Am, An,
                                B, Bm, Bn,
                                T, Tm, Tn );
}
void INSERT_TASK_ctradd( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex32_t beta,  const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_ctradd( options,
                                uplo, trans, m, n, nb,
                                alpha, A, Am, An,
                                beta,  B, Bm, Bn );
}
void INSERT_TASK_ctrasm( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_ctrasm( options,
                                storev, uplo, diag, M, N,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_ctrmm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_ctrmm( options,
                               side, uplo, transA, diag,
                               m, n, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn );
}
void INSERT_TASK_ctrsm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        CHAMELEON_Complex32_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_ctrsm( options,
                               side, uplo, transA, diag,
                               m, n, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn );
}
void INSERT_TASK_ctrssq( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag,
                         int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) {
    return _INSERT_TASK_ctrssq( options,
                                uplo, diag,
                                m, n,
                                A, Am, An,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
}
void INSERT_TASK_ctrtri( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int iinfo ) {
    return _INSERT_TASK_ctrtri( options,
                                uplo, diag, n, nb,
                                A, Am, An,
                                iinfo );
}
void INSERT_TASK_ctsmlq_hetra1( const RUNTIME_option_t *options,
                                cham_side_t side, cham_trans_t trans,
                                int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                const CHAM_desc_t *A1, int A1m, int A1n,
                                const CHAM_desc_t *A2, int A2m, int A2n,
                                const CHAM_desc_t *V, int Vm, int Vn,
                                const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_ctsmlq_hetra1( options,
                                       side, trans,
                                       m1, n1, m2, n2, k, ib, nb,
                                       A1, A1m, A1n,
                                       A2, A2m, A2n,
                                       V, Vm, Vn,
                                       T, Tm, Tn );
}
void INSERT_TASK_ctsmqr_hetra1( const RUNTIME_option_t *options,
                                cham_side_t side, cham_trans_t trans,
                                int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                const CHAM_desc_t *A1, int A1m, int A1n,
                                const CHAM_desc_t *A2, int A2m, int A2n,
                                const CHAM_desc_t *V, int Vm, int Vn,
                                const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_ctsmqr_hetra1( options,
                                       side, trans,
                                       m1, n1, m2, n2, k, ib, nb,
                                       A1, A1m, A1n,
                                       A2, A2m, A2n,
                                       V, Vm, Vn,
                                       T, Tm, Tn );
}
void INSERT_TASK_ctstrf( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *U, int Um, int Un,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *L, int Lm, int Ln,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo ) {
    return _INSERT_TASK_ctstrf( options,
                                m, n, ib, nb,
                                U, Um, Un,
                                A, Am, An,
                                L, Lm, Ln,
                                IPIV,
                                check_info, iinfo );
}
void INSERT_TASK_cunmlq( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int ib,  int nb, int k,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_cunmlq( options,
                                side, trans,
                                m, n, ib,  nb, k,
                                A, Am, An,
                                T, Tm, Tn,
                                C, Cm, Cn );
}
void INSERT_TASK_cunmqr( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_cunmqr( options,
                                side, trans,
                                m, n, k, ib, nb,
                                A, Am, An,
                                T, Tm, Tn,
                                C, Cm, Cn );
}
void INSERT_TASK_cgram( const RUNTIME_option_t *options,
                        cham_uplo_t uplo,
                        int m, int n, int mt, int nt,
                        const CHAM_desc_t *Di, int Dim, int Din,
                        const CHAM_desc_t *Dj, int Djm, int Djn,
                        const CHAM_desc_t *D, int Dm, int Dn,
                        CHAM_desc_t *A, int Am, int An) {
    return _INSERT_TASK_cgram( options,
                               uplo,
                               m, n, mt, nt,
                               Di, Dim, Din,
                               Dj, Djm, Djn,
                               D, Dm, Dn,
                               A, Am, An);
}
void INSERT_TASK_dlag2z( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dlag2z( options, uplo, m, n, A, Am, An, B, Bm, Bn );
}
void INSERT_TASK_dasum( const RUNTIME_option_t *options,
                        cham_store_t storev, cham_uplo_t uplo, int M, int N,
                        const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dasum( options, storev, uplo, M, N, A, Am, An, B, Bm, Bn );
}
void INSERT_TASK_daxpy( const RUNTIME_option_t *options,
                        int M, double alpha,
                        const CHAM_desc_t *A, int Am, int An, int incA,
                        const CHAM_desc_t *B, int Bm, int Bn, int incB ) {
    return _INSERT_TASK_daxpy( options, M, alpha, A, Am, An, incA, B, Bm, Bn, incB );
}
void INSERT_TASK_dbuild( const RUNTIME_option_t *options,
                         const CHAM_desc_t *A, int Am, int An,
                         void *user_data, void* user_build_callback ) {
    return _INSERT_TASK_dbuild( options,
                                A, Am, An,
                                user_data, user_build_callback );
}
void INSERT_TASK_dgeadd( const RUNTIME_option_t *options,
                         cham_trans_t trans, int m, int n, int nb,
                         double alpha, const CHAM_desc_t *A, int Am, int An,
                         double beta,  const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dgeadd( options,
                                trans, m, n, nb,
                                alpha, A, Am, An,
                                beta,  B, Bm, Bn );
}
void INSERT_TASK_dgelqt( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_dgelqt( options,
                                m, n, ib, nb,
                                A, Am, An,
                                T, Tm, Tn );
}
void INSERT_TASK_dgemv( const RUNTIME_option_t *options,
                        cham_trans_t trans, int m, int n,
                        double alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *X, int Xm, int Xn, int incX,
                        double beta,  const CHAM_desc_t *Y, int Ym, int Yn, int incY ) {
    return _INSERT_TASK_dgemv( options,
                               trans, m, n,
                               alpha, A, Am, An,
                               X, Xm, Xn, incX,
                               beta,  Y, Ym, Yn, incY );
}
void INSERT_TASK_dgemm( const RUNTIME_option_t *options,
                        cham_trans_t transA, cham_trans_t transB,
                        int m, int n, int k, int nb,
                        double alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn,
                        double beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_dgemm( options,
                               transA, transB,
                               m, n, k, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn,
                               beta, C, Cm, Cn );
}
void INSERT_TASK_dgeqrt( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_dgeqrt( options,
                                m, n, ib, nb,
                                A, Am, An,
                                T, Tm, Tn );
}
void INSERT_TASK_dgessm( const RUNTIME_option_t *options,
                         int m, int n, int k, int ib, int nb,
                         int *IPIV,
                         const CHAM_desc_t *L, int Lm, int Ln,
                         const CHAM_desc_t *D, int Dm, int Dn,
                         const CHAM_desc_t *A, int Am, int An ) {
    return _INSERT_TASK_dgessm( options,
                                m, n, k, ib, nb,
                                IPIV,
                                L, Lm, Ln,
                                D, Dm, Dn,
                                A, Am, An );
}
void INSERT_TASK_dgessq( const RUNTIME_option_t *options,
                         cham_store_t storev, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) {
    return _INSERT_TASK_dgessq( options,
                                storev, m, n,
                                A, Am, An,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
}
void INSERT_TASK_dgetrf( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo ) {
    return _INSERT_TASK_dgetrf( options,
                                m, n, nb,
                                A, Am, An,
                                IPIV,
                                check_info, iinfo );
}
void INSERT_TASK_dgetrf_incpiv( const RUNTIME_option_t *options,
                                int m, int n, int ib, int nb,
                                const CHAM_desc_t *A, int Am, int An,
                                const CHAM_desc_t *L, int Lm, int Ln,
                                int *IPIV,
                                cham_bool_t check_info, int iinfo ) {
    return _INSERT_TASK_dgetrf_incpiv( options,
                                       m, n, ib, nb,
                                       A, Am, An,
                                       L, Lm, Ln,
                                       IPIV,
                                       check_info, iinfo );
}
void INSERT_TASK_dgetrf_nopiv( const RUNTIME_option_t *options,
                               int m, int n, int ib, int nb,
                               const CHAM_desc_t *A, int Am, int An, int iinfo ) {
    return _INSERT_TASK_dgetrf_nopiv( options,
                                      m, n, ib, nb,
                                      A, Am, An, iinfo );
}
void INSERT_TASK_dhe2ge( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dhe2ge( options,
                                uplo,
                                m, n, mb,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_dsymm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo,
                        int m, int n, int nb,
                        double alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn,
                        double beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_dsymm( options,
                               side, uplo,
                               m, n, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn,
                               beta, C, Cm, Cn );
}
void INSERT_TASK_dsyr2k( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int n, int k, int nb,
                         double alpha, const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         double beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_dsyr2k( options,
                                uplo, trans,
                                n, k, nb,
                                alpha, A, Am, An,
                                B, Bm, Bn,
                                beta, C, Cm, Cn );
}
void INSERT_TASK_dsyrfb( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_dsyrfb( options,
                                uplo,
                                n, k, ib, nb,
                                A, Am, An,
                                T, Tm, Tn,
                                C, Cm, Cn );
}
void INSERT_TASK_dsyrk( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_trans_t trans,
                        int n, int k, int nb,
                        double alpha, const CHAM_desc_t *A, int Am, int An,
                        double beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_dsyrk( options,
                               uplo, trans,
                               n, k, nb,
                               alpha, A, Am, An,
                               beta, C, Cm, Cn );
}
void INSERT_TASK_dsyssq( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) {
    return _INSERT_TASK_dsyssq( options,
                                storev, uplo, n,
                                A, Am, An,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
}
void INSERT_TASK_dlacpy( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dlacpy( options,
                                uplo, m, n, mb,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_dlacpyx( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n, int mb,
                          int displA, const CHAM_desc_t *A, int Am, int An,
                          int displB, const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dlacpyx( options,
                                 uplo, m, n, mb,
                                 displA, A, Am, An,
                                 displB, B, Bm, Bn );
}
void INSERT_TASK_dlange( const RUNTIME_option_t *options,
                         cham_normtype_t norm, int M, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dlange( options,
                                norm, M, N, NB,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_dlange_max( const RUNTIME_option_t *options,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dlange_max( options,
                                    A, Am, An,
                                    B, Bm, Bn );
}
void INSERT_TASK_dlansy( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dlansy( options,
                                norm, uplo, N, NB,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_dlantr( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag,
                         int M, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dlantr( options,
                                norm, uplo, diag,
                                M, N, NB,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_dlascal( const RUNTIME_option_t *options,
                          cham_uplo_t uplo,
                          int m, int n, int nb,
                          double alpha,
                          const CHAM_desc_t *A, int Am, int An )  {
    return _INSERT_TASK_dlascal( options,
                                 uplo,
                                 m, n, nb,
                                 alpha,
                                 A, Am, An );
}
void INSERT_TASK_dlaset( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n1, int n2,
                         double alpha, double beta,
                         const CHAM_desc_t *tileA, int tileAm, int tileAn ) {
    return _INSERT_TASK_dlaset( options,
                                uplo, n1, n2,
                                alpha, beta,
                                tileA, tileAm, tileAn );
}
void INSERT_TASK_dlaset2( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int n1, int n2, double alpha,
                          const CHAM_desc_t *tileA, int tileAm, int tileAn ) {
    return _INSERT_TASK_dlaset2( options,
                                 uplo, n1, n2, alpha,
                                 tileA, tileAm, tileAn );
}
void INSERT_TASK_dlatro( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dlatro( options,
                                uplo, trans, m, n, mb,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_dlauum( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An ) {
    return _INSERT_TASK_dlauum( options,
                                uplo, n, nb,
                                A, Am, An );
}
void INSERT_TASK_dplgsy( const RUNTIME_option_t *options,
                         double bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed ) {
    return _INSERT_TASK_dplgsy( options,
                                bump, m, n, A, Am, An,
                                bigM, m0, n0, seed );
}
void INSERT_TASK_dplrnt( const RUNTIME_option_t *options,
                         int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed ) {
    return _INSERT_TASK_dplrnt( options,
                                m, n, A, Am, An,
                                bigM, m0, n0, seed );
}
void INSERT_TASK_dplssq( const RUNTIME_option_t *options,
                         cham_store_t storev, int M, int N,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                         const CHAM_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn ) {
    return _INSERT_TASK_dplssq( options,
                                storev, M, N,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn,
                                SCLSSQ,     SCLSSQm,     SCLSSQn );
}
void INSERT_TASK_dplssq2( const RUNTIME_option_t *options, int N,
                          const CHAM_desc_t *RESULT, int RESULTm, int RESULTn ) {
    return _INSERT_TASK_dplssq2( options, N,
                                 RESULT, RESULTm, RESULTn );
}
void INSERT_TASK_dpotrf( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int iinfo ) {
    return _INSERT_TASK_dpotrf( options,
                                uplo, n, nb,
                                A, Am, An,
                                iinfo );
}
void INSERT_TASK_dssssm( const RUNTIME_option_t *options,
                         int m1, int n1, int m2, int n2, int k, int ib, int nb,
                         const CHAM_desc_t *A1, int A1m, int A1n,
                         const CHAM_desc_t *A2, int A2m, int A2n,
                         const CHAM_desc_t *L1, int L1m, int L1n,
                         const CHAM_desc_t *L2, int L2m, int L2n,
                         const int *IPIV ) {
    return _INSERT_TASK_dssssm( options,
                                m1, n1, m2, n2, k, ib, nb,
                                A1, A1m, A1n,
                                A2, A2m, A2n,
                                L1, L1m, L1n,
                                L2, L2m, L2n,
                                IPIV );
}
void INSERT_TASK_dsytrf_nopiv( const RUNTIME_option_t *options,
                               cham_uplo_t uplo, int n, int nb,
                               const CHAM_desc_t *A, int Am, int An,
                               int iinfo ) {
    return _INSERT_TASK_dsytrf_nopiv( options,
                                      uplo, n, nb,
                                      A, Am, An,
                                      iinfo );
}
void INSERT_TASK_dtplqt( const RUNTIME_option_t *options,
                         int m, int n, int l, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_dtplqt( options,
                                m, n, l, ib, nb,
                                A, Am, An,
                                B, Bm, Bn,
                                T, Tm, Tn );
}
void INSERT_TASK_dtpmlqt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int M, int N, int K, int L, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dtpmlqt( options,
                                 side, trans,
                                 M, N, K, L, ib, nb,
                                 V, Vm, Vn,
                                 T, Tm, Tn,
                                 A, Am, An,
                                 B, Bm, Bn );
}
void INSERT_TASK_dtpmqrt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int m, int n, int k, int l, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dtpmqrt( options,
                                 side, trans,
                                 m, n, k, l, ib, nb,
                                 V, Vm, Vn,
                                 T, Tm, Tn,
                                 A, Am, An,
                                 B, Bm, Bn );
}
void INSERT_TASK_dtpqrt( const RUNTIME_option_t *options,
                         int m, int n, int l, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_dtpqrt( options,
                                m, n, l, ib, nb,
                                A, Am, An,
                                B, Bm, Bn,
                                T, Tm, Tn );
}
void INSERT_TASK_dtradd( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                         double alpha, const CHAM_desc_t *A, int Am, int An,
                         double beta,  const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dtradd( options,
                                uplo, trans, m, n, nb,
                                alpha, A, Am, An,
                                beta,  B, Bm, Bn );
}
void INSERT_TASK_dtrasm( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dtrasm( options,
                                storev, uplo, diag, M, N,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_dtrmm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        double alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dtrmm( options,
                               side, uplo, transA, diag,
                               m, n, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn );
}
void INSERT_TASK_dtrsm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        double alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dtrsm( options,
                               side, uplo, transA, diag,
                               m, n, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn );
}
void INSERT_TASK_dtrssq( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag,
                         int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) {
    return _INSERT_TASK_dtrssq( options,
                                uplo, diag,
                                m, n,
                                A, Am, An,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
}
void INSERT_TASK_dtrtri( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int iinfo ) {
    return _INSERT_TASK_dtrtri( options,
                                uplo, diag, n, nb,
                                A, Am, An,
                                iinfo );
}
void INSERT_TASK_dtsmlq_hetra1( const RUNTIME_option_t *options,
                                cham_side_t side, cham_trans_t trans,
                                int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                const CHAM_desc_t *A1, int A1m, int A1n,
                                const CHAM_desc_t *A2, int A2m, int A2n,
                                const CHAM_desc_t *V, int Vm, int Vn,
                                const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_dtsmlq_hetra1( options,
                                       side, trans,
                                       m1, n1, m2, n2, k, ib, nb,
                                       A1, A1m, A1n,
                                       A2, A2m, A2n,
                                       V, Vm, Vn,
                                       T, Tm, Tn );
}
void INSERT_TASK_dtsmqr_hetra1( const RUNTIME_option_t *options,
                                cham_side_t side, cham_trans_t trans,
                                int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                const CHAM_desc_t *A1, int A1m, int A1n,
                                const CHAM_desc_t *A2, int A2m, int A2n,
                                const CHAM_desc_t *V, int Vm, int Vn,
                                const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_dtsmqr_hetra1( options,
                                       side, trans,
                                       m1, n1, m2, n2, k, ib, nb,
                                       A1, A1m, A1n,
                                       A2, A2m, A2n,
                                       V, Vm, Vn,
                                       T, Tm, Tn );
}
void INSERT_TASK_dtstrf( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *U, int Um, int Un,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *L, int Lm, int Ln,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo ) {
    return _INSERT_TASK_dtstrf( options,
                                m, n, ib, nb,
                                U, Um, Un,
                                A, Am, An,
                                L, Lm, Ln,
                                IPIV,
                                check_info, iinfo );
}
void INSERT_TASK_dormlq( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int ib,  int nb, int k,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_dormlq( options,
                                side, trans,
                                m, n, ib,  nb, k,
                                A, Am, An,
                                T, Tm, Tn,
                                C, Cm, Cn );
}
void INSERT_TASK_dormqr( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_dormqr( options,
                                side, trans,
                                m, n, k, ib, nb,
                                A, Am, An,
                                T, Tm, Tn,
                                C, Cm, Cn );
}
void INSERT_TASK_dgram( const RUNTIME_option_t *options,
                        cham_uplo_t uplo,
                        int m, int n, int mt, int nt,
                        const CHAM_desc_t *Di, int Dim, int Din,
                        const CHAM_desc_t *Dj, int Djm, int Djn,
                        const CHAM_desc_t *D, int Dm, int Dn,
                        CHAM_desc_t *A, int Am, int An) {
    return _INSERT_TASK_dgram( options,
                               uplo,
                               m, n, mt, nt,
                               Di, Dim, Din,
                               Dj, Djm, Djn,
                               D, Dm, Dn,
                               A, Am, An);
}
void INSERT_TASK_slag2c( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_slag2c( options, uplo, m, n, A, Am, An, B, Bm, Bn );
}
void INSERT_TASK_sasum( const RUNTIME_option_t *options,
                        cham_store_t storev, cham_uplo_t uplo, int M, int N,
                        const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_sasum( options, storev, uplo, M, N, A, Am, An, B, Bm, Bn );
}
void INSERT_TASK_saxpy( const RUNTIME_option_t *options,
                        int M, float alpha,
                        const CHAM_desc_t *A, int Am, int An, int incA,
                        const CHAM_desc_t *B, int Bm, int Bn, int incB ) {
    return _INSERT_TASK_saxpy( options, M, alpha, A, Am, An, incA, B, Bm, Bn, incB );
}
void INSERT_TASK_sbuild( const RUNTIME_option_t *options,
                         const CHAM_desc_t *A, int Am, int An,
                         void *user_data, void* user_build_callback ) {
    return _INSERT_TASK_sbuild( options,
                                A, Am, An,
                                user_data, user_build_callback );
}
void INSERT_TASK_sgeadd( const RUNTIME_option_t *options,
                         cham_trans_t trans, int m, int n, int nb,
                         float alpha, const CHAM_desc_t *A, int Am, int An,
                         float beta,  const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_sgeadd( options,
                                trans, m, n, nb,
                                alpha, A, Am, An,
                                beta,  B, Bm, Bn );
}
void INSERT_TASK_sgelqt( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_sgelqt( options,
                                m, n, ib, nb,
                                A, Am, An,
                                T, Tm, Tn );
}
void INSERT_TASK_sgemv( const RUNTIME_option_t *options,
                        cham_trans_t trans, int m, int n,
                        float alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *X, int Xm, int Xn, int incX,
                        float beta,  const CHAM_desc_t *Y, int Ym, int Yn, int incY ) {
    return _INSERT_TASK_sgemv( options,
                               trans, m, n,
                               alpha, A, Am, An,
                               X, Xm, Xn, incX,
                               beta,  Y, Ym, Yn, incY );
}
void INSERT_TASK_sgemm( const RUNTIME_option_t *options,
                        cham_trans_t transA, cham_trans_t transB,
                        int m, int n, int k, int nb,
                        float alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn,
                        float beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_sgemm( options,
                               transA, transB,
                               m, n, k, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn,
                               beta, C, Cm, Cn );
}
void INSERT_TASK_sgeqrt( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_sgeqrt( options,
                                m, n, ib, nb,
                                A, Am, An,
                                T, Tm, Tn );
}
void INSERT_TASK_sgessm( const RUNTIME_option_t *options,
                         int m, int n, int k, int ib, int nb,
                         int *IPIV,
                         const CHAM_desc_t *L, int Lm, int Ln,
                         const CHAM_desc_t *D, int Dm, int Dn,
                         const CHAM_desc_t *A, int Am, int An ) {
    return _INSERT_TASK_sgessm( options,
                                m, n, k, ib, nb,
                                IPIV,
                                L, Lm, Ln,
                                D, Dm, Dn,
                                A, Am, An );
}
void INSERT_TASK_sgessq( const RUNTIME_option_t *options,
                         cham_store_t storev, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) {
    return _INSERT_TASK_sgessq( options,
                                storev, m, n,
                                A, Am, An,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
}
void INSERT_TASK_sgetrf( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo ) {
    return _INSERT_TASK_sgetrf( options,
                                m, n, nb,
                                A, Am, An,
                                IPIV,
                                check_info, iinfo );
}
void INSERT_TASK_sgetrf_incpiv( const RUNTIME_option_t *options,
                                int m, int n, int ib, int nb,
                                const CHAM_desc_t *A, int Am, int An,
                                const CHAM_desc_t *L, int Lm, int Ln,
                                int *IPIV,
                                cham_bool_t check_info, int iinfo ) {
    return _INSERT_TASK_sgetrf_incpiv( options,
                                       m, n, ib, nb,
                                       A, Am, An,
                                       L, Lm, Ln,
                                       IPIV,
                                       check_info, iinfo );
}
void INSERT_TASK_sgetrf_nopiv( const RUNTIME_option_t *options,
                               int m, int n, int ib, int nb,
                               const CHAM_desc_t *A, int Am, int An, int iinfo ) {
    return _INSERT_TASK_sgetrf_nopiv( options,
                                      m, n, ib, nb,
                                      A, Am, An, iinfo );
}
void INSERT_TASK_she2ge( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_she2ge( options,
                                uplo,
                                m, n, mb,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_ssymm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo,
                        int m, int n, int nb,
                        float alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn,
                        float beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_ssymm( options,
                               side, uplo,
                               m, n, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn,
                               beta, C, Cm, Cn );
}
void INSERT_TASK_ssyr2k( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int n, int k, int nb,
                         float alpha, const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         float beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_ssyr2k( options,
                                uplo, trans,
                                n, k, nb,
                                alpha, A, Am, An,
                                B, Bm, Bn,
                                beta, C, Cm, Cn );
}
void INSERT_TASK_ssyrfb( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_ssyrfb( options,
                                uplo,
                                n, k, ib, nb,
                                A, Am, An,
                                T, Tm, Tn,
                                C, Cm, Cn );
}
void INSERT_TASK_ssyrk( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_trans_t trans,
                        int n, int k, int nb,
                        float alpha, const CHAM_desc_t *A, int Am, int An,
                        float beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_ssyrk( options,
                               uplo, trans,
                               n, k, nb,
                               alpha, A, Am, An,
                               beta, C, Cm, Cn );
}
void INSERT_TASK_ssyssq( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) {
    return _INSERT_TASK_ssyssq( options,
                                storev, uplo, n,
                                A, Am, An,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
}
void INSERT_TASK_slacpy( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_slacpy( options,
                                uplo, m, n, mb,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_slacpyx( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n, int mb,
                          int displA, const CHAM_desc_t *A, int Am, int An,
                          int displB, const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_slacpyx( options,
                                 uplo, m, n, mb,
                                 displA, A, Am, An,
                                 displB, B, Bm, Bn );
}
void INSERT_TASK_slange( const RUNTIME_option_t *options,
                         cham_normtype_t norm, int M, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_slange( options,
                                norm, M, N, NB,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_slange_max( const RUNTIME_option_t *options,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_slange_max( options,
                                    A, Am, An,
                                    B, Bm, Bn );
}
void INSERT_TASK_slansy( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_slansy( options,
                                norm, uplo, N, NB,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_slantr( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag,
                         int M, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_slantr( options,
                                norm, uplo, diag,
                                M, N, NB,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_slascal( const RUNTIME_option_t *options,
                          cham_uplo_t uplo,
                          int m, int n, int nb,
                          float alpha,
                          const CHAM_desc_t *A, int Am, int An )  {
    return _INSERT_TASK_slascal( options,
                                 uplo,
                                 m, n, nb,
                                 alpha,
                                 A, Am, An );
}
void INSERT_TASK_slaset( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n1, int n2,
                         float alpha, float beta,
                         const CHAM_desc_t *tileA, int tileAm, int tileAn ) {
    return _INSERT_TASK_slaset( options,
                                uplo, n1, n2,
                                alpha, beta,
                                tileA, tileAm, tileAn );
}
void INSERT_TASK_slaset2( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int n1, int n2, float alpha,
                          const CHAM_desc_t *tileA, int tileAm, int tileAn ) {
    return _INSERT_TASK_slaset2( options,
                                 uplo, n1, n2, alpha,
                                 tileA, tileAm, tileAn );
}
void INSERT_TASK_slatro( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_slatro( options,
                                uplo, trans, m, n, mb,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_slauum( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An ) {
    return _INSERT_TASK_slauum( options,
                                uplo, n, nb,
                                A, Am, An );
}
void INSERT_TASK_splgsy( const RUNTIME_option_t *options,
                         float bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed ) {
    return _INSERT_TASK_splgsy( options,
                                bump, m, n, A, Am, An,
                                bigM, m0, n0, seed );
}
void INSERT_TASK_splrnt( const RUNTIME_option_t *options,
                         int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed ) {
    return _INSERT_TASK_splrnt( options,
                                m, n, A, Am, An,
                                bigM, m0, n0, seed );
}
void INSERT_TASK_splssq( const RUNTIME_option_t *options,
                         cham_store_t storev, int M, int N,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                         const CHAM_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn ) {
    return _INSERT_TASK_splssq( options,
                                storev, M, N,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn,
                                SCLSSQ,     SCLSSQm,     SCLSSQn );
}
void INSERT_TASK_splssq2( const RUNTIME_option_t *options, int N,
                          const CHAM_desc_t *RESULT, int RESULTm, int RESULTn ) {
    return _INSERT_TASK_splssq2( options, N,
                                 RESULT, RESULTm, RESULTn );
}
void INSERT_TASK_spotrf( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int iinfo ) {
    return _INSERT_TASK_spotrf( options,
                                uplo, n, nb,
                                A, Am, An,
                                iinfo );
}
void INSERT_TASK_sssssm( const RUNTIME_option_t *options,
                         int m1, int n1, int m2, int n2, int k, int ib, int nb,
                         const CHAM_desc_t *A1, int A1m, int A1n,
                         const CHAM_desc_t *A2, int A2m, int A2n,
                         const CHAM_desc_t *L1, int L1m, int L1n,
                         const CHAM_desc_t *L2, int L2m, int L2n,
                         const int *IPIV ) {
    return _INSERT_TASK_sssssm( options,
                                m1, n1, m2, n2, k, ib, nb,
                                A1, A1m, A1n,
                                A2, A2m, A2n,
                                L1, L1m, L1n,
                                L2, L2m, L2n,
                                IPIV );
}
void INSERT_TASK_ssytrf_nopiv( const RUNTIME_option_t *options,
                               cham_uplo_t uplo, int n, int nb,
                               const CHAM_desc_t *A, int Am, int An,
                               int iinfo ) {
    return _INSERT_TASK_ssytrf_nopiv( options,
                                      uplo, n, nb,
                                      A, Am, An,
                                      iinfo );
}
void INSERT_TASK_stplqt( const RUNTIME_option_t *options,
                         int m, int n, int l, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_stplqt( options,
                                m, n, l, ib, nb,
                                A, Am, An,
                                B, Bm, Bn,
                                T, Tm, Tn );
}
void INSERT_TASK_stpmlqt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int M, int N, int K, int L, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_stpmlqt( options,
                                 side, trans,
                                 M, N, K, L, ib, nb,
                                 V, Vm, Vn,
                                 T, Tm, Tn,
                                 A, Am, An,
                                 B, Bm, Bn );
}
void INSERT_TASK_stpmqrt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int m, int n, int k, int l, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_stpmqrt( options,
                                 side, trans,
                                 m, n, k, l, ib, nb,
                                 V, Vm, Vn,
                                 T, Tm, Tn,
                                 A, Am, An,
                                 B, Bm, Bn );
}
void INSERT_TASK_stpqrt( const RUNTIME_option_t *options,
                         int m, int n, int l, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_stpqrt( options,
                                m, n, l, ib, nb,
                                A, Am, An,
                                B, Bm, Bn,
                                T, Tm, Tn );
}
void INSERT_TASK_stradd( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                         float alpha, const CHAM_desc_t *A, int Am, int An,
                         float beta,  const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_stradd( options,
                                uplo, trans, m, n, nb,
                                alpha, A, Am, An,
                                beta,  B, Bm, Bn );
}
void INSERT_TASK_strasm( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_strasm( options,
                                storev, uplo, diag, M, N,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_strmm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        float alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_strmm( options,
                               side, uplo, transA, diag,
                               m, n, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn );
}
void INSERT_TASK_strsm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        float alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_strsm( options,
                               side, uplo, transA, diag,
                               m, n, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn );
}
void INSERT_TASK_strssq( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag,
                         int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) {
    return _INSERT_TASK_strssq( options,
                                uplo, diag,
                                m, n,
                                A, Am, An,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
}
void INSERT_TASK_strtri( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int iinfo ) {
    return _INSERT_TASK_strtri( options,
                                uplo, diag, n, nb,
                                A, Am, An,
                                iinfo );
}
void INSERT_TASK_stsmlq_hetra1( const RUNTIME_option_t *options,
                                cham_side_t side, cham_trans_t trans,
                                int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                const CHAM_desc_t *A1, int A1m, int A1n,
                                const CHAM_desc_t *A2, int A2m, int A2n,
                                const CHAM_desc_t *V, int Vm, int Vn,
                                const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_stsmlq_hetra1( options,
                                       side, trans,
                                       m1, n1, m2, n2, k, ib, nb,
                                       A1, A1m, A1n,
                                       A2, A2m, A2n,
                                       V, Vm, Vn,
                                       T, Tm, Tn );
}
void INSERT_TASK_stsmqr_hetra1( const RUNTIME_option_t *options,
                                cham_side_t side, cham_trans_t trans,
                                int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                const CHAM_desc_t *A1, int A1m, int A1n,
                                const CHAM_desc_t *A2, int A2m, int A2n,
                                const CHAM_desc_t *V, int Vm, int Vn,
                                const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_stsmqr_hetra1( options,
                                       side, trans,
                                       m1, n1, m2, n2, k, ib, nb,
                                       A1, A1m, A1n,
                                       A2, A2m, A2n,
                                       V, Vm, Vn,
                                       T, Tm, Tn );
}
void INSERT_TASK_ststrf( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *U, int Um, int Un,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *L, int Lm, int Ln,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo ) {
    return _INSERT_TASK_ststrf( options,
                                m, n, ib, nb,
                                U, Um, Un,
                                A, Am, An,
                                L, Lm, Ln,
                                IPIV,
                                check_info, iinfo );
}
void INSERT_TASK_sormlq( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int ib,  int nb, int k,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_sormlq( options,
                                side, trans,
                                m, n, ib,  nb, k,
                                A, Am, An,
                                T, Tm, Tn,
                                C, Cm, Cn );
}
void INSERT_TASK_sormqr( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_sormqr( options,
                                side, trans,
                                m, n, k, ib, nb,
                                A, Am, An,
                                T, Tm, Tn,
                                C, Cm, Cn );
}
void INSERT_TASK_sgram( const RUNTIME_option_t *options,
                        cham_uplo_t uplo,
                        int m, int n, int mt, int nt,
                        const CHAM_desc_t *Di, int Dim, int Din,
                        const CHAM_desc_t *Dj, int Djm, int Djn,
                        const CHAM_desc_t *D, int Dm, int Dn,
                        CHAM_desc_t *A, int Am, int An) {
    return _INSERT_TASK_sgram( options,
                               uplo,
                               m, n, mt, nt,
                               Di, Dim, Din,
                               Dj, Djm, Djn,
                               D, Dm, Dn,
                               A, Am, An);
}
void INSERT_TASK_dzasum( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int M, int N,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_dzasum( options, storev, uplo, M, N, A, Am, An, B, Bm, Bn );
}
void INSERT_TASK_zaxpy( const RUNTIME_option_t *options,
                        int M, CHAMELEON_Complex64_t alpha,
                        const CHAM_desc_t *A, int Am, int An, int incA,
                        const CHAM_desc_t *B, int Bm, int Bn, int incB ) {
    return _INSERT_TASK_zaxpy( options, M, alpha, A, Am, An, incA, B, Bm, Bn, incB );
}
void INSERT_TASK_zbuild( const RUNTIME_option_t *options,
                         const CHAM_desc_t *A, int Am, int An,
                         void *user_data, void* user_build_callback ) {
    return _INSERT_TASK_zbuild( options,
                                A, Am, An,
                                user_data, user_build_callback );
}
void INSERT_TASK_zgeadd( const RUNTIME_option_t *options,
                         cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_zgeadd( options,
                                trans, m, n, nb,
                                alpha, A, Am, An,
                                beta,  B, Bm, Bn );
}
void INSERT_TASK_zgelqt( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_zgelqt( options,
                                m, n, ib, nb,
                                A, Am, An,
                                T, Tm, Tn );
}
void INSERT_TASK_zgemv( const RUNTIME_option_t *options,
                        cham_trans_t trans, int m, int n,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *X, int Xm, int Xn, int incX,
                        CHAMELEON_Complex64_t beta,  const CHAM_desc_t *Y, int Ym, int Yn, int incY ) {
    return _INSERT_TASK_zgemv( options,
                               trans, m, n,
                               alpha, A, Am, An,
                               X, Xm, Xn, incX,
                               beta,  Y, Ym, Yn, incY );
}
void INSERT_TASK_zgemm( const RUNTIME_option_t *options,
                        cham_trans_t transA, cham_trans_t transB,
                        int m, int n, int k, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_zgemm( options,
                               transA, transB,
                               m, n, k, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn,
                               beta, C, Cm, Cn );
}
void INSERT_TASK_zgeqrt( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_zgeqrt( options,
                                m, n, ib, nb,
                                A, Am, An,
                                T, Tm, Tn );
}
void INSERT_TASK_zgessm( const RUNTIME_option_t *options,
                         int m, int n, int k, int ib, int nb,
                         int *IPIV,
                         const CHAM_desc_t *L, int Lm, int Ln,
                         const CHAM_desc_t *D, int Dm, int Dn,
                         const CHAM_desc_t *A, int Am, int An ) {
    return _INSERT_TASK_zgessm( options,
                                m, n, k, ib, nb,
                                IPIV,
                                L, Lm, Ln,
                                D, Dm, Dn,
                                A, Am, An );
}
void INSERT_TASK_zgessq( const RUNTIME_option_t *options,
                         cham_store_t storev, int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) {
    return _INSERT_TASK_zgessq( options,
                                storev, m, n,
                                A, Am, An,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
}
void INSERT_TASK_zgetrf( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo ) {
    return _INSERT_TASK_zgetrf( options,
                                m, n, nb,
                                A, Am, An,
                                IPIV,
                                check_info, iinfo );
}
void INSERT_TASK_zgetrf_incpiv( const RUNTIME_option_t *options,
                                int m, int n, int ib, int nb,
                                const CHAM_desc_t *A, int Am, int An,
                                const CHAM_desc_t *L, int Lm, int Ln,
                                int *IPIV,
                                cham_bool_t check_info, int iinfo ) {
    return _INSERT_TASK_zgetrf_incpiv( options,
                                       m, n, ib, nb,
                                       A, Am, An,
                                       L, Lm, Ln,
                                       IPIV,
                                       check_info, iinfo );
}
void INSERT_TASK_zgetrf_nopiv( const RUNTIME_option_t *options,
                               int m, int n, int ib, int nb,
                               const CHAM_desc_t *A, int Am, int An, int iinfo ) {
    return _INSERT_TASK_zgetrf_nopiv( options,
                                      m, n, ib, nb,
                                      A, Am, An, iinfo );
}
void INSERT_TASK_zhe2ge( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_zhe2ge( options,
                                uplo,
                                m, n, mb,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_zhemm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_zhemm( options,
                               side, uplo,
                               m, n, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn,
                               beta, C, Cm, Cn );
}
void INSERT_TASK_zher2k( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int n, int k, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         double beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_zher2k( options,
                                uplo, trans,
                                n, k, nb,
                                alpha, A, Am, An,
                                B, Bm, Bn,
                                beta, C, Cm, Cn );
}
void INSERT_TASK_zherfb( const RUNTIME_option_t *options,
                         cham_uplo_t uplo,
                         int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_zherfb( options,
                                uplo,
                                n, k, ib, nb,
                                A, Am, An,
                                T, Tm, Tn,
                                C, Cm, Cn );
}
void INSERT_TASK_zherk( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_trans_t trans,
                        int n, int k, int nb,
                        double alpha, const CHAM_desc_t *A, int Am, int An,
                        double beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_zherk( options,
                               uplo, trans,
                               n, k, nb,
                               alpha, A, Am, An,
                               beta, C, Cm, Cn );
}
void INSERT_TASK_zhessq( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) {
    return _INSERT_TASK_zhessq( options,
                                storev, uplo, n,
                                A, Am, An,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
}
void INSERT_TASK_zlacpy( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_zlacpy( options,
                                uplo, m, n, mb,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_zlacpyx( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int m, int n, int mb,
                          int displA, const CHAM_desc_t *A, int Am, int An,
                          int displB, const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_zlacpyx( options,
                                 uplo, m, n, mb,
                                 displA, A, Am, An,
                                 displB, B, Bm, Bn );
}
void INSERT_TASK_zlange( const RUNTIME_option_t *options,
                         cham_normtype_t norm, int M, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_zlange( options,
                                norm, M, N, NB,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_zlange_max( const RUNTIME_option_t *options,
                             const CHAM_desc_t *A, int Am, int An,
                             const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_zlange_max( options,
                                    A, Am, An,
                                    B, Bm, Bn );
}
void INSERT_TASK_zlanhe( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_zlanhe( options,
                                norm, uplo, N, NB,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_zlansy( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_zlansy( options,
                                norm, uplo, N, NB,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_zlantr( const RUNTIME_option_t *options,
                         cham_normtype_t norm, cham_uplo_t uplo, cham_diag_t diag,
                         int M, int N, int NB,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_zlantr( options,
                                norm, uplo, diag,
                                M, N, NB,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_zlascal( const RUNTIME_option_t *options,
                          cham_uplo_t uplo,
                          int m, int n, int nb,
                          CHAMELEON_Complex64_t alpha,
                          const CHAM_desc_t *A, int Am, int An )  {
    return _INSERT_TASK_zlascal( options,
                                 uplo,
                                 m, n, nb,
                                 alpha,
                                 A, Am, An );
}
void INSERT_TASK_zlaset( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n1, int n2,
                         CHAMELEON_Complex64_t alpha, CHAMELEON_Complex64_t beta,
                         const CHAM_desc_t *tileA, int tileAm, int tileAn ) {
    return _INSERT_TASK_zlaset( options,
                                uplo, n1, n2,
                                alpha, beta,
                                tileA, tileAm, tileAn );
}
void INSERT_TASK_zlaset2( const RUNTIME_option_t *options,
                          cham_uplo_t uplo, int n1, int n2, CHAMELEON_Complex64_t alpha,
                          const CHAM_desc_t *tileA, int tileAm, int tileAn ) {
    return _INSERT_TASK_zlaset2( options,
                                 uplo, n1, n2, alpha,
                                 tileA, tileAm, tileAn );
}
void INSERT_TASK_zlatro( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int mb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_zlatro( options,
                                uplo, trans, m, n, mb,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_zlauum( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An ) {
    return _INSERT_TASK_zlauum( options,
                                uplo, n, nb,
                                A, Am, An );
}
void INSERT_TASK_zplghe( const RUNTIME_option_t *options,
                         double bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed ) {
    return _INSERT_TASK_zplghe( options,
                                bump, m, n, A, Am, An,
                                bigM, m0, n0, seed );
}
void INSERT_TASK_zplgsy( const RUNTIME_option_t *options,
                         CHAMELEON_Complex64_t bump, int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed ) {
    return _INSERT_TASK_zplgsy( options,
                                bump, m, n, A, Am, An,
                                bigM, m0, n0, seed );
}
void INSERT_TASK_zplrnt( const RUNTIME_option_t *options,
                         int m, int n, const CHAM_desc_t *A, int Am, int An,
                         int bigM, int m0, int n0, unsigned long long int seed ) {
    return _INSERT_TASK_zplrnt( options,
                                m, n, A, Am, An,
                                bigM, m0, n0, seed );
}
void INSERT_TASK_zplssq( const RUNTIME_option_t *options,
                         cham_store_t storev, int M, int N,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn,
                         const CHAM_desc_t *SCLSSQ,     int SCLSSQm,     int SCLSSQn ) {
    return _INSERT_TASK_zplssq( options,
                                storev, M, N,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn,
                                SCLSSQ,     SCLSSQm,     SCLSSQn );
}
void INSERT_TASK_zplssq2( const RUNTIME_option_t *options, int N,
                          const CHAM_desc_t *RESULT, int RESULTm, int RESULTn ) {
    return _INSERT_TASK_zplssq2( options, N,
                                 RESULT, RESULTm, RESULTn );
}
void INSERT_TASK_zpotrf( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int iinfo ) {
    return _INSERT_TASK_zpotrf( options,
                                uplo, n, nb,
                                A, Am, An,
                                iinfo );
}
void INSERT_TASK_zssssm( const RUNTIME_option_t *options,
                         int m1, int n1, int m2, int n2, int k, int ib, int nb,
                         const CHAM_desc_t *A1, int A1m, int A1n,
                         const CHAM_desc_t *A2, int A2m, int A2n,
                         const CHAM_desc_t *L1, int L1m, int L1n,
                         const CHAM_desc_t *L2, int L2m, int L2n,
                         const int *IPIV ) {
    return _INSERT_TASK_zssssm( options,
                                m1, n1, m2, n2, k, ib, nb,
                                A1, A1m, A1n,
                                A2, A2m, A2n,
                                L1, L1m, L1n,
                                L2, L2m, L2n,
                                IPIV );
}
void INSERT_TASK_zsymm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn,
                        CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_zsymm( options,
                               side, uplo,
                               m, n, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn,
                               beta, C, Cm, Cn );
}
void INSERT_TASK_zsyr2k( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans,
                         int n, int k, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_zsyr2k( options,
                                uplo, trans,
                                n, k, nb,
                                alpha, A, Am, An,
                                B, Bm, Bn,
                                beta, C, Cm, Cn );
}
void INSERT_TASK_zsyrk( const RUNTIME_option_t *options,
                        cham_uplo_t uplo, cham_trans_t trans,
                        int n, int k, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        CHAMELEON_Complex64_t beta, const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_zsyrk( options,
                               uplo, trans,
                               n, k, nb,
                               alpha, A, Am, An,
                               beta, C, Cm, Cn );
}
void INSERT_TASK_zsyssq( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) {
    return _INSERT_TASK_zsyssq( options,
                                storev, uplo, n,
                                A, Am, An,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
}
void INSERT_TASK_zsytrf_nopiv( const RUNTIME_option_t *options,
                               cham_uplo_t uplo, int n, int nb,
                               const CHAM_desc_t *A, int Am, int An,
                               int iinfo ) {
    return _INSERT_TASK_zsytrf_nopiv( options,
                                      uplo, n, nb,
                                      A, Am, An,
                                      iinfo );
}
void INSERT_TASK_ztplqt( const RUNTIME_option_t *options,
                         int m, int n, int l, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_ztplqt( options,
                                m, n, l, ib, nb,
                                A, Am, An,
                                B, Bm, Bn,
                                T, Tm, Tn );
}
void INSERT_TASK_ztpmlqt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int M, int N, int K, int L, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_ztpmlqt( options,
                                 side, trans,
                                 M, N, K, L, ib, nb,
                                 V, Vm, Vn,
                                 T, Tm, Tn,
                                 A, Am, An,
                                 B, Bm, Bn );
}
void INSERT_TASK_ztpmqrt( const RUNTIME_option_t *options,
                          cham_side_t side, cham_trans_t trans,
                          int m, int n, int k, int l, int ib, int nb,
                          const CHAM_desc_t *V, int Vm, int Vn,
                          const CHAM_desc_t *T, int Tm, int Tn,
                          const CHAM_desc_t *A, int Am, int An,
                          const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_ztpmqrt( options,
                                 side, trans,
                                 m, n, k, l, ib, nb,
                                 V, Vm, Vn,
                                 T, Tm, Tn,
                                 A, Am, An,
                                 B, Bm, Bn );
}
void INSERT_TASK_ztpqrt( const RUNTIME_option_t *options,
                         int m, int n, int l, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn,
                         const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_ztpqrt( options,
                                m, n, l, ib, nb,
                                A, Am, An,
                                B, Bm, Bn,
                                T, Tm, Tn );
}
void INSERT_TASK_ztradd( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_trans_t trans, int m, int n, int nb,
                         CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                         CHAMELEON_Complex64_t beta,  const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_ztradd( options,
                                uplo, trans, m, n, nb,
                                alpha, A, Am, An,
                                beta,  B, Bm, Bn );
}
void INSERT_TASK_ztrasm( const RUNTIME_option_t *options,
                         cham_store_t storev, cham_uplo_t uplo, cham_diag_t diag, int M, int N,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_ztrasm( options,
                                storev, uplo, diag, M, N,
                                A, Am, An,
                                B, Bm, Bn );
}
void INSERT_TASK_ztrmm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_ztrmm( options,
                               side, uplo, transA, diag,
                               m, n, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn );
}
void INSERT_TASK_ztrsm( const RUNTIME_option_t *options,
                        cham_side_t side, cham_uplo_t uplo, cham_trans_t transA, cham_diag_t diag,
                        int m, int n, int nb,
                        CHAMELEON_Complex64_t alpha, const CHAM_desc_t *A, int Am, int An,
                        const CHAM_desc_t *B, int Bm, int Bn ) {
    return _INSERT_TASK_ztrsm( options,
                               side, uplo, transA, diag,
                               m, n, nb,
                               alpha, A, Am, An,
                               B, Bm, Bn );
}
void INSERT_TASK_ztrssq( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag,
                         int m, int n,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *SCALESUMSQ, int SCALESUMSQm, int SCALESUMSQn ) {
    return _INSERT_TASK_ztrssq( options,
                                uplo, diag,
                                m, n,
                                A, Am, An,
                                SCALESUMSQ, SCALESUMSQm, SCALESUMSQn );
}
void INSERT_TASK_ztrtri( const RUNTIME_option_t *options,
                         cham_uplo_t uplo, cham_diag_t diag, int n, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         int iinfo ) {
    return _INSERT_TASK_ztrtri( options,
                                uplo, diag, n, nb,
                                A, Am, An,
                                iinfo );
}
void INSERT_TASK_ztsmlq_hetra1( const RUNTIME_option_t *options,
                                cham_side_t side, cham_trans_t trans,
                                int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                const CHAM_desc_t *A1, int A1m, int A1n,
                                const CHAM_desc_t *A2, int A2m, int A2n,
                                const CHAM_desc_t *V, int Vm, int Vn,
                                const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_ztsmlq_hetra1( options,
                                       side, trans,
                                       m1, n1, m2, n2, k, ib, nb,
                                       A1, A1m, A1n,
                                       A2, A2m, A2n,
                                       V, Vm, Vn,
                                       T, Tm, Tn );
}
void INSERT_TASK_ztsmqr_hetra1( const RUNTIME_option_t *options,
                                cham_side_t side, cham_trans_t trans,
                                int m1, int n1, int m2, int n2, int k, int ib, int nb,
                                const CHAM_desc_t *A1, int A1m, int A1n,
                                const CHAM_desc_t *A2, int A2m, int A2n,
                                const CHAM_desc_t *V, int Vm, int Vn,
                                const CHAM_desc_t *T, int Tm, int Tn ) {
    return _INSERT_TASK_ztsmqr_hetra1( options,
                                       side, trans,
                                       m1, n1, m2, n2, k, ib, nb,
                                       A1, A1m, A1n,
                                       A2, A2m, A2n,
                                       V, Vm, Vn,
                                       T, Tm, Tn );
}
void INSERT_TASK_ztstrf( const RUNTIME_option_t *options,
                         int m, int n, int ib, int nb,
                         const CHAM_desc_t *U, int Um, int Un,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *L, int Lm, int Ln,
                         int *IPIV,
                         cham_bool_t check_info, int iinfo ) {
    return _INSERT_TASK_ztstrf( options,
                                m, n, ib, nb,
                                U, Um, Un,
                                A, Am, An,
                                L, Lm, Ln,
                                IPIV,
                                check_info, iinfo );
}
void INSERT_TASK_zunmlq( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int ib,  int nb, int k,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_zunmlq( options,
                                side, trans,
                                m, n, ib,  nb, k,
                                A, Am, An,
                                T, Tm, Tn,
                                C, Cm, Cn );
}
void INSERT_TASK_zunmqr( const RUNTIME_option_t *options,
                         cham_side_t side, cham_trans_t trans,
                         int m, int n, int k, int ib, int nb,
                         const CHAM_desc_t *A, int Am, int An,
                         const CHAM_desc_t *T, int Tm, int Tn,
                         const CHAM_desc_t *C, int Cm, int Cn ) {
    return _INSERT_TASK_zunmqr( options,
                                side, trans,
                                m, n, k, ib, nb,
                                A, Am, An,
                                T, Tm, Tn,
                                C, Cm, Cn );
}
void INSERT_TASK_zgram( const RUNTIME_option_t *options,
                        cham_uplo_t uplo,
                        int m, int n, int mt, int nt,
                        const CHAM_desc_t *Di, int Dim, int Din,
                        const CHAM_desc_t *Dj, int Djm, int Djn,
                        const CHAM_desc_t *D, int Dm, int Dn,
                        CHAM_desc_t *A, int Am, int An) {
    return _INSERT_TASK_zgram( options,
                               uplo,
                               m, n, mt, nt,
                               Di, Dim, Din,
                               Dj, Djm, Djn,
                               D, Dm, Dn,
                               A, Am, An);
}

void INSERT_TASK_dlag2s( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         CHAM_desc_t *A, int Am, int An, int lda,
                         CHAM_desc_t *B, int Bm, int Bn, int ldb ) {
    return _INSERT_TASK_dlag2s( options, m, n, nb, A, Am, An, lda, B, Bm, Bn, ldb );
}
void INSERT_TASK_slag2d( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         CHAM_desc_t *A, int Am, int An, int lda,
                         CHAM_desc_t *B, int Bm, int Bn, int ldb ) {
    return _INSERT_TASK_slag2d( options, m, n, nb, A, Am, An, lda, B, Bm, Bn, ldb );
}
void INSERT_TASK_clag2z( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         CHAM_desc_t *A, int Am, int An, int lda,
                         CHAM_desc_t *B, int Bm, int Bn, int ldb ) {
    return _INSERT_TASK_clag2z( options, m, n, nb, A, Am, An, lda, B, Bm, Bn, ldb );
}
void INSERT_TASK_zlag2c( const RUNTIME_option_t *options,
                         int m, int n, int nb,
                         CHAM_desc_t *A, int Am, int An, int lda,
                         CHAM_desc_t *B, int Bm, int Bn, int ldb ) {
    return _INSERT_TASK_zlag2c( options, m, n, nb, A, Am, An, lda, B, Bm, Bn, ldb );
}
