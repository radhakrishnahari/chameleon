/**
 *
 * @file runtime.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief The common runtimes API
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2020-03-03
 *
 */
#ifndef _chameleon_runtime_h_
#define _chameleon_runtime_h_

#include "chameleon/config.h"
#include "chameleon/constants.h"
#include "chameleon/runtime_struct.h"
#include "chameleon/tasks.h"

BEGIN_C_DECLS

/**
 * @name RUNTIME Context functions
 * @{
 *    These functions manipulate the context data structure before the runtime is
 *    started and after the runtime is stopped.
 */

/**
 * @brief Create the runtime specific options in the context before starting it
 *
 * @param[in,out] ctxt
 *            The runtime context to initialize for the runtime.
 */
void RUNTIME_context_create( CHAM_context_t *ctxt );
void (*CHAMELEON_RUNTIME_context_create)( CHAM_context_t *ctxt );

/**
 * @brief Destroy the specific options in the context after this last one has
 * been stop.
 *
 * @param[in,out] ctxt
 *            The runtime context in which the runtime specific options must
 *            be destroyed
 */
void RUNTIME_context_destroy( CHAM_context_t *ctxt );
void (*CHAMELEON_RUNTIME_context_destroy)( CHAM_context_t *ctxt );

/**
 * @brief Enable a global option of the runtime.
 * @warning Should be called only by CHAMELEON_Enable()
 *
 * @param[in] runtime_ctxt
 *            Pointer to the runtime data structure
 *
 * @param[in] option
 *            @arg CHAMELEON_PROFILING_MODE: start the profiling mode of the runtime.
 */
void RUNTIME_enable( void *runtime_ctxt, int option );
void (*CHAMELEON_RUNTIME_enable)( void *runtime_ctxt, int option );

/**
 * @brief Disable a global option of the runtime.
 * @warning Should be called only by CHAMELEON_Disable()
 *
 * @param[in] runtime_ctxt
 *            Pointer to the runtime data structure
 *
 * @param[in] option
 *            @arg CHAMELEON_PROFILING_MODE: stop the profiling mode of the runtime.
 */
void RUNTIME_disable( void *runtime_ctxt, int option );
void (*CHAMELEON_RUNTIME_disable)( void *runtime_ctxt, int option );

/**
 * @}
 *
 * @name RUNTIME Control functions
 * @{
 *   These functions control the global behavior of the runtime.
 */

/**
 * @brief Initialize the scheduler with the given parameters
 *
 * @param[in,out] ctxt
 *            The runtime context in which to initialize the runtime support.
 *
 * @param[in] ncpus
 *            Defines the total number of cores given to the runtime per
 *            node. Including cuda and communication workers for runtimes that
 *            dedicates cores to this. ncpus > 0, or -1 to target the whole machine.
 *
 * @param[in] ncudas
 *            Defines the number of CUDA devices used by node. If the runtime
 *            uses one core dedicated to each CUDA device, they will be taken
 *            from ncpus. If ncpus > 0, ncudas < ncpus. -1 to target all the
 *            CUDA devices available.
 *
 * @param[in] nthreads_per_worker
 *            Defines the number of threads per worker when multi-threaded tasks
 *            are enabled.  This is used to exploit parallel BLAS kernels, and
 *            defines a better binding of the workers.
 *            -1 to disable, or > 0 to enable.
 *
 * @retval -1 on failure to initialize the runtime.
 * @retval >0 on success to initialize the runtime.
 *
 */
int RUNTIME_init( CHAM_context_t *ctxt,
              int ncpus,
              int ncudas,
              int nthreads_per_worker );

int (*CHAMELEON_RUNTIME_init)( CHAM_context_t *ctxt,
              int ncpus,
              int ncudas,
              int nthreads_per_worker );

/**
 * @brief Finalize the scheduler used for the computations.
 *
 * @param[in,out] ctxt
 *            The runtime context for which the runtime system must be shut down.
 */
void RUNTIME_finalize( CHAM_context_t *ctxt );
void (*CHAMELEON_RUNTIME_finalize)( CHAM_context_t *ctxt );

/**
 * @brief Suspend the processing of new tasks submitted to the runtime system.
 *
 * @param[in] ctxt
 *            The runtime context for which the suspension must be made.
 */
void RUNTIME_pause( CHAM_context_t *ctxt );
void (*CHAMELEON_RUNTIME_pause)( CHAM_context_t *ctxt );

/**
 * @brief Resume the processing of new tasks submitted to the runtime system.
 *
 * @param[in] ctxt
 *            The runtime context for which the execution must be resumed.
 */
void RUNTIME_resume( CHAM_context_t *ctxt );
void (*CHAMELEON_RUNTIME_resume)( CHAM_context_t *ctxt );

/**
 * @brief Wait for completion of all tasks submitted to the runtime.
 *
 * @param[in] ctxt
 *            The runtime context in which the task completion is performed.
 */
void RUNTIME_barrier( CHAM_context_t *ctxt );
void (*CHAMELEON_RUNTIME_barrier)( CHAM_context_t *ctxt );

/**
 * @brief Show the progress of the computations when enabled.
 *
 * @param[in] ctxt
 *            The runtime context for which the context needs to be printed.
 */
void RUNTIME_progress( CHAM_context_t *ctxt );
void (*CHAMELEON_RUNTIME_progress)( CHAM_context_t *ctxt );

/**
 * @brief Get the rank of the current worker for the runtime.
 *
 * @param[in] ctxt
 *            The runtime context for which the thread rank is asked.
 *
 * @retval The rank of the current thread in the runtime.
 */
int RUNTIME_thread_rank( CHAM_context_t *ctxt );
int (*CHAMELEON_RUNTIME_thread_rank)( CHAM_context_t *ctxt );

/**
 * @brief Get the number of CPU workers of the runtime.
 *
 * @param[in] ctxt
 *            The runtime context for which the number of workers is requested
 *
 * @retval The number of threads currently used by the runtime.
 */
int (*CHAMELEON_RUNTIME_thread_size)( CHAM_context_t *ctxt );
int RUNTIME_thread_size( CHAM_context_t *ctxt );

/**
 * @brief Get the MPI comm rank of the current process related to the runtime.
 *
 * @param[in] ctxt
 *            The runtime context for which the rank is asked.
 *
 * @retval The rank of the process in the communicator known by the runtime.
 */
int RUNTIME_comm_rank( CHAM_context_t *ctxt );
int (*CHAMELEON_RUNTIME_comm_rank)( CHAM_context_t *ctxt );

/**
 * @brief Get the MPI comm size related to the runtime.
 *
 * @param[in] ctxt
 *            The runtime context for which the communicator size is asked.
 *
 * @retval The size of the communicator known by the runtime.
 */
int RUNTIME_comm_size( CHAM_context_t *ctxt );
int (*CHAMELEON_RUNTIME_comm_size)( CHAM_context_t *ctxt );

/**
 * @brief Set the data tags bounds based on runtime tags limit.
 *
 * @param[in] user_tag_width
 *            Total number of bits used to defined the data tags (31 by default).
 *            Must be larger than 20.
 *
 * @param[in] user_tag_sep
 *            The number of bits dedicated to identify the pieces of data per
 *            descriptor.  (user_tag_width - user_tag_sep) defines the number of
 *            bits used to enumerate the descriptors.
 */
void RUNTIME_comm_set_tag_sizes( int user_tag_width, int user_tag_sep );
void (*CHAMELEON_RUNTIME_comm_set_tag_sizes)( int user_tag_width, int user_tag_sep );

/**
 * @}
 *
 * @name RUNTIME Asynchonous functionalities
 * @{
 *    These functions manage the sequences of tasks. A sequence is a subset of
 *    related tasks belonging to a same problem.
 */

/**
 * @brief Create a sequence structure associated to a given context.
 *
 * @param[in] ctxt
 *            The runtime context in which the sequence is created.
 *
 * @param[in,out] sequence
 *            On entry the allocated runtime sequence structure, and on exit
 *            the scheduler specifics for the sequence have been initialized.
 *
 * @retval CHAMELEON_SUCCESS on success.
 * @retval CHAMELEON_ERR_OUT_OF_RESOURCES, if the sequence could not be created.
 */
int RUNTIME_sequence_create( CHAM_context_t *ctxt, RUNTIME_sequence_t *sequence );
int (*CHAMELEON_RUNTIME_sequence_create)( CHAM_context_t *ctxt, RUNTIME_sequence_t *sequence );

/**
 * @brief Destroy the sequence structure.
 *
 * @param[in] ctxt
 *            The runtime context to which the sequence belongs.
 *
 * @param[in,out] sequence
 *            On entry the sequence structure.
 *            On exit, the scheduler specifics of the sequence have been
 *            destroyed.
 *
 * @retval CHAMELEON_SUCCESS on success.
 */
int RUNTIME_sequence_destroy( CHAM_context_t *ctxt, RUNTIME_sequence_t *sequence);
int (*CHAMELEON_RUNTIME_sequence_destroy)( CHAM_context_t *ctxt, RUNTIME_sequence_t *sequence);

/**
 * @brief Wait for completion of all tasks in the given sequence.
 *
 * @param[in] ctxt
 *            The runtime context to which the sequence belongs.
 *
 * @param[in] sequence
 *            The sequence that gathers a set of tasks to complete.
 *
 * @retval CHAMELEON_SUCCESS on success.
 * @retval The algorithm return code on failure.
 */
int RUNTIME_sequence_wait( CHAM_context_t *ctxt, RUNTIME_sequence_t *sequence );
int (*CHAMELEON_RUNTIME_sequence_wait)( CHAM_context_t *ctxt, RUNTIME_sequence_t *sequence );

/**
 * @brief Cancel all remaining tasks in the given sequence.
 *
 * @param[in] ctxt
 *            The runtime context to which the sequence belongs.
 *
 * @param[in,out] sequence
 *            The sequence that gathers the set of tasks to cancel.
 *            On exit, the sequence holds the error code of the algorithm and
 *            the request that failed in the sequence.
 *
 * @param[in,out] request
 *            The request that failed in the sequence.
 *            On exit, the return status is set to status.
 *
 * @param[in] status
 *            The return status of the task that failed in the request given
 *            above.
 */
void RUNTIME_sequence_flush( CHAM_context_t *ctxt,
                             RUNTIME_sequence_t *sequence,
                             RUNTIME_request_t  *request,
                             int                 status );
void (*CHAMELEON_RUNTIME_sequence_flush)( CHAM_context_t *ctxt,
                                          RUNTIME_sequence_t *sequence,
                                          RUNTIME_request_t  *request,
                                          int                 status );

/**
 * @}
 *
 * @name RUNTIME Asynchonous functionalities
 * @{
 *    These functions manage the request of tasks. A request is a subset of
 *    related tasks belonging to a same procedure.
 */

/**
 * @brief Create a request structure associated to a given context.
 *
 * @param[in] ctxt
 *            The runtime context in which the request is created.
 *
 * @param[in,out] request
 *            On entry the allocated runtime request structure, and on exit
 *            the scheduler specifics for the request have been initialized.
 *
 * @retval CHAMELEON_SUCCESS on success.
 * @retval CHAMELEON_ERR_OUT_OF_RESOURCES, if the request could not be created.
 */
int RUNTIME_request_create( CHAM_context_t *ctxt, RUNTIME_request_t *request );
int (*CHAMELEON_RUNTIME_request_create)( CHAM_context_t *ctxt, RUNTIME_request_t *request );

/**
 * @brief Destroy the request structure.
 *
 * @param[in] ctxt
 *            The runtime context to which the request belongs.
 *
 * @param[in,out] request
 *            On entry the request structure.
 *            On exit, the scheduler specifics of the request have been
 *            destroyed.
 *
 * @retval CHAMELEON_SUCCESS on success.
 */
int RUNTIME_request_destroy( CHAM_context_t *ctxt, RUNTIME_request_t *request);
int (*CHAMELEON_RUNTIME_request_destroy)( CHAM_context_t *ctxt, RUNTIME_request_t *request);

/**
 *  @brief RUNTIME_request_set - Set RUNTIME parameter for a request.
 *
 *******************************************************************************
 *
 * @param[in] ctxt
 *            The runtime context to which the request belongs.
 *
 * @param[in,out] request
 *            On entry the request structure.
 *            On exit, the request structure modified.
 *
 * @param[in] param
 *          Feature to be enabled:
 *          @arg CHAMELEON_REQUEST_WORKERID: force tasks execution on a specific worker id
 *
 * @param[in] value
 *          Value of the parameter.
 *
 *******************************************************************************
 *
 * @retval CHAMELEON_SUCCESS successful exit
 *
 */
int RUNTIME_request_set( CHAM_context_t  *chamctxt,
                         RUNTIME_request_t *request,
                         int param, int value );
int (*CHAMELEON_RUNTIME_request_set)( CHAM_context_t  *chamctxt,
                                      RUNTIME_request_t *request,
                                      int param, int value );

/**
 * @}
 *
 * @name RUNTIME Memory management
 * @{
 ***
 * @brief Allocate size bytes through the runtime memory management system if any, or with malloc otherwise.
 *
 * This function allows to allocate pinned memory if needed, or eventually not
 * perform the allocation in simulation mode.

 * @param[in] size
 *            The number of bytes to allocate.
 *
 * @return The pointer to allocated area of size bytes on success, NULL otherwise.
 */
void *RUNTIME_malloc( size_t size );
void *(*CHAMELEON_RUNTIME_malloc)( size_t size );

/**
 * @brief Free allocated memory through RUNTIME_malloc() function call
 *
 * @param[in,out] ptr
 *            The ptr to free
 *
 * @param[in] size
 *            The size in bytes of the allocated area associated to ptr.
 */
void RUNTIME_free( void *ptr, size_t size );
void (*CHAMELEON_RUNTIME_free)( void *ptr, size_t size );

/**
 * @}
 *
 * @name RUNTIME Descriptor functions
 * @{
 ***
 * @brief Initialize runtime specific data structure to a given descriptor.
 *
 * @param[in,out] desc
 *            The descriptor to initialize.
 */
void RUNTIME_desc_create( CHAM_desc_t *desc );
void (*CHAMELEON_RUNTIME_desc_create)( CHAM_desc_t *desc );

/**
 * @brief Finalize runtime specific data structure of a given descriptor.
 *
 * @param[in,out] desc
 *            The descriptor to finalize.
 */
void RUNTIME_desc_destroy( CHAM_desc_t *desc );
void (*CHAMELEON_RUNTIME_desc_destroy)( CHAM_desc_t *desc );

/**
 * @brief Acquire in main memory an up-to-date copy of the data described by the
 * descriptor for read-write access.
 *
 * The application must call this function prior to accessing registered data
 * from main memory outside tasks. RUNTIME_desc_release() must be called once
 * the application does not need to access the data anymore.
 *
 * @param[in] desc
 *            The descriptor to acquire.
 *
 * @retval CHAMELEON_SUCCESS on success
 */
int RUNTIME_desc_acquire( const CHAM_desc_t *desc );
int (*CHAMELEON_RUNTIME_desc_acquire)( const CHAM_desc_t *desc );

/**
 * @brief Release the data described by the descriptor to be used by the runtime
 * tasks again.
 *
 * This function releases the data acquired by the application either by
 * RUNTIME_desc_acquire() or by RUNTIME_desc_acquire_async() to the runtime.
 *
 * @param[in] desc
 *            The descriptor to release.
 *
 * @retval CHAMELEON_SUCCESS on success
 */
int RUNTIME_desc_release( const CHAM_desc_t *desc );
int (*CHAMELEON_RUNTIME_desc_release)( const CHAM_desc_t *desc );

/**
 * @brief Flush all pieces of data from a descriptor.
 *
 * This function marks all pieces of data from the descriptor as unused in the future to:
 *   - cleanup the data from the distributed cache if present.
 *   - mark for eviction from the GPU if space is needed
 *   - move the data back to the main memory
 *
 * @param[in] desc
 *            The descriptor to flush.
 *
 * @param[in] sequence
 *            The sequence in which the data is used.
 */
void RUNTIME_desc_flush( const CHAM_desc_t *desc, const RUNTIME_sequence_t *sequence );
void (*CHAMELEON_RUNTIME_desc_flush)( const CHAM_desc_t *desc, const RUNTIME_sequence_t *sequence );

/**
 * @brief Flush all data submitted to the runtime systems from the distributed
 * cache, and/or mark for eviction from the GPU memory.
 *
 * This function flushes all data from the distributed cache of the runtime system.
 */
void RUNTIME_flush( );
void (*CHAMELEON_RUNTIME_flush)( );

/**
 * @brief Flush a single piece of data.
 *
 * This function marks a piece of data as unused in the future to:
 *   - cleanup the data from the distributed cache if present.
 *   - mark for eviction from the GPU if space is needed
 *   - move the data back to the main memory
 *
 * @param[in] sequence
 *            The sequence in which the data is used.
 *
 * @param[in] A
 *            The descriptor to which the piece of data belongs.
 *
 * @param[in] Am
 *            The row coordinate of the piece of data in the matrix
 *
 * @param[in] An
 *            The column coordinate of the piece of data in the matrix
 */
void RUNTIME_data_flush( const RUNTIME_sequence_t *sequence,
                         const CHAM_desc_t *A, int Am, int An );
void (*CHAMELEON_RUNTIME_data_flush)( const RUNTIME_sequence_t *sequence,
                                      const CHAM_desc_t *A, int Am, int An );

/**
 * @brief Migrate a single piece of data.
 *
 * This function migrate a piece of data from its original rank to the new_rank
 * and changes its ownership.
 *
 * @param[in] sequence
 *            The sequence in which the data is used.
 *
 * @param[in] A
 *            The descriptor to which the piece of data belongs.
 *
 * @param[in] Am
 *            The row coordinate of the piece of data in the matrix
 *
 * @param[in] An
 *            The column coordinate of the piece of data in the matrix
 *
 * @param[in] new_rank
 *            The new_rank on which to migrate the data
 */
#if defined(CHAMELEON_USE_MIGRATE)
void RUNTIME_data_migrate( const RUNTIME_sequence_t *sequence,
                           const CHAM_desc_t *A, int Am, int An, int new_rank );
void (*CHAMELEON_RUNTIME_data_migrate)( const RUNTIME_sequence_t *sequence,
                                        const CHAM_desc_t *A, int Am, int An, int new_rank );
#else
static inline void CHAMELEON_RUNTIME_data_migrate( const RUNTIME_sequence_t *sequence,
                                         const CHAM_desc_t *A, int Am, int An, int new_rank )
{ (void)sequence; (void)A; (void)Am; (void)An; (void)new_rank; }
#endif

/**
 * @brief Get the pointer to the data or the runtime handler associated to the
 * piece of data (m, n) in desc.
 *
 * @param[in] A
 *            The descriptor to which belongs the piece of data.
 *
 * @param[in] Am
 *            The row coordinate of the piece of data in the matrix
 *
 * @param[in] An
 *            The column coordinate of the piece of data in the matrix
 *
 * @retval The runtime handler address of the piece of data.
 */
void *RUNTIME_data_getaddr( const CHAM_desc_t *A, int Am, int An );

/**
 * @}
 *
 * @name RUNTIME Insert task options management
 * @{
 */

/**
 * @brief Initialize the runtime option structure
 *
 * @param[in,out] options
 *            The options structure that needs to be initialized
 *
 * @param[in] context
 *            The runtime context in which to initialize the runtime support.
 *
 * @param[in] sequence
 *            The sequence structure to associate in the options.
 *
 * @param[in] request
 *            The request structure to associate in the options.
 *
 */
void RUNTIME_options_init( RUNTIME_option_t   *options,
                           CHAM_context_t  *context,
                           RUNTIME_sequence_t *sequence,
                           RUNTIME_request_t  *request );
void (*CHAMELEON_RUNTIME_options_init)( RUNTIME_option_t   *options,
                                        CHAM_context_t  *context,
                                        RUNTIME_sequence_t *sequence,
                                        RUNTIME_request_t  *request );

/**
 * @brief Finalize the runtime options structure
 *
 * @param[in,out] options
 *            The options structure to finalize
 *
 * @param[in] context
 *            The runtime context in which to options structure has been
 *            initialized.
 *
 */
void RUNTIME_options_finalize( RUNTIME_option_t  *options,
                               CHAM_context_t *context );
void (*CHAMELEON_RUNTIME_options_finalize)( RUNTIME_option_t  *options,
                                            CHAM_context_t *context );

/**
 * @brief Declare and initialize the workspaces for the given options structure
 *
 * @param[in,out] options
 *            The options structure in which to initialize the workspaces
 *
 * @param[in] wsize
 *            The worker workspace size required. This workspace will be
 *            allocated on the device that runs the kernel.
 *
 * @param[in] hsize
 *            The host workspace size required. This workspace will always be
 *            allocated on the host and will be used by hybrid CPU/GPU kernels.
 *
 */
int RUNTIME_options_ws_alloc( RUNTIME_option_t *options,
                              size_t wsize, size_t hsize );
int (*CHAMELEON_RUNTIME_options_ws_alloc)( RUNTIME_option_t *options,
                                           size_t wsize, size_t hsize );

/**
 * @brief Submit the release of the workspaces associated to the options structure.
 *
 * @warning With some runtime, this call might have to call a synchronization.
 *
 * @param[in,out] options
 *            The options structure for which to workspcaes will be released
 *
 */
int RUNTIME_options_ws_free( RUNTIME_option_t *options );
int (*CHAMELEON_RUNTIME_options_ws_free)( RUNTIME_option_t *options );

/**
 * @}
 *
 * @name RUNTIME Kernel locality management
 * @{
 */
void RUNTIME_zlocality_allrestore ();
void RUNTIME_clocality_allrestore ();
void RUNTIME_dlocality_allrestore ();
void RUNTIME_slocality_allrestore ();

void RUNTIME_zlocality_allrestrict(uint32_t device);
void RUNTIME_zlocality_onerestrict(cham_tasktype_t task, uint32_t device);
void RUNTIME_zlocality_onerestore (cham_tasktype_t task);

void RUNTIME_clocality_allrestrict(uint32_t device);
void RUNTIME_clocality_onerestrict(cham_tasktype_t task, uint32_t device);
void RUNTIME_clocality_onerestore (cham_tasktype_t task);

void RUNTIME_dlocality_allrestrict(uint32_t device);
void RUNTIME_dlocality_onerestrict(cham_tasktype_t task, uint32_t device);
void RUNTIME_dlocality_onerestore (cham_tasktype_t task);

void RUNTIME_slocality_allrestrict(uint32_t device);
void RUNTIME_slocality_onerestrict(cham_tasktype_t task, uint32_t device);
void RUNTIME_slocality_onerestore (cham_tasktype_t task);

void (*CHAMELEON_RUNTIME_zlocality_allrestore) ();
void (*CHAMELEON_RUNTIME_clocality_allrestore) ();
void (*CHAMELEON_RUNTIME_dlocality_allrestore) ();
void (*CHAMELEON_RUNTIME_slocality_allrestore) ();

void (*CHAMELEON_RUNTIME_zlocality_allrestrict)(uint32_t device);
void (*CHAMELEON_RUNTIME_zlocality_onerestrict)(cham_tasktype_t task, uint32_t device);
void (*CHAMELEON_RUNTIME_zlocality_onerestore) (cham_tasktype_t task);

void (*CHAMELEON_RUNTIME_clocality_allrestrict)(uint32_t device);
void (*CHAMELEON_RUNTIME_clocality_onerestrict)(cham_tasktype_t task, uint32_t device);
void (*CHAMELEON_RUNTIME_clocality_onerestore) (cham_tasktype_t task);

void (*CHAMELEON_RUNTIME_dlocality_allrestrict)(uint32_t device);
void (*CHAMELEON_RUNTIME_dlocality_onerestrict)(cham_tasktype_t task, uint32_t device);
void (*CHAMELEON_RUNTIME_dlocality_onerestore) (cham_tasktype_t task);

void (*CHAMELEON_RUNTIME_slocality_allrestrict)(uint32_t device);
void (*CHAMELEON_RUNTIME_slocality_onerestrict)(cham_tasktype_t task, uint32_t device);
void (*CHAMELEON_RUNTIME_slocality_onerestore) (cham_tasktype_t task);

/**
 * @}
 *
 * @name RUNTIME Profiling
 * @{
 */
void   RUNTIME_schedprofile_display ();
void   RUNTIME_kernelprofile_display();
double RUNTIME_get_time();

void  RUNTIME_iteration_push     (CHAM_context_t*, unsigned long iteration);
void  RUNTIME_iteration_pop      (CHAM_context_t*);

void RUNTIME_start_profiling();
void RUNTIME_stop_profiling();

void RUNTIME_start_stats();
void RUNTIME_stop_stats();

void RUNTIME_zdisplay_allprofile ();
void RUNTIME_zdisplay_oneprofile (cham_tasktype_t task);
void RUNTIME_cdisplay_allprofile ();
void RUNTIME_cdisplay_oneprofile (cham_tasktype_t task);
void RUNTIME_ddisplay_allprofile ();
void RUNTIME_ddisplay_oneprofile (cham_tasktype_t task);
void RUNTIME_sdisplay_allprofile ();
void RUNTIME_sdisplay_oneprofile (cham_tasktype_t task);

void   (*CHAMELEON_RUNTIME_schedprofile_display) ();
void   (*CHAMELEON_RUNTIME_kernelprofile_display)();
double (*CHAMELEON_RUNTIME_get_time)();

void  (*CHAMELEON_RUNTIME_iteration_push)     (CHAM_context_t*, unsigned long iteration);
void  (*CHAMELEON_RUNTIME_iteration_pop)      (CHAM_context_t*);

void (*CHAMELEON_RUNTIME_start_profiling)();
void (*CHAMELEON_RUNTIME_stop_profiling)();

void (*CHAMELEON_RUNTIME_start_stats)();
void (*CHAMELEON_RUNTIME_stop_stats)();

void (*CHAMELEON_RUNTIME_zdisplay_allprofile) ();
void (*CHAMELEON_RUNTIME_zdisplay_oneprofile) (cham_tasktype_t task);
void (*CHAMELEON_RUNTIME_cdisplay_allprofile) ();
void (*CHAMELEON_RUNTIME_cdisplay_oneprofile) (cham_tasktype_t task);
void (*CHAMELEON_RUNTIME_ddisplay_allprofile) ();
void (*CHAMELEON_RUNTIME_ddisplay_oneprofile) (cham_tasktype_t task);
void (*CHAMELEON_RUNTIME_sdisplay_allprofile) ();
void (*CHAMELEON_RUNTIME_sdisplay_oneprofile) (cham_tasktype_t task);

/**
 * @}
 */

END_C_DECLS

#endif /* _chameleon_runtime_h_ */
