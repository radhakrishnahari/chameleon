/**
 *
 * @file core.c
 *
 * @copyright 2025-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon gpucublas EZtrace 2.0 module
 *
 * @version 1.3.0
 * @author Brieuc Nicolas
 * @date 2025-07-09
 *
 */
#include <eztrace-core/eztrace_config.h>
#include <eztrace-instrumentation/pptrace.h>
#include <eztrace-lib/eztrace.h>
#include <eztrace-lib/eztrace_module.h>

/* set to 1 when all the hooks are set.
 * This is usefull in order to avoid recursive calls to mutex_lock for example
 */
static volatile int _coreblas_core_initialized = 0;

struct ezt_instrumented_function *PPTRACE_SYMBOL_LIST(coreblas_core);

extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(coreblas_core_z);
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(coreblas_core_s);
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(coreblas_core_d);
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(coreblas_core_c);
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(coreblas_core_zc);
extern struct ezt_instrumented_function PPTRACE_SYMBOL_LIST(coreblas_core_ds);

static void init_coreblas_core_()
{

    if (eztrace_autostart_enabled())
        eztrace_start();

    _coreblas_core_initialized = 1;
}

static void finalize_coreblas_core()
{
    _coreblas_core_initialized = 0;

    eztrace_stop();
}

static void _coreblas_core_init(void) __attribute__((constructor));
static void _coreblas_core_init(void)
{


    // TODO: protection pour chaque precision #if defined
    size_t sizez = sizeof(PPTRACE_SYMBOL_LIST(coreblas_core_z)) / sizeof(struct ezt_instrumented_function);
    size_t sizec = sizeof(PPTRACE_SYMBOL_LIST(coreblas_core_c)) / sizeof(struct ezt_instrumented_function);
    size_t sized = sizeof(PPTRACE_SYMBOL_LIST(coreblas_core_d)) / sizeof(struct ezt_instrumented_function);
    size_t sizes = sizeof(PPTRACE_SYMBOL_LIST(coreblas_core_s)) / sizeof(struct ezt_instrumented_function);
    size_t sizezc = sizeof(PPTRACE_SYMBOL_LIST(coreblas_core_zc)) / sizeof(struct ezt_instrumented_function);
    size_t sizeds = sizeof(PPTRACE_SYMBOL_LIST(coreblas_core_ds)) / sizeof(struct ezt_instrumented_function);

    /** All the -1 are because of FUNCTION_NONE applied when calling
     *  PPTRACE_END_INTERCEPT_FUNCTIONS macro
     */
    PPTRACE_SYMBOL_LIST(coreblas_core) = malloc(sizeof(struct ezt_instrumented_function) *
                                                ((sizez - 1) + (sizec - 1) + (sized - 1) + (sizes - 1) + (sizezc - 1) + sizeds));

    struct ezt_instrumented_function *current_list_start = PPTRACE_SYMBOL_LIST(coreblas_core);
    memcpy((void *)current_list_start,
           (void *)(&PPTRACE_SYMBOL_LIST(coreblas_core_z)),
           (sizez - 1) * sizeof(struct ezt_instrumented_function));
    current_list_start += sizez - 1;

    memcpy((void *)current_list_start,
           (void *)(&PPTRACE_SYMBOL_LIST(coreblas_core_c)),
           (sizec - 1) * sizeof(struct ezt_instrumented_function));
    current_list_start += sizec - 1;

    memcpy((void *)current_list_start,
           (void *)(&PPTRACE_SYMBOL_LIST(coreblas_core_d)),
           (sized - 1) * sizeof(struct ezt_instrumented_function));
    current_list_start += sized - 1;

    memcpy((void *)current_list_start,
           (void *)(&PPTRACE_SYMBOL_LIST(coreblas_core_s)),
           (sizes - 1) * sizeof(struct ezt_instrumented_function));
    current_list_start += sizes - 1;

    memcpy((void *)current_list_start,
           (void *)(&PPTRACE_SYMBOL_LIST(coreblas_core_zc)),
           (sizezc - 1) * sizeof(struct ezt_instrumented_function));
    current_list_start += sizezc - 1;

    memcpy((void *)current_list_start,
           (void *)(&PPTRACE_SYMBOL_LIST(coreblas_core_ds)),
           sizeds * sizeof(struct ezt_instrumented_function));

    eztrace_log(dbg_lvl_debug, "eztrace_coreblas_core constructor starts\n");
    EZT_REGISTER_MODULE(coreblas_core, "Module for the coreblas_core library",
                        init_coreblas_core_, finalize_coreblas_core);
    eztrace_log(dbg_lvl_debug, "eztrace_coreblas_core constructor ends\n");
}

static void _coreblas_core_destroy(void) __attribute__((destructor));
static void _coreblas_core_destroy(void)
{
    free(PPTRACE_SYMBOL_LIST(coreblas_core));
}
