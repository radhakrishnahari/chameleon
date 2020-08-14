/**
 *
 * @file chameleon_dlopen.h
 *
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon dlopen header
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2020-08-10
 *
 */
#ifndef _chameleon_dlopen_h_
#define _chameleon_dlopen_h_

/**
 *  Routines to handle dlopen
 */
#ifdef __cplusplus
extern "C" {
#endif

void chameleon_dlopen();
void chameleon_dlclose();

#ifdef __cplusplus
}
#endif

#endif /* _chameleon_dlopen_h_ */
