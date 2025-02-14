/**
 *
 * @file getenv.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2025 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon auxiliary file to manage environment variables. Must be included after chameleon.h.
 *
 * @version 1.3.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @author Loris Lucido
 * @date 2025-01-29
 *
 */
#ifndef _chameleon_getenv_h_
#define _chameleon_getenv_h_

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Get environment variable
 */
#if defined(CHAMELEON_OS_WINDOWS)

static inline int
chameleon_setenv( const char *var, const char *value, int overwrite ) {
    return !(SetEnvironmentVariable( var, value ));
}

static inline char *
chameleon_getenv( const char *var ) {
    char *str;
    int len = 512;
    int rc;
    str = (char*)malloc( sizeof(char) * len );
    rc = GetEnvironmentVariable(var, str, len);
    if (rc == 0) {
        free(str);
        str = NULL;
    }
    return str;
}

static inline void
chameleon_cleanenv( char *str ) {
    if (str != NULL) free(str);
}

#else /* Other OS systems */

static inline int
chameleon_setenv( const char *var, const char *value, int overwrite ) {
    return setenv( var, value, overwrite );
}

static inline char *
chameleon_getenv( const char *var ) {
    return getenv( var );
}

static inline void
chameleon_cleanenv( char *str ) {
    (void)str;
}

#endif


static inline int
chameleon_env_is_set_to(char * str, char * value) {
    char * val;
    if ( (val = chameleon_getenv(str)) &&
         !strcmp(val, value))
        return 1;
    return 0;
}

static inline int
chameleon_env_on_off( char * str, int default_value ) {
    if ( chameleon_env_is_set_to(str, "1") ) {
        return CHAMELEON_TRUE;
    }
    if ( chameleon_env_is_set_to(str, "0") ) {
        return CHAMELEON_FALSE;
    }
    return default_value;
}

static inline int
chameleon_getenv_get_value_int( char * string, int default_value ) {
    extern int errno;
    long int ret;
    int      rc;
    char    *str = chameleon_getenv(string);

    if ( str == NULL ) {
        return default_value;
    }

    rc = sscanf( str, "%ld", &ret );
    if ( rc != 1 ) {
        if ( ( rc == EOF ) && ( errno != 0 ) ) {
            perror( "chameleon_getenv_get_value_int(sscanf)" );
        }
        else {
            fprintf( stderr, "%s: env variable %s expects an int value\n", __func__, string );
        }
        return default_value;
    }

    return (int)ret;
}

#ifdef __cplusplus
}
#endif

#endif /* _chameleon_getenv_h_ */
