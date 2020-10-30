/**
 *
 * @file chameleon_ztesting.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CHAMELEON_Complex64_t auxiliary testings routines
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cédric Castagnède
 * @author Lucas Barros de Assis
 * @date 2020-03-03
 * @precisions normal z -> c d s
 *
 */
#include "testings.h"
#if defined(CHAMELEON_HAVE_GETOPT_LONG)
#include <getopt.h>
#else
struct option;
#endif

/**
 * @brief Defines all the parameters of the testings
 */
static parameter_t parameters[] = {
    /* Name, helper, shname, flags, has_arg, psize, valtype, value, vallist, read, sprint */
    { "id", "Id of the run", 0, PARAM_OUTPUT, 0, 3, TestValInt, {0}, NULL, NULL, sprint_int },

    { NULL, "Options", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "help",     "Show this help",                           'h', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "check",    "Enable checking of the result",            'c', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "human",    "Enable human readable mode",               'H', PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "niter",    "Perform multiple iteration per test",      'l', PARAM_OPTION, 1, 0, TestValInt, {1}, NULL, pread_int, sprint_int },
    { "trace",    "Enable the trace generation",              -30, PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "nowarmup", "Disable the warmup run to load libraries", -31, PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "mtxfmt",   "Change the way the matrix is stored (0: global, 1: tiles, 2: OOC)", -32, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 1, 6, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "profile",  "Display the kernel profiling",             -33, PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "forcegpu", "Force kernels on GPU",                     -34, PARAM_OPTION, 0, 0, TestValInt, {0}, NULL, pread_int, sprint_int },

    { NULL, "Machine parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "threads", "Number of CPU workers per node",      't', PARAM_OPTION | PARAM_OUTPUT, 1, 7, TestValInt, {-1}, NULL, pread_int, sprint_int },
    { "gpus",    "Number of GPU workers per node",      'g', PARAM_OPTION | PARAM_OUTPUT, 1, 4, TestValInt, { 0}, NULL, pread_int, sprint_int },
    { "P",       "Rows (P) in the PxQ process grid",    'P', PARAM_OPTION | PARAM_OUTPUT, 1, 2, TestValInt, { 1}, NULL, pread_int, sprint_int },
    { "Q",       "Columns (Q) in the PxQ process grid", 'Q', PARAM_OUTPUT,                1, 2, TestValInt, { 1}, NULL, pread_int, sprint_int },

    { NULL, "Main input parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "op",   "Operation to test/time ('-o help' to get the list of functions)", 'o', PARAM_OPTION | PARAM_OUTPUT, 1, 1, TestString, {0}, NULL, pread_string, sprint_string },
    { "file", "Input file",                                                      'f', PARAM_OPTION,                1, 1, TestString, {0}, NULL, pread_string, sprint_string },

    { NULL, "Matrix definition parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "m",    "Dimension M of the operation",    'm', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "n",    "Dimension N of the operation",    'n', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "k",    "Dimension K of the operation",    'k', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "nrhs", "Dimension NRHS of the operation", 'r', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },

    { "nb", "Tile size nb",       'b', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "ib", "Inner tile size ib", 'i', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 2, TestValInt, {0}, NULL, pread_int, sprint_int },

    { "lda", "Leading dimension of the matrix A", 'A', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "ldb", "Leading dimension of the matrix B", 'B', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "ldc", "Leading dimension of the matrix C", 'C', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestValInt, {0}, NULL, pread_int, sprint_int },

    { "seedA", "Seed for the matrix A random generation", 'X', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 11, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "seedB", "Seed for the matrix B random generation", 'Y', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 11, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "seedC", "Seed for the matrix C random generation", 'Z', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 11, TestValInt, {0}, NULL, pread_int, sprint_int },

    { NULL, "Matrix generation numerical parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "bump",  "Bump value to make a matrix diagonal dominant",           'z', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 13, TestValComplex64, {0}, NULL, pread_complex64, sprint_complex64 },
    { "mode",  "Mode that specifies the eigen/singular values in xlatms", -30, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2,  4, TestValInt,       {0}, NULL, pread_int,       sprint_int    },
    { "cond",  "Conditional number of the matrix used by xlatms",         -31, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 13, TestValDouble,    {0}, NULL, pread_double,    sprint_double    },

    { NULL, "Operation specific parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "trans",  "Value of the trans parameter",  -11, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 9, TestTrans,    {0}, NULL, pread_trans, sprint_trans },
    { "transA", "Value of the transA parameter", -12, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 9, TestTrans,    {0}, NULL, pread_trans, sprint_trans },
    { "transB", "Value of the transB parameter", -13, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 9, TestTrans,    {0}, NULL, pread_trans, sprint_trans },
    { "uplo",   "Value of the uplo parameter",   -14, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 7, TestUplo,     {0}, NULL, pread_uplo,  sprint_uplo  },
    { "diag",   "Value of the diag parameter",   -15, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 7, TestDiag,     {0}, NULL, pread_diag,  sprint_diag  },
    { "side",   "Value of the side parameter",   -16, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 5, TestSide,     {0}, NULL, pread_side,  sprint_side  },
    { "norm",   "Value of the norm parameter",   -17, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 4, TestNormtype, {0}, NULL, pread_norm,  sprint_norm  },

    { NULL, "Operation specific scalar", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "alpha", "Value of the scalar alpha",                       'x', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 13, TestValComplex64, {0}, NULL, pread_complex64, sprint_complex64 },
    { "beta",  "Value of the scalar beta",                        'y', PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 13, TestValComplex64, {0}, NULL, pread_complex64, sprint_complex64 },

    { NULL, "QR/LQ parameters", 0, PARAM_OPTION, 0, 0, 0, {0}, NULL, NULL, NULL },
    { "qra",    "Size of TS domain (=RH for householder trees)",           -20, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "qrp",    "Size of high level tree for distributed",                 -21, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 3, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "llvl",   "Tree used for low level reduction insides nodes",         -22, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 4, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "hlvl",   "Tree used for high level reduction between nodes",        -23, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 4, TestValInt, {0}, NULL, pread_int, sprint_int },
    { "domino", "Enable/Disable the domino between upper and lower trees", -24, PARAM_OPTION | PARAM_INPUT | PARAM_OUTPUT, 2, 6, TestValInt, {0}, NULL, pread_int, sprint_int },

    { "time",          "Time in s",                             1000, PARAM_OUTPUT, 2, 13, TestValFixdbl, {0}, NULL, pread_fixdbl, sprint_fixdbl },
    { "gflops",        "GFlop/s",                               1001, PARAM_OUTPUT, 2, 13, TestValFixdbl, {0}, NULL, pread_fixdbl, sprint_fixdbl },
    { "RETURN",        "Result of the testing: SUCCESS/FAILED", 1002, PARAM_OUTPUT, 2,  7, TestValInt,    {0}, NULL, pread_int,    sprint_check  },
    { "||Ax-b||",      "Norm of the residual",                  1003, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||A-fact(A)||", "Norm of the residual",                  1004, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||A||",         "Norm of the matrix A",                  1005, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||B||",         "Norm of the matrix B",                  1006, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||C||",         "Norm of the matrix C",                  1007, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||R||",         "Residual norm",                         1008, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||b||",         "Norm of the vector b",                  1009, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||x||",         "Norm of the vector x",                  1010, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||Ax-b||/N/eps/(||A||||x||+||b||", "",                   1011, PARAM_OUTPUT, 2, 22, TestValDouble, {0}, NULL, pread_double, sprint_double },
    { "||I-QQ'||",     "Orthonormality of Q",                   1012, PARAM_OUTPUT, 2, 13, TestValDouble, {0}, NULL, pread_double, sprint_double },
};

#define STR_MAX_LENGTH 256

void print_usage( const char* prog_name )
{
    int rank;
    if ( CHAMELEON_Initialized() ) {
        rank = CHAMELEON_Comm_rank();
    }
    else {
        rank = 0;
    }

    if (rank == 0) {
        parameter_t *param = parameters;
        int i, nbparams = sizeof( parameters ) / sizeof( parameter_t );
        printf( "Usage:\n"
                "  %s -o|--op operation_name [options]\n"
                "  %s -f|--file input_file [options]\n",
                prog_name, prog_name );

        for (i=0; i<nbparams; i++, param++) {
            char str[STR_MAX_LENGTH];

            /* This is not an option, we skip it */
            if ( !(param->flags & PARAM_OPTION) ) {
                continue;
            }

            /* This is an option header */
            if ( param->name == NULL ) {
                printf( "\n  %s:\n", param->helper );
                continue;
            }

            if ( param->shname > 0 ) {
                snprintf( str, STR_MAX_LENGTH, "-%c, --%s",
                          param->shname, param->name );
            }
            else {
                snprintf( str, STR_MAX_LENGTH, "    --%s",
                          param->name );
            }

            /* If an argument is needed, add " x" */
            if ( param->has_arg > 0 ) {
                int len = strlen(str);
                assert( len < (STR_MAX_LENGTH-3) );

                str[ len   ] = ' ';
                str[ len+1 ] = 'x';
                str[ len+2 ] = '\0';
            }
            printf( "    %-23s %s\n",
                    str, param->helper );
        }

        printf( "\n"
                "For example: %s -H -o gemm -t 2 -m 2000 -n 2000 -k 2000 -b 200\n"
                "  will run one gemm with three matrices of size 2000x2000 each and a tile size of 200.\n"
                "  The output will be in the human readable format\n"
                "\n", prog_name );
    }
}

/**
 * @brief List of all the testings available.
 */
static testing_t *testings = NULL;

/**
 * @brief Function to register a new testing
 */
void
testing_register( testing_t *test )
{
    assert( test->next == NULL );
    test->next = testings;
    testings = test;
}

/**
 * @brief Get the testing structure associated to a test
 */
testing_t *
testing_gettest( const char *prog_name,
                 const char *func_name )
{
    testing_t *test = testings;
    int rank = CHAMELEON_Comm_rank();

    if ( func_name == NULL ) {
        print_usage( prog_name );
        exit(1);
    }

    while( test != NULL ) {
        /* Check the name with and without the precision */
        if ( (strcasecmp( func_name, test->name     ) == 0) ||
             (strcasecmp( func_name, test->name + 1 ) == 0) )
        {
            break;
        }
        test = test->next;
    }

    if ( test == NULL ) {
        if ( rank == 0 ) {
            printf( "Operation %s not found\n", func_name );
            printf( "The available operations are:\n" );
            test = testings;
            while( test != NULL ) {
                printf( "  %-10s %s\n", test->name, test->helper );
                test = test->next;
            }
        }
        exit(1);
    }

    return test;
}

int
parameters_compute_q( int p )
{
    parameter_t *param;
    int np = CHAMELEON_Comm_size();

    if ( (np % p) != 0 ) {
        fprintf( stderr, "ERROR: The number of processes (%d) must be a multiple of P (%d)\n", np, p );
        exit(1);
    }

    param = parameters_get( 'Q' );
    param->value.ival = np / p;
    return param->value.ival;
}

void
parameters_getopt_init( char           *optstring,
                        struct option **longopts )
{
    parameter_t *param = parameters;
    int i, nbparams = sizeof( parameters ) / sizeof( parameter_t );
    int nboptions = 0;
    int strpos = 0;

    for ( i=0; i<nbparams; i++, param++ ) {
        /* This is not an option, we skip it */
        if ( !(param->flags & PARAM_OPTION) ||
             (param->name == NULL) )
        {
            continue;
        }

        nboptions++;

        if ( param->shname < 0 ) {
            continue;
        }

        optstring[strpos] = param->shname;
        strpos++;
        assert( strpos < STR_MAX_LENGTH );

        if ( param->has_arg > 0 ) {
            optstring[strpos] = ':';
            strpos++;
            assert( strpos < STR_MAX_LENGTH );
        }
    }
    optstring[strpos] = '\0';

    /* Now, let's generate the long opt if needed */
#if defined(CHAMELEON_HAVE_GETOPT_LONG)
    if ( longopts != NULL ) {
        struct option *opt;
        *longopts = calloc( nboptions+1, sizeof( struct option ) );

        opt = *longopts;
        param = parameters;

        for ( i=0; i<nboptions; i++, opt++, param++ ) {

            /* Look for a valid option */
            while ( !(param->flags & PARAM_OPTION) ||
                    (param->name == NULL) )
            {
                param++;
            }

            opt->name    = param->name;
            opt->has_arg = ( param->has_arg > 0 ) ? 1 : 0;
            opt->flag    = NULL;
            opt->val     = param->shname;
        }
    }
#endif
}

parameter_t *
parameters_get( int shname )
{
    parameter_t *param = parameters;
    int i, nbparams = sizeof( parameters ) / sizeof( parameter_t );

    for ( i=0; i<nbparams; i++, param++ ) {
        /* This is not an option, we skip it */
        if ( param->name == NULL ) {
            continue;
        }

        if ( shname == param->shname ) {
            return param;
        }
    }

    fprintf( stderr, "parameters_get could not find parameter %d(%c)\n", shname, shname );
    return NULL;
}

int
parameters_getvalue_int( const char *name )
{
    parameter_t *param = parameters;
    int i, nbparams = sizeof( parameters ) / sizeof( parameter_t );

    for ( i=0; i<nbparams; i++, param++ ) {
        /* This is not an option, we skip it */
        if ( param->name == NULL ) {
            continue;
        }

        if ( strcasecmp( name, param->name ) != 0 ) {
            continue;
        }

        if ( param->has_arg > 1 ) {
            fprintf( stderr, "parameters_getvalue_int should not be called with parameter %s\n", name );
            return -1;
        }

        if ( param->valtype != TestValInt ) {
            fprintf( stderr, "parameters_getvalue_int has been called with a non integer parameter (%s)\n", name );
            return -1;
        }

        return param->value.ival;
    }

    fprintf( stderr, "parameters_getvalue_int could not find parameter %s\n", name );
    return -1;
}

char *
parameters_getvalue_str( const char *name )
{
    parameter_t *param = parameters;
    int i, nbparams = sizeof( parameters ) / sizeof( parameter_t );

    for ( i=0; i<nbparams; i++, param++ ) {
        /* This is not an option, we skip it */
        if ( param->name == NULL ) {
            continue;
        }

        if ( strcasecmp( name, param->name ) != 0 ) {
            continue;
        }

        if ( param->has_arg > 1 ) {
            fprintf( stderr, "parameters_getvalue_str should not be called with parameter %s\n", name );
            return NULL;
        }

        if ( param->valtype != TestString ) {
            fprintf( stderr, "parameters_getvalue_str has been called with a non string parameter (%s)\n", name );
            return NULL;
        }

        return param->value.str;
    }

    fprintf( stderr, "parameters_getvalue_str could not find parameter %s\n", name );
    return NULL;
}

parameter_t *
parameters_getbyname( const char *name )
{
    parameter_t *param = parameters;
    int i, nbparams = sizeof( parameters ) / sizeof( parameter_t );

    for ( i=0; i<nbparams; i++, param++ ) {
        /* This is not an option, we skip it */
        if ( param->name == NULL ) {
            continue;
        }

        if ( strcasecmp( name, param->name ) != 0 ) {
            continue;
        }

        /* if ( param->has_arg < 2 ) { */
        /*     fprintf( stderr, "parameters_getbyname should not be called with parameter %s\n", name ); */
        /*     return NULL; */
        /* } */

        return param;
    }

    fprintf( stderr, "parameters_getbyname could not find parameter %s\n", name );
    return NULL;
}

void parameters_parser( int argc, char **argv )
{
    int opt;
    char optstring[STR_MAX_LENGTH];
    struct option *longopts = NULL;
    parameter_t *param;

    parameters_getopt_init( optstring, &longopts );

#if defined(CHAMELEON_HAVE_GETOPT_LONG)
    while ((opt = getopt_long(argc, argv, optstring, longopts, NULL)) != -1)
#else
    while ((opt = getopt(argc, argv, optstring)) != -1)
#endif
    {
        switch(opt) {
        case 'h':
            print_usage(argv[0]);
            exit(0);

        case '?': /* error from getopt[_long] */
            exit(1);
            break;

        default:
            param = parameters_get( opt );
            if ( param == NULL ) {
                print_usage(argv[0]);
                exit(1);
            }
            parameters_addvalues( param, optarg );
        }
    }

    if ( longopts != NULL ) {
        free( longopts );
    }
}

void
parameters_destroy()
{
    parameter_t *param = parameters;
    int i, nbparams = sizeof( parameters ) / sizeof( parameter_t );
    vallist_t *current, *next;

    for ( i=0; i<nbparams; i++, param++ ) {
        /* This is not an option, we skip it */
        if ( param->has_arg < 2 ) {
            continue;
        }

        current = param->vallist;
        while ( current != NULL )
        {
            next = current->next;
            free( current );
            current = next;
        }
    }
    return;
}

int main (int argc, char **argv) {

    int ncores, ngpus, human, check, i, niter;
    int trace, nowarmup, profile, forcegpu;
    int rc, info = 0;
    int run_id = 0;
    char *func_name;
    char *input_file;
    run_list_t *runlist;
    testing_t * test;
    run_list_elt_t *run, *next;

    /* Reads the arguments from command line */
    parameters_parser( argc, argv );
    input_file = parameters_getvalue_str( "file" );
    if ( input_file != NULL ) {
        parameters_read_file( input_file );
        free(input_file);
    }
    ncores    = parameters_getvalue_int( "threads"  );
    ngpus     = parameters_getvalue_int( "gpus"     );
    check     = parameters_getvalue_int( "check"    );
    human     = parameters_getvalue_int( "human"    );
    func_name = parameters_getvalue_str( "op"       );
    niter     = parameters_getvalue_int( "niter"    );
    trace     = parameters_getvalue_int( "trace"    );
    nowarmup  = parameters_getvalue_int( "nowarmup" );
    profile   = parameters_getvalue_int( "profile" );
    forcegpu  = parameters_getvalue_int( "forcegpu" );

    CHAMELEON_Init( ncores, ngpus );

    /* Set ncores to the right value */
    if ( ncores == -1 ) {
        parameter_t *param;
        param = parameters_get( 't' );
        param->value.ival = CHAMELEON_GetThreadNbr();
    }

    /* Binds the right function to be called and builds the parameters combinations */
    test = testing_gettest( argv[0], func_name );
    free(func_name);

    /* Generate the cartesian product of the parameters */
    runlist = run_list_generate( test->params );

    /* Executes the tests */
    run_print_header( test, check, human );
    run = runlist->head;

    /* Force all possible kernels on GPU */
    if ( forcegpu ) {
        if ( ngpus == 0 ) {
            fprintf( stderr,
                     "--forcegpu can't be enable without GPU (-g 0).\n"
                     "  Please specify a larger number of GPU or disable this option\n" );
            return EXIT_FAILURE;
        }
        RUNTIME_zlocality_allrestrict( RUNTIME_CUDA );
    }

    /* Warmup */
    if ( !nowarmup ) {
        run_arg_list_t copy = run_arg_list_copy( &(run->args) );
        test->fptr( &copy, check );
        run_arg_list_destroy( &copy );
    }

    /* Start kernel statistics */
    if ( profile ) {
        CHAMELEON_Enable( CHAMELEON_KERNELPROFILE_MODE );
    }

    /* Start tracing */
    if ( trace ) {
        CHAMELEON_Enable( CHAMELEON_PROFILING_MODE );
    }

    /* Perform all runs */
    while ( run != NULL ) {
        for(i=0; i<niter; i++) {
            run_arg_list_t copy = run_arg_list_copy( &(run->args) );
            rc = test->fptr( &copy, check );

            /* If rc < 0, we skipped the test */
            if ( rc >= 0 ) {
                run_arg_add_int( &copy, "RETURN", rc );
                run_print_line( test, &copy, check, human, run_id );
                run_id++;
                info += rc;
            }
            run_arg_list_destroy( &copy );
        }

        /* Move to next run */
        next = run->next;
        run_list_destroy( run );
        run = next;
    }

    /* Stop tracing */
    if ( trace ) {
        CHAMELEON_Disable( CHAMELEON_PROFILING_MODE );
    }

    /* Stop kernel statistics and display results */
    if ( profile ) {
        CHAMELEON_Disable( CHAMELEON_KERNELPROFILE_MODE );
        RUNTIME_kernelprofile_display();
    }
    free( runlist );

    CHAMELEON_Finalize();
    parameters_destroy();

    return info;
}
