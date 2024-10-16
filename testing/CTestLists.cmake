#
# Check testing/
#
set(NP 2) # Amount of MPI processes
set(THREADS 2) # Amount of threads
set(N_GPUS 0) # Amount of graphic cards
set(TEST_CATEGORIES shm)
if (CHAMELEON_USE_MPI AND MPI_C_FOUND)
  set( TEST_CATEGORIES ${TEST_CATEGORIES} mpi )
endif()
if (CHAMELEON_USE_CUDA AND CUDA_FOUND)
  set(N_GPUS 0 1)
endif()
if (CHAMELEON_USE_HIP AND HIP_FOUND)
  set(N_GPUS 0 1)
endif()
if (CHAMELEON_SIMULATION)
  set(TEST_CATEGORIES simushm)
  if (CHAMELEON_USE_CUDA)
    set( TEST_CATEGORIES ${TEST_CATEGORIES} simugpu )
  endif()
endif()

set( SINGLE_PRECISIONS s c )
# list all tests that have a specific input file for single precision computations
set( SINGLE_TESTS gepdf_qdwh genm2 )

if (NOT CHAMELEON_SIMULATION)

  foreach(prec ${CHAMELEON_PRECISION})
    if ( ${prec} STREQUAL ds OR ${prec} STREQUAL zc )
      continue()
    endif()

    set (CMD ./chameleon_${prec}testing)

    #
    # Create the list of test based on precision and runtime
    #
    set( TESTS print lacpy lange lantr lansy plrnk )
    if ( ${prec} STREQUAL c OR ${prec} STREQUAL z )
      set( TESTS ${TESTS} lanhe )
    endif()
    set( TESTS ${TESTS} geadd tradd lascal gemm symm syrk syr2k trmm trsm )
    if ( ${prec} STREQUAL c OR ${prec} STREQUAL z )
      set( TESTS ${TESTS} hemm herk her2k )
    endif()
    set( TESTS ${TESTS} potrf potrs posv trtri lauum )
    if ( NOT CHAMELEON_SCHED_PARSEC )
      set( TESTS ${TESTS} potri poinv)
    endif()
    if ( ${prec} STREQUAL c OR ${prec} STREQUAL z )
      set( TESTS ${TESTS} sytrf sytrs sysv )
    endif()
    set( TESTS ${TESTS} getrf_nopiv getrs_nopiv gesv_nopiv geqrf gelqf geqrf_hqr gelqf_hqr )
    if ( ${prec} STREQUAL c OR ${prec} STREQUAL z )
      set( TESTS ${TESTS}
        ungqr     unglq     unmqr     unmlq
        ungqr_hqr unglq_hqr unmqr_hqr unmlq_hqr)
    else()
      set( TESTS ${TESTS}
        orgqr     orglq     ormqr     ormlq
        orgqr_hqr orglq_hqr ormqr_hqr ormlq_hqr)
    endif()
    #set( TESTS ${TESTS} geqrs     gelqs     )
    #set( TESTS ${TESTS} geqrs_hqr gelqs_hqr )
    set( TESTS ${TESTS} gels gels_hqr )
    set( TESTS ${TESTS} genm2 gepdf_qr gepdf_qdwh gesvd )
    set( TESTS ${TESTS} cesca gram )

    foreach( cat ${TEST_CATEGORIES} )
      foreach( gpus ${N_GPUS} )

        set( TESTSTMP ${TESTS} )

        if ( NOT ( ${gpus} EQUAL 0 ) )
          set( cat ${cat}_gpu )
          list( REMOVE_ITEM TESTSTMP gram lacpy lanhe lange lansy lantr lascal plrnk print )
        endif()

        if ( ${cat} STREQUAL "mpi" )
          set ( PREFIX mpiexec --bind-to none -n ${NP} )
          list( REMOVE_ITEM TESTSTMP gesvd )
        else()
          set ( PREFIX "" )
        endif()

        foreach( test ${TESTSTMP} )
          if ( ${test} IN_LIST SINGLE_TESTS AND ${prec} IN_LIST SINGLE_PRECISIONS )
            add_test( test_${cat}_${prec}${test} ${PREFIX} ${CMD} -c -t ${THREADS} -g ${gpus} -P 1 -f input/${test}_32.in )
          else()
            add_test( test_${cat}_${prec}${test} ${PREFIX} ${CMD} -c -t ${THREADS} -g ${gpus} -P 1 -f input/${test}.in )
          endif()
        endforeach()

        if ( CHAMELEON_SCHED_STARPU )
            add_test( test_${cat}_${prec}getrf_nopivpercol ${PREFIX} ${CMD} -c -t ${THREADS} -g ${gpus} -P 1 --diag=ChamUnit -f input/getrf_nopiv.in )
            set_tests_properties( test_${cat}_${prec}getrf_nopivpercol
                                  PROPERTIES ENVIRONMENT "CHAMELEON_GETRF_ALGO=nopivpercolumn;CHAMELEON_GETRF_BATCH_SIZE=0" )

            if ( HAVE_STARPU_NONE_NONZERO )
                add_test( test_${cat}_${prec}getrf_ppivpercol ${PREFIX} ${CMD} -c -t ${THREADS} -g ${gpus} -P 1 -f input/getrf_nopiv.in )
                set_tests_properties( test_${cat}_${prec}getrf_ppivpercol
                                    PROPERTIES ENVIRONMENT "CHAMELEON_GETRF_ALGO=ppivpercolumn;CHAMELEON_GETRF_BATCH_SIZE=0" )

                add_test( test_${cat}_${prec}getrf_ppivpercol_batch ${PREFIX} ${CMD} -c -t ${THREADS} -g ${gpus} -P 1 -f input/getrf_nopiv.in )
                set_tests_properties( test_${cat}_${prec}getrf_ppivpercol_batch
                                    PROPERTIES ENVIRONMENT "CHAMELEON_GETRF_ALGO=ppivpercolumn;CHAMELEON_GETRF_BATCH_SIZE=6" )

                add_test( test_${cat}_${prec}getrf_ppiv ${PREFIX} ${CMD} -c -t ${THREADS} -g ${gpus} -P 1 -f input/getrf.in )
                set_tests_properties( test_${cat}_${prec}getrf_ppiv
                                    PROPERTIES ENVIRONMENT "CHAMELEON_GETRF_ALGO=ppiv;CHAMELEON_GETRF_BATCH_SIZE=0" )

                add_test( test_${cat}_${prec}getrf_ppiv_batch ${PREFIX} ${CMD} -c -t ${THREADS} -g ${gpus} -P 1 -f input/getrf.in )
                set_tests_properties( test_${cat}_${prec}getrf_ppiv_batch
                                    PROPERTIES ENVIRONMENT "CHAMELEON_GETRF_ALGO=ppiv;CHAMELEON_GETRF_BATCH_SIZE=6" )
            endif()
        endif()

        list( REMOVE_ITEM TESTSTMP print gepdf_qr )

        if ( NOT (${cat} STREQUAL "mpi"))
          foreach( test ${TESTSTMP} )
            if ( ${test} IN_LIST SINGLE_TESTS AND ${prec} IN_LIST SINGLE_PRECISIONS )
              add_test( test_${cat}_${prec}${test}_std ${CMD} -c -t ${THREADS} -g ${gpus} -P 1 -f input/${test}_32.in --api=1 )
            else()
              add_test( test_${cat}_${prec}${test}_std ${CMD} -c -t ${THREADS} -g ${gpus} -P 1 -f input/${test}.in --api=1 )
            endif()
          endforeach()
        endif()
      endforeach()
    endforeach()
  endforeach()

else (NOT CHAMELEON_SIMULATION)

  # constraints for which we have perfmodels in simucore/perfmodels/
  set( TESTS "potrf")
  set(CHAMELEON_PRECISIONS_SIMU "s;d")
  set(TEST_CMD_simushm -t ${THREADS} -g 0)
  set(TEST_CMD_simugpu -t ${THREADS} -g 1)
  set(PLATFORMS "mirage;sirocco")
  set(BSIZE_mirage "320;960")
  set(BSIZE_sirocco "80;440;960;1440;1920")

  # loop over constraints
  foreach(cat ${TEST_CATEGORIES})
    foreach(prec ${CHAMELEON_PRECISIONS_SIMU})
      string(TOUPPER ${prec} PREC)
      foreach(test ${TESTS})
        if (CHAMELEON_PREC_${PREC})
          foreach(plat ${PLATFORMS})
            foreach(bsize ${BSIZE_${plat}})
              math(EXPR size "10 * ${bsize}")
              add_test(test_${cat}_${prec}${test}_${plat}_${bsize} ./chameleon_${prec}testing -o ${test} ${TEST_CMD_${cat}} -n ${size} -b ${bsize})
              set_tests_properties(test_${cat}_${prec}${test}_${plat}_${bsize} PROPERTIES
                                  ENVIRONMENT "STARPU_HOME=${CMAKE_SOURCE_DIR}/simucore/perfmodels;STARPU_HOSTNAME=${plat}"
                                  )
            endforeach()
          endforeach()
        endif()
      endforeach()
    endforeach()
  endforeach()

endif (NOT CHAMELEON_SIMULATION)
