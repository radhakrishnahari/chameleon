###
#
# @file CMakeLists.txt
#
# @copyright 2009-2014 The University of Tennessee and The University of
#                      Tennessee Research Foundation. All rights reserved.
# @copyright 2012-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
###
#
#  @project CHAMELEON
#  CHAMELEON is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
# @version 1.0.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 2020-03-03
#
###
cmake_minimum_required(VERSION 3.3)
project(CHAMELEON Fortran C CXX)

# directly make an error if in-source build
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(FATAL_ERROR "In-source builds are not allowed.\n"
   "Please create a build directory first and execute cmake configuration from "
   "this directory. Example: mkdir build && cd build && cmake ..")
endif()

# set project version number
set(CHAMELEON_VERSION_MAJOR 1)
set(CHAMELEON_VERSION_MINOR 0)
set(CHAMELEON_VERSION_MICRO 0)

set(CHAMELEON_CMAKE_DIR "" CACHE PATH "Directory of CHAMELEON CMake modules, can be external to the project")

# Add extra cmake module path and initialize morse cmake modules
# --------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)
if(CHAMELEON_CMAKE_DIR)
    list(APPEND CMAKE_MODULE_PATH "${CHAMELEON_CMAKE_DIR}/cmake_modules/morse_cmake/modules")
    set(CHAMELEON_CMAKE_MODULE_PATH ${CHAMELEON_CMAKE_DIR}/cmake_modules/morse_cmake/modules )
elseif(EXISTS "${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake")
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake/modules")
    set(CHAMELEON_CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake/modules )
else()
  message(FATAL_ERROR "Chameleon CMake system relies on morse_cmake modules developed here: "
    "https://gitlab.inria.fr/solverstack/morse_cmake. Please set CHAMELEON_CMAKE_DIR to this source "
    "directory.")
endif()

include(MorseInit)
include(GenPkgConfig)

# Ensure that the fortran compiler and c compiler specified are compatible
# ------------------------------------------------------------------------
include(FortranCInterface)
FortranCInterface_VERIFY()
FortranCInterface_HEADER(${CMAKE_CURRENT_BINARY_DIR}/include/chameleon/mangling.h
                         MACRO_NAMESPACE "CHAMELEON_")

#############################################
#                                           #
#        Compilation of CHAMELEON           #
#                                           #
#############################################


###############################################################################
# Parameters/Options #
######################

# Set the RPATH config
# --------------------

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# the RPATH to be used when installing
list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Misc options
# ------------
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Define precision supported by CHAMELEON
# -----------------------------------------
set( RP_CHAMELEON_DICTIONNARY ${CMAKE_SOURCE_DIR}/cmake_modules/local_subs.py )
set( RP_CHAMELEON_PRECISIONS  "s;d;c;z" )
include(RulesPrecisions)

# Options to select the runtime
# -----------------------------

# Create a list of possible runtime
set(CHAMELEON_SCHED_list OPENMP PARSEC QUARK STARPU
  CACHE INTERNAL "List of available runtimes" )
set(CHAMELEON_SCHED "" CACHE STRING "Runtime systems enabled" FORCE)
set(CHAMELEON_SCHED_OPENMP ON  CACHE BOOL "Enable OpenMP scheduler")
set(CHAMELEON_SCHED_PARSEC OFF CACHE BOOL "Enable PaRSEC scheduler")
set(CHAMELEON_SCHED_QUARK  OFF CACHE BOOL "Enable Quark scheduler")
set(CHAMELEON_SCHED_STARPU OFF CACHE BOOL "Enable StarPU scheduler")

# Count number of runtime set to ON
math(EXPR CHAMELEON_SCHED_NUMBER 0)
foreach (runtime ${CHAMELEON_SCHED_list})
    if (CHAMELEON_SCHED_${runtime})
        math(EXPR CHAMELEON_SCHED_NUMBER "${CHAMELEON_SCHED_NUMBER}+1")
        message("-- ${BoldGreen}CHAMELEON_SCHED_${runtime} is ON${ColourReset}")
    endif()
endforeach()
message("-- ${BoldGreen}CHAMELEON_SCHED_NUMBER = ${CHAMELEON_SCHED_NUMBER}${ColourReset}")
set(CHAMELEON_SCHED_DLOPEN OFF)
if (${CHAMELEON_SCHED_NUMBER} GREATER 1)
    if (NOT BUILD_SHARED_LIBS)
        message(WARNING "Number of activated scheduler is greater than 1 but BUILD_SHARED_LIBS is OFF. "
        "BUILD_SHARED_LIBS is set to ON in order to use schedulers as plugins to be loaded at runtime (dlopen).")
    endif()
    set(BUILD_SHARED_LIBS ON)
    set(CHAMELEON_SCHED_DLOPEN ON)
endif()

# Additional options
# ------------------
option(CHAMELEON_USE_MPI "Enable distributed memory through MPI" OFF)
if (NOT CHAMELEON_USE_MPI)
    message("-- ${BoldGreen}CHAMELEON_USE_MPI is set to OFF, turn it ON to use MPI (unsupported by OpenMP and Quark)${ColourReset}")
endif()

cmake_dependent_option(CHAMELEON_USE_MIGRATE
  "This options enables the data migration in QR algorithms" ON
  "CHAMELEON_USE_MPI" OFF)

cmake_dependent_option(CHAMELEON_USE_MPI_DATATYPES
  "Enable MPI datatypes whenever supported by the runtime" ON
  "CHAMELEON_USE_MPI" OFF)

# Enable Cuda kernels if StarPU (only if StarPU or PaRSEC runtime is enabled)
# TODO: Default should be changed to ON/OFF when it will be ok
option(CHAMELEON_USE_CUDA "Enable CUDA kernels" OFF)
if (NOT CHAMELEON_USE_CUDA)
    message("-- ${BoldGreen}CHAMELEON_USE_CUDA is set to OFF, turn it ON to use CUDA (unsupported by OpenMP and Quark)${ColourReset}")
endif()

option(CHAMELEON_RUNTIME_SYNC "Enable synchronous task submission when available to debug the code without parallelism" OFF)
if (CHAMELEON_RUNTIME_SYNC)
    message("-- ${BoldGreen}CHAMELEON_RUNTIME_SYNC is set to ON, turn it OFF to avoid synchronisation in the tasks submission${ColourReset}")
endif()

# Options to enable/disable doc, examples, and testings
# -----------------------------------------------------
option(CHAMELEON_ENABLE_DOC      "Enable documentation build"  OFF)
if (CHAMELEON_ENABLE_DOC)
    message("-- ${BoldGreen}CHAMELEON_ENABLE_DOC is set to ON, turn it OFF to avoid building docs${ColourReset}")
endif()
option(CHAMELEON_ENABLE_EXAMPLE  "Enable examples build"       ON)
if (CHAMELEON_ENABLE_EXAMPLE)
    message("-- ${BoldGreen}CHAMELEON_ENABLE_EXAMPLE is set to ON, turn it OFF to avoid building examples${ColourReset}")
endif()
option(CHAMELEON_ENABLE_TESTING  "Enable testings build"       ON)
if (CHAMELEON_ENABLE_TESTING)
    message("-- ${BoldGreen}CHAMELEON_ENABLE_TESTING is set to ON, turn it OFF to avoid building testing${ColourReset}")
endif()

# Option to activate or not simulation mode (use Simgrid through StarPU)
# ----------------------------------------------------------------------
option(CHAMELEON_SIMULATION "Enable simulation mode using Simgrid through StarPU" OFF)
if(CHAMELEON_SIMULATION)
    set(CHAMELEON_SCHED_OPENMP OFF  CACHE BOOL "Enable OpenMP scheduler" FORCE)
    set(CHAMELEON_SCHED_PARSEC OFF CACHE BOOL "Enable PaRSEC scheduler" FORCE)
    set(CHAMELEON_SCHED_QUARK  OFF CACHE BOOL "Enable Quark scheduler" FORCE)
    set(CHAMELEON_SCHED_STARPU ON CACHE BOOL "Enable StarPU scheduler" FORCE)
    set(CHAMELEON_SCHED_DLOPEN OFF)
    message("-- ${BoldGreen}CHAMELEON_SIMULATION is set to ON"
    " StarPU (compiled with SimGrid) is the only scheduler compatible, we turn OFF others${ColourReset}")
else()
    message("-- ${BoldGreen}CHAMELEON_SIMULATION is set to OFF, turn it ON to use"
            " SIMULATION mode (only with StarPU compiled with SimGrid)${ColourReset}")
endif()
if(CHAMELEON_SCHED_STARPU)
    option(CHAMELEON_ENABLE_PRUNING_STATS "Enable pruning statistics" OFF)
    if (NOT CHAMELEON_ENABLE_PRUNING_STATS)
        message("-- ${BoldGreen}CHAMELEON_ENABLE_PRUNING_STATS is set to OFF, turn it ON to build pruning statistics${ColourReset}")
    endif()
    if(CHAMELEON_ENABLE_PRUNING_STATS)
        add_definitions(-DCHAMELEON_ENABLE_PRUNING_STATS)
    endif(CHAMELEON_ENABLE_PRUNING_STATS)

    cmake_dependent_option(CHAMELEON_SIMULATION_EXTENDED
      "Enable non GPU kernels to be simulated on GPU" OFF
      "CHAMELEON_SIMULATION" OFF)
endif()

# Initially we need to generate files for different precisions
# TODO: use this option to avoid generate multiple precisions each time we launch cmake
#option(CHAMELEON_GEN_PREC "Generate source files precisions" ON)
#------------------------------------------------------------------------------

###############################################################################
# Build dependency HQR library #
################################
add_subdirectory(hqr)
include_directories(hqr/include)

###############################################################################
# Look for dependencies #
#########################
set(CHAMELEON_DEP "")

# Check for Thread library
# ------------------------
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)
if( THREADS_FOUND AND CMAKE_THREAD_LIBS_INIT )
  libraries_absolute_path(CMAKE_THREAD_LIBS_INIT "")
  list(APPEND EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
endif ()

# Add math library to the list of extra
# it normally exists on all common systems provided with a C compiler
if(UNIX OR WIN32)
    find_library(
      M_LIBRARY
      NAMES m
      )
    mark_as_advanced(M_LIBRARY)
    if (M_LIBRARY)
      list(APPEND EXTRA_LIBRARIES "${M_LIBRARY}")
    else()
      message(FATAL_ERROR "Could NOT find libm on your system."
        " Are you sure to a have a C compiler installed?")
    endif()
endif()

# Try to find librt (libposix4 - POSIX.1b Realtime Extensions library)
# on Unix systems except Apple ones because it does not exist on it
if(UNIX AND NOT APPLE)
    find_library(
        RT_LIBRARY
        NAMES rt
        )
    mark_as_advanced(RT_LIBRARY)
    if (RT_LIBRARY)
        list(APPEND EXTRA_LIBRARIES "${RT_LIBRARY}")
    else()
        message(FATAL_ERROR "Could NOT find librt on your system")
    endif()
endif()

if(CHAMELEON_SCHED_DLOPEN AND CMAKE_DL_LIBS)
    list(APPEND EXTRA_LIBRARIES "${CMAKE_DL_LIBS}")
endif()

# If simulation we don't enter in kernel functions so that we don't need to
# link with concerned libraries
if(NOT CHAMELEON_SIMULATION)
    # CHAMELEON depends on CBLAS
    #---------------------------
    find_package(CBLAS REQUIRED)
    if(CBLAS_FOUND)
        if (CBLAS_INCLUDE_DIRS_DEP)
            include_directories(${CBLAS_INCLUDE_DIRS_DEP})
        endif()
        if(CBLAS_CFLAGS_OTHER_DEP)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CBLAS_CFLAGS_OTHER_DEP}")
            string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
        endif()
        if(CBLAS_LDFLAGS_OTHER_DEP)
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CBLAS_LDFLAGS_OTHER_DEP}")
            string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
        endif()
        if(CBLAS_LIBRARY_DIRS_DEP)
            list(APPEND CMAKE_INSTALL_RPATH "${CBLAS_LIBRARY_DIRS_DEP}")
        endif()
        if (CBLAS_LIBRARIES_DEP)
            list(INSERT CHAMELEON_DEP 0 ${CBLAS_LIBRARIES_DEP})
        endif()
    else()
        if(CHAMELEON_VERBOSE_FIND_PACKAGE)
            if(CBLAS_STANDALONE OR NOT CBLAS_WORKS)
                if (NOT CBLAS_cblas.h_DIRS)
                    Print_Find_Header_Status(cblas cblas.h)
                endif ()
                if (NOT CBLAS_cblas_LIBRARY)
                    Print_Find_Library_Status(cblas libcblas)
                endif ()
            endif()
        else()
            message(WARNING "CBLAS library has not been found and CHAMELEON_VERBOSE_FIND_PACKAGE is set to OFF."
            " Try to activate CHAMELEON_VERBOSE_FIND_PACKAGE option (-DCHAMELEON_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
        endif()
        message(FATAL_ERROR "A CBLAS library is required but has not been found")
    endif()

    # CHAMELEON depends on LAPACKE
    #-----------------------------
    find_package(LAPACKE REQUIRED)
    if(LAPACKE_FOUND)
        if (LAPACKE_INCLUDE_DIRS_DEP)
            include_directories(${LAPACKE_INCLUDE_DIRS_DEP})
        endif()
        if(LAPACKE_CFLAGS_OTHER_DEP)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${LAPACKE_CFLAGS_OTHER_DEP}")
            string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
        endif()
        if(LAPACKE_LDFLAGS_OTHER_DEP)
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACKE_LDFLAGS_OTHER_DEP}")
            string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
        endif()
        if(LAPACKE_LIBRARY_DIRS_DEP)
            list(APPEND CMAKE_INSTALL_RPATH "${LAPACKE_LIBRARY_DIRS_DEP}")
        endif()
        if (LAPACKE_LIBRARIES_DEP)
            list(INSERT CHAMELEON_DEP 0 ${LAPACKE_LIBRARIES_DEP})
        endif()
    else()
        if(CHAMELEON_VERBOSE_FIND_PACKAGE)
            if (LAPACKE_STANDALONE OR NOT LAPACKE_WORKS)
                if (NOT LAPACKE_lapacke.h_DIRS)
                    Print_Find_Header_Status(lapacke lapacke.h)
                endif ()
                if (NOT LAPACKE_lapacke_LIBRARY)
                    Print_Find_Library_Status(lapacke liblapacke)
                endif ()
            endif()
        else()
            message(WARNING "LAPACKE library has not been found and CHAMELEON_VERBOSE_FIND_PACKAGE is set to OFF."
            " Try to activate CHAMELEON_VERBOSE_FIND_PACKAGE option (-DCHAMELEON_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
        endif()
        message(FATAL_ERROR "A LAPACKE library is required but has not been found")
    endif()

    # CHAMELEON depends on CUDA/CUBLAS
    #----------------------------------
    if (CHAMELEON_USE_CUDA)

        find_package(CUDA REQUIRED)
        if (CUDA_FOUND)
            message("-- ${Blue}Add definition CHAMELEON_USE_CUDA"
            " - Activate CUDA in Chameleon${ColourReset}")
            libraries_absolute_path(CUDA_LIBRARIES "$ENV{CUDA_ROOT}/lib64")
            set(EXTRA_LIBRARIES_CUDA)
            if(CUDA_VERSION VERSION_LESS "3.0")
                set(CUDA_HOST_COMPILATION_CPP OFF)
            endif(CUDA_VERSION VERSION_LESS "3.0")
            set(CUDA_BUILD_EMULATION OFF)
            if (CUDA_INCLUDE_DIRS)
                include_directories(${CUDA_INCLUDE_DIRS})
            else()
                message(WARNING "CHAMELEON_USE_CUDA requires"
                "\n   CUDA_INCLUDE_DIRS to be found. Be sure you have"
                "\n   cuda headers with your distribution of CUDA.")
            endif()
            if (CUDA_LIBRARIES)
                set(saved_include "${CMAKE_REQUIRED_INCLUDES}")
                set(saved_libs "${CMAKE_REQUIRED_LIBRARIES}")
                set(CMAKE_REQUIRED_INCLUDES  "${CMAKE_REQUIRED_INCLUDES};${CUDA_INCLUDE_DIRS}")
                set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES};${CUDA_LIBRARIES}")
                if(CUDA_VERSION VERSION_LESS "4.0")
                    set(CUDA_HAVE_PEER_DEVICE_MEMORY_ACCESS 0)
                else()
                    check_function_exists(cuDeviceCanAccessPeer CUDA_HAVE_PEER_DEVICE_MEMORY_ACCESS)
                endif()
                set(CMAKE_REQUIRED_INCLUDES  "${saved_include}")
                set(CMAKE_REQUIRED_LIBRARIES "${saved_libs}")
                # Add cublas if found
                if (CUDA_CUBLAS_LIBRARIES)
                    message("-- ${Blue}Add definition CHAMELEON_USE_CUBLAS"
                    " - Use GPU kernels from cuBLAS${ColourReset}")
                    list(APPEND EXTRA_LIBRARIES_CUDA ${CUDA_CUBLAS_LIBRARIES})
                    list(APPEND EXTRA_LIBRARIES_CUDA ${CUDA_LIBRARIES}) # need CUDA libs to link.
                else()
                    message(FATAL_ERROR "CHAMELEON_USE_CUDA requires"
                    "\n   CUDA_CUBLAS_LIBRARIES to be found. Be sure you have"
                    "\n   libcublas with your distribution of CUDA.")
                endif()
            else()
                message(FATAL_ERROR "CHAMELEON_USE_CUDA requires"
                "\n   CUDA_LIBRARIES to be found. Be sure you have"
                "\n   libcuda with your distribution of CUDA.")
            endif()
        endif (CUDA_FOUND)

        list(APPEND EXTRA_LIBRARIES ${EXTRA_LIBRARIES_CUDA})
        foreach (_lib ${EXTRA_LIBRARIES_CUDA})
            if (EXISTS ${_lib})
                get_filename_component(libdir ${_lib} PATH)
                if (EXISTS ${libdir})
                    list(APPEND CMAKE_INSTALL_RPATH "${libdir}")
                endif()
            endif()
        endforeach()

    endif(CHAMELEON_USE_CUDA)

    # CHAMELEON depends on MPI
    #-------------------------
    if (CHAMELEON_USE_MPI)

        # allows to use an external mpi compilation by setting compilers with
        # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
        # at cmake configure
        if(NOT MPI_C_COMPILER)
            set(MPI_C_COMPILER mpicc)
        endif()
        find_package(MPI REQUIRED)

        if (MPI_C_FOUND)
            message("-- ${Blue}Add definition CHAMELEON_USE_MPI"
            " - Activate MPI in Chameleon${ColourReset}")
            list(APPEND EXTRA_LIBRARIES ${MPI_C_LIBRARIES} )
            include_directories( ${MPI_C_INCLUDE_PATH} )
            # tests for intel mpi
            #list(APPEND MPI_C_COMPILE_FLAGS "-mt_mpi")
            #list(APPEND MPI_COMPILE_FLAGS "-mt_mpi")
            if(MPI_C_COMPILE_FLAGS)
                set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
                string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
            endif()
            if(MPI_C_LINK_FLAGS)
                set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")
                string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
            endif()
        endif (MPI_C_FOUND)

        foreach (_lib ${MPI_C_LIBRARIES})
            if (EXISTS ${_lib})
                get_filename_component(libdir ${_lib} PATH)
                if (EXISTS ${libdir})
                    list(APPEND CMAKE_INSTALL_RPATH "${libdir}")
                endif()
            endif()
        endforeach()

    endif (CHAMELEON_USE_MPI)

else (NOT CHAMELEON_SIMULATION)

    message("-- ${BoldBlue}Simulation mode of CHAMELEON is activated (CHAMELEON_SIMULATION=ON)."
    "\n   This mode allows you to simulate execution of algorithms with StarPU compiled with SimGrid."
    "\n   To do so, we provide some perfmodels in the simucore/perfmodels/ directory of CHAMELEON sources."
    "\n   The resulting directory where to find perfmodels is: ${CMAKE_CURRENT_SOURCE_DIR}/simucore/perfmodels."
    "\n   To use these perfmodels, please set your STARPU_HOME environment variable to ${CMAKE_CURRENT_SOURCE_DIR}/simucore/perfmodels."
    "\n   Finally, you need to set your STARPU_HOSTNAME environment variable to the name of the machine to simulate."
    "\n   For example: STARPU_HOSTNAME=mirage."
    "\n   Note that only POTRF kernels (simple and double precision) on mirage machine (block sizes of 320 or 960) and sirocco (block sizes of 80, 440, 960, 1440, 1920) are available for now."
    "\n   Database of models is subject to change, it should be enrich in a near future.${ColourReset}")

    message("-- ${Blue}Add definition CHAMELEON_SIMULATION"
    " - Activate simulation mode in Chameleon (to use with StarPU+SimGrid)${ColourReset}")

    # Guard against mis-used simulation mode
    if(NOT DEFINED CHAMELEON_USE_CUDA)
        message(FATAL_ERROR "${BoldBlue}"
        "In simulation mode CHAMELEON_USE_CUDA should be set to"
        "   ensure that the user is aware of the version to be used. If Chameleon's"
        "   kernels are available for NVIDIA CUDA GPUs and if the according"
        "   perfmodels are available in STARPU_HOME then use CHAMELEON_USE_CUDA=ON"
        "   else set CHAMELEON_USE_CUDA=OFF."
        "${ColourReset}")
    endif()

    if (NOT CHAMELEON_SCHED_STARPU)
        message(FATAL_ERROR "Simulation mode of Chameleon only works with"
        "\n   the StarPU runtime. Please use CHAMELEON_SCHED_STARPU=ON.")
    endif()

    if (CHAMELEON_ENABLE_EXAMPLE)
        set(CHAMELEON_ENABLE_EXAMPLE OFF)
        message("-- ${BoldBlue}CHAMELEON_ENABLE_EXAMPLE is set to ON but we turn it OFF.")
    endif ()

    # Simulation mode: we depend on SimGrid
    find_package(SIMGRID REQUIRED)

    if (SIMGRID_INCLUDE_DIRS)
        include_directories(${SIMGRID_INCLUDE_DIRS})
    endif()

    # CHAMELEON depends on MPI
    #-------------------------
    if (CHAMELEON_USE_MPI)

        # allows to use an external mpi compilation by setting compilers with
        # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
        # at cmake configure
        if(NOT MPI_C_COMPILER)
            set(MPI_C_COMPILER "${SIMGRID_DIR_FOUND}/bin/smpicc")
        endif()
        if(NOT MPI_C_LIBRARIES)
            set(MPI_C_LIBRARIES "${SIMGRID_LIBRARIES}")
        endif()
        if(NOT MPI_C_INCLUDE_PATH)
            set(MPI_C_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}")
            list(APPEND MPI_C_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}/smpi")
        endif()
        if(NOT MPI_CXX_COMPILER)
            set(MPI_CXX_COMPILER "${SIMGRID_DIR_FOUND}/bin/smpicxx")
        endif()
        if(NOT MPI_CXX_LIBRARIES)
            set(MPI_CXX_LIBRARIES "${SIMGRID_LIBRARIES}")
        endif()
        if(NOT MPI_CXX_INCLUDE_PATH)
            set(MPI_CXX_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}")
            list(APPEND MPI_CXX_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}/smpi")
        endif()
        if(NOT MPI_Fortran_COMPILER)
            set(MPI_Fortran_COMPILER "${SIMGRID_DIR_FOUND}/bin/smpif90")
        endif()
        if(NOT MPI_Fortran_LIBRARIES)
            set(MPI_Fortran_LIBRARIES "${SIMGRID_LIBRARIES}")
        endif()
        if(NOT MPI_Fortran_INCLUDE_PATH)
            set(MPI_Fortran_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}")
            list(APPEND MPI_Fortran_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}/smpi")
        endif()
        find_package(MPI REQUIRED)

        if (MPI_C_FOUND)
            message("-- ${Blue}Add definition CHAMELEON_USE_MPI"
            " - Activate MPI in Chameleon${ColourReset}")
            list(APPEND EXTRA_LIBRARIES ${MPI_C_LIBRARIES} )
            include_directories( ${MPI_C_INCLUDE_PATH} )
            if(MPI_C_COMPILE_FLAGS)
                set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
                string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
            endif()
            if(MPI_C_LINK_FLAGS)
                set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")
                string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
            endif()
        endif (MPI_C_FOUND)

    endif (CHAMELEON_USE_MPI)

    foreach (_lib ${SIMGRID_LIBRARIES} ${MPI_C_LIBRARIES})
        if (EXISTS ${_lib})
            get_filename_component(libdir ${_lib} PATH)
            if (EXISTS ${libdir})
                list(APPEND CMAKE_INSTALL_RPATH "${libdir}")
            endif()
        endif()
    endforeach()

endif(NOT CHAMELEON_SIMULATION)

# getopt
unset(CMAKE_REQUIRED_LIBRARIES)
check_include_files(getopt.h CHAMELEON_HAVE_GETOPT_H)
if (CHAMELEON_HAVE_GETOPT_H)
  check_function_exists(getopt_long CHAMELEON_HAVE_GETOPT_LONG)
endif()

# Add option to exploit cublas API v2
# -----------------------------------
cmake_dependent_option(CHAMELEON_USE_CUBLAS_V2
  "Enable cublas API v2" ON
  "CHAMELEON_USE_CUDA;CHAMELEON_SCHED_STARPU" OFF)

list(REMOVE_DUPLICATES CMAKE_C_FLAGS)
string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
list(REMOVE_DUPLICATES CMAKE_EXE_LINKER_FLAGS)
string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
# Fix a problem on Mac OS X when building shared libraries
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_SHARED_LINKER_FLAGS "-undefined dynamic_lookup")
endif()

# Check for the subdirectories
# ----------------------------

set(CHAMELEON_SOURCES_TARGETS "" CACHE INTERNAL "List of targets of sources")

# include headers
# ---------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)

# subdirectory for simulation mode with starpu+simgrid
if (CHAMELEON_SIMULATION)
    add_subdirectory(simucore)
else (CHAMELEON_SIMULATION)
  ###############################################################################
  # Coreblas library (kernels for CPUs) #
  #######################################
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/coreblas/include)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/coreblas/include)
  add_subdirectory(coreblas)

  ###############################################################################
  # Cudablas library (kernels for CUDAs) #
  ########################################
  if(CHAMELEON_USE_CUDA)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/cudablas/include)
    include_directories(${CMAKE_CURRENT_BINARY_DIR}/cudablas/include)
    add_subdirectory(cudablas)
  endif()
endif(CHAMELEON_SIMULATION)

###############################################################################
# Runtimes libraries #
######################
add_subdirectory(runtime)

###############################################################################
# Main library #
################
add_subdirectory(control)
add_subdirectory(compute)
add_subdirectory(include)

###############################################################################
# Executables and Tests #
#########################
enable_testing()
include(CTest)

# Necessary to compile C executables linked with ifort
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    message(STATUS "Add -nofor_main to the Fortran linker (Intel compiler)")
    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -nofor_main")
endif()

# Testing executables
if(CHAMELEON_ENABLE_EXAMPLE AND NOT CHAMELEON_SIMULATION)
  add_subdirectory(example)
endif(CHAMELEON_ENABLE_EXAMPLE AND NOT CHAMELEON_SIMULATION)

# Testing executables
if(CHAMELEON_ENABLE_TESTING AND NOT CHAMELEON_SIMULATION)
  add_subdirectory(testing)
endif(CHAMELEON_ENABLE_TESTING)
#------------------------------------------------------------------------------

# Define a target which gathers all targets of sources
# useful for generating the doc without compiling the binaries (objects, libraries)
add_custom_target(chameleon_all_sources ALL DEPENDS ${CHAMELEON_SOURCES_TARGETS})

###############################################################################
# Documentation #
#################

if(CHAMELEON_ENABLE_DOC)
    add_subdirectory(doc)
    if (DOXYGEN_EXECUTABLE AND EMACS_COMPILER AND TEX_COMPILER)
        add_custom_target(doc ALL DEPENDS doxygen-out doc-html-users_guide doc-pdf-users_guide)
    endif()
endif()
#------------------------------------------------------------------------------


###############################################################################
# Make cmake include starpu_simgrid_wrap.h (which #define main and thus makes #
# link tests fail) only after tests.                                          #
###############################################################################
if( CHAMELEON_SCHED_STARPU )
    if ( STARPU_FOUND )
        if (CHAMELEON_SIMULATION)
            list(APPEND CMAKE_REQUIRED_FLAGS "-include" "starpu_simgrid_wrap.h")
        endif()
    endif()
endif()

###############################################################################
# Config files (pkg_config) #
#############################

# Create .pc file
# ---------------
generate_pkgconfig_file()
#------------------------------------------------------------------------------

###############################################################################
# Release source #
##################
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_NAME "chameleon")
set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README.org)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Chameleon is a C library providing parallel algorithms to perform BLAS/LAPACK "
                                      " operations exploiting fully modern architectures.")
set(CHAMELEON_VERSION_PATCH ${CHAMELEON_VERSION_MICRO})
set(CPACK_PACKAGE_VERSION "${CHAMELEON_VERSION_MAJOR}.${CHAMELEON_VERSION_MINOR}.${CHAMELEON_VERSION_PATCH}")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")
set(CPACK_PACKAGE_VENDOR "Inria")
set(CPACK_PACKAGE_CONTACT "morse-devel@lists.gforge.inria.fr")
set(CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_SOURCE_DIR}/README.org)
set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/LICENCE.txt)
set(CPACK_SOURCE_IGNORE_FILES "/plasma-conversion/;/build*;.cproject;.settings;.dir-locals.el;.project;.pydevproject;/.git*")
include(CPack)

###############################################################################
# Print Options #
#################

# Save extra dependencies (all required links)
list(APPEND CHAMELEON_DEP ${EXTRA_LIBRARIES})
# We remove duplicate dependencies ensuring that redundant libraries are placed last
list(REVERSE CHAMELEON_DEP)
list(REMOVE_DUPLICATES CHAMELEON_DEP)
list(REVERSE CHAMELEON_DEP)
# we add chameleon's own library
set(CHAMELEON_LIBRARIES "chameleon")
if (NOT CHAMELEON_SIMULATION)
  list(APPEND CHAMELEON_LIBRARIES coreblas)
  if (CHAMELEON_USE_CUDA)
    list(APPEND CHAMELEON_LIBRARIES cudablas)
  endif()
endif()
if (NOT CHAMELEON_SCHED_DLOPEN)
  if (CHAMELEON_SCHED_OPENMP)
    list(APPEND CHAMELEON_LIBRARIES chameleon_openmp)
  endif()
  if (CHAMELEON_SCHED_PARSEC)
    list(APPEND CHAMELEON_LIBRARIES chameleon_parsec)
  endif()
  if (CHAMELEON_SCHED_QUARK)
    list(APPEND CHAMELEON_LIBRARIES chameleon_quark)
  endif()
  if (CHAMELEON_SCHED_STARPU)
    list(APPEND CHAMELEON_LIBRARIES chameleon_starpu)
  endif()
endif()
list(APPEND CHAMELEON_LIBRARIES hqr)
set(CHAMELEON_LIBRARIES_DEP ${CHAMELEON_LIBRARIES} ${CHAMELEON_DEP})

include(PrintOpts)

if (CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" OR CMAKE_INSTALL_PREFIX STREQUAL "C:/Program Files")
    message("-- ${Yellow}Your CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX} which is a default system path."
    " You may want to change it: set the CMAKE_INSTALL_PREFIX variable to do so${ColourReset}")
else()
    message("-- ${Yellow}Your CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX}${ColourReset}")
endif()

###
### END CMakeLists.txt
###
