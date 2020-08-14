###
#
# @file starpu/CMakeLists.txt
#
# @copyright 2009-2015 The University of Tennessee and The University of
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
cmake_minimum_required(VERSION 2.8)

set(CHAMELEON_STARPU_VERSION "1.3" CACHE STRING "oldest STARPU version desired")

# create list of components in order to make a single call to find_package(starpu...)
if(NOT CHAMELEON_SIMULATION)
  set(STARPU_COMPONENT_LIST "HWLOC")
  if(CHAMELEON_USE_CUDA)
    list(APPEND STARPU_COMPONENT_LIST "CUDA")
  endif()
else()
  set(STARPU_COMPONENT_LIST "SIMGRID")
endif()
if(CHAMELEON_USE_MPI)
  list(APPEND STARPU_COMPONENT_LIST "MPI")
endif()

find_package(STARPU ${CHAMELEON_STARPU_VERSION} REQUIRED
  COMPONENTS ${STARPU_COMPONENT_LIST})

# Add definition and include_dir if found
if ( STARPU_FOUND )
  message("-- ${Blue}Add definition CHAMELEON_SCHED_STARPU"
    " - Activate StarPU in Chameleon${ColourReset}")
  set(CHAMELEON_SCHED "${CHAMELEON_SCHED};STARPU" CACHE STRING "Runtime systems available" FORCE)
  if (STARPU_INCLUDE_DIRS_DEP)
    include_directories(${STARPU_INCLUDE_DIRS_DEP})
    set(CMAKE_REQUIRED_INCLUDES "${STARPU_INCLUDE_DIRS_DEP}")
  endif()
  if(STARPU_CFLAGS_OTHER_DEP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${STARPU_CFLAGS_OTHER_DEP}")
    string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${STARPU_CFLAGS_OTHER_DEP}")
  endif()
  if(STARPU_LDFLAGS_OTHER_DEP)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${STARPU_LDFLAGS_OTHER_DEP}")
    string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_REQUIRED_LDFLAGS "${STARPU_LDFLAGS_OTHER_DEP}")
  endif()
  if(STARPU_LIBRARY_DIRS_DEP)
    list(APPEND CMAKE_INSTALL_RPATH "${STARPU_LIBRARY_DIRS_DEP}")
  endif()
  if (STARPU_LIBRARIES_DEP)
    set(CHAMELEON_DEP  "${STARPU_LIBRARIES_DEP};${CHAMELEON_DEP}" PARENT_SCOPE)
    set(CMAKE_REQUIRED_LIBRARIES "${STARPU_LIBRARIES_DEP}")
  endif()
  string(REPLACE ";" " " CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
  check_function_exists(starpu_data_idle_prefetch_on_node HAVE_STARPU_IDLE_PREFETCH)
  if ( HAVE_STARPU_IDLE_PREFETCH )
    message("-- ${Blue}Add definition HAVE_STARPU_IDLE_PREFETCH${ColourReset}")
  endif()
  check_function_exists(starpu_iteration_push HAVE_STARPU_ITERATION_PUSH)
  if ( HAVE_STARPU_ITERATION_PUSH )
    message("-- ${Blue}Add definition HAVE_STARPU_ITERATION_PUSH${ColourReset}")
  endif()
  check_function_exists(starpu_data_wont_use HAVE_STARPU_DATA_WONT_USE)
  if ( HAVE_STARPU_DATA_WONT_USE )
    message("-- ${Blue}Add definition HAVE_STARPU_DATA_WONT_USE${ColourReset}")
  endif()
  check_function_exists(starpu_data_set_coordinates HAVE_STARPU_DATA_SET_COORDINATES)
  if ( HAVE_STARPU_DATA_SET_COORDINATES )
    message("-- ${Blue}Add definition HAVE_STARPU_DATA_SET_COORDINATES${ColourReset}")
  endif()
  check_function_exists(starpu_malloc_on_node_set_default_flags HAVE_STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS)
  if ( HAVE_STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS )
    message("-- ${Blue}Add definition HAVE_STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS${ColourReset}")
  endif()
  check_function_exists(starpu_data_set_ooc_flag HAVE_STARPU_DATA_SET_OOC_FLAG)
  if ( HAVE_STARPU_DATA_SET_OOC_FLAG )
    message("-- ${Blue}Add definition HAVE_STARPU_DATA_SET_OOC_FLAG${ColourReset}")
  endif()
  check_function_exists(starpu_mpi_init_conf HAVE_STARPU_MPI_INIT_CONF)
  if ( HAVE_STARPU_MPI_INIT_CONF )
    message("-- ${Blue}Add definition HAVE_STARPU_MPI_INIT_CONF${ColourReset}")
  endif()
  check_function_exists(starpu_mpi_wait_for_all HAVE_STARPU_MPI_WAIT_FOR_ALL)
  if ( HAVE_STARPU_MPI_WAIT_FOR_ALL )
    message("-- ${Blue}Add definition HAVE_STARPU_MPI_WAIT_FOR_ALL${ColourReset}")
  endif()
  check_function_exists(starpu_mpi_data_migrate HAVE_STARPU_MPI_DATA_MIGRATE)
  if ( HAVE_STARPU_MPI_DATA_MIGRATE )
    message("-- ${Blue}Add definition HAVE_STARPU_MPI_DATA_MIGRATE${ColourReset}")
  elseif(CHAMELEON_USE_MIGRATE)
    set(CHAMELEON_USE_MIGRATE "OFF")
    message("-- ${Blue}CHAMELEON_USE_MIGRATE is turned OFF because starpu_mpi_data_migrate not found${ColourReset}")
  endif()
  check_function_exists(starpu_mpi_interface_datatype_register HAVE_STARPU_MPI_INTERFACE_DATATYPE_REGISTER)
  if ( HAVE_STARPU_MPI_INTERFACE_DATATYPE_REGISTER )
    message("-- ${Blue}Add definition HAVE_STARPU_MPI_INTERFACE_DATATYPE_REGISTER${ColourReset}")
  else()
    if( CHAMELEON_USE_MPI_DATATYPES )
      set( CHAMELEON_USE_MPI_DATATYPES OFF CACHE BOOL "Enable MPI datatypes whenever supported by the runtime" FORCE)
    endif()
  endif()
  check_function_exists(starpu_interface_copy2d HAVE_STARPU_INTERFACE_COPY2D)
  if ( HAVE_STARPU_INTERFACE_COPY2D )
    message("-- ${Blue}Add definition HAVE_STARPU_INTERFACE_COPY2D${ColourReset}")
  endif()
  if (CHAMELEON_USE_MPI)
    # Check if a specific function exist
    check_function_exists(starpu_mpi_data_register_comm HAVE_STARPU_MPI_DATA_REGISTER)
    if ( HAVE_STARPU_MPI_DATA_REGISTER )
      message("-- ${Blue}Add definition HAVE_STARPU_MPI_DATA_REGISTER - Activate"
        " use of starpu_mpi_data_register() in Chameleon with StarPU${ColourReset}")
    endif()
    check_function_exists(starpu_mpi_comm_rank HAVE_STARPU_MPI_COMM_RANK)
    if ( HAVE_STARPU_MPI_COMM_RANK )
      message("-- ${Blue}Add definition HAVE_STARPU_MPI_COMM_RANK - Activate"
        " use of starpu_mpi_comm_rank() in Chameleon with StarPU${ColourReset}")
    endif()
    check_function_exists(starpu_mpi_cached_receive HAVE_STARPU_MPI_CACHED_RECEIVE)
    if ( HAVE_STARPU_MPI_CACHED_RECEIVE )
      message("-- ${Blue}Add definition HAVE_STARPU_MPI_CACHED_RECEIVE${ColourReset}")
    endif()
    check_function_exists(starpu_mpi_comm_get_attr HAVE_STARPU_MPI_COMM_GET_ATTR)
    if ( HAVE_STARPU_MPI_COMM_GET_ATTR )
      message("-- ${Blue}Add definition HAVE_STARPU_MPI_COMM_GET_ATTR${ColourReset}")
    endif()
  endif()
else ( STARPU_FOUND )
  if(CHAMELEON_VERBOSE_FIND_PACKAGE)
    if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
      if (NOT HWLOC_hwloc.h_DIRS)
        Print_Find_Header_Status(hwloc hwloc.h)
      endif ()
      if (NOT HWLOC_hwloc_LIBRARY)
        Print_Find_Library_Status(hwloc libhwloc)
      endif ()
    endif()
    if(CHAMELEON_SIMULATION AND (NOT SIMGRID_FOUND OR NOT SIMGRID_LIBRARIES))
      if (NOT SIMGRID_simgrid.h_DIRS)
        Print_Find_Header_Status(simgrid simgrid.h)
      endif ()
      if (NOT SIMGRID_simgrid_LIBRARY)
        Print_Find_Library_Status(simgrid libsimgrid)
      endif ()
    endif()
    if( (NOT STARPU_SHM_FOUND) OR (NOT STARPU_SHM_LIBRARIES) OR
        ( STARPU_LOOK_FOR_MPI AND (NOT STARPU_MPI_FOUND OR NOT STARPU_MPI_LIBRARIES) )
        )
      foreach(starpu_hdr ${STARPU_hdrs_to_find})
        if (NOT STARPU_${starpu_hdr}_INCLUDE_DIRS)
          Print_Find_Header_Status(starpu ${starpu_hdr})
        endif ()
      endforeach()
      if(STARPU_VERSION_STRING)
        foreach(starpu_lib ${STARPU_libs_to_find})
          if (NOT STARPU_${starpu_lib}_LIBRARY)
            Print_Find_Library_Status(starpu ${starpu_lib})
          endif ()
        endforeach()
      endif ()
    endif ()
  else(CHAMELEON_VERBOSE_FIND_PACKAGE)
    message(WARNING "StarPU library has not been found and CHAMELEON_VERBOSE_FIND_PACKAGE is set to OFF."
      " Try to activate CHAMELEON_VERBOSE_FIND_PACKAGE option (-DCHAMELEON_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
  endif(CHAMELEON_VERBOSE_FIND_PACKAGE)
  if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
    message(FATAL_ERROR "hwloc library is required but has not been found")
  endif()
  if(CHAMELEON_SIMULATION AND (NOT SIMGRID_FOUND OR NOT SIMGRID_LIBRARIES))
    message(FATAL_ERROR "SimGrid library is required but has not been found")
  endif()
  if( (NOT STARPU_SHM_FOUND) OR (NOT STARPU_SHM_LIBRARIES) OR
      ( STARPU_LOOK_FOR_MPI AND (NOT STARPU_MPI_FOUND OR NOT STARPU_MPI_LIBRARIES) )
      )
    message(FATAL_ERROR "StarPU library is required but has not been found")
  endif()
  return()
endif ( STARPU_FOUND )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )
include_directories( ${CMAKE_CURRENT_BINARY_DIR}/include )

configure_file("include/chameleon_starpu.h.in"
               "include/chameleon_starpu.h"
               @ONLY)

# Generate headers for all possible precisions
# --------------------------------------------
set(RUNTIME_HDRS_GENERATED "")
set(ZHDR
  include/runtime_codelet_z.h
  )

precisions_rules_py(RUNTIME_HDRS_GENERATED "${ZHDR}"
  PRECISIONS "s;d;c;z;ds;zc"
  TARGETDIR  "include")

# Define the list of headers
# --------------------------
set(RUNTIME_HDRS
  ${CMAKE_CURRENT_BINARY_DIR}/include/chameleon_starpu.h
  include/runtime_codelet_profile.h
  include/runtime_codelets.h
  include/runtime_profiling.h
  include/runtime_workspace.h
  )

# Add generated headers
# ---------------------
foreach( hdr_file ${RUNTIME_HDRS_GENERATED} )
  list(APPEND RUNTIME_HDRS ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
endforeach()

# Force generation of headers
# ---------------------------
add_custom_target(
  runtime_starpu_include
  ALL SOURCES ${RUNTIME_HDRS})

# Installation
# ------------
install(
  FILES ${RUNTIME_HDRS}
  DESTINATION include/runtime/starpu )

# Generate the Chameleon common for all possible precisions
# ---------------------------------------------------------
set(RUNTIME_COMMON_GENERATED "")
set(ZSRC
  control/runtime_zprofiling.c
  control/runtime_zlocality.c
  )

precisions_rules_py(
  RUNTIME_COMMON_GENERATED "${ZSRC}"
  PRECISIONS "${CHAMELEON_PRECISION}"
  TARGETDIR "control")

set(RUNTIME_COMMON
  control/runtime_async.c
  control/runtime_context.c
  control/runtime_control.c
  control/runtime_descriptor.c
  control/runtime_options.c
  control/runtime_profiling.c
  control/runtime_workspace.c
  interface/cham_tile_interface.c
  ${RUNTIME_COMMON_GENERATED}
  )

set(flags_to_add "")
foreach(_prec ${CHAMELEON_PRECISION})
  set(flags_to_add "${flags_to_add} -DPRECISION_${_prec}")
endforeach()
set_source_files_properties(control/runtime_profiling.c PROPERTIES COMPILE_FLAGS "${flags_to_add}")

# Generate the Chameleon sources for all possible precisions
# ----------------------------------------------------------
set(RUNTIME_SRCS_GENERATED "")
set(ZSRC
  codelets/codelet_zcallback.c
  ${CODELETS_ZSRC}
  )

precisions_rules_py(RUNTIME_SRCS_GENERATED "${ZSRC}"
  PRECISIONS "${CHAMELEON_PRECISION}"
  TARGETDIR "codelets")

set(RUNTIME_SRCS
  ${RUNTIME_COMMON}
  ${RUNTIME_SRCS_GENERATED}
  ${CODELETS_SRC}
  )

# Force generation of sources
# ---------------------------
add_custom_target(starpu_sources ALL SOURCES ${RUNTIME_SRCS})
set(CHAMELEON_SOURCES_TARGETS "${CHAMELEON_SOURCES_TARGETS};runtime_starpu_include;starpu_sources" CACHE INTERNAL "List of targets of sources")

# Add library
# -----------
add_library(chameleon_starpu ${RUNTIME_SRCS})
set_property(TARGET chameleon_starpu PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET chameleon_starpu PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

target_link_libraries(chameleon_starpu
  ${STARPU_LIBRARIES_DEP})
if (NOT CHAMELEON_SIMULATION)
  target_link_libraries(chameleon_starpu
    coreblas)
  if(CHAMELEON_USE_CUDA)
    target_link_libraries(chameleon_starpu
      cudablas)
  endif(CHAMELEON_USE_CUDA)
endif(NOT CHAMELEON_SIMULATION)

add_dependencies(chameleon_starpu
  chameleon_include
  control_include
  runtime_starpu_include
  starpu_sources
  )
if (NOT CHAMELEON_SIMULATION)
  add_dependencies(chameleon_starpu coreblas_include)
  if (CHAMELEON_USE_CUDA)
    add_dependencies(chameleon_starpu cudablas_include)
  endif()
endif()

# installation
# ------------
install(TARGETS chameleon_starpu
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

###
### END CMakeLists.txt
###
