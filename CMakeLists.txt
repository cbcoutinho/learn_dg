project(learn_dg Fortran)
cmake_minimum_required(VERSION 3.7)

#-----------------------------------------------------------------------------------
#                      CMake macros for adding source files
#-----------------------------------------------------------------------------------

# Includes the cmake_dependent_option macro used in profiling setting below
include( CMakeDependentOption )

include(cmake/macros/add_exe_sources.txt)
include(cmake/macros/add_lib_sources.txt)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
# set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/lib)


#-----------------------------------------------------------------------------------
#                               Set default options
#-----------------------------------------------------------------------------------

# Set default options

# Default release build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build. Options are: Debug Release" FORCE)
endif()

if(NOT PROFILE)
    set(PROFILE OFF CACHE BOOL "Choose to enable profiling. Options are: ON OFF" FORCE)
endif()

# If NOT CMAKE_CROSSCOMPILING, create a variable called PROFILE with default off
cmake_dependent_option(PROFILE
    "Use Profiling" OFF
    "NOT CMAKE_CROSSCOMPILING" OFF
)

message(STATUS " ")
message(STATUS "    Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "    Profiling: ${PROFILE}")

#-----------------------------------------------------------------------------------
#                               Set compiler flags
#-----------------------------------------------------------------------------------

SET(FFLAGS "-std=f2008 -fPIC -fmax-errors=1")

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set (USER_COMPILE_FLAGS_DEBUG
        "-O0 -fpp -assume realloc_lhs -debug -g \
        -traceback -check all -mcmodel=medium \
        -fp-stack-check -heap-arrays 1600"
    )
    set (USER_COMPILE_FLAGS_RELEASE "-O3 -fpp -assume realloc_lhs -mcmodel=medium \
        -align array32byte -traceback"
    )
    set (CMAKE_Fortran_LINK_FLAGS_DEBUG "")
    set (CMAKE_Fortran_LINK_FLAGS_RELEASE "")

    set (PROFILE_FLAGS "-p")

    #set (CMAKE_AR xiar)        # for -ipo
    #set (CMAKE_LINKER xild)    # for -ipo

elseif (CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 6.3)
        message (FATAL_ERROR "Minimum required gfortran version is 6.3")
    endif()

    set (USER_COMPILE_FLAGS_DEBUG   "-O0 -cpp --coverage -ffree-form"
        "-ffree-line-length-none -fbacktrace"
        "-fstack-protector -fcheck=bounds"
        "-g"
    )
    string(REPLACE ";" " " USER_COMPILE_FLAGS_DEBUG "${USER_COMPILE_FLAGS_DEBUG}")

    set (USER_COMPILE_FLAGS_RELEASE "-O3 -cpp -ffree-form "
        "-ffree-line-length-none -fexternal-blas "
        "-march=native -funroll-loops "
    )
    string(REPLACE ";" " " USER_COMPILE_FLAGS_RELEASE "${USER_COMPILE_FLAGS_RELEASE}")

    set (CMAKE_Fortran_LINK_FLAGS_DEBUG "-O0 --coverage ")
    set (CMAKE_Fortran_LINK_FLAGS_RELEASE "-O3 ")

    set (PROFILE_FLAGS "-p")

endif()

# Append user flags to debug and release variables
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${FFLAGS} ${USER_COMPILE_FLAGS_DEBUG}")
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${FFLAGS} ${USER_COMPILE_FLAGS_RELEASE}")

# Append profile flags

if (PROFILE)
    set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${PROFILE_FLAGS}")
    set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${PROFILE_FLAGS}")
endif()



#----------------------------------------------------------------------------------------------
#                                         Find Lapack
#----------------------------------------------------------------------------------------------

message(STATUS " ")
message(STATUS "BLAS/LAPACK Summary:")
find_package(LAPACK REQUIRED)
message(STATUS "    LAPACK Libraries: ")
foreach (lib ${LAPACK_LIBRARIES})
    message(STATUS "        " ${lib})
endforeach()

#----------------------------------------------------------------------------------------------
#                                         Find OpenMP
#----------------------------------------------------------------------------------------------

message(STATUS " ")
message(STATUS "Looking for OpenMP Libraries")
find_package(OpenMP)
if (OPENMP_FOUND)
    message(STATUS "    OpenMP Libraries Found")
    SET (USE_OPENMP 1)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()



#----------------------------------------------------------------------------------------------
#                                   Add source tree and includes
#----------------------------------------------------------------------------------------------

add_subdirectory(src)

#----------------------------------------------------------------------------------------------
#                                         Build libraries
#----------------------------------------------------------------------------------------------

set(LIB_NAME "core")
message(STATUS " ")

add_library(core-static STATIC ${LIB_SRCS}) # Static library
add_library(core-shared SHARED ${LIB_SRCS}) # Shared library

foreach( mylib core-shared core-static )
    target_link_libraries( ${mylib}
        ${LAPACK_LIBRARIES}
        ${BLAS_LIBRARIES}
    )

    set_target_properties( ${mylib}
        PROPERTIES
        OUTPUT_NAME     ${LIB_NAME}
        COMPILE_FLAGS   "${CMAKE_Fortran_FLAGS}"
        LINK_FLAGS      "${CMAKE_Fortran_LINK_FLAGS}"
        Fortran_MODULE_DIRECTORY ${mylib}
    )
endforeach()

#----------------------------------------------------------------------------------------------
#                                       Build Interfaces
#----------------------------------------------------------------------------------------------
# Add driver executable linked to core library
# add_executable(main ${EXE_SRCS})

# set (lib_name core-static)
set (lib_name core-shared)
foreach( SOURCE_FILE ${EXE_SRCS} )
    string( REPLACE ".f90" "" SOURCE_NAME ${SOURCE_FILE} )
    string( REPLACE "src/" "" SOURCE_NAME ${SOURCE_NAME} )
    add_executable( ${SOURCE_NAME} ${SOURCE_FILE} )
    target_link_libraries( ${SOURCE_NAME} ${lib_name} )
    target_include_directories(${SOURCE_NAME} PUBLIC ${CMAKE_BINARY_DIR}/${lib_name})
    # include_directories(${SOURCE_NAME} ${CMAKE_BINARY_DIR}/${lib_name})
    set_target_properties( ${SOURCE_NAME}
        PROPERTIES
        COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
        LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS}"
        # Fortran_MODULE_DIRECTORY ${lib_name}
    )
endforeach( SOURCE_FILE ${EXE_SRCS} )
