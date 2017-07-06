project(learn_dg Fortran)
cmake_minimum_required(VERSION 3.7)


#-----------------------------------------------------------------------------------
#                      CMake macros for adding source files
#-----------------------------------------------------------------------------------
include(cmake/macros/add_exe_sources.txt)
include(cmake/macros/add_lib_sources.txt)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/lib)


#-----------------------------------------------------------------------------------
#                               Set default options
#-----------------------------------------------------------------------------------

# Set default options

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build. Options are: Debug Release" FORCE)
endif()

message(STATUS "        Build Type: ${CMAKE_BUILD_TYPE}")

# if(NOT BUILD_SHARED_LIBRARIES)
#     set(BUILD_SHARED_LIBRARIES OFF CACHE BOOL "Build shared libraries - Default OFF. Options are: ON OFF" FORCE)
# endif()




#-----------------------------------------------------------------------------------
#                               Set compiler flags
#-----------------------------------------------------------------------------------
SET(FFLAGS "-std=f2008 -fPIC -fmax-errors=1")

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set (USER_COMPILE_FLAGS_DEBUG   "-O0 -fpp -assume realloc_lhs -debug -g \
                                     -traceback -check all -mcmodel=medium \
                                     -fp-stack-check -heap-arrays 1600")
    set (USER_COMPILE_FLAGS_RELEASE "-O3 -fpp -assume realloc_lhs -mcmodel=medium \
                                     -align array32byte -traceback")
    set (CMAKE_Fortran_LINK_FLAGS_DEBUG "")
    set (CMAKE_Fortran_LINK_FLAGS_RELEASE "")

    set (PROFILE_FLAGS "-p")

    #set (CMAKE_AR xiar)        # for -ipo
    #set (CMAKE_LINKER xild)    # for -ipo

elseif (CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 6.3)
        message (FATAL_ERROR "Minimum required gfortran version is 6.3")
    endif()

    set (USER_COMPILE_FLAGS_DEBUG   "${FFLAGS} -O0 -cpp --coverage -ffree-form \
                                     -ffree-line-length-none -fbacktrace \
                                     -fstack-protector -fcheck=bounds \
                                     -fmax-stack-var-size=1600")
                                    # -Wall -Werror")
    set (USER_COMPILE_FLAGS_RELEASE "${FFLAGS} -O3 -cpp -ffree-form \
                                     -ffree-line-length-none -fexternal-blas \
                                     -march=native -funroll-loops")

    set (CMAKE_Fortran_LINK_FLAGS_DEBUG "-O0 --coverage")
    set (CMAKE_Fortran_LINK_FLAGS_RELEASE "-O3")

    set (PROFILE_FLAGS "-p -g")

endif()

# Append user flags to debug and release variables
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${USER_COMPILE_FLAGS_DEBUG}")
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${USER_COMPILE_FLAGS_RELEASE}")

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
message(STATUS "        LAPACK Libraries: ${LAPACK_LIBRARIES}")





#----------------------------------------------------------------------------------------------
#                                   Add source tree and includes
#----------------------------------------------------------------------------------------------
add_subdirectory(src)


#----------------------------------------------------------------------------------------------
#                                         Build library
#----------------------------------------------------------------------------------------------

set(LIB_NAME "core")

message(STATUS " ")

# Static library
add_library(core-static STATIC ${LIB_SRCS})
target_link_libraries(core-static
                    ${LAPACK_LIBRARIES}
                    ${BLAS_LIBRARIES}
                    )
set_target_properties(core-static PROPERTIES
                    COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                    LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS}"
                    OUTPUT_NAME ${LIB_NAME}
                    )

# Shared library
add_library(core-shared SHARED ${LIB_SRCS})
target_link_libraries(core-shared
                    ${LAPACK_LIBRARIES}
                    ${BLAS_LIBRARIES}
                    )
set_target_properties(core-shared PROPERTIES
                    COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                    LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS}"
                    OUTPUT_NAME ${LIB_NAME}
                    )



#----------------------------------------------------------------------------------------------
#                                       Build Interfaces
#----------------------------------------------------------------------------------------------
# Add driver executable linked to core library
# add_executable(main ${EXE_SRCS})

foreach( SOURCE_FILE ${EXE_SRCS})
    string( REPLACE ".f90" "" SOURCE_NAME ${SOURCE_FILE})
    string( REPLACE "src/" "" SOURCE_NAME ${SOURCE_NAME})
    add_executable(${SOURCE_NAME} ${SOURCE_FILE})
    target_link_libraries(${SOURCE_NAME}  core-static)
    set_target_properties(${SOURCE_NAME}  PROPERTIES
                        COMPILE_FLAGS "${CMAKE_Fortran_FLAGS}"
                        LINK_FLAGS    "${CMAKE_Fortran_LINK_FLAGS}")
endforeach( SOURCE_FILE ${EXE_SRCS} )
