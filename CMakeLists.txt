cmake_minimum_required(VERSION 3.7)

project(main Fortran)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# FFLAGS
set (CMAKE_Fortran_FLAGS "-std=f2008 -fPIC -fmax-errors=1")
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_DEBUG   "-Og -g -fcheck=all -fbacktrace -Wimplicit-interface -Wall -Wextra")
  set (CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native -ffast-math -funroll-loops")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

# Default to Debug mode if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Debug' as none was specified.")
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()
# Put message
message( STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} )

find_library(CMAKE_BLAS_LIBRARY blas DOC "Path to the BLAS implementation")
find_library(CMAKE_LAPACK_LIBRARY lapack DOC "Path to the LAPACK implementation")

# find_package(CMAKE_BLAS_LIBRARY blas DOC "Path to the BLAS implementation")
# find_package(CMAKE_LAPACK_LIBRARY lapack DOC "Path to the LAPACK implementation")


message ( STATUS "Blas: ${CMAKE_BLAS_LIBRARY}")
message ( STATUS "Lapack: ${CMAKE_LAPACK_LIBRARY}")

# find_package( LAPACK REQUIRED ${CMAKE_LAPACK_LIBRARY} )
# find_package( BLAS REQUIRED ${CMAKE_BLAS_LIBRARY})

#
# Compile.
#
file(GLOB sources  "src/*.f90" "src/fortranlib/src/lib_array.f90")
list(REMOVE_ITEM sources ${CMAKE_CURRENT_SOURCE_DIR}/src/doubleint.f90)
add_executable("main" ${sources})
target_link_libraries("main" lapack blas)

# LINK_DIRECTORIES(${CMAKE_LAPACK_LIBRARY} ${CMAKE_BLAS_LIBRARY})
