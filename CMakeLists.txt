cmake_minimum_required(VERSION 3.7)

project(main)
enable_language(Fortran)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# FFLAGS
set (CMAKE_Fortran_FLAGS "-std=f2008 -fPIC -fmax-errors=1 -Wimplicit-interface -Wall -Wextra")
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_DEBUG   "-Og -g -fcheck=all -fbacktrace")
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


### Output Directories

# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/ar_cmake)
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib_cmake)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin_cmake)




#
# Compile.
#
file(GLOB sources  "src/*.f90" "src/fortranlib/src/lib_array.f90")
list(REMOVE_ITEM sources ${CMAKE_CURRENT_SOURCE_DIR}/src/doubleint.f90)
add_executable("main" ${sources})
TARGET_LINK_LIBRARIES("main" blas lapack)
