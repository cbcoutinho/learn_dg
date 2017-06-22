SET(CMAKE_CROSSCOMPILING TRUE)

# this one is important
SET(CMAKE_SYSTEM_NAME Windows)
#this one not so much
SET(CMAKE_SYSTEM_VERSION 1)

# specify the cross compiler
SET(CMAKE_C_COMPILER        x86_64-w64-mingw32-gcc)
SET(CMAKE_CXX_COMPILER      x86_64-w64-mingw32-g++)
SET(CMAKE_Fortran_COMPILER  x86_64-w64-mingw32-gfortran)

# set(BLAS_LIBRARIES /usr/x86_64-w64-mingw32/sys-root/mingw/lib/libblas.dll.a)
# set(LAPACK_LIBRARIES /usr/x86_64-w64-mingw32/sys-root/mingw/lib/liblapack.dll.a;/usr/x86_64-w64-mingw32/sys-root/mingw/lib/libblas.dll.a)

# add_library(${LAPACK_LIBRARIES})

# SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
# SET(BUILD_SHARED_LIBRARIES OFF)
# SET(CMAKE_EXE_LINKER_FLAGS "-static")

# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-libgcc")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++")
# SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -static-libgfortran")

# SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --sysroot=${CMAKE_FIND_ROOT_PATH}")
# SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --sysroot=${CMAKE_FIND_ROOT_PATH}")
# SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} --sysroot=${CMAKE_FIND_ROOT_PATH}")
# SET(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

# BLAS and LAPACK libraries


# where is the target environment
SET(CMAKE_FIND_ROOT_PATH /usr/x86_64-w64-mingw32 /usr/x86_64-w64-mingw32/sys-root/mingw )

# search for programs in the build host directories
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
# for libraries and headers in the target directories
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
