SET(CMAKE_CROSSCOMPILING TRUE)

# this one is important
SET(CMAKE_SYSTEM_NAME Windows)
#this one not so much
SET(CMAKE_SYSTEM_VERSION 1)

# specify the cross compiler
SET(CMAKE_C_COMPILER        i686-w64-mingw32-gcc)
SET(CMAKE_CXX_COMPILER      i686-w64-mingw32-g++)
SET(CMAKE_Fortran_COMPILER  i686-w64-mingw32-gfortran)
SET(OpenMP_gomp_LIBRARY     libgomp-1.dll)


# where is the target environment
SET(CMAKE_FIND_ROOT_PATH /usr/i686-w64-mingw32 /usr/i686-w64-mingw32/sys-root/mingw )

# search for programs in the build host directories
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
# for libraries and headers in the target directories
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
