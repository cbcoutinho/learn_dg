[![Build Status](https://travis-ci.org/cbcoutinho/learn_dg.svg)](https://travis-ci.org/cbcoutinho/learn_dg)

[![forthebadge](http://forthebadge.com/images/badges/fuck-it-ship-it.svg)](http://forthebadge.com)

## Documentation
Documentation about the project can be found here:

[Learn_dg Documentation](https://cbcoutinho.github.io/learn_dg/)

## Dependencies
This project has been refactored to be built on Linux using CMake. It's still possible to create a windows executable, but that is all done by cross-compilers on the Linux side. Further, the project also requires BLAS/Lapack.

To build the project, execute `make cmake` at the project root directory (cmake is a make recipe right now). If that is strange for you, make a build subdirectory in this directory and execute the following:

    cmake ..

To create a windows executable/library of the project, either execute `make cmake_win` at the project root, or create a build subdicretory and execute the following:

    cmake -DCMAKE_TOOLCHAIN_FILE:STRING=../cmake/Toolchain-x64-mingw32.cmake ..

Cross compiling the sources for Windows requires a few exotic packages. On openSUSE Leap 42.2, I had to first add the 'windows_mingw_win64' repo and install the following packages. These should cover the current state of the project.

```
mingw64-cross-gcc
mingw64-cross-g++
mingw64-cross-gfortran
mingw64-lapack-devel
mingw64-blas-devel
```
