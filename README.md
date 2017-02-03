[![Build Status](https://travis-ci.org/cbcoutinho/learn_dg.svg)](https://travis-ci.org/cbcoutinho/learn_dg)

[![forthebadge](http://forthebadge.com/images/badges/fuck-it-ship-it.svg)](http://forthebadge.com)

## Dependencies
This project uses lapack/blas as the numerical backend. Linking to these
libraries was no problem on Linux, but on Windows you need to build lapack/blas
from scratch.

First download lapack/blas from [Netlib](http://www.netlib.org/lapack/#_lapack_for_windows)

Install using `cmake`, and then move the `libblas.dll.a` and `liblapack.dll.a`
libraries into wherever `mingw` expects to find them:

On my machine, I use the following target directory:

    C:\Development\MinGW\lib

## Documentation
The documentation is built using [FORD](https://github.com/cmacmackin/ford), a python package that produces documentation from source code similar to Doxygen, but specifically designed for Fortran.
FORD takes inline documentation, prefixed by `!!`, and the markdown files located in `docs/user_guide`, and wraps them all up into an attractive website, available [here](https://cbcoutinho.github.io/learn_dg).
The documentation website is automatically updated whenever a commit is added to the master branch. This is done through Travis-CI.

To build the documentation locally run

    make docs

This will build the documentation website in the `./docs/html` directory. The organization of this project was heavily influenced by [FIDASIM](https://github.com/D3DEnergetic/FIDASIM)
