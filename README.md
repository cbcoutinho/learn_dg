[![Build Status](https://travis-ci.org/cbcoutinho/learn_dg.svg?branch=master)](https://travis-ci.org/cbcoutinho/learn_dg)


## Dependencies
This project uses lapack/blas as the numerical backend. Linking to these
libraries was no problem on Linux, but on Windows you need to build lapack/blas
from scratch.

First download lapack/blas from here:
'http://www.netlib.org/lapack/#_lapack_for_windows'



Install using CMAKE, and then move the libblas.dll.a and liblapack.dll.a
libraries into wherever mingw expects to find libraries:

Example:
'C:\Development\MinGW\lib'

## Documentation Website
The documentation is built using [FORD](https://github.com/cmacmackin/ford).
FORD takes inline documentation, prefixed by `!+`, and the markdown files located in `docs/user_guide` and wraps them all up in a pretty website.
The documentation website is automatically updated whenever a commit is added to the master branch. This is done through Travis-CI.

To build the documentation locally run
```
make docs
```
This will build the documentation website in the `docs/html` directory.

This project organization was heavily influenced by [FIDASIM](https://github.com/D3DEnergetic/FIDASIM)
