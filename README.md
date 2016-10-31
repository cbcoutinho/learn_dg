This project uses lapack/blas as the numerical backend. Linking to these
libraries was no problem on Linux, but on Windows you need to build lapack/blas
from scratch.

First download lapack/blas from here:
'http://www.netlib.org/lapack/#_lapack_for_windows'



Install using CMAKE, and then move the libblas.dll.a and liblapack.dll.a
libraries into wherever mingw expects to find libraries:

Example:
'C:\Development\MinGW\lib'
