ERmod - a tool for obtaining histograms of solute-solvent interaction energy
=========================

## Typical Installation
This package is built with autotools. To compile the program,

1. configure the package with "configure",
2. then compile the package with "make".

If you have Intel MKL, configure program with:

    $ ./configure --with-mkl

If your computer have BLAS & LAPACK library, and [FFTW](http://www.fftw.org/) version 3, try:

    $ ./configure --with-fft=fftw

If configuration finishes successfully, type

    $ make

to start compilation.

## Non-typical Installation
If you want to use the specific version of Intel MKL, try

    ./configure --with-mkl=(version number) --with-fft=mkl

If you want to use the specific BLAS library, try

    ./configure --with-blas=(BLAS library specification)

Current configuration script only supports gfortran and ifort as a compiler.
For other compilers, please specify F77 / FC / FFLAGS / FCFLAGS environment variables to appropriate ones.
