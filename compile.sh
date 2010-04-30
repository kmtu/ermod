#!/bin/sh

F90FLAGS_noOPT="-g -traceback -xP -mcmodel=large -autodouble -extend-source -shared-intel -fnoomit-frame-pointer"
F90FLAGS_noOPT="-g -march=native -fdefault-real-8 -ffixed-line-length-132 -fno-omit-frame-pointer -fno-math-errno"
F90FLAGS="$F90FLAGS_noOPT -O3"
LAPACKLIBS="-lmkl_em64t -lguide"
LAPACKLIBS="-lblas -llapack"
#FFTFLAGS="-DMKL"
#FFTLIBS="/opt/intel/mkl/10.0.011/include/mkl_dfti.f90"
FFTFLAGS="-DFFTW -I$HOME/include"
FFTLIBS="-L$HOME/lib -lfftw3"
LIBS="$FFTLIBS $LAPACKLIBS"
#mpif90 -Wall -o energromacs engmain.f90 realcal_blk.f90 spline.f90 fft_iface.f90 enganal.f $F90FLAGS -cpp -Dtrjctry -DGROMACS ${FFTFLAGS} ${LIBS} || exit 1
mpif90 -Wall -o energromacs_p engmain.f90 realcal_blk.f90 spline.f90 fft_iface.f90 enganal.f $F90FLAGS -cpp -Dtrjctry -DGROMACS ${FFTFLAGS} ${LIBS} -L$HOME/lib -lprofiler || exit 1
#mpif90 -Wall -o energromacs_noMPI engmain.f90 realcal_blk.f90 spline.f90 fft_iface.f90 enganal.f $F90FLAGS_noOPT -cpp -Dtrjctry -DGROMACS -DnoMPI ${FFTFLAGS} ${LIBS} -L$HOME/lib -lprofiler || exit 1
#mpif90 -Wall -o enernamd engmain.f90 realcal_blk.f90 spline.f90 fft_iface.f90 enganal.f $F90FLAGS -cpp -Dtrjctry -DNAMD ${FFTFLAGS} ${LIBS} || exit 1
#mpiifort -g -o enganal.o -c enganal.f $F90FLAGS -cpp -Dtrjctry -DMARBLE -DMKL 
#mpiifort -o insertion.o -c insertion.f $F90FLAGS -cpp -Dtrjctry -DMARBLE -DMKL 
#mpiifort -o setconf.o -c setconf.f $F90FLAGS -cpp -Dtrjctry -DMARBLE -DMKL 
#mpiifort -g -o enermarble enganal.o /opt/intel/mkl/10.0.011/include/mkl_dfti.f90 -lmkl_em64t -lguide -L$HOME/lib -lprofiler
#mpiifort -g -o enganal.o -c enganal.f $F90FLAGS -cpp -Dtrjctry -DGROMACS -DMKL 
#mpiifort -g -o energromacs enganal.o /opt/intel/mkl/10.0.011/include/mkl_dfti.f90 -lmkl_em64t -lguide -L$HOME/lib -lprofiler
#mpiifort -g -o enganal.o -c enganal.f $F90FLAGS -cpp -Dtrjctry -DNAMD -DMKL 
#mpiifort -g -o enernamd enganal.o /opt/intel/mkl/10.0.011/include/mkl_dfti.f90 -lmkl_em64t -lguide -L$HOME/lib -lprofiler
#mpiifort -o energromacs -O3 -xP -fpic -shared-intel -mcmodel=large –autodouble -extend-source -cpp -Dtrjctry -DGROMACS -DMKL /opt/intel/mkl/10.0.011/include/mkl_dfti.f90 enganal.f insertion.f setconf.f -lmkl_em64t -lguide -L/usr/lib –lpthread
#mpiifort -o enernamd -O3 -xP -fpic -shared-intel -mcmodel=large -autodouble -extend-source -cpp -Dtrjctry -DNAMD -DMKL /opt/intel/mkl/10.0.011/include/mkl_dfti.f90 enganal.f insertion.f setconf.f -lmkl_em64t -lguide -L/usr/lib –lpthread
