#!/bin/sh

F90FLAGS_noOPT="-g -traceback -xP -mcmodel=large -autodouble -extend-source -shared-intel -fomit-frame-pointer"
F90FLAGS="$F90FLAGS_noOPT -O3"
FFTFLAGS="-DMKL"
FFTLIBS="/opt/intel/mkl/10.0.011/include/mkl_dfti.f90 -lmkl_em64t -lguide"
#FFTFLAGS="-DMKL"
#FFTLIBS="/opt/intel/mkl/10.0.011/include/mkl_dfti.f90 -lmkl_em64t -lguide"
mpiifort -warn all -static -o energromacs spline.f90 fft_iface.f90 enganal.f $F90FLAGS -cpp -Dtrjctry -DGROMACS ${FFTFLAGS} ${FFTLIBS} || exit 1
mpiifort -warn all -o energromacs_p spline.f90 fft_iface.f90 enganal.f $F90FLAGS -cpp -Dtrjctry -DGROMACS ${FFTFLAGS} ${FFTLIBS} -L$HOME/lib -lprofiler || exit 1
mpiifort -warn all -o energromacs_noMPI spline.f90 fft_iface.f90 enganal.f $F90FLAGS_noOPT -cpp -O0 -Dtrjctry -DGROMACS -DMKL -DnoMPI /opt/intel/mkl/10.0.011/include/mkl_dfti.f90 -lmkl_em64t -lguide -L$HOME/lib -lprofiler || exit 1
mpiifort -warn all -static -o enernamd spline.f90 fft_iface.f90 enganal.f $F90FLAGS -cpp -Dtrjctry -DNAMD ${FFTFLAGS} ${FFTLIBS} || exit 1
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
mpiifort --version
