#!/bin/sh

F90FLAGS="-O3 -xP -mcmodel=large -autodouble -extend-source -shared-intel"
mpiifort -g -o enganal.o -c enganal.f $F90FLAGS -cpp -Dtrjctry -DMARBLE -DMKL 
#mpiifort -o insertion.o -c insertion.f $F90FLAGS -cpp -Dtrjctry -DMARBLE -DMKL 
#mpiifort -o setconf.o -c setconf.f $F90FLAGS -cpp -Dtrjctry -DMARBLE -DMKL 
mpiifort -g -o enermarble enganal.o /opt/intel/mkl/10.0.011/include/mkl_dfti.f90 -lmkl_em64t -lguide
mpiifort -g -o enganal.o -c enganal.f $F90FLAGS -cpp -Dtrjctry -DGROMACS -DMKL 
mpiifort -g -o energromacs enganal.o /opt/intel/mkl/10.0.011/include/mkl_dfti.f90 -lmkl_em64t -lguide
mpiifort -g -o enganal.o -c enganal.f $F90FLAGS -cpp -Dtrjctry -DNAMD -DMKL 
mpiifort -g -o enernamd enganal.o /opt/intel/mkl/10.0.011/include/mkl_dfti.f90 -lmkl_em64t -lguide
#mpiifort -o energromacs -O3 -xP -fpic -shared-intel -mcmodel=large –autodouble -extend-source -cpp -Dtrjctry -DGROMACS -DMKL /opt/intel/mkl/10.0.011/include/mkl_dfti.f90 enganal.f insertion.f setconf.f -lmkl_em64t -lguide -L/usr/lib –lpthread
#mpiifort -o enernamd -O3 -xP -fpic -shared-intel -mcmodel=large -autodouble -extend-source -cpp -Dtrjctry -DNAMD -DMKL /opt/intel/mkl/10.0.011/include/mkl_dfti.f90 enganal.f insertion.f setconf.f -lmkl_em64t -lguide -L/usr/lib –lpthread
