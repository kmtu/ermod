
! mpi module

module mpiproc                                                   ! MPI
  implicit none
#ifndef noMPI
  ! MPI
  include "mpif.h"
#endif
  integer ierror,myrank,nprocs                                     ! MPI
contains
  subroutine mpi_setup(type)                                       ! MPI
    character*4 type                                                 ! MPI
#ifndef noMPI
    real(4) :: cputime                                               ! MPI
    real(8), save :: walltime                                        ! MPI
    if(type.eq.'init') then                                          ! MPI
       call mpi_init(ierror)                                          ! MPI
#ifdef PERF
       walltime = MPI_WTIME()
#endif
    endif
    if(type.eq.'stop') then                                          ! MPI
#ifdef PERF
       call CPU_TIME(cputime)
       print *, "rank = ", myrank, ", CPUtime = ", cputime
       if(myrank == 1) then
          print *, "Wall-clock time: ", MPI_WTIME() - walltime
       endif
#endif
       call mpi_finalize(ierror)                                      ! MPI
    endif
#endif
    return
  end subroutine mpi_setup                                                   ! MPI

  subroutine mpi_info                                              ! MPI
    nprocs=1
    myrank=0
#ifndef noMPI
    call mpi_comm_size(mpi_comm_world,nprocs,ierror)                 ! MPI
    call mpi_comm_rank(mpi_comm_world,myrank,ierror)                 ! MPI
#endif
    return
  end subroutine mpi_info                                                   ! MPI

  subroutine mpi_abend()
    integer :: ierror
#ifndef noMPI
    call mpi_abort(mpi_comm_world, 1, ierror)                        ! MPI
#endif
  end subroutine mpi_abend

end module mpiproc                                                       ! MPI
