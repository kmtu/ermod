
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

  ! FIXME
  subroutine halt_with_error(type)
    use engmain, only: stdout
    character*3 type
    if(type.eq.'typ') write(stdout,991)
    if(type.eq.'num') write(stdout,992)
    if(type.eq.'ins') write(stdout,993)
    if(type.eq.'par') write(stdout,994)
    if(type.eq.'slt') write(stdout,995)
    if(type.eq.'crd') write(stdout,996)
    if(type.eq.'eng') write(stdout,997)
    if(type.eq.'siz') write(stdout,998)
    if(type.eq.'min') write(stdout,999)
    if(type.eq.'ecd') write(stdout,981)
    if(type.eq.'fst') write(stdout,982)
991 format(' The number of solute types is incorrectly set')
992 format(' The number of solute molecules is incorrectly set')
993 format(' The solute numbering is incorrect for insertion')
994 format(' The input parameter is incorrectly set')
995 format(' The input parameter is incorrect for solute')
996 format(' The coordinate system is incorrectly set')
997 format(' Inconsistency is present in the program')
998 format(' The number of energy-coordinate meshes is too large')
999 format(' The minimum of the energy coordinate is too large')
981 format(' The energy-coordinate system is inconsistent')
982 format(' The first particle needs to be the solute')
#ifndef noMPI
    call mpi_abend()                                                ! MPI
#endif
    stop
  end subroutine halt_with_error

end module mpiproc                                                       ! MPI
