
! mpi module

module mpiproc                                                   ! MPI
#ifndef noMPI
  ! MPI
  use mpi
#endif
  implicit none

  ! mpi common variables
  integer ierror, mpistatus(mpi_status_size), myrank, nprocs         ! MPI
  integer, parameter :: tag_cell = 11, tag_coord = 12
  integer :: mpi_comm_activeprocs

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
       call mpi_info
    endif
    if(type.eq.'stop') then                                          ! MPI
#ifdef PERF
       call CPU_TIME(cputime)
       print *, "rank = ", myrank, ", CPUtime = ", cputime
       if(myrank == 0) then
          print *, "Wall-clock time: ", MPI_WTIME() - walltime
       endif
#endif
       call mpi_finalize(ierror)                                      ! MPI
    endif
#endif
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

  subroutine mpi_init_active_group(nproc)
    implicit none
    integer, intent(in) :: nproc
    integer, allocatable :: grpprocs(:)
    integer :: i
    integer :: allgrp, activegrp

    allocate(grpprocs(nproc))
    do i = 1, nproc
       grpprocs(i) = i - 1
    end do
#ifndef noMPI
    call mpi_comm_group(mpi_comm_world, allgrp, ierror)
    call mpi_group_incl(allgrp, nproc, grpprocs, activegrp, ierror)
    call mpi_comm_create(mpi_comm_world, activegrp, mpi_comm_activeprocs, ierror)
    call mpi_group_free(activegrp, ierror)
    call mpi_group_free(allgrp, ierror)

    if(myrank >= nproc .and. mpi_comm_activeprocs /= mpi_comm_null) then
       stop "failed @ mpi_init_active_group"
    endif
#endif
    
    deallocate(grpprocs)
  end subroutine mpi_init_active_group

  subroutine mpi_finish_active_group()
    implicit none
#ifndef noMPI
    if(mpi_comm_activeprocs /= mpi_comm_null) then
       call mpi_comm_free(mpi_comm_activeprocs, ierror)
    end if
#endif
  end subroutine mpi_finish_active_group

  ! helper library for reduce variables
  ! Note: mpi_reduce(mpi_in_place, ...) seems to be allowed only on MPI 2.2+
  subroutine mympi_reduce_real(data, data_size, operation, rootrank)
    implicit none
    integer, intent(in) :: data_size, operation, rootrank
    real, intent(inout) :: data(data_size)
    real, allocatable :: buf(:)
    integer :: mpitype

#ifndef noMPI
    allocate(buf(data_size))
    select case(kind(data))
    case(4)
       mpitype = mpi_real
    case(8)
       mpitype = mpi_double_precision
    case default
       stop "invalid kind(real) value"
    end select

    call mpi_reduce(data, buf, data_size, mpitype, operation, rootrank, mpi_comm_world, ierror)
    data(:) = buf(:)
    deallocate(buf)
#endif
  end subroutine mympi_reduce_real



  ! Stop calculation with error message
  subroutine halt_with_error(type)
    use engmain, only: stdout
    implicit none
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
    if(type.eq.'slb') write(stdout,"(A)"), " Slab condition can only used in periodic system"
    if(type.eq.'bug') write(stdout,"(A)"), " Critical failure in the program detected"
    if(type.eq.'trj') write(stdout,"(A)"), " Trajectory is shorter than specified in MDinfo"
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
    call mpi_abend()                                                     ! MPI
    stop
  end subroutine halt_with_error

  subroutine warning(typ)
    use engmain, only: stdout, force_calculation
    implicit none
    character(len=4), intent(in) :: typ
    if(typ == 'mbin') write(stdout, '(A)'), " Warning: the maximum bin coordinate is too low for this species"
    if(force_calculation) return
    write(stdout, '(A)'), "The program aborts becaue there is a warning"
    write(stdout, '(A)'), "If you wish to force program running, specify 'force_calculation = .true.' in parameters_er."
    call mpi_abend()
    stop
  end subroutine warning
end module mpiproc                                                       ! MPI
