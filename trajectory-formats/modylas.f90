! -*- F90 -*-

! Trajectory module for Modylas
module trajectory
  type handle
     integer :: iohandle
  end type handle
  
contains

  subroutine init_trajectory()
    implicit none
  end subroutine init_trajectory

  subroutine finish_trajectory()
    implicit none
  end subroutine finish_trajectory

  ! Open trajectory and returns handle as htraj. 
  ! Should open fail, the program abends.
  subroutine open_trajectory(htraj, fname)
    use utility, only: newunit
    implicit none
    type(handle), intent(inout) :: htraj
    character(len=*), intent(in) :: fname

    open(unit=newunit(htraj%iohandle), file=fname, action="READ", form="UNFORMATTED")

  end subroutine open_trajectory

  ! Close trajectory specified by handle
  subroutine close_trajectory(htraj)
    implicit none
    type(handle), intent(inout) :: htraj

    close(htraj%iohandle)

    ! release extra resources if necessary

  end subroutine close_trajectory

  ! Read trajectory and returns [crd] as a coordinates, and [cell] as a periodic cell, represented in Angstrom.
  ! [status] is non-zero if any error occurs. In such a case, [crd] and [cell] can be an arbitrary value.
  ! [cell] may be an arbitrary value if the trajectory does not contain cell information.
  ! The coordinate is not guaranteed to be within a unit cell.
  subroutine read_trajectory(htraj, natom, is_periodic, crd, cell, status)
    use utility, only: angles_to_cell_vector
    implicit none
    type(handle), intent(in) :: htraj
    integer, intent(in) :: natom
    logical, intent(in) :: is_periodic
    real(8), intent(out) :: crd(3,natom)
    real(8), intent(out) :: cell(3,3)
    integer, intent(out) :: status

    integer :: n
    real(8) :: xcell, ycell, zcell, alpha, beta, gamma
    real(8) :: lengths(3), angles(3)

    read(htraj%iohandle, err = 999, end = 999) ! step, time
    read(htraj%iohandle, err = 999, end = 999) n
    if(n /= natom) then
       write(5, "(A,I8,A,I8,A)"), "modylas.f90, read_trajectory: natom = ", natom, ", but trajectory has ", n, " atoms"
       goto 999
    end if
    read(htraj%iohandle, err = 999, end = 999) crd
    read(htraj%iohandle, err = 999, end = 999) ! # of Nose-Hoover chains (temperature)
    read(htraj%iohandle, err = 999, end = 999) ! Nose-Hoover coordinates / velocities
    read(htraj%iohandle, err = 999, end = 999) ! # of Nose-Hoover chains (pressure)
    read(htraj%iohandle, err = 999, end = 999) ! Nose-Hoover coordinates / velocities
    read(htraj%iohandle, err = 999, end = 999) xcell, ycell, zcell, alpha, beta, gamma

    lengths = (/ xcell, ycell, zcell /)
    angles = (/ alpha, beta, gamma /)
    call angles_to_cell_vector(lengths, angles, cell)

    status = 0
    return
    
999 status = 1
    return
  end subroutine read_trajectory

end module trajectory
