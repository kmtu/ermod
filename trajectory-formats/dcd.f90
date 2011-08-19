! module for DCD I/O, as an example of fortran I/O module
module trajectory
  type handle
     integer :: iohandle
     logical :: have_cell_info
  end type handle
  
contains

  ! Open trajectory and returns handle as htraj. 
  ! Should open fails, the program abends.
  subroutine open_trajectory(htraj, fname)
    use utility, only: newunit
    implicit none
    type(handle), intent(inout) :: htraj
    character(len=*), intent(in) :: fname
    integer(4) :: dcd_header(21)
    integer(4), parameter :: dcd_magic = X'434f5244'

    open(unit=newunit(htraj%iohandle), file=fname, action="READ", form="UNFORMATTED")
    ! Check dcd magic
    read(htraj%iohandle) dcd_header(:)
    if(dcd_header(1) /= dcd_magic) stop "incorrect dcd format (maybe different endian?)"
    htraj%have_cell_info = (dcd_header(12) == 1)
    read(htraj%iohandle)
    read(htraj%iohandle)
  end subroutine open_trajectory

  ! Close trajectory specified by handle
  subroutine close_trajectory(htraj)
    implicit none
    type(handle), intent(inout) :: htraj

    close(htraj%iohandle)
  end subroutine close_trajectory

  ! Read trajectory and returns [crd] as a coordinates, and [cell] as a periodic cell, represented in Angstrom.
  ! [status] is non-zero if any error occurs. In such a case, [crd] and [cell] can be an arbitrary value.
  ! [cell] may be an arbitrary value if the trajectory does not contain cell information.
  ! The coordinate is not guaranteed to be within a unit cell.
  subroutine read_trajectory(htraj, natom, crd, cell, status)
    implicit none
    type(handle), intent(in) :: htraj
    integer, intent(in) :: natom
    real(8), intent(out) :: crd(3,natom)
    real(8), intent(out) :: cell(3,3)
    integer, intent(out) :: status
    real(4), allocatable :: buffer(:)

    cell(:, :) = 0.
    if(htraj%have_cell_info) then
       read(htraj%iohandle) cell(1,1), cell(1,2), cell(2,2), cell(1,3), cell(2,3), cell(3,3)
    end if

    allocate(buffer(natom))
    read(htraj%iohandle) buffer(:)
    crd(1, :) = buffer(:)
    read(htraj%iohandle) buffer
    crd(2, :) = buffer
    read(htraj%iohandle) buffer
    crd(3, :) = buffer
    deallocate(buffer)
  end subroutine read_trajectory

end module trajectory
