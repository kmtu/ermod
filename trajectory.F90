! module that governs trajectory I/O

module trajectory
  type handle
     integer(8) :: vmdhandle
  end type handle
  
contains

  ! Open trajectory and returns handle as htraj. 
  ! Should open fails, the program abends.
  subroutine open_trajectory(htraj, fname)
    implicit none
    type(handle), intent(inout) :: htraj
    character(len=*), intent(in) :: fname

    integer :: status
    external vmdfio_open_traj

    call vmdfio_open_traj(htraj%vmdhandle, fname, len_trim(fname), status)
    if(status /= 0) then
       stop "vmdfio_open_traj: unable to open trajectory. HISTORY must be a symlink"
    endif
  end subroutine open_trajectory

  ! Close trajectory specified by handle
  subroutine close_trajectory(htraj)
    implicit none
    type(handle), intent(inout) :: htraj

    external vmdfio_close_traj
    call vmdfio_close_traj(htraj%vmdhandle)
  end subroutine close_trajectory

  ! Read trajectory and returns [crd] as a coordinates, and [cell] as a periodic cell, represented in Angstrom.
  ! [status] is non-zero if any error occurs. In such a case, [crd] and [cell] can be an arbitrary value.
  ! [cell] may be an arbitrary value if the trajectory does not contain cell information.
  ! The coordinate is not guaranteed to be within a unit cell.
  subroutine read_trajectory(htraj, natom, is_periodic, crd, cell, status)
    implicit none
    type(handle), intent(in) :: htraj
    integer, intent(in) :: natom
    logical, intent(in) :: is_periodic
    real, intent(out) :: crd(3,natom)
    real, intent(out) :: cell(3,3)
    integer, intent(out) :: status
    
    external vmdfio_read_traj_step

    if(kind(crd) /= 8) stop "vmdfio: write interfacing wrapper"

    call vmdfio_read_traj_step(htraj%vmdhandle, crd, cell, natom, status)
  end subroutine read_trajectory

end module trajectory
