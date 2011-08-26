

module utility
  implicit none

contains
  integer(8) function hash(arr, size) 
    implicit none
    integer, intent(in) :: size
    real, intent(in) :: arr(size)
    integer(8) :: ret
    external hash_double, hash_float
    select case(kind(arr))
    case(4)
       call hash_float(arr, size, ret)
    case(8)
       call hash_double(arr, size, ret)
    case default
       stop "Error: hash(): unknown real type"
    end select
    hash = ret
  end function hash

  ! The following function is a snippet from Fortran wiki and in public domain.
  ! 
  ! This is a simple function to search for an available unit.
  ! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
  ! The UNIT value is returned by the function, and also by the optional
  ! argument. This allows the function to be used directly in an OPEN
  ! statement, and optionally save the result in a local variable.
  ! If no units are available, -1 is returned.
  integer function newunit(unit)
    implicit none
    integer, intent(out), optional :: unit
    ! local
    integer, parameter :: LUN_MIN=100, LUN_MAX=110
    logical :: opened
    integer :: lun
    ! begin
    newunit=-1
    do lun=LUN_MIN,LUN_MAX
       inquire(unit=lun,opened=opened)
       if (.not. opened) then
          newunit=lun
          exit
       end if
    end do
    if (present(unit)) unit=newunit
  end function newunit

  ! convert cell-length & (alpha, beta, gamma) to cell vectors
  subroutine angles_to_cell_vector(cell_len, angles, out_cell_vectors)
    use engmain, only: pi
    implicit none
    real, intent(in) :: cell_len(3)
    real, intent(in) :: angles(3)
    real, intent(out) :: out_cell_vectors(3, 3)
    real :: alpha, beta, gamma, x, y, u, v, w

    alpha = angles(1) * pi / 180.0 ! for b-c axes
    beta  = angles(2) * pi / 180.0 ! for a-c axes
    gamma = angles(3) * pi / 180.0 ! for a-b axes

    ! ~a = (1, 0, 0)
    ! ~b = (x, y, 0)
    ! ~c = (u, v, w)
    ! ~a.~b = x = cos gamma
    ! |~a*~b| = y = sin gamma
    ! ~a.~c = u = cos beta
    ! ~b.~c = xu + yv = cos alpha

    x = cos(gamma)
    y = sin(gamma)
    u = cos(beta)
    v = (cos(alpha) - x * u) / y ! FIXME: potential underflow risk
    w = sqrt(1 - u * u - v * v)  ! FIXME: same above
    
    out_cell_vectors(1, 1) = cell_len(1)
    out_cell_vectors(2, 1) = 0.0
    out_cell_vectors(3, 1) = 0.0

    out_cell_vectors(1, 2) = cell_len(2) * x
    out_cell_vectors(2, 2) = cell_len(2) * y
    out_cell_vectors(3, 2) = 0.0

    out_cell_vectors(1, 3) = cell_len(3) * u
    out_cell_vectors(2, 3) = cell_len(3) * v
    out_cell_vectors(3, 3) = cell_len(3) * w

  end subroutine angles_to_cell_vector
end module utility
  
