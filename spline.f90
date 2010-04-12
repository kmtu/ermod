module spline
  implicit none
  real(8), allocatable :: coeff(:)
  integer :: order
contains

  subroutine spline_init(spline_order)
    integer, intent(in) :: spline_order
    integer :: i, k
    real(8) :: factor
    order = spline_order
    allocate( coeff(0:order) )
    do i = 0,order
       factor=1.0e0
       do k = 1, i ! pass thru when i == 0
          factor = factor * dble(k)
       end do
       do k = 1, order - i ! pass thru when i == order
          factor = factor * dble(k)
       end do
       factor = order / factor
       if(mod(i,2) == 1) factor = -factor
       coeff(i) = factor
    end do
  end subroutine spline_init

  ! FIXME: speed it up
  real(8) function spline_value(rst)
    real(8), intent(in) :: rst
    integer :: i, k
    real(8) :: f
    f = 0.0e0
    if((rst > 0.0) .and. (rst < order)) then
       k = int(rst)
       do i = 0, k
          f = f + coeff(i) * ((rst-i)**(order-1))
       end do
    endif
    spline_value = f
  end function spline_value

  ! never called in usual case
  subroutine spline_cleanup()
    deallocate(coeff)
  end subroutine spline_cleanup
end module spline
