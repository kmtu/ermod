! -*- F90 -*-
!----------------
! best fit RMSD minimization module
! Based on reviews of:
!   Horn B. K. P. "Closed-form solution of 
!     absolute orientation using unit quaternions." 
!     Opt. Soc. Am. A, 4(4), 629 (1987).
! see also:
!   Theobald D. L. "Rapid calculation of RMSDs"
!     Acta Cryst., A61, 478 (2005).
!----------------

module quaternion
  implicit none
contains
  subroutine array_of_quaternion(q, a) 
    real(8), intent(in)  :: q(0:3)
    real(8), intent(out) :: a(3)
    a = q(1:3)
  end subroutine array_of_quaternion
  
  subroutine quaternion_of_array(a, q)
    real(8), intent(in)  :: a(3)
    real(8), intent(out) :: q(0:3)
    q(0) = 0.
    q(1:3) = a
  end subroutine quaternion_of_array

  subroutine prod(q1, q2, p)
    real(8), intent(in)  :: q1(0:3), q2(0:3)
    real(8), intent(out) :: p(0:3)
    real(8) :: temp(3)
    p(0) = q1(0) * q2(0) - dot_product(q1(1:3), q2(1:3))
    call cross_product(q1(1:3), q2(1:3), temp)
    p(1:3) = temp + q1(0) * q2(1:3) + q2(0) * q1(1:3)

  contains
    subroutine cross_product(u, v, r)
      implicit none
      real(8), intent(in) :: u(3), v(3)
      real(8), intent(out) :: r(3)
      r = cshift(u, 1) * cshift(v, -1) - cshift(u, -1) * cshift(v, 1)
    end subroutine cross_product
  end subroutine prod
  
  subroutine conjugate(q, r)
    real(8), intent(in)  :: q(0:3)
    real(8), intent(out) :: r(0:3)
    r(0)   =  q(0)
    r(1:3) = -q(1:3)
  end subroutine conjugate

  subroutine rotate(q, r, res)
    real(8), intent(in)  :: q(1:3), r(0:3)
    real(8), intent(out) :: res(1:3)
    real(8) :: t(0:3), cj(0:3)
    real(8) :: q2(0:3), res2(0:3)
    q2(0) = 0.
    q2(1:3) = q(:)
    call conjugate(r, cj)
    call prod(r, q2,  t)
    call prod(t, cj, res2)
    res(:) = res2(1:3)
  end subroutine rotate
end module quaternion

module bestfit
  implicit none
contains  
  ! find a rotation that satisfies 
  subroutine find_rotation_quaternion(n, refPt, movedPt, masses, rotation)
    implicit none
    external dsyev
    integer :: n, info
    real(8), intent(in) :: refPt(3, n), movedPt(3, n), masses(n)
    real(8), intent(out) :: rotation(0:3)
    real(8) :: inner_prod(3, 3)
    real(8) :: matmax(4, 4)
    real(8) :: eigenvalue(4)
    integer, parameter :: lwork = 256
    real(8) :: work(lwork)
    integer :: i, j, k
    
    do i = 1, 3
       do j = 1, 3
          inner_prod(j, i) = inner_prod(j, i) + sum(refPt(i, :) * movedPt(j, :) * masses(:))
       end do
    end do
    
    !-- Calculate a linear operator that describes sums of products.
    !-- Regrettably, there is no succinct expression to write; we just write it "as-is" in page 635.
    !-- We need upper half of the matrix only, therefore omitting lower.
    !-- (Be aware that LAPACK's matrix is
    !-- a(1, 1) a(1, 2) a(1, 3) ...
    !-- a(2, 1) a(2, 2) a(2, 3) ...) 
    matmax(1, 1) = + inner_prod(1, 1) + inner_prod(2, 2) + inner_prod(3, 3)
    matmax(2, 2) = + inner_prod(1, 1) - inner_prod(2, 2) - inner_prod(3, 3)
    matmax(3, 3) = - inner_prod(1, 1) + inner_prod(2, 2) - inner_prod(3, 3)
    matmax(4, 4) = - inner_prod(1, 1) - inner_prod(2, 2) + inner_prod(3, 3)

    matmax(1, 2) = + inner_prod(2, 3) - inner_prod(3, 2)
    matmax(1, 3) = + inner_prod(3, 1) - inner_prod(1, 3)
    matmax(1, 4) = + inner_prod(1, 2) - inner_prod(2, 1)

    matmax(2, 3) = + inner_prod(1, 2) + inner_prod(2, 1)
    matmax(3, 4) = + inner_prod(2, 3) + inner_prod(3, 2)
    matmax(2, 4) = + inner_prod(3, 1) + inner_prod(1, 3)

    !-- solve eigenvalue and eigenvector
    call dsyev('V', 'U', 4, matmax, 4, eigenvalue, work, lwork, info)

    if (info /= 0) then
       print *, "error on LAPACK/DSYEV"
       print *, "info =", info
    end if

    ! value with maximum eigenvalue is the required vector
    rotation(0:3) = matmax(1:4, 4)
  end subroutine find_rotation_quaternion

  subroutine center_of_mass(n, points, masses, center)
    integer, intent(in) :: n
    real(8), intent(in) :: points(3, n), masses(n)
    real(8), intent(out) :: center(3)
    real(8) :: sumOfMasses
    integer :: i
    
    sumOfMasses = sum(masses)
    do i = 1, 3
       center(i) = dot_product(points(i, :), masses) / sumOfMasses
    end do
  end subroutine center_of_mass

  subroutine fit(n, refcoord, coord, masses, outcoord)
    use quaternion
    integer, intent(in) :: n
    real(8), intent(in) :: refcoord(3, n), masses(n)
    real(8), intent(in) :: coord(3, n)
    real(8), intent(out) :: outcoord(3, n)
    real(8) :: translatedRef(3, n), translatedMoved(3, n)
    real(8) :: centerOfMassRef(3)
    real(8) :: centerOfMassMoved(3)
    real(8) :: rotation(0:3), temp(3)
    integer :: i

    call center_of_mass(n, refcoord, masses, centerOfMassRef)
    call center_of_mass(n, coord,    masses, centerOfMassMoved)
    
    do i = 1, n
       translatedRef  (:, i) = refcoord(:, i) - centerOfMassRef  (:)
       translatedMoved(:, i) = coord   (:, i) - centerOfMassMoved(:)
    end do

    call find_rotation_quaternion(n, translatedRef, translatedMoved, masses, rotation)

    do i = 1,n
       call rotate(translatedMoved(:, i), rotation, temp)
       outcoord(:, i) = temp(:) + centerOfMassRef(:)
    end do
  end subroutine fit
end module bestfit


