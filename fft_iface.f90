#ifdef MKL
#  include "mkl_dfti.f90"
#endif

module fft_iface
#ifdef MKL
  use MKL_DFTI
#endif
  implicit none
  
  integer :: fftsize(3)
#ifdef MKL
  type(Dfti_Descriptor), pointer :: desc_ctc, desc_c, desc_rtc
#endif
#ifdef FFTW
  ! unsupported
  integer(8) :: plan_ctc, plan_ctc_backward, plan_inplace, plan_r2c_inplace
#include "fftw3.f"
#endif
  
contains 

  subroutine fft_set_size(fftsize_in)
    integer, intent(in) :: fftsize_in(3)  
    fftsize(:) = fftsize_in(:)
  end subroutine fft_set_size
#ifdef MKL
  ! 3D-FFT, MKL version
  
  ! initialize FFT, complex to complex, out-of-place
  subroutine fft_init_ctc(in, out)
    use MKL_DFTI
    integer :: stat
    complex :: in(fftsize(1), fftsize(2), fftsize(3)), out(fftsize(1), fftsize(2), fftsize(3))
    real :: dummy
    if(kind(dummy) == 8) then
       stat = DftiCreateDescriptor(desc_ctc, DFTI_DOUBLE, DFTI_COMPLEX, 3, fftsize)
    else
       stat = DftiCreateDescriptor(desc_ctc, DFTI_FLOAT, DFTI_COMPLEX, 3, fftsize)
    endif
    if(stat /= 0) stop "MKL-FFT: failed to execute DftiCreateDescriptor"
    stat = DftiSetValue(desc_ctc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    if(stat /= 0) stop "MKL-FFT: failed to set DFTI_NOT_INPLACE"
    stat = DftiCommitDescriptor(desc_ctc)
    if(stat /= 0) stop "MKL-FFT: failed to execute DftiCommitDescriptor"
  end subroutine fft_init_ctc

  subroutine fft_init_ctc_backward(in, out)
    complex :: in(fftsize(1), fftsize(2), fftsize(3)), out(fftsize(1), fftsize(2), fftsize(3))
    call fft_init_ctc(in, out)
  end subroutine fft_init_ctc_backward
  
  subroutine fft_init_inplace(inout)
    use MKL_DFTI
    integer :: stat
    complex :: inout(fftsize(1), fftsize(2), fftsize(3))
    real :: dummy
    if(kind(dummy) == 8) then
       stat = DftiCreateDescriptor(desc_c, DFTI_DOUBLE, DFTI_COMPLEX, 3, fftsize)
    else
       stat = DftiCreateDescriptor(desc_c, DFTI_FLOAT, DFTI_COMPLEX, 3, fftsize)
    endif
    if(stat /= 0) stop "MKL-FFT: failed to execute DftiCreateDescriptor"
    stat = DftiCommitDescriptor(desc_c)
    if(stat /= 0) stop "MKL-FFT: failed to execute DftiCommitDescriptor"
  end subroutine fft_init_inplace

  subroutine fft_init_rtc(in, out)
    use MKL_DFTI
    integer :: stat
    real :: in(fftsize(1), fftsize(2), fftsize(3))
    complex :: out(fftsize(1) / 2 + 1, fftsize(2) / 2 + 1, fftsize(3) / 2 + 1)
    real :: dummy
    if(kind(dummy) == 8) then
       stat = DftiCreateDescriptor(desc_rtc, DFTI_DOUBLE, DFTI_REAL, 3, fftsize)
    else
       stat = DftiCreateDescriptor(desc_rtc, DFTI_FLOAT, DFTI_REAL, 3, fftsize)
    end if
    if(stat /= 0) stop "MKL-FFT: failed to execute DftiCreateDescriptor"
    stat = DftiSetValue(desc_rtc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    if(stat /= 0) stop "MKL-FFT: failed to set DFTI_NOT_INPLACE"
    stat = DftiCommitDescriptor(desc_rtc)
    if(stat /= 0) stop "MKL-FFT: failed to execute DftiCreateDescriptor"
  end subroutine fft_init_rtc

  ! Should be BACKWARD
!   subroutine fft_init_ctr(fftsize, in, out)
!     use MKL_DFTI
!     integer :: stat, fftsize(3)
!     complex(8) :: in(:)
!     real(8) :: out(:)
!     stat = DftiCreateDescriptor(desc_ctr, DFTI_DOUBLE, DFTI_COMPLEX_REAL, 3, fftsize)
!     if(stat /= 0) stop "MKL-FFT: failed to execute DftiCreateDescriptor"
!     stat = DftiSetValue(desc_ctr, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
!     if(stat /= 0) stop "MKL-FFT: failed to set DFTI_NOT_INPLACE"
!     stat = DftiCommitDescriptor(desc_ctr)
!     if(stat /= 0) stop "MKL-FFT: failed to execute DftiCreateDescriptor"
!   end subroutine fft_init_ctr

  ! complex-to-complex, out-of-place FFT
  subroutine fft_ctc(in, out)
    use MKL_DFTI
    integer :: stat
    ! use fortran77-style size to bypass type-check. Note that this invalidates type-check!
    complex, intent(in) :: in(*)
    complex, intent(out) :: out(*)
    stat = DftiComputeForward(desc_ctc, in, out)
    if(stat /= 0) stop "MKL-FFT: failed to execute compute_forward_z_out (oops!)"
  end subroutine fft_ctc

  ! complex-to-complex, out-of-place FFT, backward
  subroutine fft_ctc_backward(in, out)
    use MKL_DFTI
    integer :: stat
    ! use fortran77-style size to bypass type-check. Note that this invalidates type-check!
    complex, intent(in) :: in(*)
    complex, intent(out) :: out(*)
    stat = DftiComputeBackward(desc_ctc, in, out)
    if(stat /= 0) stop "MKL-FFT: failed to execute compute_forward_z_out (oops!)"
  end subroutine fft_ctc_backward

  ! complex-to-complex, in-place FFT
  subroutine fft_inplace(inout)
    use MKL_DFTI
    integer :: stat
    complex, intent(inout) :: inout(*)
    stat = DftiComputeForward(desc_c, inout)
    if(stat /= 0) stop "MKL-FFT: failed to execute compute_forward_z (oops!)"
  end subroutine fft_inplace

  ! real-to-complex, out-of-place FFT
  subroutine fft_rtc(in, out)
    use MKL_DFTI
    integer :: stat
    real, intent(in) :: in(*)
    complex, intent(out) :: out(*)
    stat = DftiComputeForward(desc_rtc, in, out)
    if(stat /= 0) stop "MKL-FFT: failed to execute compute_forward_dz (oops!)"
  end subroutine fft_rtc

!   ! complex-to-real, out-of-place FFT
!   subroutine fft_ctr(in, out)
!     use MKL_DFTI
!     integer :: stat
!     complex(8), intent(in) :: in(:, :, :)
!     real(8), intent(out) :: out(:, :, :)
!     stat = dfti_compute_forward_z_out_d(desc_ctc, in, out)
!   end subroutine fft_ctr

  ! clean-up fft handle
  subroutine fft_cleanup_ctc()
    integer :: stat
    stat = DftiFreeDescriptor(desc_ctc)
    if(stat /= 0) stop "MKL-FFT: failed to execute DftiFreeDescriptor (oops!)"
  end subroutine fft_cleanup_ctc

  ! clean-up fft handle
  subroutine fft_cleanup_inplace()
    integer :: stat
    stat = DftiFreeDescriptor(desc_c)
    if(stat /= 0) stop "MKL-FFT: failed to execute DftiFreeDescriptor (oops!)"
  end subroutine fft_cleanup_inplace

  ! clean-up fft handle
  subroutine fft_cleanup_rtc()
    integer :: stat
    stat = DftiFreeDescriptor(desc_rtc)
    if(stat /= 0) stop "MKL-FFT: failed to execute DftiFreeDescriptor (oops!)"
  end subroutine fft_cleanup_rtc

!   ! clean-up fft handle
!   subroutine fft_cleanup_ctr()
!     stat = DftiFreeDescriptor(desc_ctr)
!     if(stat /= 0) stop "MKL-FFT: failed to execute DftiFreeDescriptor (oops!)"
!   end subroutine fft_cleanup_ctr

#endif
#ifdef FFTW
  ! 3D-FFT, FFTW version
  ! This one is unsupported; use at your own risk.

  subroutine fft_init_ctc(in, out)
    integer :: stat
    complex, intent(in) :: in(fftsize(1), fftsize(2), fftsize(3))
    complex, intent(out) :: out(fftsize(1), fftsize(2), fftsize(3))

    real :: dummy
    if(kind(dummy) == 8) then
       call dfftw_import_system_wisdom(stat)
       ! FIXME: change to FFTW_PATIENT if it becomes a REAL bottleneck
       call dfftw_plan_dft_3d(plan_ctc, fftsize(1), fftsize(2), fftsize(3), &
            in, out, &
            FFTW_FORWARD, FFTW_MEASURE)
    else
       ! FFTW need separate compilation for single precision,
       ! and many people incorrecly report this as a bug.
       stop "fftw for single precision not supported"
    end if
       
  end subroutine fft_init_ctc

  subroutine fft_init_ctc_backward(in, out)
    integer :: stat
    complex, intent(in) :: in(fftsize(1), fftsize(2), fftsize(3))
    complex, intent(out) :: out(fftsize(1), fftsize(2), fftsize(3))

    real :: dummy
    if(kind(dummy) == 8) then
       call dfftw_import_system_wisdom(stat)
       ! FIXME: change to FFTW_PATIENT if it becomes a REAL bottleneck
       call dfftw_plan_dft_3d(plan_ctc_backward, fftsize(1), fftsize(2), fftsize(3), &
            in, out, &
            FFTW_BACKWARD, FFTW_MEASURE)
    else
       ! FFTW need separate compilation for single precision,
       ! and many people incorrecly report this as a bug.
       stop "fftw for single precision not supported"
    end if
       
  end subroutine fft_init_ctc_backward

  subroutine fft_init_inplace(inout)
    integer :: stat
    complex, intent(inout) :: inout(fftsize(1), fftsize(2), fftsize(3))
    real :: dummy
    if(kind(dummy) == 8) then
       call dfftw_import_system_wisdom(stat)
       call dfftw_plan_dft_3d(plan_inplace, fftsize(1), fftsize(2), fftsize(3), &
            inout, inout, &
            FFTW_FORWARD, FFTW_MEASURE)
    else
       stop "fftw for single precision not supported"
    endif

  end subroutine fft_init_inplace

  subroutine fft_init_rtc(in, out)
    integer :: stat
    real(8), intent(in) :: in(fftsize(1), fftsize(2), fftsize(3))
    complex(8), intent(out) :: out(fftsize(1), fftsize(2), fftsize(3))

    stop "implement fft_init_rtc"
  end subroutine fft_init_rtc

  subroutine fft_init_r2c_inplace(inout)
    integer :: stat
    complex, intent(inout) :: inout(fftsize(1)/2+1, fftsize(2), fftsize(3))
    real :: dummy

    if(kind(dummy) == 8) then
       call dfftw_import_system_wisdom(stat)
       call dfftw_plan_dft_r2c_3d(plan_r2c_inplace, fftsize(1), fftsize(2), fftsize(3), &
            inout, inout, &
            FFTW_MEASURE)
    else
       stop "fftw for single precision not supported"
    endif
  end subroutine fft_init_r2c_inplace

  subroutine fft_init_c2r_inplace(inout)
    integer :: stat
    complex, intent(inout) :: inout(fftsize(1)/2+1, fftsize(2), fftsize(3))
    real :: dummy

    if(kind(dummy) == 8) then
       call dfftw_import_system_wisdom(stat)
       call dfftw_plan_dft_c2r_3d(plan_r2c_inplace, fftsize(1), fftsize(2), fftsize(3), &
            inout, inout, &
            FFTW_MEASURE)
    else
       stop "fftw for single precision not supported"
    endif
  end subroutine fft_init_c2r_inplace

  ! complex-to-complex, out-of-place FFT
  subroutine fft_ctc(in, out)
    complex, intent(in) :: in(fftsize(1), fftsize(2), fftsize(3))
    complex, intent(out) :: out(fftsize(1), fftsize(2), fftsize(3))
    real :: dummy
    if(kind(dummy) == 8) then
       call dfftw_execute(plan_ctc)
    else
       stop
    endif
  end subroutine fft_ctc

  subroutine fft_ctc_backward(in, out)
    complex, intent(in) :: in(fftsize(1), fftsize(2), fftsize(3))
    complex, intent(out) :: out(fftsize(1), fftsize(2), fftsize(3))
    real :: dummy
    if(kind(dummy) == 8) then
       call dfftw_execute(plan_ctc_backward)
    else
       stop
    endif
  end subroutine fft_ctc_backward

  subroutine fft_r2c_inplace(inout)
    complex, intent(inout) :: inout(fftsize(1)/2+1, fftsize(2), fftsize(3))
    real :: dummy
    if(kind(dummy) == 8) then
       call dfftw_execute(plan_r2c_inplace)
    else
       stop
    endif
  end subroutine fft_r2c_inplace

  ! complex-to-complex, in-place FFT
  subroutine fft_inplace(inout)
    complex, intent(inout) :: inout(fftsize(1), fftsize(2), fftsize(3))
    real :: dummy
    if(kind(dummy) == 8) then
       call dfftw_execute(plan_inplace)
    else
       stop
    endif

  end subroutine fft_inplace

  ! real-to-complex, out-of-place FFT
  subroutine fft_rtc(in, out)
    real(8), intent(in) :: in(fftsize(1), fftsize(2), fftsize(3))
    complex(8), intent(out) :: out(fftsize(1) / 2 + 1, fftsize(2) / 2 + 1, fftsize(3) / 2 + 1)

    stop "implement fft_rtc"
  end subroutine fft_rtc

  ! clean-up fft handle
  subroutine fft_cleanup_ctc()
    real :: dummy
    if(kind(dummy) == 8) then
       call dfftw_destroy_plan(plan_ctc)
    else
       stop
    endif
  end subroutine fft_cleanup_ctc

  ! clean-up fft handle
  subroutine fft_cleanup_inplace()
    real :: dummy
    if(kind(dummy) == 8) then
       call dfftw_destroy_plan(plan_inplace)
    else
       stop
    endif

  end subroutine fft_cleanup_inplace

  ! clean-up fft handle
  subroutine fft_cleanup_rtc()
    stop "implement fft_cleanup_rtc"
  end subroutine fft_cleanup_rtc
#endif

end module fft_iface
