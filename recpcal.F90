! -*- F90 -*-
! ERmod - Eneregy Representation Module
! Copyright (C) 2000-2015 Nobuyuki Matubayasi
! Copyright (C) 2010-2015 Shun Sakuraba
! 
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

module reciprocal
  use fft_iface, only: fft_handle
  implicit none
  integer :: rc1min, rc1max, rc2min, rc2max, rc3min, rc3max
  integer :: ccesize, ccemax
  integer, allocatable :: slvtag(:)
  real,    allocatable :: engfac(:,:,:)
  complex, allocatable :: rcpslv(:,:,:,:)
  complex, allocatable :: rcpslt(:,:,:)
  real,    allocatable :: splslv(:,:,:)
  integer, allocatable :: grdslv(:,:)
  real,    allocatable :: cnvslt(:,:,:)
  real,    allocatable :: splfc1(:), splfc2(:), splfc3(:)
  complex, allocatable :: fft_buf(:, :, :)

  real :: solute_self_energy

  type(fft_handle) :: handle_c2r, handle_r2c

contains
  subroutine recpcal_init(slvmax, tagpt)
    use engmain, only:  nummol, numsite, splodr, ms1max, ms2max, ms3max
    use spline, only: spline_init
    use fft_iface, only: fft_init_ctr, fft_init_rtc, fft_set_size
    implicit none
    integer, intent(in) :: slvmax, tagpt(slvmax)
    integer :: m, k
    integer :: gridsize(3), ptrnk

    allocate( slvtag(nummol) )
    slvtag(:) = -1
    ptrnk = 0
    do k = 1, slvmax
       m = tagpt(k)
       slvtag(m) = ptrnk + 1
       ptrnk = ptrnk + numsite(m)
    enddo

    rc1min = 0 ; rc1max = ms1max - 1
    rc2min = 0 ; rc2max = ms2max - 1
    rc3min = 0 ; rc3max = ms3max - 1
    ccesize = ms1max / 2 + 1; ccemax = ccesize - 1
    call spline_init(splodr)
    allocate( splslv(0:splodr-1, 3, ptrnk), grdslv(3, ptrnk) )
    allocate( cnvslt(rc1min:rc1max, rc2min:rc2max, rc3min:rc3max) )
    ! initialize spline table for all axes
    allocate( splfc1(rc1min: rc1max) )
    allocate( splfc2(rc2min: rc2max) )
    allocate( splfc3(rc3min: rc3max) )
    call init_spline_axis(rc1min, rc1max, splfc1(rc1min:rc1max))
    call init_spline_axis(rc2min, rc2max, splfc2(rc2min:rc2max))
    call init_spline_axis(rc3min, rc3max, splfc3(rc3min:rc3max))
    gridsize(1) = ms1max
    gridsize(2) = ms2max
    gridsize(3) = ms3max
    call fft_set_size(gridsize)
    allocate( engfac(rc1min:ccemax, rc2min:rc2max, rc3min:rc3max) )
    allocate( rcpslt(rc1min:ccemax, rc2min:rc2max, rc3min:rc3max) )
    ! init fft
    call fft_init_rtc(handle_r2c, cnvslt, rcpslt)
    call fft_init_ctr(handle_c2r, rcpslt, cnvslt)
  end subroutine recpcal_init

  subroutine init_spline_axis(imin, imax, splfc)
    use engmain, only: splodr, PI
    use spline, only: spline_value
    implicit none
    integer, intent(in) :: imin, imax
    real, intent(out) :: splfc(imin:imax)
    real :: chr, factor, rtp2
    real :: cosk, sink
    complex :: rcpi
    integer :: rci, spi
    do rci = imin, imax
       rcpi = (0.0, 0.0)
       do spi = 0, splodr - 2
          chr = spline_value(real(spi + 1))
          rtp2 = 2.0 * PI * real(spi * rci) / real(imax + 1)
          cosk = chr * cos(rtp2)
          sink = chr * sin(rtp2)
          rcpi = rcpi + cmplx(cosk, sink)
       end do
       factor = real(rcpi * conjg(rcpi))
       splfc(rci) = factor
    end do
  end subroutine init_spline_axis

  subroutine recpcal_spline_greenfunc()
    use engmain, only: invcl, ms1max, ms2max, ms3max, splodr, volume, screen, PI
    implicit none
    integer :: rc1, rc2, rc3, rci, m, rcimax
    real :: factor, rtp2, chr
    real :: inm(3), xst(3)
    do rc3 = rc3min, rc3max
       do rc2 = rc2min, rc2max
          do rc1 = rc1min, ccemax
             factor = 0.0
             if(rc1 == 0 .and. rc2 == 0 .and. rc3 == 0) cycle
             do m = 1, 3
                if(m == 1) rci = rc1
                if(m == 2) rci = rc2
                if(m == 3) rci = rc3

                if(m == 1) rcimax = ms1max
                if(m == 2) rcimax = ms2max
                if(m == 3) rcimax = ms3max

                if((mod(splodr, 2) == 1) .and. (2*abs(rci) == rcimax)) goto 3219
                if(rci <= rcimax / 2) then
                   inm(m) = real(rci)
                else
                   inm(m) = real(rci - rcimax)
                endif
             end do
             do m = 1, 3
                xst(m) = dot_product(invcl(:, m), inm(:))
             end do
             rtp2 = sum(xst(1:3) ** 2)
             chr = (PI ** 2) * rtp2 / (screen ** 2 )
             factor = exp(-chr) / rtp2 / PI / volume
             rtp2 = splfc1(rc1) * splfc2(rc2) * splfc3(rc3)
             factor = factor / rtp2
3219         continue
             engfac(rc1, rc2, rc3) = factor
          end do
       end do
    end do
    engfac(0, 0, 0) = 0.0
  end subroutine recpcal_spline_greenfunc

  ! note: this routine is named as "solvent", but may include solute molecule, when mutiple solute is used.
  subroutine recpcal_prepare_solvent(i)
    use engmain, only: numsite
    use mpiproc, only: halt_with_error
    implicit none
    integer, intent(in) :: i
    integer :: svi, stmax

    svi = slvtag(i)
    if(svi <= 0) call halt_with_error('rcp_cns')

    stmax = numsite(i)
    call calc_spline_molecule(i, stmax, splslv(:,:,svi:svi+stmax-1), grdslv(:,svi:svi+stmax-1))
  end subroutine recpcal_prepare_solvent

  subroutine recpcal_prepare_solute(tagslt)
    use engmain, only: ms1max, ms2max, ms3max, sitepos, invcl, numsite, splodr, specatm, charge
    use fft_iface, only: fft_ctr, fft_rtc
    use mpiproc, only: perf_time
    implicit none
    integer, intent(in) :: tagslt
    integer :: rc1, rc2, rc3, sid, ati, cg1, cg2, cg3, stmax
    real :: factor, chr
    real, allocatable :: splval(:,:,:)
    integer, allocatable :: grdval(:,:)

    stmax = numsite(tagslt)
    allocate( splval(0:splodr-1, 3, stmax), grdval(3, stmax) )
    call perf_time("kchg")
    call calc_spline_molecule(tagslt, stmax, splval(:,:,1:stmax), grdval(:,1:stmax))
    cnvslt(:,:,:) = 0.0
    do sid = 1, stmax
       ati = specatm(sid, tagslt)
       chr = charge(ati)
       do cg3 = 0, splodr - 1
          do cg2 = 0, splodr - 1
             do cg1 = 0, splodr - 1
                rc1 = modulo(grdval(1, sid) - cg1, ms1max)
                rc2 = modulo(grdval(2, sid) - cg2, ms2max)
                rc3 = modulo(grdval(3, sid) - cg3, ms3max)
                factor = chr * splval(cg1, 1, sid) * splval(cg2, 2, sid) &
                                                   * splval(cg3, 3, sid)
                cnvslt(rc1, rc2, rc3) = cnvslt(rc1, rc2, rc3) + factor
             end do
          end do
       end do
    end do
    call perf_time()

    call perf_time("kfft")
    call fft_rtc(handle_r2c, cnvslt, rcpslt)                         ! 3D-FFT
    call perf_time()

    call perf_time("kslf")
    ! original form is:
    ! 0.5 * sum(engfac(:, :, :) * real(rcpslt_c(:, :, :)) * conjg(rcpslt_c(:, :, :)))
    ! where rcpslt_c(rc1, rc2, rc3) = conjg(rcpslt_buf(ms1max - rc1, ms2max - rc2, ms3max - rc3))
    ! Here we use symmetry of engfac to calculate efficiently
    if(mod(ms1max, 2) == 0) then
       solute_self_energy = &
            sum(engfac(1:(ccemax-1), :, :) * real(rcpslt(1:(ccemax-1), :, :) * conjg(rcpslt(1:(ccemax-1), :, :)))) + &
            0.5 * sum(engfac(0,      :, :) * real(rcpslt(0,      :, :) * conjg(rcpslt(0,      :, :)))) + &
            0.5 * sum(engfac(ccemax, :, :) * real(rcpslt(ccemax, :, :) * conjg(rcpslt(ccemax, :, :))))
    else
       solute_self_energy = &
            sum(engfac(1:ccemax, :, :) * real(rcpslt(1:ccemax, :, :) * conjg(rcpslt(1:ccemax, :, :)))) + &
            0.5 * sum(engfac(0,      :, :) * real(rcpslt(0,      :, :) * conjg(rcpslt(0,      :, :))))
    endif
    call perf_time()

    call perf_time("kuve")
    rcpslt(:, :, :) = engfac(:, :, :) * rcpslt(:, :, :)
    call perf_time()

    call perf_time("kfft")
    call fft_ctr(handle_c2r, rcpslt, cnvslt)                    ! 3D-FFT
    call perf_time()

    deallocate( splval,grdval )
  end subroutine recpcal_prepare_solute

  subroutine calc_spline_molecule(imol, stmax, store_spline, store_grid)
    use engmain, only: ms1max, ms2max, ms3max, splodr, specatm, sitepos, invcl
    use spline, only: spline_value
    implicit none

    integer, intent(in) :: imol, stmax
    real, intent(out) :: store_spline(0:splodr-1, 3, 1:stmax)
    integer, intent(out) :: store_grid(3, 1:stmax)
    
    integer :: sid, ati, rcimax, m, k, rci, spi
    real :: xst(3), inm(3)
    real :: factor, rtp2

    do sid = 1, stmax
       ati = specatm(sid, imol)
       xst(:) = sitepos(:, ati)
       do k = 1, 3
          factor = dot_product(invcl(k,:), xst(:))
          factor = factor - floor(factor)
          inm(k) = factor
       end do
       do m = 1, 3
          if(m == 1) rcimax = ms1max
          if(m == 2) rcimax = ms2max
          if(m == 3) rcimax = ms3max
          factor = inm(m) * real(rcimax)
          rci = int(factor)
          do spi = 0, splodr - 1
             rtp2 = factor - real(rci - spi)
             store_spline(spi, m, sid) = spline_value(rtp2)
          end do
          store_grid(m, sid) = rci
       end do
    end do
  end subroutine calc_spline_molecule

  subroutine recpcal_self_energy(pairep)
    implicit none
    real, intent(out) :: pairep

    pairep = solute_self_energy
  end subroutine recpcal_self_energy

  subroutine recpcal_energy(tagslt, i, pairep)
    use engmain, only: ms1max, ms2max, ms3max, splodr, numsite, specatm, sluvid, charge
    use mpiproc, only: halt_with_error, perf_time
    implicit none
    integer, intent(in) :: tagslt, i
    real, intent(inout) :: pairep
    integer :: cg1, cg2, cg3, k
    integer :: rc1, rc2, rc3, ptrnk, sid, ati, svi, stmax
    real :: fac1, fac2, fac3, chr
    integer :: grid1
    complex :: rcpt

    call perf_time("kuve")
    pairep = 0.0
    k = sluvid(tagslt)
    if(k == 0) call halt_with_error('rcp_fst')
    if(tagslt == i) then              ! solute self-energy
       call recpcal_self_energy(pairep)
    else                              ! solute-solvent pair interaction
       svi = slvtag(i)
       if(svi <= 0) call halt_with_error('rcp_cns')
       stmax = numsite(i)
       do sid = 1, stmax
          ptrnk = svi + sid - 1
          ati = specatm(sid, i)
          chr = charge(ati)
          do cg3 = 0, splodr - 1
             fac1 = chr * splslv(cg3, 3, ptrnk)
             rc3 = modulo(grdslv(3, ptrnk) - cg3, ms3max)
             do cg2 = 0, splodr - 1
                fac2 = fac1 * splslv(cg2, 2, ptrnk)
                rc2 = modulo(grdslv(2, ptrnk) - cg2, ms2max)
                grid1 = grdslv(1, ptrnk)
                if(grid1 >= splodr-1 .and. grid1 < ms1max) then
                   do cg1 = 0, splodr - 1
                      fac3 = fac2 * splslv(cg1, 1, ptrnk)
                      rc1 = grid1 - cg1
                      pairep = pairep + fac3 * real(cnvslt(rc1, rc2, rc3))
                   enddo
                else
                   do cg1 = 0, splodr - 1
                      fac3 = fac2 * splslv(cg1, 1, ptrnk)
                      rc1 = mod(grid1 + ms1max - cg1, ms1max) ! speedhack
                      pairep = pairep + fac3 * real(cnvslt(rc1, rc2, rc3))
                   end do
                endif
             end do
          end do
       end do
    endif
    call perf_time()
  end subroutine recpcal_energy

end module reciprocal
