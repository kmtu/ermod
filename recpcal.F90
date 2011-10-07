! -*- F90 -*-
module reciprocal
  use fft_iface, only: fft_handle
  implicit none
  integer :: rc1min,rc1max,rc2min,rc2max,rc3min,rc3max
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
    use engmain, only:  nummol,maxsite,numatm,numsite,sluvid,&
         cltype,screen,splodr,charge,&
         ms1max,ms2max,ms3max,&
         sitepos,invcl,volume,&
         pi
    use spline, only: spline_init
    use fft_iface, only: fft_init_ctr, fft_init_rtc, &
         fft_set_size
    implicit none
    integer, intent(in) :: slvmax, tagpt(slvmax)
    integer :: m, k, rci, rcimax, spi
    complex :: rcpi
    integer :: gridsize(3), ptrnk

    allocate( slvtag(nummol) )
    do m=1,nummol
       slvtag(m)=-1
    enddo
    ptrnk=0
    do k=1,slvmax
       m=tagpt(k)
       slvtag(m)=ptrnk+1
       ptrnk=ptrnk+numsite(m)
    enddo

    rc1min=0 ; rc1max=ms1max-1
    rc2min=0 ; rc2max=ms2max-1
    rc3min=0 ; rc3max=ms3max-1
    ccesize = ms1max / 2 + 1; ccemax = ccesize - 1
    call spline_init(splodr)
    allocate( splslv(0:splodr-1,3,ptrnk),grdslv(3,ptrnk) )
    allocate( cnvslt(rc1min:rc1max,rc2min:rc2max,rc3min:rc3max) )
    ! initialize spline table for all axes
    allocate( splfc1(rc1min:rc1max) )
    allocate( splfc2(rc2min:rc2max) )
    allocate( splfc3(rc3min:rc3max) )
    call init_spline_axis(rc1min, rc1max, splfc1(rc1min:rc1max))
    call init_spline_axis(rc2min, rc2max, splfc2(rc2min:rc2max))
    call init_spline_axis(rc3min, rc3max, splfc3(rc3min:rc3max))
    gridsize(1) = ms1max
    gridsize(2) = ms2max
    gridsize(3) = ms3max
    call fft_set_size(gridsize)
    allocate( engfac(rc1min:ccemax,rc2min:rc2max,rc3min:rc3max) )
    allocate( rcpslt(rc1min:ccemax,rc2min:rc2max,rc3min:rc3max) )
    ! init fft
    call fft_init_rtc(handle_r2c, cnvslt, rcpslt)
    call fft_init_ctr(handle_c2r, rcpslt, cnvslt)
  end subroutine recpcal_init

  subroutine init_spline_axis(imin, imax, splfc)
    use engmain, only: splodr, pi
    use spline, only: spline_value
    implicit none
    integer, intent(in) :: imin, imax
    real, intent(out) :: splfc(imin:imax)
    real :: chr, factor, rtp2
    real :: cosk, sink
    complex :: rcpi
    integer :: rci, spi
    do rci=imin,imax
       rcpi=(0.0e0,0.0e0)
       do spi=0,splodr-2
          chr=spline_value(real(spi+1))
          rtp2=2.0e0*pi*real(spi*rci)/real(imax+1)
          cosk=chr*cos(rtp2)
          sink=chr*sin(rtp2)
          rcpi=rcpi+cmplx(cosk,sink)
       end do
       factor=real(rcpi*conjg(rcpi))
       splfc(rci)=factor
    end do
  end subroutine init_spline_axis

  subroutine recpcal_spline_greenfunc()
    use engmain, only: invcl, ms1max, ms2max, ms3max, splodr, volume, screen, pi
    implicit none
    integer :: rc1, rc2, rc3, rci, m, rcimax
    real :: factor, rtp2, chr
    real :: inm(3), xst(3)
    do rc3 = rc3min, rc3max
       do rc2 = rc2min, rc2max
          do rc1 = rc1min, ccemax
             factor=0.0e0
             if(rc1 == 0 .and. rc2 == 0 .and. rc3 == 0) cycle
             do m=1,3
                if(m.eq.1) rci = rc1
                if(m.eq.2) rci = rc2
                if(m.eq.3) rci = rc3

                if(m.eq.1) rcimax=ms1max
                if(m.eq.2) rcimax=ms2max
                if(m.eq.3) rcimax=ms3max
                if((mod(splodr,2).eq.1).and.(2*abs(rci).eq.rcimax)) then
                   go to 3219
                endif
                if(rci.le.rcimax/2) inm(m)=real(rci)
                if(rci.gt.rcimax/2) inm(m)=real(rci-rcimax)
             end do
             do m=1,3
                xst(m)=dot_product(invcl(:, m), inm(:))
             end do
             rtp2=xst(1)*xst(1)+xst(2)*xst(2)+xst(3)*xst(3)
             chr=pi*pi*rtp2/screen/screen
             factor=exp(-chr)/rtp2/pi/volume
             rtp2=splfc1(rc1)*splfc2(rc2)*splfc3(rc3)
             factor=factor/rtp2
3219         continue
             engfac(rc1,rc2,rc3)=factor
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

    svi=slvtag(i)
    if(svi <= 0) call halt_with_error('eng')

    stmax=numsite(i)
    call calc_spline_molecule(i, stmax, splslv(:,:,svi:svi+stmax-1), grdslv(:,svi:svi+stmax-1))
  end subroutine recpcal_prepare_solvent

  subroutine recpcal_prepare_solute(tagslt)
    use engmain, only: ms1max, ms2max, ms3max, sitepos, invcl, numsite, splodr, specatm, charge
    use fft_iface, only: fft_ctr, fft_rtc
    implicit none
    integer, intent(in) :: tagslt
    real :: xst(3), inm(3)
    integer :: rc1, rc2, rc3, rci, sid, m, k, ati, cg1, cg2, cg3, &
         stmax, ptrnk, rcimax, svi, uvi, spi, ccenodup
    real :: factor, rtp2, chr
    complex :: rcpi, rcptemp
    real, allocatable :: splval(:,:,:)
    integer, allocatable :: grdval(:,:)

    stmax=numsite(tagslt)
    allocate( splval(0:splodr-1,3,stmax),grdval(3,stmax) )
    call calc_spline_molecule(tagslt, stmax, splval(:,:,1:stmax), grdval(:,1:stmax))
    cnvslt(:, :, :)=0.0e0
    do sid=1,stmax
       ati=specatm(sid,tagslt)
       chr=charge(ati)
       do cg3=0,splodr-1
          do cg2=0,splodr-1
             do cg1=0,splodr-1
                rc1=modulo(grdval(1,sid)-cg1,ms1max)
                rc2=modulo(grdval(2,sid)-cg2,ms2max)
                rc3=modulo(grdval(3,sid)-cg3,ms3max)
                factor=chr*splval(cg1,1,sid)*splval(cg2,2,sid)&
                     *splval(cg3,3,sid)
                cnvslt(rc1, rc2, rc3) = cnvslt(rc1, rc2, rc3) + factor
             end do
          end do
       end do
    end do

    call fft_rtc(handle_r2c, cnvslt, rcpslt)                         ! 3D-FFT

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

    rcpslt(:, :, :) = engfac(:, :, :) * rcpslt(:, :, :)
    call fft_ctr(handle_c2r, rcpslt, cnvslt)                    ! 3D-FFT

    deallocate( splval,grdval )
  end subroutine recpcal_prepare_solute

  subroutine calc_spline_molecule(imol, stmax, store_spline, store_grid)
    use engmain, only: ms1max, ms2max, ms3max, splodr, specatm, sitepos, invcl
    use mpiproc, only: halt_with_error
    use spline, only: spline_value
    implicit none

    integer, intent(in) :: imol, stmax
    real, intent(out) :: store_spline(0:splodr-1, 3, 1:stmax)
    integer, intent(out) :: store_grid(3, 1:stmax)
    
    integer :: sid, ati, rcimax, m, k, rci, spi
    real :: xst(3), inm(3)
    real :: factor, rtp2

    do sid=1,stmax
       ati=specatm(sid,imol)
       xst(:) = sitepos(:,ati)
       do k=1,3
          factor = dot_product(invcl(k,:),xst(:))
          factor = factor - floor(factor)
          inm(k)=factor
       end do
       do m=1,3
          if(m.eq.1) rcimax=ms1max
          if(m.eq.2) rcimax=ms2max
          if(m.eq.3) rcimax=ms3max
          factor=inm(m)*real(rcimax)
          rci=int(factor)
          do spi=0,splodr-1
             rtp2=factor-real(rci-spi)
             store_spline(spi,m,sid)=spline_value(rtp2)
          end do
          store_grid(m,sid)=rci
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
    use mpiproc, only: halt_with_error
    implicit none
    integer, intent(in) :: tagslt, i
    real, intent(inout) :: pairep
    integer :: cg1, cg2, cg3, k
    integer :: rc1, rc2, rc3, ptrnk, sid, ati, svi, stmax
    real :: fac1, fac2, fac3, chr
    integer :: grid1
    complex :: rcpt

    pairep=0.0e0
    k=sluvid(tagslt)
    if(k.eq.0) call halt_with_error('fst')
    if(tagslt.eq.i) then              ! solute self-energy
       call recpcal_self_energy(pairep)
    endif

    if(tagslt.ne.i) then              ! solute-solvent pair interaction
       svi=slvtag(i)
       if(svi.le.0) call halt_with_error('eng')
       stmax=numsite(i)
       do sid=1,stmax
          ptrnk=svi+sid-1
          ati=specatm(sid,i)
          chr=charge(ati)
          do cg3=0,splodr-1
             fac1 = chr * splslv(cg3,3,ptrnk)
             rc3=modulo(grdslv(3,ptrnk)-cg3,ms3max)
             do cg2=0,splodr-1
                fac2 = fac1 * splslv(cg2,2,ptrnk)
                rc2=modulo(grdslv(2,ptrnk)-cg2,ms2max)
                grid1=grdslv(1,ptrnk)
                if(grid1 >= splodr-1 .and. grid1 < ms1max) then
                   do cg1=0,splodr-1
                      fac3 = fac2 * splslv(cg1,1,ptrnk)
                      rc1=grid1-cg1
                      pairep=pairep+fac3*real(cnvslt(rc1,rc2,rc3))
                   enddo
                else
                   do cg1=0,splodr-1
                      fac3 = fac2 * splslv(cg1,1,ptrnk)
                      rc1=mod(grid1+ms1max-cg1,ms1max) ! speedhack
                      pairep=pairep+fac3*real(cnvslt(rc1,rc2,rc3))
                   end do
                endif
             end do
          end do
       end do
    endif
  end subroutine recpcal_energy

end module reciprocal
