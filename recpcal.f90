module reciprocal
  implicit none
  integer :: rc1min,rc1max,rc2min,rc2max,rc3min,rc3max
  integer, allocatable :: slvtag(:)
  real,    allocatable :: engfac(:,:,:)
  complex, allocatable :: rcpslv(:,:,:,:)
  complex, allocatable :: rcpslt(:,:,:)
  real,    allocatable :: splslv(:,:,:)
  integer, allocatable :: grdslv(:,:)
  complex, allocatable :: cnvslt(:,:,:)
  real,    allocatable :: splfc1(:), splfc2(:), splfc3(:)
  complex, allocatable :: fft_buf(:, :, :)

contains
  subroutine recpcal(tagslt,i,pairep,slvmax,tagpt,scheme)
    !
    use engmain, only:  nummol,maxsite,numatm,numsite,sluvid,&
         cltype,screen,splodr,charge,&
         ms1max,ms2max,ms3max,&
         specatm,sitepos,invcl,volume
    use spline, only: spline_init, spline_value
    use fft_iface, only: fft_init_ctc, fft_init_inplace, &
         fft_ctc, fft_inplace,&
         fft_set_size
    implicit none
    integer, intent(in) :: tagslt, i, slvmax, tagpt(slvmax)
    real, intent(inout) :: pairep
    character(len=6), intent(in) :: scheme
    integer ptrnk
    integer svi,uvi,ati,sid,stmax,m,k
    integer rc1,rc2,rc3,rci,rcimax,spi,cg1,cg2,cg3,grid1
    real pi,chr,xst(3),inm(3),rtp2,cosk,sink,factor,fac1,fac2,fac3
    complex rcpi,rcpt
    !
    real, dimension(:,:,:),      allocatable :: splval
    integer, dimension(:,:),     allocatable :: grdval
    integer :: si, gridsize(3)
    !
    if(cltype.eq.0) then                                   ! bare Coulomb
       pairep=0.0e0
       return
    endif
    if(cltype == 1) stop "recpcal.f90: Ewald force is no longer supported"
    !
    pi=real(4)*atan(real(1))
    !
    if(scheme.eq.'alloct') then
       call recpcal_init(slvmax, tagpt)
    endif
    !
    if(scheme.eq.'preeng') then
       call recpcal_spline_greenfunc()
    endif
    !
    if((scheme.eq.'slvenv').or.(scheme.eq.'sltsys')) then
       if(tagslt.ne.i) call eng_stop('eng')
       if(scheme.eq.'slvenv') uvi=0                         ! solvent
       if(scheme.eq.'sltsys') uvi=1                         ! solute
       if(uvi.eq.0) svi=slvtag(i)
       if((uvi.eq.0).and.(svi.le.0)) call eng_stop('eng')
       stmax=numsite(i)
       if(cltype.eq.2) then                                 ! PME
          allocate( splval(0:splodr-1,3,stmax),grdval(3,stmax) )
          do sid=1,stmax
             ati=specatm(sid,i)
             do m=1,3
                xst(m)=sitepos(m,ati)
             end do
             do k=1,3
                factor=0.0e0
                do m=1,3
                   factor=factor+invcl(k,m)*xst(m)
                end do
                if(factor.lt.0.0e0) factor=factor+1.0e0
                if(factor.gt.1.0e0) factor=factor-1.0e0
                if((factor.lt.0.0e0).or.(factor.gt.1.0e0)) then
                   call eng_stop('crd')
                endif
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
                   splval(spi,m,sid)=spline_value(rtp2)
                end do
                grdval(m,sid)=rci
             end do
          end do
          if(uvi.eq.0) then
             do sid=1,stmax
                ptrnk=svi+sid-1
                do m=1,3
                   do spi=0,splodr-1
                      splslv(spi,m,ptrnk)=splval(spi,m,sid)
                   end do
                   grdslv(m,ptrnk)=grdval(m,sid)
                end do
             end do
          endif
          if(uvi.gt.0) then
             rcpslt(:, :, :)=(0.0e0,0.0e0)
             do sid=1,stmax
                ati=specatm(sid,i)
                chr=charge(ati)
                do cg3=0,splodr-1
                   do cg2=0,splodr-1
                      do cg1=0,splodr-1
                         rc1=modulo(grdval(1,sid)-cg1,ms1max)
                         rc2=modulo(grdval(2,sid)-cg2,ms2max)
                         rc3=modulo(grdval(3,sid)-cg3,ms3max)
                         factor=chr*splval(cg1,1,sid)*splval(cg2,2,sid)&
                              *splval(cg3,3,sid)
                         rcpi=cmplx(factor,0.0e0)
                         rcpslt(rc1,rc2,rc3)=rcpslt(rc1,rc2,rc3)+rcpi
                      end do
                   end do
                end do
             end do
             ! FIXME: rewrite to real-to-complex transform
             call fft_inplace(rcpslt)                         ! 3D-FFT
             do rc3=rc3min,rc3max
                do rc2=rc2min,rc2max
                   do rc1=rc1min,rc1max
                      rcpi=cmplx(engfac(rc1,rc2,rc3),0.0e0)
                      fft_buf(rc1,rc2,rc3)=rcpi*conjg(rcpslt(rc1,rc2,rc3))
                   end do
                end do
             end do
             call fft_ctc(fft_buf, cnvslt)                    ! 3D-FFT
          endif
          deallocate( splval,grdval )
       endif
    endif
    !
    if(scheme.eq.'energy') then
       pairep=0.0e0
       k=sluvid(tagslt)
       if(k.eq.0) call eng_stop('fst')
       if(tagslt.eq.i) then              ! solute self-energy
          do rc3=rc3min,rc3max
             do rc2=rc2min,rc2max
                do rc1=rc1min,rc1max
                   rcpt=rcpslt(rc1,rc2,rc3)
                   pairep=pairep+engfac(rc1,rc2,rc3)*real(rcpt*conjg(rcpt))
                end do
             end do
          end do
          pairep=pairep/2.0e0
       endif
       if(tagslt.ne.i) then              ! solute-solvent pair interaction
          svi=slvtag(i)
          if(svi.le.0) call eng_stop('eng')
          if(cltype.eq.2) then                               ! PME
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
                      if(grid1 >= splodr-1) then
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
       endif
    endif
    !
    return
  end subroutine recpcal

  subroutine recpcal_init(slvmax, tagpt)
    use engmain, only:  nummol,maxsite,numatm,numsite,sluvid,&
         cltype,screen,splodr,charge,&
         ms1max,ms2max,ms3max,&
         specatm,sitepos,invcl,volume,&
         pi
    use spline, only: spline_init
    use fft_iface, only: fft_init_ctc, fft_init_inplace, &
         fft_ctc, fft_inplace,&
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
    ! allocate fft-buffers
    allocate( fft_buf(rc1min:rc1max,rc2min:rc2max,rc3min:rc3max) )
    gridsize(1) = ms1max
    gridsize(2) = ms2max
    gridsize(3) = ms3max
    call fft_set_size(gridsize)
    allocate( engfac(rc1min:rc1max,rc2min:rc2max,rc3min:rc3max) )
    allocate( rcpslt(rc1min:rc1max,rc2min:rc2max,rc3min:rc3max) )
    ! init fft
    call fft_init_inplace(rcpslt)
    call fft_init_ctc(fft_buf, cnvslt)
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
          do rc1 = rc1min, rc1max
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

  ! FIXME
  subroutine eng_stop(type)
    use engmain, only: io6
    use mpiproc                                                      ! MPI
    character*3 type
    if(type.eq.'typ') write(io6,991)
    if(type.eq.'num') write(io6,992)
    if(type.eq.'ins') write(io6,993)
    if(type.eq.'par') write(io6,994)
    if(type.eq.'slt') write(io6,995)
    if(type.eq.'crd') write(io6,996)
    if(type.eq.'eng') write(io6,997)
    if(type.eq.'siz') write(io6,998)
    if(type.eq.'min') write(io6,999)
    if(type.eq.'ecd') write(io6,981)
    if(type.eq.'fst') write(io6,982)
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
    call mpi_abend()                                                ! MPI
    stop
  end subroutine eng_stop

end module reciprocal
