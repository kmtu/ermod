! -*- F90 -*-
! ERMOD - Eneregy Representation Module
! Copyright (C) 2000-2012 Nobuyuki Matubayasi
! Copyright (C) 2010-2012 Shun Sakuraba
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

module sysvars
  implicit none

  character*5 :: clcond='merge'

  integer :: pecore=200,numprm=15
  integer :: numsln=10,numref=5,numdiv=-1

  character*3 :: peread='not', uvread='yes'
  character*3 :: slfslt='yes', infchk='not'
  character*4 :: zerosft='orig',wgtfnform='harm'
  character*3 :: refmerge='yes',extsln='lin'
  character*3 :: wgtf2smpl='yes',slncor='not'
  character*3 :: normalize='yes',showdst='not',cumuint='not'
  character*3 :: wrtzrsft='not',readwgtfl='yes'

  real :: inptemp=300.0e0                                     ! Kelvin
  real, parameter :: zero=0.0e0
  real :: error=1.0e-8, tiny=1.0e-10
  integer :: pickgr=3,msemin=1,msemax=5
  real :: mesherr=0.10e0                                      ! kcal/mol
  integer :: maxmesh=30000, large=500000, itrmax=100
  integer :: extthres_soln=1, extthres_refs=1
  
  character(len=1024) :: solndirec='soln'
  character(len=1024) :: refsdirec='refs'
  character(len=1024) :: wgtslnfl='weight_soln'
  character(len=1024) :: wgtreffl='weight_refs'
  character(len=1024) :: slndnspf='engsln'
  character(len=1024) :: slncorpf='corsln'
  character(len=1024) :: refdnspf='engref'
  character(len=1024) :: refcorpf='corref'
  character(len=1024) :: aveuvfile='aveuv.tt'
  character(len=1024) :: ecdinfofl='EcdInfo'
  character(len=1024) :: cumuintfl='cumsfe'
  character(*), parameter :: numbers='0123456789'
  
  integer prmmax,maxsln,maxref,numrun
  integer numslv,ermax
  real temp,kT,slfeng

  real, dimension(:),     allocatable :: nummol
  integer, dimension(:),  allocatable :: rduvmax,rduvcore
  real, dimension(:),     allocatable :: rdcrd,rddst,rddns
  real, dimension(:,:),   allocatable :: rdslc,rdcor
  integer, dimension(:),  allocatable :: rdspec
  real, dimension(:,:,:), allocatable :: chmpt
  real, dimension(:),     allocatable :: aveuv
  real, dimension(:,:),   allocatable :: uvene,blkuv
  integer, dimension(:),  allocatable :: svgrp,svinf
  real, dimension(:),     allocatable :: wgtsln,wgtref
  
  namelist /fevars/ clcond, pecore, numprm, numsln, numref, numdiv, &
       peread, uvread, slfslt, infchk, zerosft, wgtfnform, &
       refmerge, extsln, extthres_soln, extthres_refs, &
       wgtf2smpl, slncor, normalize, showdst, wrtzrsft, readwgtfl, &
       inptemp, pickgr, msemin, msemax, mesherr, &
       maxmesh, large, itrmax, error, tiny, &
       solndirec, refsdirec, wgtslnfl, wgtreffl, &
       slndnspf, slncorpf, refdnspf, refcorpf, &
       cumuint, cumuintfl

contains

  subroutine init_sysvars
    character(len=*), parameter :: parmfname = 'parameters_fe'
    integer, parameter :: iounit = 191
    integer :: iostat
    
    open(unit=iounit, file= parmfname, action='read', status='old', iostat=iostat)
    if(iostat /= 0) goto 99
    read(iounit, nml=fevars)
    close(iounit)
    
99  if(numdiv == -1) numdiv = numsln
    
  end subroutine init_sysvars
end module sysvars

module sysread
  use sysvars, only: clcond,uvread,slfslt,slncor,&
       maxsln,maxref,numrun,&
       numslv,ermax,slfeng,nummol,&
       rdcrd,rddst,rddns,rdslc,rdcor,rdspec,&
       aveuv,uvene,blkuv,wgtsln,wgtref
  character*85 engfile(5)
contains

  subroutine defcond

    use sysvars, only: peread,infchk,readwgtfl,&
         solndirec,refsdirec,wgtslnfl,wgtreffl,slndnspf,aveuvfile,ecdinfofl,&
         numprm,prmmax,numsln,numref,numdiv,&
         inptemp,temp,kT,&
         pecore,maxmesh,large,&
         rduvmax,rduvcore,&
         chmpt,svgrp,svinf

    integer group,inft,prmcnt,iduv,i,j,k,m,pti,cnt
    real factor
    character*85 opnfile

    if((clcond.ne.'basic').and.(clcond.ne.'range') .and. &
       (clcond.ne.'merge')) stop ' The clcond parameter is incorrect'
    !
    if(clcond.ne.'merge') then
       write(6,781) ; read(5,*) engfile(1)
       if(slncor.eq.'yes') then
          write(6,782) ; read(5,*) engfile(2)
       endif
       write(6,783) ; read(5,*) engfile(3)
       write(6,784) ; read(5,*) engfile(4)
       if(peread.eq.'yes') then
          write(6,785) ; read(5,*) engfile(5)
       endif
781    format(' What is the energy correlation in solution?')
782    format(' What is the energy correlation in solution?')
783    format(' What is the energy density for insertion?')
784    format(' What is the energy density for insertion?')
785    format(' Which file describes energy-coordinate parameters?')
       maxsln=1 ; maxref=1 ; numrun=1
    endif
    if(clcond.eq.'merge') then
       maxsln=numsln ; maxref=numref ; numrun=numdiv
    endif
    if(clcond.eq.'basic') prmmax=1
    if(clcond.ne.'basic') prmmax=numprm
    !
    if(clcond.ne.'merge') opnfile=engfile(1)
    if(clcond.eq.'merge') opnfile=trim(solndirec)//'/'//trim(slndnspf)//'.01'
    open(unit=71,file=opnfile,status='old')
    ermax=0 ; numslv=0 ; k=0
    do iduv=1,large
       read(71,*,end=7809) factor,pti,factor
       if(pti.ne.k) then
          numslv=numslv+1 ; k=pti
       endif
       ermax=ermax+1
    end do
7809 continue
    close(71)
    !
    allocate( nummol(numslv) )
    allocate( rduvmax(numslv),rduvcore(numslv) )
    allocate( rdcrd(ermax),rddst(ermax),rddns(ermax) )
    if(slncor.eq.'yes') allocate( rdslc(ermax,ermax) )
    allocate( rdcor(ermax,ermax) )
    allocate( rdspec(ermax) )
    allocate( chmpt(0:numslv,prmmax,numrun),aveuv(numslv) )
    if((uvread.eq.'yes').and.(clcond.eq.'merge')) then
       allocate( uvene(numslv,maxsln),blkuv(0:numslv,numrun) )
    endif
    allocate( svgrp(prmmax),svinf(prmmax) )
    allocate( wgtsln(maxsln),wgtref(maxref) )
    !
    if(peread.ne.'yes') then
       cnt=ermax/numslv
       rduvmax(:)=cnt ; rduvcore(:)=pecore

       ! check consistency
       if(clcond.ne.'merge') opnfile=engfile(5)
       if(clcond.eq.'merge') opnfile=trim(solndirec)//'/'//trim(ecdinfofl)
       open(unit=71, file=opnfile, status='old', err=7899)

       ! Something is wrong ... 
       close(71)
       print *, 'Warning: EcdInfo file exists, but peread is not "yes"'
       print *, "Perhaps you forgot to set peread?"

7899   continue
    endif
    if(peread.eq.'yes') then
       open(unit=71,file=opnfile,status='old')
       k=0
       do iduv=1,ermax
          read(71,*) factor,pti,factor
          if(pti.ne.k) then
             rduvmax(pti)=1 ; k=pti
          else
             rduvmax(pti)=rduvmax(pti)+1
          endif
       end do
       close(71)
       do pti=1,numslv
          rduvcore(pti)=pecore
          if(clcond.ne.'merge') opnfile=engfile(5)
          if(clcond.eq.'merge') opnfile=trim(solndirec)//'/'//trim(ecdinfofl)
          open(unit=71,file=opnfile,status='old')
          read(71,*)        ! comment line
          read(71,*)        ! line for solute-self energy
          do i=1,large
             read(71,*,end=7829) k
             if(k.eq.pti) then
                backspace(71)
                read(71,*) m,(factor,j=1,9),rduvcore(pti)
                exit
             endif
          end do
7829      continue
          close(71)
       end do
       k=0
       do pti=1,numslv
          k=k+rduvmax(pti)
       end do
       if(k.ne.ermax) stop ' The file format is incorrect'
    endif
    !
    if(ermax.gt.maxmesh) stop ' The number of meshes is too large'
    !
    if(clcond.eq.'basic') then
       write(6,*) ' How many data are grouped into one?'
       read(5,*) group
       inft = 0
       if(infchk == 'yes') then
          write(6,*) ' How many large-energy meshes are merged ? (in %)'
          read(5,*) inft
       endif
       svgrp(1)=group ; svinf(1)=inft
       write(6,*) ' What is the temperature in Kelvin?'
       read(5,*) temp
    endif
    if(clcond.ne.'basic') then
       do prmcnt=1,prmmax                 ! inft in % (percentage)
          if(infchk == 'yes') then
             if(prmcnt.eq.1) then  ; group=1  ; inft=0   ; endif
             if(prmcnt.eq.2) then  ; group=1  ; inft=60  ; endif
             if(prmcnt.eq.3) then  ; group=1  ; inft=80  ; endif
             if(prmcnt.eq.4) then  ; group=1  ; inft=100 ; endif
             if(prmcnt.eq.5) then  ; group=2  ; inft=0   ; endif
             if(prmcnt.eq.6) then  ; group=3  ; inft=0   ; endif
             if(prmcnt.eq.7) then  ; group=4  ; inft=0   ; endif
             if(prmcnt.eq.8) then  ; group=5  ; inft=0   ; endif
             if(prmcnt.eq.9) then  ; group=5  ; inft=60  ; endif
             if(prmcnt.eq.10) then ; group=5  ; inft=80  ; endif
             if(prmcnt.eq.11) then ; group=5  ; inft=100 ; endif
             if(prmcnt.eq.12) then ; group=8  ; inft=0   ; endif
             if(prmcnt.eq.13) then ; group=10 ; inft=0   ; endif
             if(prmcnt.eq.14) then ; group=15 ; inft=0   ; endif
             if(prmcnt.eq.15) then ; group=20 ; inft=0   ; endif
          else
             inft = 0
             if(prmcnt <= 10) group = prmcnt
             if(prmcnt > 10) group = 10 + (prmcnt - 10) * 2
          endif
          svgrp(prmcnt)=group ; svinf(prmcnt)=inft
       end do
       temp=inptemp
    endif
    kT=temp*8.314510e-3/4.184e0               ! kcal/mol
    !
    if(uvread.eq.'yes') then
       if(clcond.ne.'merge') then
          write(6,891)
891       format(' What is average solute-solvent energy in solution?')
          read(5,*) (aveuv(pti), pti=1,numslv)
       endif
       if(clcond.eq.'merge') then
          opnfile=trim(solndirec)//'/'//trim(aveuvfile)
          open(unit=82,file=opnfile,status='old')
          do k=1,maxsln
             read(82,*) m,(uvene(pti,k), pti=1,numslv)
          end do
          close(82)
       endif
    endif

    do cnt=1,2
       if(cnt.eq.1) j=maxsln
       if(cnt.eq.2) j=maxref
       if(cnt.eq.1) wgtsln(1:j)=1.0e0
       if(cnt.eq.2) wgtref(1:j)=1.0e0
       if((clcond.eq.'merge').and.(readwgtfl.eq.'yes')) then
          if(cnt.eq.1) opnfile=trim(solndirec)//'/'//trim(wgtslnfl)
          if(cnt.eq.2) opnfile=trim(refsdirec)//'/'//trim(wgtreffl)
          open(unit=81,file=opnfile,status='old')
          do i=1,j
             read(81,*) k,factor
             if(cnt.eq.1) wgtsln(i)=factor
             if(cnt.eq.2) wgtref(i)=factor
          end do
          close(81)
          if(cnt.eq.1) then
            factor=sum(wgtsln(1:j))
            wgtsln(1:j)=wgtsln(1:j)/factor
          endif
          if(cnt.eq.2) then
            factor=sum(wgtref(1:j))
            wgtref(1:j)=wgtref(1:j)/factor
          endif
       endif
    end do

    if(slfslt.eq.'yes') then
       if(clcond.ne.'merge') then
          write(6,881) ; read(5,*) slfeng
881       format(' What is the solute self-energy?')
       endif
       if(clcond.eq.'merge') then
          if(readwgtfl.eq.'not') then
             write(6,882) ; stop
882          format(' readwgtfl needs to be yes when slfslt is yes')
          endif
          slfeng=0.0e0
          opnfile=trim(refsdirec)//'/'//trim(wgtreffl)
          open(unit=81,file=opnfile,status='old')
          do i=1,maxref
             read(81,*) k,factor,factor
             slfeng=slfeng+wgtref(i)*factor
          end do
          close(81)
       endif
    endif

    return
  end subroutine

  subroutine datread(cntrun)
    use sysvars, only: refmerge,zero,tiny,&
         solndirec,refsdirec,slndnspf,slncorpf,refdnspf,refcorpf,numbers
    integer cntrun,slnini,slnfin,refini,reffin,ecmin,ecmax
    integer iduv,iduvp,i,k,m,pti,cnt
    real factor,ampl
    real, allocatable :: corref_temp(:, :)
    character*85 opnfile
    character*3 suffnum

    if(clcond.ne.'merge') then
       slnini=1 ; slnfin=1 ; refini=1 ; reffin=1
    endif
    if(clcond.eq.'merge') then
       if(maxsln.ge.numrun) then
          k=maxsln/numrun
          slnini=(cntrun-1)*k+1 ; slnfin=cntrun*k
       endif
       if(maxsln.lt.numrun) then
          slnini=mod(cntrun-1,maxsln)+1 ; slnfin=slnini
       endif
       if(refmerge.eq.'not') then
          if(maxref.ge.numrun) then
             m=maxref/numrun
             refini=(cntrun-1)*m+1 ; reffin=cntrun*m
          endif
          if(maxref.lt.numrun) then
             refini=mod(cntrun-1,maxref)+1 ; reffin=refini
          endif
       endif
       if(refmerge.eq.'yes') then
          refini=1 ; reffin=maxref
       endif
    endif
    
    rddst(:)=0.0e0
    if(slncor.eq.'yes') rdslc(:,:)=0.0e0
    if((cntrun.eq.1).or.(refmerge.eq.'not')) then
       rddns(:)=0.0e0 ; rdcor(:,:)=0.0e0
    endif
    
    ! FIXME: this part is kinda spaghetti and should be rewritten WITHOUT looping by cnt!
    do cnt=1,4
       if((cnt.eq.2).and.(slncor.ne.'yes')) cycle
       if((cnt.ge.3).and.(cntrun.gt.1) .and. &
            (refmerge.eq.'yes')) cycle
       
       if(cnt.le.2) ecmin=slnini
       if(cnt.le.2) ecmax=slnfin
       if(cnt.ge.3) ecmin=refini
       if(cnt.ge.3) ecmax=reffin
       
       if(cnt.le.2) then
         factor=sum(wgtsln(ecmin:ecmax))
         wgtsln(ecmin:ecmax)=wgtsln(ecmin:ecmax)/factor
       endif
       if(cnt.ge.3) then
         factor=sum(wgtref(ecmin:ecmax))
         wgtref(ecmin:ecmax)=wgtref(ecmin:ecmax)/factor
       endif

       do i=ecmin,ecmax
          if(clcond.ne.'merge') opnfile=engfile(cnt)
          if(clcond.eq.'merge') then
             m=i/10 ; k=i-10*m
             suffnum='.'//numbers(m+1:m+1)//numbers(k+1:k+1)
             if(cnt.eq.1) opnfile=trim(solndirec)//'/'//trim(slndnspf)//suffnum
             if(cnt.eq.2) opnfile=trim(solndirec)//'/'//trim(slncorpf)//suffnum
             if(cnt.eq.3) opnfile=trim(refsdirec)//'/'//trim(refdnspf)//suffnum
             if(cnt.eq.4) opnfile=trim(refsdirec)//'/'//trim(refcorpf)//suffnum
          endif
          if(cnt == 4) open(unit=71,file=opnfile,status='old', form="UNFORMATTED")
          if(cnt /= 4) open(unit=71,file=opnfile,status='old')
          if((cnt.eq.1).or.(cnt.eq.3)) then
             k=0 ; m=0
          endif
          if(cnt == 4) then
             allocate(corref_temp(ermax, ermax))
             read(71) corref_temp
             rdcor(:, :) = rdcor(:, :) + wgtref(i) * corref_temp(:, :)
             deallocate(corref_temp)
          else
             do iduv=1,ermax
                if((cnt.eq.1).or.(cnt.eq.3)) then
                   read(71,*) rdcrd(iduv),pti,factor
                   if(pti.ne.k) then
                      if(factor.gt.zero) then
                         write(6,*) ' Incorrect energy range with species ',pti
                         stop
                      endif
                      k=pti ; m=m+1
                   endif
                   if(cnt.eq.1) rddst(iduv)=rddst(iduv)+wgtsln(i)*factor
                   if(cnt.eq.3) rddns(iduv)=rddns(iduv)+wgtref(i)*factor
                   rdspec(iduv)=m
                endif
                if(cnt.eq.2) then
                   do iduvp=1,ermax
                      read(71,*) factor
                      rdslc(iduvp,iduv)=rdslc(iduvp,iduv)+wgtsln(i)*factor
                   end do
                endif
             end do
          endif
          close(71)
       end do
    end do
!
    if(cntrun.eq.1) write(6,*)
    do pti=1,numslv
       factor=0.0e0 ; ampl=0.0e0
       do iduv=1,ermax
          if(rdspec(iduv).eq.pti) then
             factor=factor+rddst(iduv) ; ampl=ampl+rddns(iduv)
          endif
       end do
       k=0
       if(abs(factor-ampl).gt.tiny) k=1
       if(cntrun.gt.1) then
          if(nint(factor).ne.nint(nummol(pti))) k=1
          if(nint(ampl).ne.nint(nummol(pti))) k=1
       endif
       if(k.eq.1) then
          write(6,718) pti ; stop
       endif
       if(cntrun.eq.1) then
          nummol(pti)=real(nint(factor))
          write(6,719) pti,nint(nummol(pti))
       endif
718    format('  Incorrect normalization at ',i4)
719    format('  Number of the ',i3,'-th solvent  = ',i12)
    end do
    if(cntrun.eq.1) write(6,*)
!
    if((uvread.eq.'yes').and.(clcond.eq.'merge')) then
       do pti=1,numslv
         aveuv(pti)=sum(wgtsln(slnini:slnfin)*uvene(pti,slnini:slnfin))
       end do
       blkuv(1:numslv,cntrun)=aveuv(1:numslv)
       blkuv(0,cntrun)=sum(blkuv(1:numslv,cntrun))
       if(slfslt.eq.'yes') blkuv(0,cntrun)=blkuv(0,cntrun)+slfeng
    endif

    return
  end subroutine datread
end module

module sfecalc
  use sysvars, only: zerosft,wgtfnform,slncor, &
                     numslv,ermax,nummol,kT,itrmax,zero,error, &
                     rduvmax,rduvcore, &
                     rdcrd,rddst,rddns,rdslc,rdcor,rdspec
  integer, dimension(:), allocatable :: idrduv,uvmax
  real, dimension(:),    allocatable :: uvcrd,edist,edens
  real, dimension(:,:),  allocatable :: edscr,ecorr
  integer, dimension(:), allocatable :: uvspec
  real, dimension(:),    allocatable :: slncv,inscv,sdrcv
  real, dimension(:),    allocatable :: zrsln,zrref,zrsdr
  integer gemax
contains
  !
  subroutine syevr_wrap(n, mat, eigval, info)
    integer, intent(in) :: n
    real, intent(inout) :: mat(n, n)
    real, intent(out) :: eigval(n)
    integer, intent(out) :: info
    real, allocatable :: z(:, :)
    real, allocatable :: work(:)
    real :: worksize
    integer :: lwork, liwork
    integer, allocatable :: iwork(:)
    integer, allocatable :: isuppz(:)
    real :: dummyr, abstol
    integer :: dummyi

    allocate(isuppz(2 * n))
    allocate(z(n, n))

    abstol = 0.0
    lwork = -1
    liwork = 10 * n
    allocate(iwork(liwork))
    call DSYEVR('V', 'A', 'U', n, mat, n, dummyr, dummyr, &
         dummyi, dummyi, abstol, dummyi, eigval, &
         z, n, isuppz, worktmp, lwork, iwork, liwork, info)
    if (info /= 0) then
       deallocate(isuppz)
       deallocate(z)
       deallocate(iwork)
       return
    endif

    lwork = worktmp
    allocate(work(lwork))
    call DSYEVR('V', 'A', 'U', n, mat, n, dummyr, dummyr, &
         dummyi, dummyi, abstol, dummyi, eigval, &
         z, n, isuppz, work, lwork, iwork, liwork, info)

    mat(:, :) = z(:, :)

    deallocate(isuppz)
    deallocate(z)
    deallocate(iwork)
    deallocate(work)
  end subroutine syevr_wrap

  subroutine chmpot(prmcnt,cntrun)
    !
    use sysvars, only: uvread,slfslt,normalize,showdst,wrtzrsft, &
                       slfeng,chmpt,aveuv,svgrp,svinf, &
                       pickgr,cumuint,cumuintfl
    !
    integer prmcnt,cntrun,group,inft
    integer iduv,iduvp,pti,cnt,j,k,m
    real factor,ampl,slvfe,uvpot,lcent,lcsln,lcref
    integer, dimension(:), allocatable :: gpnum
    real, dimension(:,:), allocatable :: cumsfe
    integer, parameter :: cumu_io = 51
    !
    group=svgrp(prmcnt) ; inft=svinf(prmcnt)
    !
    allocate( idrduv(ermax),uvmax(numslv) )
    !
    gemax=0
    do pti=1,numslv
       k=rduvcore(pti)*inft/100
       m=(rduvmax(pti)-k)/group
       uvmax(pti)=m ; gemax=gemax+m
       cnt=1 ; j=0
       if(pti.gt.1) then
          do k=1,pti-1
             cnt=cnt+rduvmax(k) ; j=j+uvmax(k)
          end do
       endif
       do iduv=cnt,cnt+rduvmax(pti)-1
          k=(iduv-cnt)/group+1
          if(k.ge.m) k=m
          idrduv(iduv)=k+j
       end do
       if((pti .eq. numslv) .and. (cnt + rduvmax(pti) - 1 .ne. ermax)) then
          print *, "Error: The total no. of meshes does not match with input"
          print *, "(Sum should be", ermax, " but was", cnt + rduvmax(pti) - 1, ")"
          stop
       endif
    end do
    !
    allocate( uvcrd(gemax),edist(gemax),edens(gemax),uvspec(gemax) )
    uvcrd(:)=0.0e0 ; edist(:)=0.0e0 ; edens(:)=0.0e0
    !
    if(slncor.eq.'yes') then
       allocate( edscr(gemax,gemax) )
       edscr(:,:)=0.0e0
    endif
    allocate( ecorr(gemax,gemax) )
    ecorr(:,:)=0.0e0
    !
    allocate( slncv(gemax),inscv(gemax),zrsln(numslv),zrref(numslv) )
    slncv(:)=0.0e0 ; inscv(:)=0.0e0 ; zrsln(:)=0.0e0 ; zrref(:)=0.0e0
    if(slncor.eq.'yes') then
       allocate( sdrcv(gemax),zrsdr(numslv) )
       sdrcv(:)=0.0e0 ; zrsdr(:)=0.0e0
    endif
    !
    allocate( gpnum(gemax) )
    gpnum(:)=0
    do iduv=1,ermax
       k=idrduv(iduv)
       uvcrd(k)=uvcrd(k)+rdcrd(iduv) ; gpnum(k)=gpnum(k)+1
       uvspec(k)=rdspec(iduv)
    end do
    do k=1,gemax
       if(gpnum(k).gt.0) uvcrd(k)=uvcrd(k)/real(gpnum(k))
    end do
    cnt=rduvmax(1) ; k=uvmax(1)
    do pti=1,numslv
       if(pti.gt.1) cnt=cnt+rduvmax(pti)
       if(pti.gt.1) k=k+uvmax(pti)
       uvcrd(k)=rdcrd(cnt)
    end do
    deallocate( gpnum )
    !
    do cnt=1,2
       do iduv=1,ermax
          k=idrduv(iduv)
          if(cnt.eq.1) edist(k)=edist(k)+rddst(iduv)
          if(cnt.eq.2) edens(k)=edens(k)+rddns(iduv)
       end do
       if((cnt.eq.1).and.(slncor.ne.'yes')) goto 1115
       do iduv=1,ermax
          do iduvp=1,ermax
             k=idrduv(iduv) ; m=idrduv(iduvp)
             if(cnt.eq.1) edscr(m,k)=edscr(m,k)+rdslc(iduvp,iduv)
             if(cnt.eq.2) ecorr(m,k)=ecorr(m,k)+rdcor(iduvp,iduv)
          end do
       end do
1115   continue
    end do
    !
    if(normalize.eq.'yes') call distnorm
    if(showdst.eq.'yes')   call distshow
    !
    call getslncv
    call getinscv
    !
    if((uvread.ne.'yes').and.(prmcnt.eq.1)) then
      do pti=1,numslv
        uvpot=0.0e0
        do iduv=1,gemax
          if(uvspec(iduv).eq.pti) uvpot=uvpot+uvcrd(iduv)*edist(iduv)
        enddo
        aveuv(pti)=uvpot
      enddo
    endif
    !
    ! if some correction such as LJ long-range is added, it should be here
    !
    if((cumuint.eq.'yes').and.(group.eq.pickgr).and.(inft.eq.0)) then
      j=gemax/numslv
      if(any(uvmax(1:numslv).ne.j)) stop ' Incorrect file format for storing the running integral'
      allocate( cumsfe(numslv,j) )
    endif
    !
    do pti=1,numslv
       uvpot=0.0e0 ; slvfe=0.0e0
       do iduv=1,gemax
          if(uvspec(iduv).eq.pti) then
             if((edist(iduv).le.zero).and.(edens(iduv).le.zero)) goto 5009
             uvpot=uvpot+uvcrd(iduv)*edist(iduv)
             slvfe=slvfe-kT*(edist(iduv)-edens(iduv))
             lcent=-(slncv(iduv)+zrsln(pti)+uvcrd(iduv))  ! kT*log(edist/edens)
             if((slncor.eq.'yes') .and. &
                  (edist(iduv).gt.zero).and.(edens(iduv).le.zero)) then
                ampl=lcent*edens(iduv)/edist(iduv)
                lcent=ampl-(zrsln(pti)+uvcrd(iduv))&
                     *(1.0e0-edens(iduv)/edist(iduv))
             endif
             slvfe=slvfe+lcent*edist(iduv)
             do cnt=1,2
                if(cnt.eq.1) lcsln=pyhnc(slncv(iduv),cnt)
                if(cnt.eq.2) lcref=pyhnc(inscv(iduv),cnt)
                if((slncor.eq.'yes') .and. (cnt.eq.1) .and. &
                     (edist(iduv).gt.zero).and.(edens(iduv).le.zero)) then
                   lcsln=pyhnc(sdrcv(iduv)+zrsdr(pti),3)
                endif
             end do
             ampl=sfewgt(edist(iduv),edens(iduv))
             factor=ampl*lcsln+(1.0e0-ampl)*lcref
             slvfe=slvfe+kT*factor*(edist(iduv)-edens(iduv))
5009         continue
             if((cumuint.eq.'yes').and.(group.eq.pickgr).and.(inft.eq.0)) then
               m=mod(iduv-1,j)+1
               cumsfe(pti,m)=uvpot+slvfe
             endif
          endif
       end do
       chmpt(pti,prmcnt,cntrun)=slvfe+aveuv(pti)
    end do
    !
    if((cumuint.eq.'yes').and.(group.eq.pickgr).and.(inft.eq.0)) then
      open(cumu_io,file=cumuintfl,status='replace')
      do iduv=1,j
        factor=sum(cumsfe(1:numslv,iduv))
        if(numslv.eq.1) write(cumu_io,511) iduv,uvcrd(iduv),factor
        if(numslv.gt.1) write(cumu_io,511) iduv,uvcrd(iduv),factor,cumsfe(1:numslv,iduv)
      enddo
511   format(i6,g15.5,9999f12.5)
      endfile(cumu_io) ; close(cumu_io) ; deallocate( cumsfe )
    endif
    !
    chmpt(0,prmcnt,cntrun)=sum(chmpt(1:numslv,prmcnt,cntrun))
    if(slfslt.eq.'yes') chmpt(0,prmcnt,cntrun)=chmpt(0,prmcnt,cntrun)+slfeng
    !
    if(wrtzrsft.eq.'yes') then
       do cnt=1,3
          if(cnt.eq.1) write(6,381) (zrsln(pti),pti=1,numslv)
          if(cnt.eq.2) write(6,382) (zrref(pti),pti=1,numslv)
          if((cnt.eq.3).and.(slncor.eq.'yes')) then
             write(6,383) (zrsdr(pti),pti=1,numslv)
          endif
       end do
    endif
381 format('  Zero shift for solution             = ', 9999f12.4)
382 format('  Zero shift for reference solvent    = ', 9999f12.4)
383 format('  Zero shift for solution correlation = ', 9999f12.4)
    !
    deallocate( slncv,inscv,zrsln,zrref )
    if(slncor.eq.'yes') deallocate( sdrcv,zrsdr,edscr )
    deallocate( uvcrd,edist,edens,ecorr,uvspec,idrduv,uvmax )
    !
    return
  end subroutine chmpot

  subroutine getslncv
    use sysvars, only: extsln,extthres_soln,extthres_refs
    integer iduv,iduvp,pti,j,k,m
    real factor,ampl,lcsln,lcref,min_rddst,min_rddns
    real, dimension(:), allocatable :: work
    integer, dimension(:), allocatable :: ext_target

    min_rddst=minval(rddst, mask=rddst.gt.zero)
    min_rddns=minval(rddns, mask=rddns.gt.zero)
    allocate( ext_target(gemax) ) ; ext_target(:)=0
    !
    do iduv=1,gemax
      j=nint(edist(iduv)/min_rddst)
      k=nint(edens(iduv)/min_rddns)
      if((j.lt.extthres_soln).or.(k.lt.extthres_refs)) ext_target(iduv)=1
    enddo
    !
    do iduv=1,gemax
       if(ext_target(iduv).eq.0) then
          factor=edist(iduv)/edens(iduv)
          slncv(iduv)=-kT*log(factor)-uvcrd(iduv)
       endif
    end do
    !
    do iduv=1,gemax
       if(ext_target(iduv).eq.1) then
          if(edist(iduv).le.zero) then
             slncv(iduv)=0.0e0
             cycle
          endif
          pti=uvspec(iduv)
          m=1 ; k=gemax
          do iduvp=1,iduv-1
             if((uvspec(iduvp).eq.pti).and.(m.lt.iduvp).and.&
                (ext_target(iduvp).eq.0)) m=iduvp
          end do
          do iduvp=gemax,iduv+1,-1
             if((uvspec(iduvp).eq.pti).and.(k.gt.iduvp).and.&
                (ext_target(iduvp).eq.0)) k=iduvp
          end do
          !
          if(extsln.eq.'sim') then
             if(abs(m-iduv).lt.abs(k-iduv)) factor=slncv(m)
             if(abs(m-iduv).gt.abs(k-iduv)) factor=slncv(k)
             if(abs(m-iduv).eq.abs(k-iduv)) then
                factor=(slncv(m)+slncv(k))/2.0e0
             endif
          else
             j=k ; if(abs(m-iduv).lt.abs(k-iduv)) j=m
             allocate( work(gemax) ) ; work(:)=0.0e0
             do iduvp=1,gemax
                if((uvspec(iduvp).eq.pti).and.(ext_target(iduvp).eq.0)) then
                   if(iduvp.eq.iduv) stop ' A bug in program or data'
                   factor=uvcrd(iduvp)-uvcrd(j)
                   if(iduvp.lt.iduv) then
                      factor=-factor-2.0e0*(uvcrd(j)-uvcrd(iduv))
                   endif
                   ampl=wgtdst(iduvp,1,'extsl',wgtfnform)
                   work(iduvp)=exp(-factor/kT)*ampl
                endif
             end do
             factor=sum(work, mask=work.gt.zero)
             do iduvp=1,gemax
                work(iduvp)=work(iduvp)/factor
             end do
             factor=0.0e0 ; ampl=0.0e0 ; lcsln=0.0e0 ; lcref=0.0e0
             do iduvp=1,gemax
                if(work(iduvp).gt.zero) then
                   factor=factor+work(iduvp)*uvcrd(iduvp)
                   ampl=ampl+work(iduvp)*uvcrd(iduvp)*uvcrd(iduvp)
                   lcsln=lcsln+work(iduvp)*slncv(iduvp)
                   lcref=lcref+work(iduvp)*uvcrd(iduvp)*slncv(iduvp)
                endif
             end do
             work(1)=(ampl*lcsln-factor*lcref)/(ampl-factor*factor)
             work(2)=(lcref-factor*lcsln)/(ampl-factor*factor)
             factor=work(1)+work(2)*uvcrd(iduv)
             deallocate( work )
          endif
          slncv(iduv)=factor
       endif
    end do
    !
    deallocate( ext_target )
    !
    do pti=1,numslv
       select case(zerosft)
       case('orig')
          factor=0.0e0
       case('mxco')
          lcref=-cvfcen(pti,2,'uvcrd','smpl','yes')
          lcsln=cvfcen(pti,1,'slncv',wgtfnform,'not')
          ampl=wgtmxco(pti)
          factor=ampl*lcsln+(1.0e0-ampl)*lcref
       case('zero')
          factor=cvfcen(pti,1,'slncv',wgtfnform,'yes')
       case('cntr')
          factor=cvfcen(pti,1,'slncv',wgtfnform,'not')
       case default
          stop ' zerosft not properly set '
       end select
       do iduv=1,gemax
          if(uvspec(iduv).eq.pti) slncv(iduv)=slncv(iduv)-factor
       end do
       zrsln(pti)=factor
    end do
    !
    return
  end subroutine getslncv
  !
  !
  subroutine getinscv
    !
    integer iduv,iduvp,pti,cnt,wrksz,k
    real factor,ampl,lcsln,lcref
    real, dimension(:),   allocatable :: work,egnvl,zerouv
    real, dimension(:,:), allocatable :: edmcr
    !
    do cnt=1,2
       if((cnt.eq.1).and.(slncor.ne.'yes')) cycle
       wrksz=gemax*gemax
       allocate( work(wrksz),egnvl(gemax),edmcr(gemax,gemax) )
       if(cnt.eq.1) edmcr(:,:)=edscr(:,:)
       if(cnt.eq.2) edmcr(:,:)=ecorr(:,:)
       do iduv=1,gemax
          do iduvp=1,gemax
             if(cnt.eq.1) then
                factor=edist(iduvp) ; ampl=edist(iduv)
             endif
             if(cnt.eq.2) then
                factor=edens(iduvp) ; ampl=edens(iduv)
             endif
             lcref=edmcr(iduvp,iduv)-factor*ampl
             if((factor.le.zero).or.(ampl.le.zero)) then
                if(iduv.eq.iduvp) lcref=1.0e0
                if(iduv.ne.iduvp) lcref=0.0e0
             endif
             edmcr(iduvp,iduv)=lcref
          end do
       end do
       !
       call syevr_wrap(gemax, edmcr, egnvl, k)
       ! call DSYEV('V','U',gemax,edmcr,gemax,egnvl,work,wrksz,k)
       pti=numslv+1
       do iduv=pti,gemax
          factor=0.0e0
          do iduvp=1,gemax
             if(cnt.eq.1) ampl=edist(iduvp)
             if(cnt.eq.2) ampl=edens(iduvp)
             if(ampl.gt.zero) factor=factor+(edist(iduvp)-edens(iduvp))&
                  *edmcr(iduvp,iduv)
          end do
          work(iduv)=factor/egnvl(iduv)
       end do
       do iduv=1,gemax
          factor=0.0e0
          do iduvp=pti,gemax
             factor=factor+edmcr(iduv,iduvp)*work(iduvp)
          end do
          if(cnt.eq.1) sdrcv(iduv)=-kT*factor
          if(cnt.eq.2) inscv(iduv)=-kT*factor
       end do
       deallocate( work,egnvl,edmcr )
       !
       allocate( zerouv(numslv) )
       do pti=1,numslv
          select case(zerosft)
          case('orig')
             k=zeroec(pti,cnt)
             if(cnt.eq.1) factor=sdrcv(k)
             if(cnt.eq.2) factor=inscv(k)
          case('mxco')
             factor=cvfcen(pti,cnt,'inscv','smpl','yes')
          case default
             factor=0.0e0
          end select
          zerouv(pti)=factor
       end do
       !
       do iduv=1,gemax
          if(cnt.eq.1) ampl=edist(iduv)
          if(cnt.eq.2) ampl=edens(iduv)
          if(ampl.gt.zero) then
             factor=-kT*(edist(iduv)-edens(iduv))/ampl
             if(cnt.eq.1) lcref=factor-sdrcv(iduv)
             if(cnt.eq.2) lcref=factor-inscv(iduv)
          endif
          if(ampl.le.zero) lcref=0.0e0
          if(cnt.eq.1) sdrcv(iduv)=lcref
          if(cnt.eq.2) inscv(iduv)=lcref
       end do
       !
       do pti=1,numslv
          select case(zerosft)
          case('orig')
             factor=-zerouv(pti)
          case('mxco')
             lcref=-zerouv(pti)
             lcsln=cvfcen(pti,cnt,'inscv',wgtfnform,'not')
             ampl=wgtmxco(pti)
             factor=ampl*lcsln+(1.0e0-ampl)*lcref
          case('zero')
             factor=cvfcen(pti,cnt,'inscv',wgtfnform,'yes')
          case('cntr')
             factor=cvfcen(pti,cnt,'inscv',wgtfnform,'not')
          case default
             stop ' zerosft not properly set '
          end select
          do iduv=1,gemax
             if(uvspec(iduv).eq.pti) then
                if(cnt.eq.1) sdrcv(iduv)=sdrcv(iduv)-factor
                if(cnt.eq.2) inscv(iduv)=inscv(iduv)-factor
             endif
          end do
          if(cnt.eq.1) zrsdr(pti)=factor
          if(cnt.eq.2) zrref(pti)=factor
       end do
       !
       deallocate( zerouv )
    end do
    !
    return
  end subroutine getinscv
  !
  !
  real function wgtmxco(pti)
    integer pti
    real numpt
    numpt=nummol(pti)
    wgtmxco=1.0e0/numpt
    return
  end function wgtmxco
  !
  !
  real function cvfcen(pti,cnt,systype,wgttype,engtype)
    integer pti,cnt,iduv,errtag
    real factor,cvfnc
    real, dimension(:), allocatable :: weight
    character*5 systype
    character*4 wgttype
    character*3 engtype
    allocate( weight(gemax) )
    call getwght(weight,pti,cnt,systype,wgttype,engtype)
    factor=0.0e0 ; errtag=0
    do iduv=1,gemax
       if(uvspec(iduv).eq.pti) then
          select case(systype)
          case('slncv')
             if(cnt.eq.1) cvfnc=slncv(iduv)
             if(cnt.eq.2) errtag=1
          case('inscv')
             if(cnt.eq.1) cvfnc=sdrcv(iduv)
             if(cnt.eq.2) cvfnc=inscv(iduv)
          case('uvcrd')
             cvfnc=uvcrd(iduv)
          case default
             errtag=1
          end select
          if(errtag.ne.0) stop ' Bug in the program'
          factor=factor+cvfnc*weight(iduv)
       endif
    end do
    cvfcen=factor
    deallocate( weight )
    return
  end function cvfcen
  !
  !
  subroutine getwght(weight,pti,cnt,systype,wgttype,engtype)
    integer pti,cnt,iduv
    real weight(gemax),minuv,ampl
    character*5 systype
    character*4 wgttype
    character*3 engtype
    weight(:)=0.0e0
    do iduv=1,gemax
       if(uvspec(iduv).eq.pti) then
          weight(iduv)=wgtdst(iduv,cnt,systype,wgttype)
       endif
    end do
    if(engtype.eq.'yes') then
       minuv=abs(uvcrd(1))
       do iduv=1,gemax
          if((uvspec(iduv).eq.pti).and.(weight(iduv).gt.zero)) then
             if(abs(uvcrd(iduv)).lt.minuv) minuv=abs(uvcrd(iduv))
          endif
       end do
       do iduv=1,gemax
          if(uvspec(iduv).eq.pti) then
             !         ampl=abs(uvcrd(iduv))-minuv
             ampl=nummol(pti)*(abs(uvcrd(iduv))-minuv)
             weight(iduv)=exp(-ampl/kT)*weight(iduv)
          endif
       end do
    endif
    ampl=0.0e0
    do iduv=1,gemax
       if(uvspec(iduv).eq.pti) ampl=ampl+weight(iduv)
    end do
    if(ampl.gt.zero) weight(:)=weight(:)/ampl
    if(ampl.le.zero) then
       write(6,*) ' Zero weight at ',pti ; stop
    endif
  end subroutine getwght
  !
  !
  integer function zeroec(pti,cnt)
    integer iduv,cnt,pti,k
    real factor,ampl,lcsln,lcref
    do iduv=1,gemax-1
       if(uvspec(iduv).eq.pti) then
          if((uvcrd(iduv).le.0).and.(uvcrd(iduv+1).ge.0)) then
             factor=abs(uvcrd(iduv)) ; ampl=uvcrd(iduv+1)
             if(ampl.gt.factor+zero) k=iduv
             if(factor.gt.ampl+zero) k=iduv+1
             if(abs(factor-ampl).le.zero) then
                if(cnt.eq.1) then
                   lcsln=edist(iduv) ; lcref=edist(iduv+1)
                endif
                if(cnt.eq.2) then
                   lcsln=edens(iduv) ; lcref=edens(iduv+1)
                endif
                if(lcsln.ge.lcref) k=iduv
                if(lcsln.lt.lcref) k=iduv+1
             endif
          endif
       endif
    end do
    zeroec=k
    return
  end function zeroec
  !
  !
  real function wgtdst(iduv,cnt,systype,wgttype)
    use sysvars, only: wgtf2smpl
    real fsln,fref,wght,factor
    integer iduv,cnt,jdg,errtag
    character*5 systype
    character*4 wgttype
    fsln=edist(iduv) ; fref=edens(iduv) ; errtag=0
    jdg=1
    if(cnt.eq.2) then
       if(systype.eq.'slncv')  jdg=9
       if(systype.eq.'extsl')  jdg=9
       if(wgttype.eq.'smpl')   jdg=2
       if(wgtf2smpl.eq.'yes')  jdg=2
    endif
    if((jdg.ne.1).and.(jdg.ne.2)) errtag=1
    if(jdg.eq.2) wght=fref
    if(jdg.eq.1) then
       wght=0.0e0
       select case(wgttype)
       case('smpl')
          wght=fsln
          if(cnt.eq.2) errtag=1
       case('geom')
          factor=fsln*fref
          if(factor.gt.zero) wght=sqrt(factor)
       case default              ! corresponding to wgttype = 'harm'
          factor=fsln+fref
          if(factor.gt.zero) wght=fsln*fref/factor
       end select
    endif
    if(errtag.ne.0) stop ' Bug in the program'
    wgtdst=wght
    return
  end function wgtdst
  !
  !
  real function sfewgt(fsln,fref)
    real fsln,fref,wght,factor
    if(fsln.ge.fref) wght=1.0e0
    if(fsln.lt.fref) then
       factor=(fsln-fref)/(fsln+fref)
       wght=1.0e0-factor*factor
    endif
    sfewgt=wght
    return
  end function sfewgt
  !
  !
  real function pyhnc(indpmf,cnt)
    real intg,indpmf,factor
    integer cnt
    factor=indpmf/kT
    if(factor.lt.-zero) then
       if(cnt.eq.1) intg=factor+factor/(exp(-factor)-1.0e0)
       if(cnt.eq.2) intg=(log(1.0e0-factor))*(1.0e0/factor-1.0e0)
       intg=intg+1.0e0
    endif
    if(factor.ge.-zero) intg=factor/2.0e0
    if(cnt.eq.3) then
       if(factor.ge.zero) then
          intg=1.0e0-(log(1.0e0+factor))*(1.0e0/factor+1.0e0)
       endif
       if(factor.lt.zero) intg=-factor/2.0e0
    endif
    pyhnc=intg
    return
  end function pyhnc
  !
  !
  subroutine distnorm
    real, dimension(:), allocatable :: correc
    integer iduv,iduvp,pti,cnt,itrcnt
    real factor,ampl,lcsln,lcref,errtmp
    allocate( correc(gemax) )
    do cnt=1,2
       do pti=1,numslv
          factor=0.0e0
          do iduv=1,gemax
             if(uvspec(iduv).eq.pti) then
                if(cnt.eq.1) factor=factor+edist(iduv)
                if(cnt.eq.2) factor=factor+edens(iduv)
             endif
          end do
          if(factor.gt.zero) factor=nummol(pti)/factor
          if(factor.le.zero) factor=0.0e0
          do iduv=1,gemax
             if(uvspec(iduv).eq.pti) then
                if(cnt.eq.1) edist(iduv)=factor*edist(iduv)
                if(cnt.eq.2) edens(iduv)=factor*edens(iduv)
             endif
          end do
       end do
       !
       if((cnt.eq.1).and.(slncor.ne.'yes')) cycle
       errtmp=error+1.0e0 ; itrcnt=0
       do iduv=1,gemax
          correc(iduv)=1.0e0
       end do
       do while((errtmp.gt.error).and.(itrcnt.le.itrmax))
          do iduv=1,gemax
             lcsln=0.0e0
             do pti=1,numslv
                ampl=0.0e0
                do iduvp=1,gemax
                   if(uvspec(iduvp).eq.pti) then
                      if(cnt.eq.1) lcref=edscr(iduvp,iduv)
                      if(cnt.eq.2) lcref=ecorr(iduvp,iduv)
                      ampl=ampl+correc(iduvp)*lcref
                   endif
                end do
                if(ampl.gt.zero) lcsln=lcsln+nummol(pti)/ampl
             end do
             lcsln=lcsln/real(numslv)
             if(cnt.eq.1) correc(iduv)=lcsln*edist(iduv)
             if(cnt.eq.2) correc(iduv)=lcsln*edens(iduv)
          end do
          do iduv=1,gemax
             do iduvp=1,gemax
                ampl=correc(iduv)*correc(iduvp)
                if(cnt.eq.1) edscr(iduvp,iduv)=ampl*edscr(iduvp,iduv)
                if(cnt.eq.2) ecorr(iduvp,iduv)=ampl*ecorr(iduvp,iduv)
             end do
          end do
          errtmp=0.0e0 ; itrcnt=itrcnt+1
          do iduv=1,gemax
             if(cnt.eq.1) ampl=edist(iduv)
             if(cnt.eq.2) ampl=edens(iduv)
             if(ampl.gt.zero) then
                factor=abs(correc(iduv)-1.0e0)
                if(factor.gt.errtmp) errtmp=factor
             endif
          end do
          if(itrcnt.eq.itrmax) then
             write(6,*) ' The optimzation of the correlation matrix'
             write(6,*) '  did not converge with an error of ',errtmp
             stop
          endif
       end do
    end do
    deallocate( correc )
    return
  end subroutine distnorm
  !
  !
  subroutine distshow
    integer iduv,iduvp,pti,cnt,ecmin,ecmax,i,k
    real factor
    write(6,*)
    do pti=1,numslv
       do cnt=1,2
          if((cnt.eq.1).and.(numslv.gt.1)) write(6,211) pti
          if((cnt.eq.2).and.(numslv.gt.1)) write(6,212) pti
          if((cnt.eq.1).and.(numslv.eq.1)) write(6,*) 'SOLUTION'
          if((cnt.eq.2).and.(numslv.eq.1)) write(6,*) 'INSERTION'
211       format('SOLUTION for',i4,'-th species')
212       format('INSERTION for',i4,'-th species')
          ecmin=gemax ; ecmax=1
          do iduv=1,gemax
             if(uvspec(iduv).eq.pti) then
                i=0
                if((cnt.eq.1).and.(edist(iduv).gt.zero)) i=1
                if((cnt.eq.2).and.(edens(iduv).gt.zero)) i=1
                if(i.eq.1) then
                   if(ecmin.gt.iduv) ecmin=iduv
                   if(ecmax.lt.iduv) ecmax=iduv
                endif
             endif
          end do
          k=0 ; factor=0.0e0
          do iduv=ecmin,ecmax
             i=0
             if((cnt.eq.1).and.(edist(iduv).gt.zero)) i=1
             if((cnt.eq.2).and.(edens(iduv).gt.zero)) i=1
             if(i.eq.1) then
                k=k+1
                if(cnt.eq.1) factor=factor+edist(iduv)
                if(cnt.eq.2) factor=factor+edens(iduv)
             endif
             if((cnt.eq.2).and.(edens(iduv).le.zero)) then
                write(6,215) iduv,uvcrd(iduv)
215             format('     No sampling at ',i5,' with energy ',g14.6)
             endif
          end do
          if(cnt.eq.1) write(6,216) real(k)/real(ecmax-ecmin+1)
          if(cnt.eq.2) write(6,217) real(k)/real(ecmax-ecmin+1)
          write(6,218) factor
216       format('     Nonzero component ratio in solution  = ',g12.4)
217       format('     Nonzero component ratio at insertion = ',g12.4)
218       format('          Number of interacting molecules = ',g12.4)
          if(cnt.eq.1) write(6,221) uvcrd(ecmin)
          if(cnt.eq.1) write(6,222) uvcrd(ecmax)
          if(cnt.eq.2) write(6,223) uvcrd(ecmin)
          if(cnt.eq.2) write(6,224) uvcrd(ecmax)
221       format('     Minimum energy in solution   =',g15.7)
222       format('     Maximum energy in solution   =',g15.7)
223       format('     Minimum energy at insertion  =',g15.7)
224       format('     Maximum energy at insertion  =',g15.7)
       end do
    end do
    return
  end subroutine distshow
  !
end module sfecalc
!
!
!
module opwrite
  use sysvars, only: clcond,uvread,slfslt,infchk,prmmax,numrun,numslv,&
       pickgr,msemin,msemax,mesherr,&
       slfeng,chmpt,aveuv,blkuv,svgrp,svinf
  integer grref
  real, dimension(:), allocatable :: mshdif
contains
  !
  subroutine wrtresl
    !
    integer prmcnt,pti,k,group,inft
    real factor,valcp
    !
    if(slfslt.eq.'yes') write(6,321) slfeng
321 format('  Self-energy of the solute   =   ',f12.4, '  kcal/mol')
    !
    if(clcond.ne.'merge') then
       write(6,*)
       if(numslv.gt.1) then
          if(uvread.ne.'yes') write(6,331) aveuv(1:numslv)
       endif
       factor=sum(aveuv(1:numslv))
       if(slfslt.eq.'yes') factor=factor+slfeng
       write(6,332) factor
331    format('  Solute-solvent energy       =   ', 9999f12.4)
332    format('  Total solvation energy      =   ',f12.4, '  kcal/mol')
    endif
    !
    if(clcond.eq.'basic') then
       if(numslv.gt.1) write(6,351) chmpt(1:numslv,1,1)
       write(6,352) chmpt(0,1,1)
351    format('  Solvation free energy       =   ', 9999f12.4)
352    format('  Total solvation free energy =   ',f12.4, '  kcal/mol')
    endif
    !
    if(clcond.ne.'basic') then
       do prmcnt=1,prmmax
          group=svgrp(prmcnt)
          if(group.eq.pickgr) then
             grref=prmcnt
             exit
          endif
       end do
       allocate( mshdif(msemin:msemax) ) ; mshdif(:)=-1.0e0
    endif
    !
    if(clcond.eq.'range') then
       if(numslv.eq.1) k=0
       if(numslv.gt.1) k=numslv
       do pti=0,k
          write(6,*)
          if(infchk == 'yes') then
            if(pti.eq.0) write(6,591)
            if(pti.ne.0) write(6,592) pti
          else
            if(pti.eq.0) write(6,593)
            if(pti.ne.0) write(6,594) pti
          endif
591       format(' group  inft  solvation free energy   difference')
592       format('               ',i3,'-th component       difference')
593       format(' group    solvation free energy   difference')
594       format('           ',i3,'-th component       difference')
          do prmcnt=1,prmmax
             group=svgrp(prmcnt) ; inft=svinf(prmcnt)
             valcp=chmpt(pti,prmcnt,1)
             factor=valcp-chmpt(pti,grref,1)
             if(infchk == 'yes') then
               write(6,661) group,inft,valcp,factor
             else
               write(6,662) group,valcp,factor
             endif
             if((pti.eq.0).and.(inft.eq.0) .and.&
                  (group >= msemin) .and. (group <= msemax)) mshdif(group)=abs(factor)
          end do
       end do
661    format(i4,i7,f17.5,f18.5)
662    format(i4,f20.5,f18.5)
    endif
    !
    if(clcond.eq.'merge') call wrtmerge
    !
    if(clcond.ne.'basic') then
       factor=0.0e0
       do group=msemin,msemax
          valcp=mshdif(group)
          if((valcp.gt.zero).and.(factor.lt.valcp)) factor=valcp
       end do
       if(factor.gt.mesherr) then
          write(6,*)
          write(6,571) factor,mesherr
       endif
571    format(' Warning: mesh error is ',f8.3,' kcal/mol and is larger',&
            ' than the recommended value of ',g12.3,' kcal/mol')
       deallocate( mshdif )
    endif
    !
    return
  end subroutine wrtresl
  !
  !
  subroutine wrtmerge
    !
    use sysvars, only: large
    integer prmcnt,cntrun,group,inft,pti,i,j,k,m
    real avecp,stdcp,avcp0,factor,slvfe,shcp(large)
    real, dimension(:,:), allocatable :: wrtdata
    !
    allocate( wrtdata(0:numslv,numrun) )
    if(uvread.eq.'yes') then
       wrtdata(0:numslv,1:numrun)=blkuv(0:numslv,1:numrun)
       write(6,*) ; write(6,*) ; write(6,773)
       call wrtcumu(wrtdata) ; write(6,*)
    endif
773 format(' cumulative average & 95% error for solvation energy')
    !
    do pti=0,numslv
       if((numslv.eq.1).and.(pti.ne.0)) cycle
       avcp0=sum(chmpt(pti,grref,1:numrun))/real(numrun)
       do prmcnt=1,prmmax
          group=svgrp(prmcnt)
          inft=svinf(prmcnt)
          avecp=0.0e0
          stdcp=0.0e0
          do cntrun=1,numrun
             slvfe=chmpt(pti,prmcnt,cntrun)
             avecp=avecp+slvfe
             stdcp=stdcp+slvfe*slvfe
          end do
          factor=real(numrun)
          avecp=avecp/factor
          stdcp=sqrt(factor/(factor-1.0e0))&
               *sqrt(stdcp/factor-avecp*avecp)
          stdcp=2.0e0*stdcp/sqrt(factor)
          if(prmcnt.eq.1) then
             write(6,*)
             if(pti.eq.0) then
               if(infchk == 'yes') then
                 write(6,671)
               else
                 write(6,672)
               endif
             endif
671          format(' group  inft  solvation free energy     error',&
                  '          difference')
672          format(' group    solvation free energy     error',&
                  '          difference')
             if(numslv.gt.1) then
                if(pti.eq.0) write(6,661)
                if(pti.ge.1) write(6,662) pti
661             format('  total solvation free energy')
662             format('  contribution from ',i2,'-th solvent component')
             endif
          endif
          if(infchk == 'yes') then
            write(6,673) group,inft,avecp,stdcp,avecp-avcp0
          else
            write(6,674) group,avecp,stdcp,avecp-avcp0
          endif
673       format(i4,i7,f17.5,2f18.5)
674       format(i4,f20.5,2f18.5)

          if((pti.eq.0).and.(inft.eq.0) .and. &
             (group >= msemin) .and. (group <= msemax)) mshdif(group)=abs(avecp-avcp0)
       end do
    end do
    !
    write(6,*) ; write(6,*)
    do pti=0,numslv
       if((numslv.eq.1).and.(pti.ne.0)) cycle
       do prmcnt=1,prmmax
          group=svgrp(prmcnt)
          inft=svinf(prmcnt)
          shcp(1:numrun)=chmpt(pti,prmcnt,1:numrun)
          if(prmcnt.eq.1) then
            if(infchk == 'yes') then
              if(numslv.eq.1) write(6,681)
              if((numslv.gt.1).and.(pti.eq.0)) write(6,682)
              if((numslv.gt.1).and.(pti.ge.1)) then
                write(6,*) ; write(6,683) pti
              endif
681           format(' group  inft   Estimated free energy (kcal/mol)')
682           format(' group  inft   Estimated free energy:',&
                     ' total (kcal/mol)')
683           format(' group  inft   Estimated free energy:',&
                     i2,'-th solvent contribution (kcal/mol)')
            else
              if(numslv.eq.1) write(6,686)
              if((numslv.gt.1).and.(pti.eq.0)) write(6,687)
              if((numslv.gt.1).and.(pti.ge.1)) then
                write(6,*) ; write(6,688) pti
              endif
686           format(' group   Estimated free energy (kcal/mol)')
687           format(' group   Estimated free energy:',&
                     ' total (kcal/mol)')
688           format(' group   Estimated free energy:',&
                     i2,'-th solvent contribution (kcal/mol)')
            endif
          endif
          k=(numrun-1)/5
          if(infchk == 'yes') then
            if(k.eq.0) write(6,121) group,inft,shcp(1:numrun)
            if(k.ge.1) then
              write(6,122) group,inft,shcp(1:5)
              if(k.gt.1) then
                do i=1,k-1 ; write(6,123) (shcp(5*i+m), m=1,5) ; enddo
              endif
            endif
            j=numrun-5*k
            if(j.gt.0) write(6,124) (shcp(m), m=5*k+1,numrun)
121         format(i4,i7,9999f13.4)
122         format(i4,i7,5f13.4)
123         format('           ',5f13.4)
124         format('           ',9999f13.4)
          else
            if(k.eq.0) write(6,126) group,shcp(1:numrun)
            if(k.ge.1) then
              write(6,127) group,shcp(1:5)
              if(k.gt.1) then
                do i=1,k-1 ; write(6,128) (shcp(5*i+m), m=1,5) ; enddo
              endif
            endif
            j=numrun-5*k
            if(j.gt.0) write(6,129) (shcp(m), m=5*k+1,numrun)
126         format(i4,'  ',9999f13.4)
127         format(i4,'  ',5f13.4)
128         format('      ',5f13.4)
129         format('      ',9999f13.4)
          endif
       end do
    end do
    !
    wrtdata(0:numslv,1:numrun)=chmpt(0:numslv,grref,1:numrun)
    write(6,*) ; write(6,*) ; write(6,770)
770 format(' cumulative average & 95% error for solvation free energy')
    call wrtcumu(wrtdata)
    deallocate( wrtdata )

    return
  end subroutine wrtmerge


  subroutine wrtcumu(wrtdata)
    use sysvars, only: large,zero
    integer cntrun,pti
    real avecp,factor,slvfe,shcp(large),wrtdata(0:numslv,numrun)
    real, dimension(:),   allocatable :: runcp,runer
    allocate( runcp(0:numslv),runer(0:numslv) )
    if(numslv.eq.2) write(6,771)
771 format('              total             1st component',&
         '         2nd component')
    do pti=0,numslv
       runcp(pti)=0.0e0 ; runer(pti)=0.0e0
    end do
    do cntrun=1,numrun
       factor=real(cntrun)
       do pti=0,numslv
          slvfe=wrtdata(pti,cntrun)
          runcp(pti)=runcp(pti)+slvfe
          runer(pti)=runer(pti)+slvfe*slvfe
          avecp=runcp(pti)/factor ; shcp(2*pti+1)=avecp
          if(cntrun.gt.1) then
             slvfe=runer(pti)/factor-avecp*avecp
             if(slvfe.le.zero) shcp(2*pti+2)=0.0e0
             if(slvfe.gt.zero) shcp(2*pti+2)=(2.0e0/sqrt(factor))&
                  *sqrt(factor/(factor-1.0e0))*sqrt(slvfe)
          endif
       end do
       if(cntrun.eq.1) then
          do pti=0,numslv
             shcp(pti+1)=shcp(2*pti+1)
          end do
          if(numslv.eq.1) write(6,772) cntrun,shcp(1)
          if(numslv.gt.1) write(6,773) cntrun,shcp(1), &
               (shcp(pti), pti=2,numslv+1)
       endif
       if(cntrun.gt.1) then
          if(numslv.eq.1) write(6,774) cntrun,(shcp(pti), pti=1,2)
          if(numslv.gt.1) write(6,775) cntrun,&
               (shcp(pti), pti=1,2*numslv+2)
       endif
    end do
772 format(i3,f11.4)
773 format(i3,f11.4,9999f22.4)
774 format(i3,2f11.4)
775 format(i3,9999f11.4)
    deallocate( runcp,runer )
    return
  end subroutine wrtcumu
  !
end module opwrite

program sfemain
  use sysvars, only: numrun,prmmax,init_sysvars
  use sysread, only: defcond,datread
  use sfecalc, only: chmpot
  use opwrite, only: wrtresl
  integer cntrun,prmcnt
  call init_sysvars
  call defcond
  do cntrun=1,numrun
     call datread(cntrun)
     do prmcnt=1,prmmax
        call chmpot(prmcnt,cntrun)
     end do
  end do
  call wrtresl
  stop
end program sfemain
