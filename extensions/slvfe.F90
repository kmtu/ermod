! -*- F90 -*-
! ERmod - Eneregy Representation Module
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

  character(len=5) :: clcond = 'merge'

  integer :: pecore = 200, numprm = 15
  integer :: numsln = 10,  numref = 5,  numdiv = -1

  character(len=3) :: peread = 'not',    uvread = 'yes'
  character(len=3) :: slfslt = 'yes',    infchk = 'not',   ljlrc = 'not'
  character(len=4) :: zerosft = 'orig',  wgtfnform = 'harm'
  character(len=3) :: refmerge = 'yes',  extsln = 'lin'
  character(len=3) :: wgtf2smpl = 'yes', slncor = 'not'
  character(len=3) :: normalize = 'yes', showdst= 'not',   cumuint = 'not'
  character(len=3) :: wrtzrsft = 'not',  readwgtfl = 'yes'

  real :: inptemp = 300.0                                     ! Kelvin
  real, parameter :: zero = 0.0
  real :: error = 1.0e-8, tiny = 1.0e-10
  integer :: pickgr = 3, msemin = 1, msemax = 5
  real :: mesherr = 0.10                                      ! kcal/mol
  integer :: maxmesh = 30000, large = 500000, itrmax = 100
  integer :: extthres_soln = 1, extthres_refs = 1
  
  character(len=1024) :: solndirec = 'soln'
  character(len=1024) :: refsdirec = 'refs'
  character(len=1024) :: wgtslnfl  = 'weight_soln'
  character(len=1024) :: wgtreffl  = 'weight_refs'
  character(len=1024) :: slndnspf  = 'engsln'
  character(len=1024) :: slncorpf  = 'corsln'
  character(len=1024) :: refdnspf  = 'engref'
  character(len=1024) :: refcorpf  = 'corref'
  character(len=1024) :: aveuvfile = 'aveuv.tt'
  character(len=1024) :: ecdinfofl = 'EcdInfo'
  character(len=1024) :: cumuintfl = 'cumsfe'
  character(len=10), parameter :: numbers='0123456789'
  
  integer prmmax, maxsln, maxref, numrun
  integer numslv, ermax
  real temp, kT, slfeng

  real, dimension(:),     allocatable :: nummol
  integer, dimension(:),  allocatable :: rduvmax, rduvcore
  real, dimension(:),     allocatable :: rdcrd, rddst, rddns
  real, dimension(:,:),   allocatable :: rdslc, rdcor
  integer, dimension(:),  allocatable :: rdspec
  real, dimension(:,:,:), allocatable :: chmpt
  real, dimension(:),     allocatable :: aveuv
  real, dimension(:,:),   allocatable :: uvene, blkuv
  integer, dimension(:),  allocatable :: svgrp, svinf
  real, dimension(:),     allocatable :: wgtsln, wgtref
  
  namelist /fevars/ clcond, pecore, numprm, numsln, numref, numdiv, &
       peread, uvread, slfslt, infchk, zerosft, wgtfnform, &
       ljlrc, &
       refmerge, extsln, extthres_soln, extthres_refs, &
       wgtf2smpl, slncor, normalize, showdst, wrtzrsft, readwgtfl, &
       inptemp, pickgr, msemin, msemax, mesherr, &
       maxmesh, large, itrmax, error, tiny, &
       solndirec, refsdirec, wgtslnfl, wgtreffl, &
       slndnspf, slncorpf, refdnspf, refcorpf, &
       cumuint, cumuintfl

contains

  subroutine init_sysvars
    implicit none
    character(len=*), parameter :: parmfname = 'parameters_fe'
    integer, parameter :: iounit = 191
    integer :: ioerr
    
    open(unit = iounit, file = parmfname, action = 'read', status = 'old', iostat = ioerr)
    if(ioerr /= 0) goto 99
    read(iounit, nml = fevars)
    close(iounit)
    
99  if(numdiv == -1) numdiv = numsln
    
  end subroutine init_sysvars
end module sysvars

module sysread
  use sysvars, only: clcond, uvread, slfslt, slncor, &
       maxsln, maxref, numrun, &
       numslv, ermax, slfeng, nummol, &
       rdcrd, rddst, rddns, rdslc, rdcor, rdspec, &
       aveuv, uvene, blkuv, wgtsln, wgtref
  implicit none
  character(len=1024) :: engfile(5)
contains

  subroutine defcond
    use sysvars, only: peread, infchk, readwgtfl, &
         solndirec, refsdirec, wgtslnfl, wgtreffl, &
         slndnspf, aveuvfile, ecdinfofl, &
         numprm, prmmax, numsln, numref, numdiv, &
         inptemp, temp, kT, &
         pecore, maxmesh, large, &
         rduvmax, rduvcore, &
         chmpt, svgrp, svinf
    implicit none
    integer :: group, inft, prmcnt, iduv, i, j, k, m, pti, cnt
    real :: factor
    character(len=1024) :: opnfile

    select case(clcond)
    case('basic', 'range', 'merge')
    case default
       stop ' The clcond parameter is incorrect'
    end select
    !
    select case(clcond)
    case('basic', 'range')
       write(6, "(A)") " What is the energy distribution in solution?"
       read(5, *) engfile(1)
       if(slncor == 'yes') then
          write(6, "(A)") " What is the energy correlation in solution?"
          read(5, *) engfile(2)
       endif
       write(6, "(A)") " What is the energy density for insertion?"
       read(5, *) engfile(3)
       write(6, "(A)") " What is the energy correlation for insertion?"
       read(5, *) engfile(4)
       if(peread == 'yes') then
          write(6, "(A)") " Which file describes energy-coordinate parameters?"
          read(5, *) engfile(5)
       endif
       maxsln=1
       maxref=1
       numrun=1
       if(clcond == 'basic') prmmax=1
       if(clcond == 'range') prmmax=numprm
    case('merge')
       maxsln=numsln
       maxref=numref
       numrun=numdiv
       prmmax=numprm
    end select
    !
    if(clcond == 'merge') then
       opnfile = trim(solndirec) // '/' // trim(slndnspf) // '.01'
    else
       opnfile = engfile(1)
    endif
    open(unit = 71, file = opnfile, status = 'old')
    ermax = 0
    numslv = 0
    k = 0
    do iduv = 1, large
       read(71, *, end = 781) factor, pti, factor
       if(pti /= k) then
          numslv = numslv + 1
          k = pti
       endif
       ermax = ermax + 1
    end do
781 continue
    close(71)
    !
    allocate( nummol(numslv) )
    allocate( rduvmax(numslv), rduvcore(numslv) )
    allocate( rdcrd(ermax), rddst(ermax), rddns(ermax) )
    if(slncor == 'yes') allocate( rdslc(ermax, ermax) )
    allocate( rdcor(ermax,ermax) )
    allocate( rdspec(ermax) )
    allocate( chmpt(0:numslv, prmmax, numrun), aveuv(numslv) )
    if((uvread == 'yes') .and. (clcond == 'merge')) then
       allocate( uvene(numslv, maxsln), blkuv(0:numslv, numrun) )
    endif
    allocate( svgrp(prmmax), svinf(prmmax) )
    allocate( wgtsln(maxsln), wgtref(maxref) )
    !
    if(peread /= 'yes') then
       cnt = ermax / numslv
       rduvmax(:) = cnt
       rduvcore(:) = pecore

       ! check consistency
       if(clcond == 'merge') then
          opnfile = trim(solndirec) // '/' // trim(ecdinfofl)
       else
          opnfile = engfile(5)
       endif
       open(unit = 71, file = opnfile, status = 'old', err = 7899)
       ! Something is wrong ... 
       close(71)
       print *, 'Warning: EcdInfo file exists, but peread is not "yes"'
       print *, "Perhaps you forgot to set peread?"
7899   continue
    else
       ! opnfile is still engfile(1) or soln/engsln.01
       open(unit = 71, file = opnfile, status = 'old')
       k = 0
       do iduv = 1, ermax
          read(71, *) factor, pti, factor
          if(pti /= k) then
             rduvmax(pti) = 1
             k = pti
          else
             rduvmax(pti) = rduvmax(pti) + 1
          endif
       end do
       close(71)
       do pti = 1, numslv
          rduvcore(pti) = pecore
          if(clcond == 'merge') then
             opnfile = trim(solndirec) // '/' // trim(ecdinfofl)
          else
             opnfile = engfile(5)
          endif
          open(unit = 71, file = opnfile, status = 'old')
          read(71, *)        ! comment line
          read(71, *)        ! line for solute-self energy
          do i = 1, large
             read(71, *, end = 7829) k
             if(k == pti) then
                backspace(71)
                read(71, *) m, (factor,j = 1, 9), rduvcore(pti)
                exit
             endif
          end do
7829      continue
          close(71)
       end do
       k = 0
       do pti = 1, numslv
          k = k + rduvmax(pti)
       end do
       if(k /= ermax) stop ' The file format is incorrect'
    endif
    !
    if(ermax > maxmesh) stop ' The number of meshes is too large'
    !
    select case(clcond)
    case('basic')
       write(6, "(A)") " How many data are grouped into one?"
       read(5, *) group
       inft = 0
       if(infchk == 'yes') then
          write(6, "(A)") " How many large-energy meshes are merged ? (in %)"
          read(5, *) inft
       endif
       svgrp(1) = group
       svinf(1) = inft
       write(6, "(A)") " What is the temperature in Kelvin?"
       read(5, *) temp
    case('range', 'merge')
       do prmcnt = 1, prmmax                 ! inft in % (percentage)
          if(infchk == 'yes') then
             if(prmcnt == 1) then  ; group=1  ; inft=0   ; endif
             if(prmcnt == 2) then  ; group=1  ; inft=60  ; endif
             if(prmcnt == 3) then  ; group=1  ; inft=80  ; endif
             if(prmcnt == 4) then  ; group=1  ; inft=100 ; endif
             if(prmcnt == 5) then  ; group=2  ; inft=0   ; endif
             if(prmcnt == 6) then  ; group=3  ; inft=0   ; endif
             if(prmcnt == 7) then  ; group=4  ; inft=0   ; endif
             if(prmcnt == 8) then  ; group=5  ; inft=0   ; endif
             if(prmcnt == 9) then  ; group=5  ; inft=60  ; endif
             if(prmcnt == 10) then ; group=5  ; inft=80  ; endif
             if(prmcnt == 11) then ; group=5  ; inft=100 ; endif
             if(prmcnt == 12) then ; group=8  ; inft=0   ; endif
             if(prmcnt == 13) then ; group=10 ; inft=0   ; endif
             if(prmcnt == 14) then ; group=15 ; inft=0   ; endif
             if(prmcnt == 15) then ; group=20 ; inft=0   ; endif
          else
             inft = 0
             if(prmcnt <= 10) then
                group = prmcnt
             else
                group = 10 + (prmcnt - 10) * 2
             endif
          endif
          svgrp(prmcnt) = group
          svinf(prmcnt) = inft
       end do
       temp = inptemp
    end select
    kT = temp * 8.314510e-3 / 4.184              ! kcal/mol
    !
    if(uvread == 'yes') then
       select case(clcond)
       case('basic', 'range')
          write(6, "(A)") " What is average solute-solvent energy in solution?"
          read(5, *) aveuv(1:numslv)
       case('merge')
          opnfile = trim(solndirec) // '/' // trim(aveuvfile)
          open(unit = 82, file = opnfile, status = 'old')
          do k = 1, maxsln
             read(82, *) m, uvene(1:numslv, k)
          end do
          close(82)
       end select
    endif

    ! solution system
    wgtsln(1:maxsln) = 1.0
    if((clcond == 'merge') .and. (readwgtfl == 'yes')) then
       opnfile = trim(solndirec) // '/' // trim(wgtslnfl)
       open(unit = 81, file = opnfile, status = 'old')
       do i = 1, maxsln
          read(81, *) k, wgtsln(i)
       end do
       close(81)
    endif
    factor = sum( wgtsln(1:maxsln) )
    wgtsln(1:maxsln) = wgtsln(1:maxsln) / factor

    ! reference solvent system
    wgtref(1:maxref) = 1.0
    if((clcond == 'merge') .and. (readwgtfl == 'yes')) then
       opnfile = trim(refsdirec) // '/' // trim(wgtreffl)
       open(unit = 81, file = opnfile, status = 'old')
       do i = 1, maxref
          read(81, *) k, wgtref(i)
       end do
       close(81)
    endif
    factor = sum( wgtref(1:maxref) )
    wgtref(1:maxref) = wgtref(1:maxref) / factor

    if(slfslt == 'yes') then
       select case(clcond)
       case('basic', 'range')
          write(6, "(A)") " What is the solute self-energy?"
          read(5, *) slfeng
       case('merge')
          if(readwgtfl == 'not') then
             stop "readwgtfl needs to be yes when slfslt is yes"
          endif
          slfeng = 0.0
          opnfile = trim(refsdirec) // '/' // trim(wgtreffl)
          open(unit = 81, file = opnfile, status = 'old')
          do i = 1, maxref
             read(81,*) k, factor, factor
             slfeng = slfeng + wgtref(i) * factor
          end do
          close(81)
       end select
    endif

    return
  end subroutine

  subroutine datread(cntrun)
    use sysvars, only: refmerge, zero, tiny,&
         solndirec, refsdirec, slndnspf, slncorpf, refdnspf, refcorpf, numbers
    implicit none
    integer, intent(in) :: cntrun
    integer :: slnini, slnfin, refini, reffin, ecmin, ecmax
    integer :: iduv, iduvp, i, k, m, pti, cnt
    real :: factor, ampl
    logical :: num_different
    real, allocatable :: cormat_temp(:, :)
    character(len=1024) :: opnfile
    character(len=3) :: suffnum

    select case(clcond)
    case('basic', 'range')
       slnini = 1
       slnfin = 1
       refini = 1
       reffin = 1
    case('merge')
       if(maxsln >= numrun) then
          k = maxsln / numrun
          slnini = (cntrun - 1) * k + 1
          slnfin = cntrun * k
       else
          slnini = mod(cntrun - 1, maxsln) + 1
          slnfin = slnini
       endif
       if(refmerge == 'not') then
          if(maxref >= numrun) then
             m = maxref / numrun
             refini = (cntrun - 1) * m + 1
             reffin = cntrun * m
          else
             refini = mod(cntrun - 1, maxref) + 1
             reffin = refini
          endif
       else
          refini = 1
          reffin = maxref
       endif
    end select
    
    rddst(:) = 0.0
    if(slncor == 'yes') rdslc(:,:)=0.0
    if((cntrun == 1) .or. (refmerge == 'not')) then
       rddns(:) = 0.0
       rdcor(:,:) = 0.0
    endif
    
    ! FIXME: this part is kinda spaghetti and should be rewritten WITHOUT looping by cnt!
    do cnt = 1, 4
       if((cnt == 2) .and. (slncor /= 'yes')) cycle
       if((cnt >= 3) .and. (cntrun > 1) .and. (refmerge == 'yes')) cycle
       
       if(cnt <= 2) then                         ! solution
          ecmin = slnini
          ecmax = slnfin
          factor = sum( wgtsln(ecmin:ecmax) )
          wgtsln(ecmin:ecmax) = wgtsln(ecmin:ecmax) / factor
       endif
       if(cnt >= 3) then                         ! reference solvent
          ecmin = refini
          ecmax = reffin
          factor = sum( wgtref(ecmin:ecmax) )
          wgtref(ecmin:ecmax) = wgtref(ecmin:ecmax) / factor
       endif
       
       do i = ecmin, ecmax
          select case(clcond)
          case('basic', 'range')
             opnfile = engfile(cnt)
          case('merge')
             m = i / 10
             k = mod(i ,10)
             suffnum = '.' // numbers(m+1:m+1) // numbers(k+1:k+1)
             if(cnt == 1) opnfile = trim(solndirec) // '/' // trim(slndnspf)
             if(cnt == 2) opnfile = trim(solndirec) // '/' // trim(slncorpf)
             if(cnt == 3) opnfile = trim(refsdirec) // '/' // trim(refdnspf)
             if(cnt == 4) opnfile = trim(refsdirec) // '/' // trim(refcorpf)
             opnfile = trim(opnfile) // suffnum
          end select
          if((cnt == 1) .or. (cnt == 3)) then    ! 1-D distribution
             open(unit = 71, file = opnfile, status = 'old')
             k = 0
             m = 0
             do iduv = 1, ermax
                read(71, *) rdcrd(iduv), pti, factor
                if(pti /= k) then
                   if(factor > zero) then
                      write(6, *) ' Incorrect energy range with species ', pti
                      stop
                   endif
                   k = pti
                   m = m + 1
                endif
                if(cnt == 1) rddst(iduv) =rddst(iduv) + wgtsln(i) * factor
                if(cnt == 3) rddns(iduv) =rddns(iduv) + wgtref(i) * factor
                rdspec(iduv) = m
             enddo
             close(71)
          endif
          if((cnt == 2) .or. (cnt == 4)) then    ! 2-D correlation matrix
             allocate(cormat_temp(ermax, ermax))
             open(unit = 72, file = opnfile, status = 'old', form = "UNFORMATTED")
             read(72) cormat_temp
             close(72)
             if(cnt == 2) rdslc(:, :) = rdslc(:, :) + wgtsln(i) * cormat_temp(:, :)
             if(cnt == 4) rdcor(:, :) = rdcor(:, :) + wgtref(i) * cormat_temp(:, :)
             deallocate(cormat_temp)
          endif
       end do
    end do
!
    if(cntrun == 1) write(6, *)
    do pti=1,numslv
       factor = sum( rddst, mask = (rdspec == pti) )
       ampl   = sum( rddns, mask = (rdspec == pti) )
       num_different = .false.
       if(abs(factor-ampl) > tiny) num_different = .true.
       if(cntrun > 1) then
          if(nint(factor) /= nint(nummol(pti))) num_different = .true.
          if(nint(ampl) /= nint(nummol(pti))) num_different = .true.
       endif
       if(num_different) then
          write(6, '(A,i4)') '  Incorrect normalization at ', pti
          stop
       endif
       if(cntrun == 1) then
          nummol(pti) = real(nint(factor))
          write(6, '(A,i3,A,i12)') '  Number of the ', pti, '-th solvent  = ', nint(nummol(pti))
       endif
    end do
    if(cntrun == 1) write(6, *)
!
    if((uvread == 'yes') .and. (clcond == 'merge')) then
       do pti = 1, numslv
         aveuv(pti) = sum( wgtsln(slnini:slnfin) * uvene(pti, slnini:slnfin) )
       end do
       blkuv(1:numslv, cntrun) = aveuv(1:numslv)
       blkuv(0, cntrun) = sum( blkuv(1:numslv, cntrun) )
       if(slfslt == 'yes') blkuv(0, cntrun) = blkuv(0, cntrun) + slfeng
    endif

    return
  end subroutine datread
end module


module uvcorrect
  character(len=*), parameter :: volmfname = 'volumedata'
  real :: volm

  character(len=*), parameter :: ermodfname = 'parameters_er'
  integer :: ljformat, ljswitch, cmbrule
  real :: lwljcut, upljcut
  integer, parameter :: LJFMT_EPS_cal_SGM_nm = 0, LJFMT_EPS_Rminh = 1, &
                        LJFMT_EPS_J_SGM_A = 2, LJFMT_A_C = 3, &
                        LJFMT_C12_C6 = 4, LJFMT_TABLE = 5
  integer, parameter :: LJSWT_POT_CHM = 0, LJSWT_POT_GMX = 1, LJSWT_FORCE = 2
  integer, parameter :: LJCMB_ARITH = 0, LJCMB_GEOM = 1

  character(len=*), parameter :: sltfile = 'SltInfo'
  character(len=*), parameter :: prmfile = 'MolPrm'
  character(len=*), parameter :: ljtablefile = 'LJTable'
  integer                              :: ljtype_max
  integer, dimension(:,:), allocatable :: ljtype
  real, dimension(:,:),    allocatable :: ljlensq_mat, ljene_mat
  integer, dimension(:),   allocatable :: ljsite
  real, dimension(:),      allocatable :: ljcorr
  integer, dimension(:), allocatable :: ptsite

  ! Only ljformat, ljswitch, cmbrule, lwljcut, and upljcut used in this program
  namelist /ene_param/ iseed, &
       skpcnf, corrcal, selfcal, &
       slttype, sltpick, wgtslf, wgtins, wgtsys, &
       estype,boxshp, &
       insorigin, insposition, insorient, insstructure, inscnd, inscfg, &
       hostspec, lwreg, upreg, lwstr, upstr, &
       ljformat, ljswitch, &
       inptemp, temp, &
       engdiv, maxins, &
       intprm, elecut, lwljcut, upljcut, &
       cmbrule, cltype, screen, ewtoler, splodr, plmode, &
       ew1max, ew2max, ew3max, ms1max, ms2max, ms3max, &
       block_threshold, force_calculation

contains
  subroutine ljcorrect(cntrun)
    use sysvars, only: uvread, clcond, numslv, aveuv, blkuv
    implicit none
    integer, intent(in) :: cntrun
    logical, save :: first_time = .true.
    integer :: pti
    if(first_time) then
       call set_keyparam
       write(6,550) lwljcut, upljcut
550    format('  Be sure that the solvent distribution is homogeneous (radial distribution function is essentially unity) when the solvent molecule is separated beyond distance of',f7.1,' or ',f7.1,' Angstrom in any direction from any atom within the solute molecule')
       write(6,*)
       call get_ljtable
       allocate( ljcorr(numslv) )
       do pti = 1, numslv
          call calc_ljlrc(pti, ljcorr(pti))
       end do
       write(6,551) ljcorr(1:numslv)
       write(6,*)
551    format('  LJ long-range correction    =   ', 9999f12.4)
       first_time = .false.
    endif
    aveuv(1:numslv) = aveuv(1:numslv) + ljcorr(1:numslv)
    if((uvread.eq.'yes').and.(clcond.eq.'merge')) then
       blkuv(1:numslv, cntrun) = blkuv(1:numslv, cntrun) + ljcorr(1:numslv)
       blkuv(0, cntrun) = blkuv(0, cntrun) + sum(ljcorr(1:numslv))
    endif
    return
  end subroutine ljcorrect

  subroutine set_keyparam
    use sysvars, only: refsdirec
    implicit none
    character(len=80) :: keyfile
    integer, parameter :: iounit = 555
    real, parameter :: volm_min = 3.0e3
    integer :: stat
    logical :: found
    ljformat = LJFMT_EPS_Rminh                    ! default setting
    ljswitch = LJSWT_POT_CHM                      ! default setting
    cmbrule = LJCMB_ARITH                         ! default setting
    upljcut = 12.0                                ! default setting
    lwljcut = upljcut - 2.0                       ! default setting
    keyfile = trim(refsdirec)//'/'//ermodfname
    open(unit = iounit, file = keyfile, action='read', status='old', iostat=stat)
    if(stat == 0) then
       read(iounit, nml = ene_param)
    else
       stop "The parameters_er file is not found in the refs directory"
    endif
    close(iounit)

    inquire(file = volmfname, exist = found)
    if(found) then
       open(unit = iounit, file = volmfname, status='old')
       read(iounit, *) volm
       close(iounit)
    else
       write(6, '(A)') "  What is the average volume of reference solvent? (in Angstrom^3)"
       read(5, *) volm
    endif
    if(volm < volm_min) then
       write(6, '(A)') "  Warning: your input volume seems too small"
       write(6, '(A, F8.1)') "           This warning appears when your input is less than ",volm_min
       write(6, '(A)') "  Re-type the volume in Angstrom^3 (NOT in nm^3)"
       read(5, *) volm
    endif
  end subroutine set_keyparam

  subroutine get_ljtable   ! taken essentially from setconf.F90
    use sysvars, only: numslv, refsdirec
    implicit none
    real, parameter :: sgmcnv = 1.7817974362806784e0 ! from Rmin/2 to sigma, 2.0**(5.0/6.0)
    real, parameter :: lencnv = 1.0e1                ! from nm to Angstrom
    real, parameter :: engcnv = 1.0e0/4.184e0        ! from kJ/mol to kcal/mol
    integer :: pti, sid, stmax, maxsite, i, m
    real :: factor, xst(3)
    integer, allocatable :: ljtype_temp(:)
    real, dimension(:), allocatable :: ljlen_temp, ljene_temp
    real, dimension(:), allocatable :: ljlen_temp_table, ljene_temp_table
    integer :: ljtype_found
    logical :: lj_is_new
    character(len=5) :: atmtype
    character(len=80) :: molfile
    integer, parameter :: molio = 71                 ! IO for molfile
    integer, parameter :: ljtableio = 70             ! IO for LJ table
    character(len=9) :: numbers = '123456789'

    allocate( ptsite(0:numslv) )
    do pti = 0, numslv
       if(pti == 0) then
          molfile = sltfile                      ! solute
       else
          molfile = prmfile//numbers(pti:pti)    ! solvent
       endif
       molfile = trim(refsdirec)//'/'//molfile
       open(unit = molio, file = molfile, status='old')
       stmax = 0
       do
          read(molio, *, end = 99) m
          stmax = stmax + 1
       end do
99     close(molio)
       ptsite(pti) = stmax
    end do

    ! large enough LJ table size
    allocate( ljlen_temp_table(1:sum(ptsite(:))), &
              ljene_temp_table(1:sum(ptsite(:))) )
    ! temporary set of LJ
    maxsite = maxval( ptsite(0:numslv) )
    allocate( ljtype_temp(maxsite), ljlen_temp(maxsite), ljene_temp(maxsite) )

    allocate( ljtype(maxsite, 0:numslv) )
    ljtype(:, :) = 0

    do pti = 0, numslv
       if(pti == 0) then
          molfile = sltfile                      ! solute
       else
          molfile = prmfile//numbers(pti:pti)    ! solvent
       endif
       molfile = trim(refsdirec)//'/'//molfile
       stmax = ptsite(pti)
       open(unit = molio, file = molfile, status = 'old')
       do sid = 1, stmax
          read(molio,*) m, atmtype, xst(1:3)
          if(ljformat == LJFMT_EPS_Rminh) xst(3) = sgmcnv * xst(3)
          if((ljformat == LJFMT_A_C) .or. (ljformat == LJFMT_C12_C6)) then
             if(xst(3) /= 0.0) then
                factor = (xst(2) / xst(3)) ** (1.0 / 6.0)
                xst(2) = xst(3) / (4.0 * (factor ** 6))
                xst(3) = factor
             else
                xst(2) = 0.0
             endif
          endif
          if((ljformat == LJFMT_EPS_J_SGM_A) .or. (ljformat == LJFMT_C12_C6)) then
             xst(2) = engcnv * xst(2)
             xst(3) = lencnv * xst(3)
          endif
          ljene_temp(sid) = xst(2)
          ljlen_temp(sid) = xst(3)
       end do
       close(molio)

      if(ljformat == LJFMT_TABLE) then
          ljtype_temp(1:stmax) = ljene_temp(1:stmax)
       else
          do sid = 1, stmax
             lj_is_new = .true.
             do i = 1, ljtype_max
                ! linear search LJ table
                if((ljlen_temp_table(i) == ljlen_temp(sid)) .and. &
                   (ljene_temp_table(i) == ljene_temp(sid))) then
                   ljtype_found = i
                   lj_is_new = .false.
                   exit
                endif
             end do
             if(lj_is_new) then
                ! new LJ type
                ljtype_max = ljtype_max + 1
                ljlen_temp_table(ljtype_max) = ljlen_temp(sid)
                ljene_temp_table(ljtype_max) = ljene_temp(sid)
                ljtype_found = ljtype_max
             endif
             ljtype_temp(sid) = ljtype_found
          end do
       endif

       ljtype(1:stmax, pti) = ljtype_temp(1:stmax)
    end do
    deallocate( ljlen_temp, ljene_temp, ljtype_temp )

    ! Fill LJ table
    if(ljformat == LJFMT_TABLE) then
       ! From table (directly)
       open(unit = ljtableio, file = trim(refsdirec)//'/'//ljtablefile, status = 'old', action = 'read')
       read(ljtableio, *) ljtype_max
       allocate( ljlensq_mat(ljtype_max, ljtype_max), &
                 ljene_mat(ljtype_max, ljtype_max) )
       do i = 1, ljtype_max
          read (ljtableio, *) ljlensq_mat(i, 1:ljtype_max)
          ljlensq_mat(i, 1:ljtype_max) = ljlensq_mat(i, 1:ljtype_max) ** 2
       end do
       do i = 1, ljtype_max
          read (ljtableio, *) ljene_mat(i, 1:ljtype_max)
       end do
       close(ljtableio)
    else
       ! From LJ data
       allocate( ljlensq_mat(ljtype_max, ljtype_max), &
                 ljene_mat(ljtype_max, ljtype_max) )
       do i = 1, ljtype_max
          select case(cmbrule)
          case(LJCMB_ARITH)    ! arithmetic mean
             ljlensq_mat(1:ljtype_max, i) = (( ljlen_temp_table(1:ljtype_max) &
                                             + ljlen_temp_table(i) ) / 2.0) ** 2
          case(LJCMB_GEOM)     ! geometric mean
             ljlensq_mat(1:ljtype_max, i) = ljlen_temp_table(1:ljtype_max) &
                                          * ljlen_temp_table(i)
          case default
             stop "Incorrect cmbrule"
          end select
          ljene_mat(1:ljtype_max, i) = sqrt( ljene_temp_table(1:ljtype_max)  &
                                           * ljene_temp_table(i) )
       end do
    endif
    deallocate(ljlen_temp_table, ljene_temp_table)

    return
  end subroutine get_ljtable

  subroutine calc_ljlrc(pti, correction)
    use sysvars, only: nummol
    implicit none
    integer, intent(in) :: pti
    real, intent(out) :: correction
    real :: dens, ljeps, ljsgm2
    integer :: ui, vi
    dens = nummol(pti) / volm
    correction = 0.0
    do ui = 1, ptsite(0)           ! sum over solute sites
       do vi = 1, ptsite(pti)      ! sum over solvent sites
          ljeps = ljene_mat( ljtype(ui, 0), ljtype(vi, pti) )
          ljsgm2 = ljlensq_mat( ljtype(ui, 0), ljtype(vi, pti) )
          correction = correction + dens * enelj(ljeps, ljsgm2)
       end do
    end do
    return
  end subroutine calc_ljlrc

  real function enelj(ljeps, ljsgm2)
    implicit none
    real, save :: rbin = 1.0e-3
    integer, save :: numbin
    real, intent(in) :: ljeps, ljsgm2
    real, parameter :: PI = 3.1415926535897932
    real :: ljint, edev, dist, r
    real :: upljcut2, lwljcut2, upljcut3, lwljcut3, upljcut6, lwljcut6
    real :: invr2, invr3, invr6, ljsgm3, ljsgm6, vdwa, vdwb, swth, swfac
    logical, save :: do_swth, first_time = .true.
    integer :: i

    if(first_time) then
       if(lwljcut > upljcut) then
          stop "Incorrect setting of lwljcut and upljcut (lwljcut > upljcut)"
       else
          numbin = nint((upljcut - lwljcut) / rbin)
          if(numbin >= 1) then
             do_swth = .true.
             rbin = (upljcut - lwljcut) / real(numbin)
          else
             do_swth = .false.
          endif
       endif
       first_time = .false.
    endif

    lwljcut2 = lwljcut ** 2
    lwljcut3 = lwljcut2 * lwljcut
    lwljcut6 = lwljcut3 * lwljcut3

    upljcut2 = upljcut ** 2
    upljcut3 = upljcut2 * upljcut
    upljcut6 = upljcut3 * upljcut3

    ljsgm6 = ljsgm2 * ljsgm2 * ljsgm2
    ljsgm3 = sqrt(ljsgm6)

    ljint = 0.0

    ! r < lwljcut
    select case(ljswitch)
    case(LJSWT_POT_CHM, LJSWT_POT_GMX)    ! potential switch
       ! do nothing
    case(LJSWT_FORCE)                     ! force switch
       vdwa = ljsgm6 * ljsgm6 / (lwljcut6 * upljcut6)
       vdwb = ljsgm6 / (lwljcut3 * upljcut3)
       edev = 4.0 * ljeps * (vdwa - vdwb)
       ljint = ljint + (4.0 * PI /3.0) * lwljcut3 * edev
    case default
       stop "Unknown ljswitch"
    end select

    ! lwljcut < r < upljcut
    if(do_swth) then
       do i = 1, numbin
          r = lwljcut + (real(i) - 0.5) * rbin
          dist = r * r
          invr2 = ljsgm2 / dist
          invr6 = invr2 * invr2 * invr2
          select case(ljswitch)
          case(LJSWT_POT_CHM)             ! potential switch (CHRAMM form)
             swth = (2.0 * dist + upljcut2 - 3.0 * lwljcut2)        &
                  * ((dist - upljcut2) ** 2) / ((upljcut2 - lwljcut2) ** 3)
             edev = 4.0 * ljeps * invr6 * (invr6 - 1.0) * (1.0 - swth)
          case(LJSWT_POT_GMX)             ! potential switch (GROMACS form)
             swfac = (r - lwljcut) / (upljcut - lwljcut)
             swth = 1.0 - 10.0 * (swfac ** 3)                      &
                        + 15.0 * (swfac ** 4) - 6.0 * (swfac ** 5) 
             edev = 4.0 * ljeps * invr6 * (invr6 - 1.0) * (1.0 - swth)
          case(LJSWT_FORCE)               ! force switch
             invr3 = sqrt(invr6)
             vdwa = upljcut6 / (upljcut6 - lwljcut6)          &
                  * ( (invr6 - ljsgm6 / upljcut6) ** 2 )
             vdwb = upljcut3 / (upljcut3 - lwljcut3)          &
                  * ( (invr3 - ljsgm3 / upljcut3) ** 2 )
             edev = 4.0 * ljeps * ( invr6 * (invr6 - 1.0) - (vdwa - vdwb) )
          end select
          ljint = ljint + 4.0 * PI * r * r * rbin * edev
       end do
    endif

    ! r > upljcut
    invr3 = ljsgm3 / upljcut3
    ljint = ljint + ljeps * ljsgm3 &
                          * (16.0 * PI / 3.0) * ( (invr3 ** 3 / 3.0) - invr3 )

    enelj = ljint
    return
  end function enelj
end module


module sfecalc
  use sysvars, only: zerosft, wgtfnform, slncor, &
                     numslv, ermax, nummol, kT, itrmax, zero, error, &
                     rduvmax, rduvcore, &
                     rdcrd, rddst, rddns, rdslc, rdcor, rdspec
  implicit none
  integer, dimension(:), allocatable :: idrduv, uvmax
  real, dimension(:),    allocatable :: uvcrd, edist, edens
  real, dimension(:,:),  allocatable :: edscr, ecorr
  integer, dimension(:), allocatable :: uvspec
  real, dimension(:),    allocatable :: slncv, inscv, sdrcv
  real, dimension(:),    allocatable :: zrsln, zrref, zrsdr
  integer gemax
contains
  ! TODO: write the cases for (kind(real) /= 8).
  subroutine syevr_wrap(n, mat, eigval, info)
    implicit none
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
         z, n, isuppz, worksize, lwork, iwork, liwork, info)
    if (info /= 0) then
       deallocate(isuppz)
       deallocate(z)
       deallocate(iwork)
       return
    endif

    lwork = worksize
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

  subroutine chmpot(prmcnt, cntrun)
    !
    use sysvars, only: uvread, slfslt, normalize, showdst, wrtzrsft, &
                       ljlrc, &
                       slfeng, chmpt, aveuv, svgrp, svinf, &
                       pickgr, cumuint, cumuintfl
    use uvcorrect, only: ljcorrect
    !
    implicit none
    integer, intent(in) :: prmcnt, cntrun
    integer :: group, inft
    integer :: iduv, iduvp, pti, cnt, j, k, m, ge_perslv
    real :: factor, ampl, slvfe, uvpot, lcent, lcsln, lcref
    integer, dimension(:), allocatable :: gpnum
    real, dimension(:,:), allocatable :: cumsfe
    integer, parameter :: cumu_io = 51
    !
    group = svgrp(prmcnt)
    inft = svinf(prmcnt)
    !
    allocate( idrduv(ermax), uvmax(numslv) )
    !
    gemax = 0
    do pti = 1, numslv
       k = rduvcore(pti) * inft / 100
       m = (rduvmax(pti) - k) / group
       uvmax(pti) = m
       gemax = gemax + m
       if(pti == 1) then
          cnt = 1
          j = 0
       else
          cnt = 1 + sum( rduvmax(1: pti-1) )
          j = sum( uvmax(1: pti-1) )
       endif
       do iduv = cnt, cnt + rduvmax(pti) - 1
          k = (iduv - cnt) / group + 1
          if(k >= m) k = m
          idrduv(iduv) = k + j
       end do
       if((pti == numslv) .and. (cnt + rduvmax(pti) - 1 /= ermax)) then
          print *, "Error: The total no. of meshes does not match with input"
          print *, "(Sum should be", ermax, " but was", cnt + rduvmax(pti) - 1, ")"
          stop
       endif
    end do
    !
    allocate( uvcrd(gemax), edist(gemax), edens(gemax), uvspec(gemax) )
    uvcrd(:) = 0.0
    edist(:) = 0.0
    edens(:) = 0.0
    !
    if(slncor == 'yes') then
       allocate( edscr(gemax, gemax) )
       edscr(:,:) = 0.0
    endif
    allocate( ecorr(gemax, gemax) )
    ecorr(:,:) = 0.0
    !
    allocate( slncv(gemax), inscv(gemax), zrsln(numslv), zrref(numslv) )
    slncv(:) = 0.0
    inscv(:) = 0.0
    zrsln(:) = 0.0
    zrref(:) = 0.0
    if(slncor == 'yes') then
       allocate( sdrcv(gemax), zrsdr(numslv) )
       sdrcv(:) = 0.0
       zrsdr(:) = 0.0
    endif
    !
    allocate( gpnum(gemax) )
    gpnum(:) = 0
    do iduv = 1, ermax
       k = idrduv(iduv)
       uvcrd(k) = uvcrd(k) + rdcrd(iduv)
       gpnum(k) = gpnum(k) + 1
       uvspec(k) = rdspec(iduv)
    end do
    do k = 1, gemax
       if(gpnum(k) > 0) uvcrd(k) = uvcrd(k) / real(gpnum(k))
    end do
    cnt = rduvmax(1)
    k = uvmax(1)
    do pti = 1, numslv
       if(pti > 1) then
          cnt = cnt + rduvmax(pti)
          k = k + uvmax(pti)
       endif
       uvcrd(k) = rdcrd(cnt)
    end do
    deallocate( gpnum )
    !
    do cnt = 1, 2
       do iduv = 1, ermax
          k = idrduv(iduv)
          if(cnt == 1) edist(k) = edist(k) + rddst(iduv)
          if(cnt == 2) edens(k) = edens(k) + rddns(iduv)
       end do
       if((cnt == 1) .and. (slncor /= 'yes')) goto 1115
       do iduv = 1, ermax
          do iduvp = 1, ermax
             k = idrduv(iduv)
             m = idrduv(iduvp)
             if(cnt == 1) edscr(m,k) = edscr(m,k) + rdslc(iduvp, iduv)
             if(cnt == 2) ecorr(m,k) = ecorr(m,k) + rdcor(iduvp, iduv)
          end do
       end do
1115   continue
    end do
    !
    if(normalize == 'yes') call distnorm
    if(showdst == 'yes')   call distshow
    !
    call getslncv
    call getinscv
    !
    if(prmcnt == 1) then
       if(uvread /= 'yes') then
          do pti=1,numslv
             aveuv(pti) = sum( uvcrd * edist, mask = (uvspec == pti) )
          end do
       endif
       if(ljlrc == 'yes') call ljcorrect(cntrun)   ! LJ long-range correction
    endif
    !
    if((cumuint == 'yes') .and. (group == pickgr) .and. (inft == 0)) then
       ge_perslv = gemax / numslv
       if(any(uvmax(1:numslv) /= ge_perslv)) stop ' Incorrect file format for storing the running integral'
       allocate( cumsfe(numslv, ge_perslv) )
    endif
    !
    do pti = 1, numslv
       uvpot = 0.0
       slvfe = 0.0
       do iduv = 1, gemax
          if(uvspec(iduv) == pti) then
             if((edist(iduv) <= zero) .and. (edens(iduv) <= zero)) goto 5009
             uvpot = uvpot + uvcrd(iduv) * edist(iduv)
             slvfe = slvfe - kT * (edist(iduv) - edens(iduv))

             ! kT*log(edist/edens)
             lcent = - (slncv(iduv) + zrsln(pti) + uvcrd(iduv))
             if((slncor == 'yes') .and. &
                     (edist(iduv) > zero) .and. (edens(iduv) <= zero)) then
                ampl = lcent * edens(iduv) / edist(iduv)
                lcent = ampl - (zrsln(pti) + uvcrd(iduv)) &
                             * (1.0 - edens(iduv) / edist(iduv))
             endif
             slvfe = slvfe + lcent * edist(iduv)

             lcsln = pyhnc(slncv(iduv), 1)
             lcref = pyhnc(inscv(iduv), 2)
             if((slncor == 'yes') .and. &
                     (edist(iduv) > zero) .and. (edens(iduv) <= zero)) then
                lcsln = pyhnc(sdrcv(iduv) + zrsdr(pti), 3)
             endif
             ampl = sfewgt(edist(iduv), edens(iduv))
             factor = ampl * lcsln + (1.0 - ampl) * lcref
             slvfe = slvfe + kT * factor * (edist(iduv) - edens(iduv))
5009         continue
             if((cumuint == 'yes') .and. (group == pickgr) .and. (inft == 0)) then
               m = mod(iduv - 1, ge_perslv) + 1
               cumsfe(pti, m) = uvpot + slvfe
             endif
          endif
       end do
       chmpt(pti, prmcnt, cntrun) = slvfe + aveuv(pti)
    end do
    !
    if((cumuint == 'yes') .and. (group == pickgr) .and. (inft == 0)) then
      open(unit = cumu_io, file = cumuintfl, status = 'replace')
      do iduv = 1, ge_perslv
        factor = sum( cumsfe(1:numslv, iduv) )
        if(numslv == 1) then
           write(cumu_io, 511) iduv, uvcrd(iduv), factor
        else
           write(cumu_io, 511) iduv, uvcrd(iduv), factor, cumsfe(1:numslv,iduv)
        endif
      enddo
511   format(i6, g15.5, 9999f12.5)
      endfile(cumu_io)
      close(cumu_io)
      deallocate( cumsfe )
    endif
    !
    chmpt(0, prmcnt, cntrun) = sum( chmpt(1:numslv, prmcnt, cntrun) )
    if(slfslt == 'yes') chmpt(0, prmcnt, cntrun) = chmpt(0, prmcnt, cntrun) + slfeng
    !
    if(wrtzrsft == 'yes') then
       write(6, 381) zrsln(1:numslv)
       write(6, 382) zrref(1:numslv)
       if(slncor == 'yes') write(6, 383) zrsdr(1:numslv)
381    format('  Zero shift for solution             = ', 9999f12.4)
382    format('  Zero shift for reference solvent    = ', 9999f12.4)
383    format('  Zero shift for solution correlation = ', 9999f12.4)
    endif
    !
    deallocate( slncv, inscv, zrsln, zrref )
    if(slncor == 'yes') deallocate( sdrcv, zrsdr, edscr )
    deallocate( uvcrd, edist, edens, ecorr, uvspec, idrduv, uvmax )
    !
    return
  end subroutine chmpot

  subroutine getslncv
    use sysvars, only: extsln, extthres_soln, extthres_refs
    implicit none
    integer :: iduv, iduvp, pti, j, k, m
    real :: factor, ampl, lcsln, lcref, min_rddst, min_rddns
    real, dimension(:), allocatable :: work
    integer, parameter :: ofdmp = 10 ! factor to suppress the integer overflow
    logical, dimension(:), allocatable :: ext_target

    min_rddst = minval( rddst, mask = (rddst > zero) )
    min_rddns = minval( rddns, mask = (rddns > zero) )
    allocate( ext_target(gemax) )
    ext_target(:) = .true.
    !
    do iduv = 1, gemax
       factor = edist(iduv) / min_rddst
       m = ofdmp * extthres_soln
       if(factor > real(m)) then
          j = m
       else
          j = nint(factor)
       endif
       factor = edens(iduv) / min_rddns
       m = ofdmp * extthres_refs
       if(factor > real(m)) then
          k = m
       else
          k = nint(factor)
       endif
       if((j < extthres_soln) .or. (k < extthres_refs)) ext_target(iduv) = .false.
    enddo
    !
    do iduv = 1, gemax
       if(ext_target(iduv)) then
          slncv(iduv) = - kT * log(edist(iduv) / edens(iduv)) - uvcrd(iduv)
       endif
    end do
    !
    do iduv = 1, gemax
       if(.not. ext_target(iduv)) then
          if(edist(iduv) <= zero) then
             slncv(iduv) = 0.0
             cycle
          endif
          pti = uvspec(iduv)
          m = 1
          k = gemax
          do iduvp = 1, iduv - 1
             if((uvspec(iduvp) == pti) .and. (ext_target(iduvp)) .and. &
                (m < iduvp)) m = iduvp
          end do
          do iduvp = gemax, iduv + 1, -1
             if((uvspec(iduvp) == pti) .and. (ext_target(iduvp)) .and. &
                (k > iduvp)) k = iduvp
          end do
          !
          if(extsln == 'sim') then
             if(abs(m - iduv) <  abs(k - iduv)) factor = slncv(m)
             if(abs(m - iduv) >  abs(k - iduv)) factor = slncv(k)
             if(abs(m - iduv) == abs(k - iduv)) then
                factor = (slncv(m) + slncv(k)) / 2.0
             endif
          else
             j = k
             if(abs(m - iduv) < abs(k - iduv)) j = m
             allocate( work(gemax) )
             work(:) = 0.0
             do iduvp = 1, gemax
                if((uvspec(iduvp) == pti) .and. (ext_target(iduvp))) then
                   if(iduvp == iduv) stop ' A bug in program or data'
                   factor = uvcrd(iduvp) - uvcrd(j)
                   if(iduvp < iduv) then
                      factor = - factor - 2.0 * (uvcrd(j) - uvcrd(iduv))
                   endif
                   ampl = wgtdst(iduvp, 1, 'extsl', wgtfnform)
                   work(iduvp) = exp(- factor / kT) * ampl
                endif
             end do
             factor = sum( work, mask = (work > zero) )
             do iduvp = 1, gemax
                work(iduvp) = work(iduvp) / factor
             end do
             factor = 0.0
             ampl = 0.0
             lcsln = 0.0
             lcref = 0.0
             do iduvp = 1, gemax
                if(work(iduvp) > zero) then
                   factor = factor + work(iduvp) * uvcrd(iduvp)
                   ampl = ampl + work(iduvp)* uvcrd(iduvp) * uvcrd(iduvp)
                   lcsln = lcsln + work(iduvp) * slncv(iduvp)
                   lcref = lcref + work(iduvp) * uvcrd(iduvp) * slncv(iduvp)
                endif
             end do
             work(1) = (ampl * lcsln - factor * lcref) / (ampl - factor ** 2)
             work(2) = (lcref - factor * lcsln) / (ampl - factor ** 2)
             factor = work(1) + work(2) * uvcrd(iduv)
             deallocate( work )
          endif
          slncv(iduv) = factor
       endif
    end do
    !
    deallocate( ext_target )
    !
    do pti = 1, numslv
       select case(zerosft)
       case('orig')
          factor = 0.0
       case('mxco')
          lcref = - cvfcen(pti, 2, 'uvcrd', 'smpl', 'yes')
          lcsln = cvfcen(pti, 1, 'slncv', wgtfnform, 'not')
          ampl = wgtmxco(pti)
          factor = ampl * lcsln + (1.0 - ampl) * lcref
       case('zero')
          factor = cvfcen(pti, 1, 'slncv', wgtfnform, 'yes')
       case('cntr')
          factor = cvfcen(pti, 1, 'slncv', wgtfnform, 'not')
       case default
          stop ' zerosft not properly set '
       end select
       do iduv = 1, gemax
          if(uvspec(iduv) == pti) slncv(iduv) = slncv(iduv) - factor
       end do
       zrsln(pti) = factor
    end do
    !
    return
  end subroutine getslncv
  !
  !
  subroutine getinscv
    !
    implicit none
    integer :: iduv, iduvp, pti, cnt, wrksz, k
    real :: factor, ampl, lcsln, lcref
    real, dimension(:),   allocatable :: work, egnvl, zerouv
    real, dimension(:,:), allocatable :: edmcr
    !
    do cnt = 1, 2     ! cnt = 1: solution   cnt = 2: reference solvent
       if((cnt == 1) .and. (slncor /= 'yes')) cycle
       wrksz = gemax ** 2
       allocate( work(wrksz), egnvl(gemax), edmcr(gemax, gemax) )
       if(cnt == 1) edmcr(:,:) = edscr(:,:)
       if(cnt == 2) edmcr(:,:) = ecorr(:,:)
       do iduv = 1, gemax
          do iduvp = 1, gemax
             if(cnt == 1) then
                factor = edist(iduvp)
                ampl = edist(iduv)
             endif
             if(cnt == 2) then
                factor = edens(iduvp)
                ampl = edens(iduv)
             endif
             lcref = edmcr(iduvp, iduv) - factor * ampl
             if((factor <= zero) .or. (ampl <= zero)) then
                if(iduv == iduvp) then
                   lcref=1.0
                else
                   lcref=0.0
                endif
             endif
             edmcr(iduvp, iduv) = lcref
          end do
       end do
       !
       call syevr_wrap(gemax, edmcr, egnvl, k)
       ! call DSYEV('V', 'U', gemax, edmcr, gemax, egnvl, work, wrksz, k)
       pti = numslv + 1
       do iduv = pti, gemax
          factor = 0.0
          do iduvp = 1, gemax
             if(cnt == 1) ampl = edist(iduvp)
             if(cnt == 2) ampl = edens(iduvp)
             if(ampl > zero) factor = factor + (edist(iduvp) - edens(iduvp)) &
                                             * edmcr(iduvp, iduv)
          end do
          work(iduv) = factor / egnvl(iduv)
       end do
       do iduv = 1, gemax
          factor = 0.0
          do iduvp = pti, gemax
             factor = factor + edmcr(iduv, iduvp) * work(iduvp)
          end do
          if(cnt == 1) sdrcv(iduv) = - kT * factor
          if(cnt == 2) inscv(iduv) = - kT * factor
       end do
       deallocate( work, egnvl, edmcr )
       !
       allocate( zerouv(numslv) )
       do pti = 1, numslv
          select case(zerosft)
          case('orig')
             k = zeroec(pti, cnt)
             if(cnt == 1) factor = sdrcv(k)
             if(cnt == 2) factor = inscv(k)
          case('mxco')
             factor = cvfcen(pti, cnt, 'inscv', 'smpl', 'yes')
          case default
             factor = 0.0
          end select
          zerouv(pti) = factor
       end do
       !
       do iduv = 1, gemax
          if(cnt == 1) ampl = edist(iduv)
          if(cnt == 2) ampl = edens(iduv)
          if(ampl > zero) then
             factor = - kT * (edist(iduv) - edens(iduv)) / ampl
             if(cnt == 1) lcref = factor - sdrcv(iduv)
             if(cnt == 2) lcref = factor - inscv(iduv)
          endif
          if(ampl <= zero) lcref=0.0
          if(cnt == 1) sdrcv(iduv) = lcref
          if(cnt == 2) inscv(iduv) = lcref
       end do
       !
       do pti = 1, numslv
          select case(zerosft)
          case('orig')
             factor = - zerouv(pti)
          case('mxco')
             lcref = - zerouv(pti)
             lcsln = cvfcen(pti, cnt, 'inscv', wgtfnform, 'not')
             ampl = wgtmxco(pti)
             factor = ampl * lcsln + (1.0 - ampl) * lcref
          case('zero')
             factor = cvfcen(pti, cnt, 'inscv', wgtfnform, 'yes')
          case('cntr')
             factor = cvfcen(pti, cnt, 'inscv', wgtfnform, 'not')
          case default
             stop ' zerosft not properly set '
          end select
          do iduv = 1, gemax
             if(uvspec(iduv) == pti) then
                if(cnt == 1) sdrcv(iduv) = sdrcv(iduv) - factor
                if(cnt == 2) inscv(iduv) = inscv(iduv) - factor
             endif
          end do
          if(cnt == 1) zrsdr(pti) = factor
          if(cnt == 2) zrref(pti) = factor
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
    implicit none
    integer, intent(in) :: pti
    real :: numpt
    numpt = nummol(pti)
    wgtmxco = 1.0 / numpt
    return
  end function wgtmxco
  !
  !
  real function cvfcen(pti, cnt, systype, wgttype, engtype)
    implicit none
    ! pti : identifier of solvent species
    ! cnt = 1: solution   cnt = 2: reference solvent
    integer, intent(in) :: pti, cnt
    character(len=5), intent(in) :: systype
    character(len=4), intent(in) :: wgttype
    character(len=3), intent(in) :: engtype
    integer :: iduv
    logical :: errtag
    real :: factor, cvfnc
    real, dimension(:), allocatable :: weight
    allocate( weight(gemax) )
    call getwght(weight, pti, cnt, systype, wgttype, engtype)
    factor = 0.0
    errtag = .false.
    do iduv = 1, gemax
       if(uvspec(iduv) == pti) then
          select case(systype)
          case('slncv')
             if(cnt == 1) cvfnc = slncv(iduv)
             if(cnt == 2) errtag = .true.
          case('inscv')
             if(cnt == 1) cvfnc = sdrcv(iduv)
             if(cnt == 2) cvfnc = inscv(iduv)
          case('uvcrd')
             cvfnc = uvcrd(iduv)
          case default
             errtag = .true.
          end select
          if(errtag) stop ' Bug in the program'
          factor = factor + cvfnc * weight(iduv)
       endif
    end do
    cvfcen = factor
    deallocate( weight )
    return
  end function cvfcen
  !
  !
  subroutine getwght(weight, pti, cnt, systype, wgttype, engtype)
    implicit none
    real, intent(out) :: weight(gemax)
    ! pti : identifier of solvent species
    ! cnt = 1: solution   cnt = 2: reference solvent
    integer, intent(in) :: pti, cnt
    character(len=5), intent(in) :: systype
    character(len=4), intent(in) :: wgttype
    character(len=3), intent(in) :: engtype
    integer :: iduv
    real :: minuv, ampl
    weight(:) = 0.0
    do iduv = 1, gemax
       if(uvspec(iduv) == pti) then
          weight(iduv) = wgtdst(iduv, cnt, systype, wgttype)
       endif
    end do
    if(engtype == 'yes') then
       minuv = minval( abs(uvcrd(:)), &
                       mask = ((uvspec == pti) .and. (weight > zero)) )
       do iduv = 1, gemax
          if(uvspec(iduv) == pti) then
             !  ampl = abs(uvcrd(iduv)) - minuv
             ampl = nummol(pti) * (abs(uvcrd(iduv)) - minuv)
             weight(iduv) = exp(- ampl / kT) * weight(iduv)
          endif
       end do
    endif
    ampl = sum( weight, mask = (uvspec == pti) )
    if(ampl > zero) then
       weight(:) = weight(:) / ampl
    else
       write(6, *) ' Zero weight at ', pti
       stop
    endif
  end subroutine getwght
  !
  !
  integer function zeroec(pti, cnt)
    implicit none
    ! pti : identifier of solvent species
    ! cnt = 1: solution   cnt = 2: reference solvent
    integer, intent(in) :: pti, cnt
    integer :: iduv, k
    real :: factor, ampl, lcsln, lcref
    do iduv = 1, gemax - 1
       if(uvspec(iduv) == pti) then
          if((uvcrd(iduv) <= 0.0) .and. (uvcrd(iduv + 1) >= 0.0)) then
             factor = abs(uvcrd(iduv))
             ampl = uvcrd(iduv + 1)
             if(ampl > factor + zero) k = iduv
             if(factor > ampl + zero) k = iduv + 1
             if(abs(factor - ampl) <= zero) then
                if(cnt == 1) then
                   lcsln = edist(iduv)
                   lcref = edist(iduv + 1)
                endif
                if(cnt == 2) then
                   lcsln = edens(iduv)
                   lcref = edens(iduv + 1)
                endif
                if(lcsln >= lcref) then
                   k = iduv
                else
                   k = iduv + 1
                endif
             endif
          endif
       endif
    end do
    zeroec = k
    return
  end function zeroec
  !
  !
  real function wgtdst(iduv, cnt, systype, wgttype)
    use sysvars, only: wgtf2smpl
    implicit none
    ! cnt = 1: solution   cnt = 2: reference solvent
    integer, intent(in) :: iduv, cnt
    character(len=5), intent(in) :: systype
    character(len=4), intent(in) :: wgttype
    real :: fsln, fref, wght, factor
    logical :: errtag
    integer :: jdg
    fsln = edist(iduv)
    fref = edens(iduv)
    errtag = .false.
    jdg = 1
    if(cnt == 2) then
       if(systype == 'slncv')  jdg = 9
       if(systype == 'extsl')  jdg = 9
       if(wgttype == 'smpl')   jdg = 2
       if(wgtf2smpl == 'yes')  jdg = 2
    endif
    if((jdg /= 1) .and. (jdg /= 2)) errtag = .true.
    if(jdg == 2) wght = fref
    if(jdg == 1) then
       wght = 0.0
       select case(wgttype)
       case('smpl')
          wght = fsln
          if(cnt == 2) errtag = .true.
       case('geom')
          factor = fsln * fref
          if(factor > zero) wght = sqrt(factor)
       case default              ! corresponding to wgttype = 'harm'
          factor = fsln + fref
          if(factor > zero) wght = fsln * fref / factor
       end select
    endif
    if(errtag) stop ' Bug in the program'
    wgtdst = wght
    return
  end function wgtdst
  !
  !
  real function sfewgt(fsln, fref)
    implicit none
    real, intent(in) :: fsln, fref
    real :: wght, factor
    if(fsln >= fref) then
       wght = 1.0
    else
       factor = (fsln - fref) / (fsln + fref)
       wght = 1.0 - factor ** 2
    endif
    sfewgt = wght
    return
  end function sfewgt
  !
  !
  real function pyhnc(indpmf, cnt)
    implicit none
    real, intent(in) :: indpmf
    integer, intent(in) :: cnt
    real :: intg, factor
    factor = indpmf / kT
    select case(cnt)
    case(1, 2)
       if(factor < - zero) then
          if(cnt == 1) intg = factor + factor / (exp(- factor) - 1.0)
          if(cnt == 2) intg = log(1.0 - factor) * (1.0 / factor - 1.0)
          intg = intg + 1.0
       else
          intg = factor / 2.0
       endif
    case(3)
       if(factor >= zero) then
          intg = 1.0 - log(1.0 + factor) * (1.0 / factor + 1.0)
       else
          intg = - factor / 2.0
       endif
    case default
       stop "Incorrct cnt argument in pyhnc"
    end select
    pyhnc = intg
    return
  end function pyhnc
  !
  !
  subroutine distnorm
    implicit none
    integer :: iduv, iduvp, pti, cnt, itrcnt
    real :: factor, ampl, lcsln, lcref, errtmp
    real, dimension(:), allocatable :: correc, edhst
    real, dimension(:,:), allocatable :: edmcr
    allocate( correc(gemax), edhst(gemax), edmcr(gemax, gemax) )
    do cnt = 1, 2     ! cnt = 1: solution   cnt = 2: reference solvent
       if(cnt == 1) then
          edhst(:) = edist(:)
          if(slncor == 'yes') edmcr(:,:) = edscr(:,:)
       endif
       if(cnt == 2) then
          edhst(:) = edens(:)
          edmcr(:,:) = ecorr(:,:)
       endif

       do pti = 1, numslv
          factor = sum( edhst, mask = (uvspec == pti) )
          if(factor > zero) then
             factor = nummol(pti) / factor
          else
             factor = 0.0
          endif
          do iduv = 1, gemax
             if(uvspec(iduv) == pti) then
                edhst(iduv) = factor * edhst(iduv)
             endif
          end do
       end do
        
       if((cnt == 1) .and. (slncor /= 'yes')) goto 5555

       errtmp = error + 1.0
       itrcnt = 0
       correc(:) = 1.0
       do while((errtmp > error) .and. (itrcnt <= itrmax))
          do iduv = 1, gemax
             lcsln = 0.0
             do pti = 1, numslv
                ampl = sum( correc(:) * edmcr(:, iduv), mask = (uvspec == pti) )
                if(ampl > zero) lcsln = lcsln + nummol(pti) / ampl
             end do
             lcsln = lcsln / real(numslv)
             correc(iduv) = lcsln * edhst(iduv)
          end do
          do iduv = 1, gemax
             do iduvp = 1, gemax
                ampl = correc(iduv) * correc(iduvp)
                edmcr(iduvp,iduv) = ampl * edmcr(iduvp,iduv)
             end do
          end do
          errtmp = maxval( abs(correc(:) - 1.0), mask = (edhst(:) > zero) )
          itrcnt = itrcnt + 1
          if(itrcnt >= itrmax) then
             write(6, *) ' The optimzation of the correlation matrix'
             write(6, *) '  did not converge with an error of ', errtmp
             stop
          endif
       end do

5555   continue
       if(cnt == 1) then
          edist(:) = edhst(:)
          if(slncor == 'yes') edscr(:,:) = edmcr(:,:)
       endif
       if(cnt == 2) then
          edens(:) = edhst(:)
          ecorr(:,:) = edmcr(:,:)
       endif
    end do

    deallocate( correc, edhst, edmcr )

    return
  end subroutine distnorm
  !
  !
  subroutine distshow
    implicit none
    integer :: iduv, pti, cnt, ecmin, ecmax, k, ilist(gemax)
    real :: factor, ratio
    real, dimension(:), allocatable :: edhst

    allocate( edhst(gemax) )
    do iduv = 1, gemax
       ilist(iduv) = iduv
    end do
    write(6, *)

    do pti = 1, numslv
       do cnt = 1, 2     ! cnt = 1: solution   cnt = 2: reference solvent
          if(cnt == 1) then
             edhst(:) = edist(:)
             if(numslv == 1) then
                write(6, "(A)") " SOLUTION"
             else
                write(6, "(A,i4,A)") " SOLUTION for", pti, "-th species"
             endif
          endif
          if(cnt == 2) then
             edhst(:) = edens(:)
             if(numslv == 1) then
                write(6, "(A)") " INSERTION"
             else
                write(6, "(A,i4,A)") " INSERTION for", pti, "-th species"
             endif
          endif

          ecmin = minval( ilist, &
                          mask = ((uvspec == pti) .and. (edhst > zero)) )
          ecmax = maxval( ilist, &
                          mask = ((uvspec == pti) .and. (edhst > zero)) )
          k = count( mask = (edhst(ecmin:ecmax) > zero) )
          ratio = real(k) / real(ecmax - ecmin + 1)
          factor = sum( edhst(ecmin:ecmax), mask = (edhst(ecmin:ecmax) > zero) )

          do iduv = ecmin, ecmax
             if((cnt == 2) .and. (edhst(iduv) <= zero)) then
                write(6, "(A,i5,A,g14.6)") "     No sampling at ", iduv, &
                                           " with energy ", uvcrd(iduv)
             endif
          end do

          if(cnt == 1) write(6, "(A,g12.4)") "     Nonzero component ratio in solution  = ", ratio
          if(cnt == 2) write(6, "(A,g12.4)") "     Nonzero component ratio at insertion = ", ratio
          write(6, "(A,g12.4)") "          Number of interacting molecules = ", factor

          if(cnt == 1) then
             write(6, "(A,g15.7)") "     Minimum energy in solution   =", uvcrd(ecmin)
             write(6, "(A,g15.7)") "     Maximum energy in solution   =", uvcrd(ecmax)
          endif
          if(cnt == 2) then
             write(6, "(A,g15.7)") "     Minimum energy at insertion  =", uvcrd(ecmin)
             write(6, "(A,g15.7)") "     Maximum energy at insertion  =", uvcrd(ecmax)
          endif
       end do
    end do

    deallocate( edhst )

    return
  end subroutine distshow
  !
end module sfecalc
!
!
!
module opwrite
  use sysvars, only: clcond, uvread, slfslt, infchk, &
       prmmax, numrun, numslv, zero, &
       pickgr, msemin, msemax, mesherr, &
       slfeng, chmpt, aveuv, blkuv, svgrp, svinf
  implicit none
  integer :: grref
  real, dimension(:), allocatable :: mshdif
contains
  !
  subroutine wrtresl
    !
    implicit none
    integer :: prmcnt, pti, k, group, inft
    real :: factor, valcp
    !
    if(slfslt == 'yes') write(6, "(A,f12.4,A)") "  Self-energy of the solute   =   ", slfeng, "  kcal/mol"
    !
    if((clcond == 'basic') .or. (clcond == 'range')) then
       write(6,*)
       if((numslv > 1) .and. (uvread /= 'yes')) write(6, 331) aveuv(1:numslv)
331    format('  Solute-solvent energy       =   ', 9999f12.4)
       factor = sum( aveuv(1:numslv) )
       if(slfslt == 'yes') factor = factor + slfeng
       write(6, "(A,f12.4,A)") "  Total solvation energy      =   ", factor, "  kcal/mol"
    endif
    !
    if(clcond == 'basic') then
       if(numslv > 1) write(6, 351) chmpt(1:numslv, 1, 1)
351    format('  Solvation free energy       =   ', 9999f12.4)
       write(6, "(A,f12.4,A)") "  Total solvation free energy =   ", chmpt(0,1,1), "  kcal/mol"
    endif
    !
    if((clcond == 'range') .or. (clcond == 'merge')) then
       do prmcnt = 1, prmmax
          group = svgrp(prmcnt)
          if(group == pickgr) then
             grref = prmcnt
             exit
          endif
       end do
       allocate( mshdif(msemin:msemax) )
       mshdif(:) = -1.0
    endif
    !
    if(clcond == 'range') then
       if(numslv == 1) then
          k = 0
       else
          k = numslv
       endif
       do pti = 0, k
          write(6, *)
          if(infchk == 'yes') then
             if(pti == 0) then
                write(6, "(A)") " group  inft  solvation free energy difference"
             else
                write(6, "(A,i3,A)") "               ", pti, "-th component difference"
             endif
          else
             if(pti == 0) then
                write(6, "(A)") " group    solvation free energy   difference"
             else
                write(6, "(A,i3,A)") "           ", pti, "-th component       difference"
             endif
          endif
          do prmcnt = 1, prmmax
             group = svgrp(prmcnt)
             inft = svinf(prmcnt)
             valcp = chmpt(pti, prmcnt, 1)
             factor = valcp - chmpt(pti, grref, 1)
             if(infchk == 'yes') then
                write(6, "(i4,i7,f17.5,f18.5)") group, inft, valcp, factor
             else
                write(6, "(i4,f20.5,f18.5)") group, valcp, factor
             endif
             if((pti == 0) .and. (inft == 0) .and. &
                (group >= msemin) .and. (group <= msemax)) then
                mshdif(group) = abs(factor)
             endif
          end do
       end do
    endif
    !
    if(clcond == 'merge') call wrtmerge
    !
    if((clcond == 'range') .or. (clcond == 'merge')) then
       factor = maxval( mshdif(msemin:msemax), &
                        mask = (mshdif(msemin:msemax) > zero) )
       if(factor > mesherr) then
          write(6, *)
          write(6, "(A,f8.3,A,g12.3)") " Warning: mesh error is ", factor, &
              " kcal/mol and is larger than the recommended value of ", &
              mesherr, " kcal/mol"
       endif
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
    implicit none
    integer :: prmcnt, cntrun, group, inft, pti, i, j, k, m
    real :: avecp, stdcp, avcp0, recnt, slvfe, shcp(large)
    real, dimension(:,:), allocatable :: wrtdata
    !
    allocate( wrtdata(0:numslv, numrun) )
    if(uvread == 'yes') then
       wrtdata(0:numslv, 1:numrun) = blkuv(0:numslv, 1:numrun)
       write(6, *)
       write(6, *)
       write(6, "(A)") " cumulative average & 95% error for solvation energy"
       call wrtcumu(wrtdata)
       write(6, *)
    endif
    !
    do pti = 0, numslv
       if((numslv == 1) .and. (pti /= 0)) cycle
       avcp0 = sum( chmpt(pti, grref, 1:numrun) ) / real(numrun)
       do prmcnt = 1, prmmax
          group = svgrp(prmcnt)
          inft = svinf(prmcnt)
          avecp = 0.0
          stdcp = 0.0
          do cntrun = 1, numrun
             slvfe = chmpt(pti, prmcnt, cntrun)
             avecp = avecp + slvfe
             stdcp = stdcp + slvfe ** 2
          end do
          recnt = real(numrun)
          avecp = avecp / recnt
          stdcp = sqrt(recnt / (recnt - 1.0)) &
                * sqrt(stdcp / recnt - avecp ** 2)
          stdcp = 2.0 * stdcp / sqrt(recnt)
          if(prmcnt == 1) then
             write(6, *)
             if(pti == 0) then
                if(infchk == 'yes') then
                   write(6, "(A)") " group  inft  solvation free energy     error          difference"
                else
                   write(6, "(A)") " group    solvation free energy     error          difference"
                endif
             endif
             if(numslv > 1) then
                if(pti == 0) then
                   write(6, "(A)") "  total solvation free energy"
                else
                   write(6, "(A,i2,A)") "  contribution from ", pti, "-th solvent component"
                endif
             endif
          endif
          if(infchk == 'yes') then
             write(6, "(i4,i7,f17.5,2f18.5)") group, inft, avecp, stdcp, (avecp - avcp0)
          else
             write(6, "(i4,f20.5,2f18.5)") group, avecp, stdcp, (avecp - avcp0)
          endif

          if((pti == 0) .and. (inft == 0) .and. &
             (group >= msemin) .and. (group <= msemax)) then
             mshdif(group) = abs(avecp - avcp0)
          endif
       end do
    end do
    !
    write(6, *)
    write(6, *)

    do pti = 0, numslv
       if((numslv == 1) .and. (pti /= 0)) cycle
       do prmcnt = 1, prmmax
          group = svgrp(prmcnt)
          inft = svinf(prmcnt)
          shcp(1:numrun) = chmpt(pti, prmcnt, 1:numrun)
          if(prmcnt == 1) then
             if(infchk == 'yes') then
                if(numslv == 1) then
                   write(6, "(A)") " group  inft   Estimated free energy (kcal/mol)"
                else
                   if(pti == 0) then
                      write(6, "(A)") " group  inft   Estimated free energy: total (kcal/mol)"
                   else
                      write(6, *)
                      write(6, "(A,i2,A)") " group  inft   Estimated free energy:", pti, "-th solvent contribution (kcal/mol)"
                   endif
                endif
             else
                if(numslv == 1) then
                   write(6, "(A)") " group   Estimated free energy (kcal/mol)"
                else
                   if(pti == 0) then
                      write(6, "(A)") " group   Estimated free energy: total (kcal/mol)"
                   else
                      write(6, *)
                      write(6, "(A,i2,A)") " group   Estimated free energy:", pti, "-th solvent contribution (kcal/mol)"
                   endif
                endif
             endif
          endif

          k = (numrun - 1) / 5
          if(infchk == 'yes') then
             if(k == 0) then
                write(6, "(i4,i7,5f13.4)") group, inft, shcp(1:numrun)
             else
                write(6, "(i4,i7,5f13.4)") group, inft, shcp(1:5)
                if(k > 1) then
                   do i = 1, k - 1
                      write(6, "('           ',5f13.4)") &
                               (shcp(5 * i + m), m = 1, 5)
                   enddo
                endif
                if(5 * k < numrun) then
                   write(6, "('           ',5f13.4)") &
                            (shcp(m), m = 5 * k + 1, numrun)
                endif
             endif
          else
             if(k == 0) then
                write(6, "(i4,'  ',5f13.4)") group, shcp(1:numrun)
             else
                write(6, "(i4,'  ',5f13.4)") group, shcp(1:5)
                if(k > 1) then
                   do i = 1, k - 1
                      write(6, "('      ',5f13.4)") &
                               (shcp(5 * i + m), m = 1, 5)
                   enddo
                endif
                if(5 * k < numrun) then
                   write(6, "('      ',5f13.4)") &
                            (shcp(m), m = 5 * k + 1, numrun)
                endif
            endif
          endif
       end do
    end do
    !
    wrtdata(0:numslv, 1:numrun) = chmpt(0:numslv, grref, 1:numrun)
    write(6, *)
    write(6, *)
    write(6, "(A)") " cumulative average & 95% error for solvation free energy"
    call wrtcumu(wrtdata)
    deallocate( wrtdata )

    return
  end subroutine wrtmerge


  subroutine wrtcumu(wrtdata)
    use sysvars, only: large, zero
    implicit none
    real, intent(in) :: wrtdata(0:numslv, numrun)
    integer :: cntrun, pti
    real :: avecp, factor, slvfe, recnt, shcp(large)
    real, dimension(:), allocatable :: runcp, runer

    allocate( runcp(0:numslv), runer(0:numslv) )
    runcp(:) = 0.0
    runer(:) = 0.0

    if(numslv == 2) write(6, "(A)") "              total             1st component         2nd component"

    do cntrun = 1, numrun
       recnt = real(cntrun)
       do pti = 0, numslv
          slvfe = wrtdata(pti, cntrun)
          runcp(pti) = runcp(pti) + slvfe
          runer(pti) = runer(pti) + slvfe ** 2
          avecp = runcp(pti) / recnt
          shcp(2 * pti + 1) = avecp
          if(cntrun > 1) then
             factor = runer(pti) / recnt - avecp ** 2
             if(factor <= zero) then
                shcp(2 * pti + 2) = 0.0
             else
                shcp(2 * pti + 2) = (2.0 / sqrt(recnt) ) &
                                  *sqrt(recnt / (recnt - 1.0)) * sqrt(factor)
             endif
          endif
       end do
       if(cntrun == 1) then
          do pti = 0, numslv
             shcp(pti + 1) = shcp(2 * pti + 1)
          end do
          if(numslv == 1) then
             write(6, "(i3,f11.4)") cntrun, shcp(1)
          else
             write(6, 771) cntrun, shcp(1), shcp(2:numslv+1)
771          format(i3,f11.4,9999f22.4)
          endif
       else
          if(numslv == 1) then
             write(6, "(i3,2f11.4)") cntrun, shcp(1:2)
          else
             write(6, 772) cntrun, shcp(1), shcp(1:2*numslv+2)
772          format(i3,9999f11.4)
          endif
       endif
    end do

    deallocate( runcp, runer )

    return
  end subroutine wrtcumu
  !
end module opwrite

program sfemain
  use sysvars, only: numrun, prmmax, init_sysvars
  use sysread, only: defcond, datread
  use sfecalc, only: chmpot
  use opwrite, only: wrtresl
  implicit none
  integer :: cntrun, prmcnt
  call init_sysvars
  call defcond
  do cntrun = 1, numrun
     call datread(cntrun)
     do prmcnt = 1, prmmax
        call chmpot(prmcnt, cntrun)
     end do
  end do
  call wrtresl
  stop
end program sfemain
