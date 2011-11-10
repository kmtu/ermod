! -*- F90 -*-
module ptinsrt
  implicit none
  !
  !  test particle insertion of the solute
  !
  real, save :: unrn
  !
  !  single-solute trajectrory file           used only when slttype = 3
  character(*), parameter :: slttrj='SltConf'    ! solute filename
  integer, parameter :: slcnf=31                 ! solute file ID
  character(*), parameter :: sltwgt='SltWght'    ! solute weight filename
  integer, parameter :: swinf=32                 ! solute weight ID
  !
  !  insertion against reference structure    used only when inscnd = 3
  !   refmlid : superposition reference among solvent species
  !             --- 0 : not reference
  !                 1 : reference solvent  2 : reference solute
  !             value set in subroutine setparam
  !   refsatm_impl : specification of the reference site
  !             --- 0 : not reference
  !                 1 : reference solvent  2 : reference solute
  !   refspos : coordiantes of interaction site for reference structure
  !   sltcen : coordinate of the solute center
  !   sltqrn : quarternion for the solute orientation
  !   movmax : number of Monte Carlo moves
  !   trmax : maximum of translational Monte Carlo move
  !   agmax : maximum of orientational Monte Carlo move
  integer, dimension(:),   allocatable :: refsatm_impl
  real, dimension(:,:),    allocatable :: refspos
  real, save :: sltcen(3),sltqrn(0:3)
  integer, parameter :: movmax=10
  real, parameter :: trmax=0.20e0
  real, parameter :: agmax=0.10e0
  !   file for reference structure
  character(*), parameter :: reffile='RefInfo'  ! reference structure
  integer, parameter :: refio=71                ! reference structure IO
  !   specifier to treat the reference as the total or as the system part only
  integer, parameter :: reftot=99, refsys=98
  !
  !
contains
  subroutine instslt(wgtslcf,caltype)
    use engmain, only: nummol,slttype,insfit,inscnd,inscfg,numslt,sltlist,iseed, &
         numsite, mol_begin_index, sitepos, bfcoord
    implicit none
    integer insml,m
    real wgtslcf,pcom(3),qrtn(0:3)
    character*4 caltype
    character*3 insyn

    integer :: nsite, molb, mole
    
    call instwgt(wgtslcf,caltype)

    if(caltype.eq.'init') call urand_init(iseed)
    if((caltype.eq.'init').or.(caltype.eq.'last')) then
       if(slttype.eq.3) call getsolute(0,caltype)       ! flexible solute
       return
    endif

    m=0
    if(numslt.ne.1) m=9
    if((numslt.eq.1).and.(sltlist(1).ne.nummol)) m=9
    if(m.ne.0) call insrt_stop('set')    ! incorrect solute specification

    insml=nummol                                       ! inserted solute
    if(slttype.eq.3) call getsolute(insml,'conf')      ! flexible solute
    if(insfit == 1) then
       call reffit(insml)
    else
       nsite = numsite(insml)
       molb = mol_begin_index(insml)
       mole = mol_begin_index(insml + 1) - 1
       
       sitepos(1:3, molb:mole) = bfcoord(1:3, 1:nsite)
    end if

    insyn='not'
    do while(insyn.eq.'not')
       call set_shift_com(insml, wgtslcf)
       call apply_orientation(insml)
       call insscheme(insml,insyn)       ! user-defined scheme to apply change / reject configuration
    end do
    
    return
  end subroutine instslt

  subroutine set_shift_com(insml, weight)
    use engmain, only: insposition, inscnd, &
         numsite, mol_begin_index, &
         lwreg, upreg, &
         sitepos, cell, invcl, celllen, &
         boxshp, SYS_NONPERIODIC, &
         pi
    implicit none
    integer, intent(in) :: insml
    real, intent(inout) :: weight
    
    integer :: i
    real :: com(3), syscen(3), r(3), norm, dir, t, s

    logical, save :: use_uniform = .false. ! used only for insposition = 4

    select case(insposition)
    case(0)
       ! fully random position within periodic box
       if(boxshp == SYS_NONPERIODIC) call insrt_stop('geo') ! system has to be periodic
       do i = 1, 3
          call urand(r(i))
       end do
       com(:) = matmul(cell(:, :), r(:))

       call set_solute_com(insml, com)
       return
    case(1) ! spherical random position
       call get_system_com(syscen)
       ! rejection method. Probability may not be good esp. lwreg ~ upreg.
       do
          do i = 1, 3
             call urand(r(i))
          end do
          r(:) = r(:) * 2.0 - 1.0
          norm = sum(r ** 2)
          if(norm < 1 .and. norm > (lwreg/upreg) ** 2) exit
       end do
       com(:) = syscen(:) + r(:) * upreg
       call set_solute_com(insml, com)
       return
    case(2) ! slab position
       if(boxshp == SYS_NONPERIODIC) call insrt_stop('geo') ! system has to be periodic
       
       call get_system_com(syscen)
       
       call urand(r(1))
       call urand(r(2))
       
       call urand(t)
       call urand(dir)
       
       t = t * (upreg - lwreg) + lwreg
       if(dir > 0.5) t = -t
       r(3) = t / celllen(3) + dot_product(invcl(3, :), syscen(:))
       
       com(:) = matmul(cell(:, :), r(:))
       call set_solute_com(insml, com)
       return
    case(3)
       ! fixed position
       ! do nothing
       return
    case(4)
       ! 50% mixture of uniform distribution and weighted distribution
       use_uniform = .not. use_uniform
       if(use_uniform) then
          ! use random position
          do i = 1, 3
             call urand(r(i))
          end do
          com(:) = matmul(cell(:, :), r(:))

          call set_solute_com(insml, com)
          return
       else
          ! for weighted insertion
          ! FIXME: this routine does not work for skewed periodic box
          ! get three N(0, 1) values       
          do i = 1, 3
             do 
                r(i) = nrand() * upreg ! to N(0, upreg)
                if(abs(r(i)) < celllen(i) / 2) exit
             end do
             ! N(0, upreg) \propto \exp [- x^2 / (2 * upreg ^ 2)]
             ! Z = \int_{-L/2}^{L/2} \exp [ - 1/2 x^2 / (2 * upreg^2)] dx
             ! weight should cancel this value ...
             ! also adjusted by celllen(i) to be equal weight to uniform distribution.
             weight = weight * exp(r(i) ** 2 / (2 * upreg ** 2)) *&
                  sqrt(2.0e0 * pi) * upreg * &
                  erf(celllen(i) / (2.0e0 * sqrt(2.0e0) * upreg)) / celllen(i)
          end do

          call shift_solute_com(insml, r)
          return
       endif
    case default
       stop "Unknown insposition"
    end select

  contains
    subroutine set_solute_com(insml, com)
      use engmain, only: numsite, &
           mol_begin_index, mol_end_index, &
           sitepos, sitemass
      use bestfit, only: com_shift, com_unshift
      implicit none
      integer, intent(in) :: insml
      real, intent(in) :: com(3)
      integer :: insb, inse, i, n
      real :: tempcom(3)
      
      n = numsite(insml)
      insb = mol_begin_index(insml)
      inse = mol_end_index(insml)
      
      call com_shift(n, sitepos(1:3, insb:inse), sitemass(insb:inse), tempcom)
      do i = 1, 3
         sitepos(i, insb:inse) = sitepos(i, insb:inse) + com(i)
      end do
    end subroutine set_solute_com

    subroutine shift_solute_com(insml, com)
      use engmain, only: numsite, &
           mol_begin_index, mol_end_index, &
           sitepos, sitemass
      implicit none
      integer, intent(in) :: insml
      real, intent(in) :: com(3)
      integer :: insb, inse, i, n
      
      n = numsite(insml)
      insb = mol_begin_index(insml)
      inse = mol_end_index(insml)
      
      do i = 1, 3
         sitepos(i, insb:inse) = sitepos(i, insb:inse) + com(i)
      end do
    end subroutine shift_solute_com

  end subroutine set_shift_com

  subroutine apply_orientation(insml)
    use engmain, only: insorient, &
         numsite, &
         mol_begin_index, mol_end_index, &
         sitepos, sitemass
    use quaternion, only: rotate_inplace
    use bestfit, only: com_shift, com_unshift
    implicit none
    integer, intent(in) :: insml
    integer :: insb, inse, i, n
    real :: com(3), randq(0:3)
    
    n = numsite(insml)
    insb = mol_begin_index(insml)
    inse = mol_end_index(insml)

    select case(insorient)
    case(0)
       ! no orientational change
       return
    case(1)
       ! random orientation
       call com_shift(n, sitepos(1:3, insb:inse), sitemass(insb:inse), com)

       ! get random quaternion (prob. success ~ 0.3)
       ! if this is really bad implement:
       ! Marsaglia, G. "Choosing a Point from the Surface of a Sphere."
       !   Ann. Math. Stat. 43, 645-646, 1972.
       do
          do i = 0, 3
             call urand(randq(i))
          end do
          randq(:) = randq(:) * 2.0 - 1.0
          if(sum(randq ** 2) < 1) exit
       end do

       randq(:) = randq(:) / sqrt(sum(randq**2)) ! set on unit hyper-sphere surface
       call rotate_inplace(numsite(insml), sitepos(1:3, insb:inse), randq)

       call com_unshift(n, sitepos(1:3, insb:inse), sitemass(insb:inse), com)
       return
    end select
  end subroutine apply_orientation

  ! FIXME: cleanup
  subroutine sltpstn(sltstat,pcom,type,tagslt)
    use engmain, only: nummol,maxsite,numatm,&
         boxshp,inscnd,inscfg,hostspec,&
         moltype,numsite,specatm,sitepos,cell,invcl,&
         lwreg,upreg
    implicit none
    integer sltstat,tagslt,stmax,sid,ati,pti,i,m,k,centag(numatm)
    real rdum,clm(3),pcom(3),qrtn(0:3),rst,dis,syscen(3),elen
    character*6 type

    if(inscfg.eq.2) then
       sltstat=1
       return
    endif

    if(type.ne.'insert') stop "BUG"
    !
    if(inscnd.eq.0) then   ! solute with random position
       if(boxshp.eq.0) call insrt_stop('geo') ! system has to be periodic
       do k=1,3
          call URAND(rdum)
          clm(k)=rdum-0.50e0
       end do
       do m=1,3
          pcom(m)=dot_product(cell(m, :), clm(:))
       end do
    endif
    !
    if((inscnd.eq.1).or.(inscnd.eq.2)) then     ! system center
       call get_system_com(syscen)
    endif
    !
    if(inscnd.eq.1) then   ! solute in spherical object or isolated droplet
       call rndmvec('p',qrtn,(lwreg/upreg))

       pcom(:)=upreg*qrtn(1:3)+syscen(:)
       dis=0.0e0
       do m=1,3
          rst=pcom(m)-syscen(m)
          dis=dis+rst*rst
       end do
       dis=sqrt(dis)
    endif
    !
    if(inscnd.eq.2) then   ! solute in slab geometry
       if(boxshp.eq.0) call insrt_stop('geo')    ! system has to be periodic
       elen=0.0e0
       do m=1,3
          elen=elen+cell(m,3)*cell(m,3)
       end do
       elen=sqrt(elen)

       do k=1,2
          call URAND(rdum)
          clm(k)=rdum-0.50e0
       end do
       do m=1,3
          rst=0.0e0
          do  k=1,2
             rst=rst+cell(m,k)*clm(k)
          end do
          pcom(m)=rst
       end do
       call URAND(rdum)
       rst=lwreg+rdum*(upreg-lwreg)
       call URAND(rdum)
       if(rdum.le.0.50e0) rst=-rst
       dis=0.0e0
       do m=1,3
          dis=dis+invcl(3,m)*syscen(m)
       end do
       dis=dis+rst/elen
       do m=1,3
          pcom(m)=pcom(m)+dis*cell(m,3)
       end do

       rst=0.0e0
       do m=1,3
          rst=rst+invcl(3,m)*(pcom(m)-syscen(m))
       end do
       dis=abs(rst)*elen
    endif
    !
    if(inscnd.eq.3) then   ! solute against reference structure
       if(type.eq.'insert') call insrt_stop('bug')
       call refsdev(dis,2,'system')
    endif
    !
    if(inscnd.eq.0) sltstat=1
    if(inscnd.gt.0) then
       if((lwreg.le.dis).and.(dis.le.upreg)) then
          sltstat=1
       else
          call insrt_stop('bug')
       endif
    endif
    !
    return
  end subroutine sltpstn
  !
  subroutine get_molecule_com(target_mol, com)
    use engmain, only: numatm, specatm, numsite

    implicit none
    integer, intent(in) :: target_mol
    real, intent(out) :: com(3)
    integer :: centag(1:numatm), ati, sid, stmax

    centag(:)=0
    stmax=numsite(target_mol)
    do sid=1,stmax
       ati=specatm(sid,target_mol)
       centag(ati)=1
    end do
    call getcen(centag,com)
  end subroutine get_molecule_com

  subroutine get_system_com(com)
    use engmain, only: nummol, numatm, specatm, numsite, hostspec, moltype
    implicit none
    real, intent(out) :: com(3)
    integer :: centag(1:numatm), ati, sid, stmax, i, m, pti

    do i=1,nummol
       pti=moltype(i)
       if(pti.eq.hostspec) m=1
       if(pti.ne.hostspec) m=0
       stmax=numsite(i)
       do sid=1,stmax
          ati=specatm(sid,i)
          centag(ati)=m
       end do
    end do
    call getcen(centag, com)
  end subroutine get_system_com

  subroutine getcen(centag,cen)             ! getting the center of mass
    use engmain, only: numatm,sitemass,sitepos
    implicit none
    integer, intent(in) :: centag(numatm)
    real, intent(out) :: cen(3)
    integer ati,m
    real wgt,sitm

    wgt=0.0e0
    cen(:)=0.0e0

    do ati=1,numatm
       if(centag(ati).eq.1) then
          sitm=sitemass(ati)
          wgt=wgt+sitm

          cen(:)=cen(:)+sitm*sitepos(:,ati)
       endif
    end do
    cen(:)=cen(:)/wgt

  end subroutine getcen
  !
  subroutine rndmvec(vectp,qrtn,lwbnd)
    implicit none
    character vectp
    integer m,inim
    real qrtn(0:3),lwbnd,rdum,factor
    if(vectp.eq.'p') inim=1
    if(vectp.eq.'q') inim=0
    factor=2.0e0
    qrtn(:) = 0.0e0
    do while((factor.gt.1.0e0).or.(factor.lt.lwbnd))
       do m=inim,3
          call URAND(rdum)
          qrtn(m)=2.0e0*rdum-1.0e0
       end do
       factor=sqrt(sum(qrtn(inim:3) ** 2))
    end do
    if(vectp.eq.'q') then
       do m=0,3
          qrtn(m)=qrtn(m)/factor
       end do
    endif
  end subroutine rndmvec
  !
  !
  subroutine coordinate(i,pcom,qrtn)
    use engmain, only: nummol,maxsite,numatm,inscfg,&
         numsite,bfcoord,specatm,sitepos
    implicit none
    integer stmax,sid,ati,m,k,i
    real pcom(3),qrtn(0:3),rotmat(3,3),rst
    stmax=numsite(i)
    if(inscfg.eq.0) call getrot(qrtn,rotmat)        ! rotation matrix
    do  sid=1,stmax
       ati=specatm(sid,i)
       do  m=1,3
          if(inscfg.eq.0) then                        ! translation + rotation
             rst=dot_product(rotmat(:, m), bfcoord(:, sid))
             rst=rst+pcom(m)
          endif
          if(inscfg.eq.1) rst=pcom(m)+bfcoord(m,sid)  ! translation only
          if(inscfg.eq.2) rst=bfcoord(m,sid)          ! no change
          sitepos(m,ati)=rst
       end do
    end do
    return
  end subroutine coordinate
!

! convert quaternion to rotation matrix
  subroutine getrot(qrtn,rotmat)
    implicit none
    real, intent(in) :: qrtn(0:3)
    real, intent(out) :: rotmat(3,3)
    rotmat(1,1)=qrtn(0)*qrtn(0)+qrtn(1)*qrtn(1)&
         -qrtn(2)*qrtn(2)-qrtn(3)*qrtn(3)
    rotmat(2,2)=qrtn(0)*qrtn(0)-qrtn(1)*qrtn(1)&
         +qrtn(2)*qrtn(2)-qrtn(3)*qrtn(3)
    rotmat(3,3)=qrtn(0)*qrtn(0)-qrtn(1)*qrtn(1)&
         -qrtn(2)*qrtn(2)+qrtn(3)*qrtn(3)
    rotmat(1,2)=2.0e0*(qrtn(1)*qrtn(2)+qrtn(0)*qrtn(3))
    rotmat(2,1)=2.0e0*(qrtn(1)*qrtn(2)-qrtn(0)*qrtn(3))
    rotmat(1,3)=2.0e0*(qrtn(1)*qrtn(3)-qrtn(0)*qrtn(2))
    rotmat(3,1)=2.0e0*(qrtn(1)*qrtn(3)+qrtn(0)*qrtn(2))
    rotmat(2,3)=2.0e0*(qrtn(2)*qrtn(3)+qrtn(0)*qrtn(1))
    rotmat(3,2)=2.0e0*(qrtn(2)*qrtn(3)-qrtn(0)*qrtn(1))
  end subroutine getrot
  ! user-defined scheme to specify inserted molecule
  ! user may reject snapshot by specifying insyn to 'not', or
  ! set coordinate in specatm
  subroutine insscheme(insml,insyn)                
    use engmain, only: nummol,maxsite,numatm,numsite,specatm,sitepos
    integer insml,stmax,sid,ati
    character*3 insyn
    insyn='yes'
    stmax=numsite(insml)
    return
  end subroutine insscheme

  ! FIXME: move to mpiproc
  subroutine insrt_stop(type)
    use engmain, only: io6
    use mpiproc                                                      ! MPI
    character(len=3), intent(in) :: type
    if(type.eq.'set') write(io6,791)
    if(type.eq.'geo') write(io6,792)
    if(type.eq.'bug') write(io6,793)
791 format(' The solute specification is incorrectly set')
792 format(' The system geometry is incorrectly set')
793 format(' There is a bug in the insertion program')
    call mpi_setup('stop')                                           ! MPI
    stop
  end subroutine insrt_stop
!
!
  subroutine getsolute(i,caltype)
    use trajectory, only: open_trajectory, close_trajectory, read_trajectory
    use engmain, only: nummol,maxsite,insposition,numsite,bfcoord
    use setconf, only: molcen
    use OUTname, only: iofmt,bxiso,toptp, OUTconfig, solute_trajectory
    use mpiproc
    character*4 caltype
    integer i,sid,stmax,m
    real xst(3),dumcl(3,3),factor
    real*4 sglfct
    character rddum
    real, dimension(:,:), allocatable :: psite
!
    if(caltype.eq.'init') then
       if(myrank /= 0) return
       call open_trajectory(solute_trajectory, slttrj)
       return
    endif
    if(caltype.eq.'last') then
       if(myrank /= 0) return
       call close_trajectory(solute_trajectory)
       return
    endif
!
    stmax=numsite(i)
    allocate( psite(3,stmax) )
    if(bxiso.eq.'not') sid=0
    if(bxiso.eq.'yes') sid=1
    if(myrank == 0) then
       call OUTconfig(psite,dumcl,stmax,sid,slcnf,'trj')
    endif
#ifndef noMPI
    ! distribute to non rank-0 nodes
    call mpi_bcast(psite, 3 * stmax, mpi_double_precision, &
         0, mpi_comm_activeprocs, ierror)
#endif

    if(insposition /= 4) call molcen(i,psite,xst,'com')
    do sid=1,stmax
       do m=1,3
          if(insposition /= 4) bfcoord(m,sid)=psite(m,sid)-xst(m)
          if(insposition == 4) bfcoord(m,sid)=psite(m,sid)
       end do
    end do
    deallocate( psite )
  end subroutine getsolute

  subroutine instwgt(wgtslcf,caltype)
    use engmain, only: slttype,wgtins
    integer m
    real wgtslcf
    character*4 caltype
    if((slttype.ne.3).or.(wgtins.ne.1)) then
       wgtslcf=1.0e0
       return
    endif
    if(caltype.eq.'init') open(unit=swinf,file=sltwgt,status='old')
    if(caltype.eq.'last') close(swinf)
    if(caltype.eq.'proc') then
       read(swinf,*) m,wgtslcf
       read(swinf,*,END=2199) m
       backspace(swinf)
       goto 2195
2199   rewind(swinf)
2195   continue
    endif
    return
  end subroutine instwgt

  ! returns random value from [0,1)
  ! Any sane compiler implements random_number
  ! (which is included in fortran 95 standards)
  subroutine urand(rndm)          ! uniform random number generator
    implicit none
    real, intent(out) :: rndm
    call random_number(rndm)
  end subroutine urand

  ! Normal random variable N(0,1)
  ! uses Box-Muller method
  real function nrand()
    use engmain, only: pi
    implicit none
    real :: r1, r2
    call urand(r1)
    call urand(r2)
    ! get (0,1] instead of [0, 1)
    r1 = 1 - r1
    nrand = sqrt(-2.0 * log(r1)) * cos(2 * pi * r2)
  end function nrand

  subroutine urand_init(seed)
    use mpiproc, only: myrank
    implicit none
    integer, intent(in) :: seed
    integer :: seedsize
    integer, allocatable :: seedarray(:)

    call random_seed(size = seedsize)
    allocate(seedarray(seedsize))
    seedarray(:) = 1

    seedarray(1) = myrank + seed
    if(seed == 0) call system_clock(count = seedarray(1))

    call random_seed(put = seedarray)
    deallocate(seedarray)
  end subroutine urand_init

  subroutine refmc
    use engmain, only: nummol,maxsite,numatm,slttype,numsite,refmlid,&
         lwreg,upreg,bfcoord,specatm,sitepos
    use bestfit
    integer i,k,m,q,ati,rfi,sid,stmax
    real xst(3),centg(3),cenrf(3)
    real factor,bfqrn(0:3),rtmbf(3,3),sltsite(3,maxsite)
    character*4 caltype
    character*8 atmtype,dump
    character eletype
    

    if(caltype.eq.'init') then
       allocate( refsatm_impl(numatm),refspos(3,numatm) )
       refsatm_impl(:) = 0
       refspos(:, :) = 0.0e0

       !  read the reference structure
       open(unit=refio,file=reffile,status='old')
       do i=1,nummol
          rfi=refmlid(i)
          if(rfi.ne.0) then
             stmax=numsite(i)
             do sid=1,stmax
                read(refio, '(12X, A4, 14X, 3F8.3)') atmtype, (xst(m), m=1,3)
                ! FIXME: this element selection looks very ugly
                atmtype = trim(atmtype)
                eletype = atmtype(1:1)
                ! atom name can be like "1HG"
                if('0' <= eletype .and. eletype <= '9') eletype = atmtype(2:2)
                ati=specatm(sid,i)
                if(eletype.eq.'H') refsatm_impl(ati) = 0      ! hydrogen atom
                if(eletype.ne.'H') refsatm_impl(ati) = rfi    ! heavy atom
                do m=1,3
                   refspos(m,ati)=xst(m)
                end do
             end do
          endif
       end do
       close(refio)
       !
       !  build the initial configuration of the solute
       if(slttype.ge.2) then
          call dispref(centg,cenrf,bfqrn,1)     ! matching reference
          call getrot(bfqrn,rtmbf)
          stmax=numsite(nummol)                 ! inserted solute molecule
          do  sid=1,stmax
             ati=specatm(sid,nummol)
             do  m=1,3
                factor=0.0e0
                do  k=1,3
                   factor=factor+rtmbf(k,m)*(refspos(k,ati)-cenrf(k))
                end do
                sitepos(m,ati)=factor+centg(m)
             end do
          end do
          call dispref(sltcen,cenrf,sltqrn,2)   ! solute center & orientation
       endif
    endif
    !
    if(caltype.eq.'inst') then                ! inserted solute molecule
       stmax=numsite(nummol)
       do sid=1,stmax
          ati=specatm(sid,nummol)
          do m=1,3
             sltsite(m,sid)=sitepos(m,ati)
          end do
          do m=1,3
             sitepos(m,ati)=bfcoord(m,sid)
          end do
       end do
       call refcen(xst,sitepos,2)      ! solute center in bfcoord
       do sid=1,stmax
          ati=specatm(sid,nummol)
          do m=1,3
             bfcoord(m,sid)=bfcoord(m,sid)-xst(m)
          end do
       end do
       call dispref(centg,cenrf,bfqrn,2)       ! matching solute
       call getrot(bfqrn,rtmbf)
       k=0
       do q=1,movmax
          call sltmove(sltcen,sltqrn,rtmbf,'frwd')
          call refsdev(factor,2,'extend')
          if((factor.lt.lwreg).or.(factor.gt.upreg)) then
             k=k+1
             call sltmove(sltcen,sltqrn,rtmbf,'back')
          endif
       end do
       if(k.eq.movmax) then                    ! when all MC are rejected
          do sid=1,stmax
             ati=specatm(sid,nummol)

             sitepos(:,ati)=sltsite(:,sid)
          end do
       endif
    endif
    !
    return
  end subroutine refmc
!
!
  subroutine refcen(cen,posatm,refmol)
    use engmain, only: nummol,maxsite,numatm,numsite,specatm
    implicit none
    integer refmol,rfyn,i,m,ati,rfi,sid,stmax
    real cen(3),posatm(3,numatm),totwgt
    totwgt=0.0e0

    cen(:)=0.0e0
    do i=1,nummol
       call getrfyn(i,rfi,rfyn,refmol)
       if(rfyn.eq.1) then
          stmax=numsite(i)
          do sid=1,stmax
             if(refsatm(sid,i).eq.rfi) then
                ati=specatm(sid,i)
                totwgt=totwgt+1.0e0

                cen(:) = cen(:) + posatm(:,ati)
             endif
          end do
       endif
    end do

    cen(:) = cen(:) / totwgt
    return
  end subroutine refcen
!
!
  subroutine sltmove(pcen,qrtn,rtmbf,mvtype)
    use engmain, only: nummol,maxsite,numatm,&
         numsite,bfcoord,specatm,sitepos
    implicit none
    integer m,k,q,ati,sid,stmax,ax1,ax2
    real pcen(3),qrtn(0:3),rtmbf(3,3)
    real rdum,factor,rotmat(3,3),rtmmov(3,3),axis(3)
    character*4 mvtype
    real, save :: odcn(3),oqrn(0:3)
    stmax=numsite(nummol)                       ! inserted solute molecule
    do m=1,3
       if(mvtype.eq.'frwd') odcn(m)=pcen(m)
       if(mvtype.eq.'back') pcen(m)=odcn(m)
    end do
    do m=0,3
       if(mvtype.eq.'frwd') oqrn(m)=qrtn(m)
       if(mvtype.eq.'back') qrtn(m)=oqrn(m)
    end do
    if(mvtype.eq.'frwd') then
       do m=1,3                             ! translation
          call URAND(rdum)
          pcen(m)=pcen(m)+trmax*(rdum+rdum-1.0e0)
       end do
       factor=2.0e0                              ! rotation
       do while(factor.gt.1.0e0)
          do m=1,3
             call URAND(rdum)
             axis(m)=rdum+rdum-1.0e0
          end do
          factor=axis(1)*axis(1)+axis(2)*axis(2)+axis(3)*axis(3)
       end do
       factor=sqrt(factor)
       do m=1,3
          axis(m)=axis(m)/factor
       end do
       call URAND(rdum)
       ax1=floor(rdum * 4)
       ax2=ax1
       do while(ax1.eq.ax2)
          call URAND(rdum)
          ax2=floor(rdum * 4)
       end do
       call URAND(rdum)
       factor=agmax*(rdum+rdum-1.0e0)
       qrtn(ax1)=cos(factor)*oqrn(ax1)-sin(factor)*oqrn(ax2)
       qrtn(ax2)=sin(factor)*oqrn(ax1)+cos(factor)*oqrn(ax2)
       factor=0.0e0
       do m=0,3
          factor=factor+qrtn(m)*qrtn(m)
       end do
       factor=sqrt(factor)
       do m=0,3
          qrtn(m)=qrtn(m)/factor
       end do
    endif
    call getrot(qrtn,rtmmov)                    ! site coordinate
    do m=1,3
       do k=1,3
          factor=0.0e0
          do q=1,3
             factor=factor+rtmmov(q,m)*rtmbf(q,k)
          end do
          rotmat(m,k)=factor
       end do
    end do
    do sid=1,stmax
       ati=specatm(sid,nummol)
       do m=1,3
          factor=0.0e0
          do k=1,3
             factor=factor+rotmat(m,k)*bfcoord(k,sid)
          end do
          sitepos(m,ati)=factor+pcen(m)
       end do
    end do
    return
  end subroutine sltmove
!
  subroutine dispref(centg,cenrf,qrtn,refmol)
    use engmain, only: nummol,maxsite,numatm,numsite,specatm,sitepos
    implicit none
    integer refmol,rfyn,i,m,k,ati,rfi,sid,stmax
    real centg(3),cenrf(3),qrtn(0:3),qrtmat(4,4),xsm(3),xrf(3)
    real totm,xp,yp,zp,xm,ym,zm,dumv(4),work(16)
    call refcen(centg,sitepos,refmol)
    call refcen(cenrf,refspos,refmol)
    totm=0.0e0
    qrtmat(:, :) = 0.0e0
    do i=1,nummol
       call getrfyn(i,rfi,rfyn,refmol)
       if(rfyn.eq.1) then
          stmax=numsite(i)
          do sid=1,stmax
             if(refsatm(sid,i).eq.rfi) then
                ati=specatm(sid,i)
                totm=totm+1.0e0
                do m=1,3
                   xsm(m)=sitepos(m,ati)-centg(m)
                end do
                do m=1,3
                   xrf(m)=refspos(m,ati)-cenrf(m)
                end do
                xp=xrf(1)+xsm(1) ; yp=xrf(2)+xsm(2) ; zp=xrf(3)+xsm(3)
                xm=xrf(1)-xsm(1) ; ym=xrf(2)-xsm(2) ; zm=xrf(3)-xsm(3)
                qrtmat(1,1)=qrtmat(1,1)+xm*xm+ym*ym+zm*zm
                qrtmat(2,1)=qrtmat(2,1)+yp*zm-ym*zp
                qrtmat(3,1)=qrtmat(3,1)+xm*zp-xp*zm
                qrtmat(4,1)=qrtmat(4,1)+xp*ym-xm*yp
                qrtmat(2,2)=qrtmat(2,2)+yp*yp+zp*zp+xm*xm
                qrtmat(3,2)=qrtmat(3,2)+xm*ym-xp*yp
                qrtmat(4,2)=qrtmat(4,2)+xm*zm-xp*zp
                qrtmat(3,3)=qrtmat(3,3)+xp*xp+zp*zp+ym*ym
                qrtmat(4,3)=qrtmat(4,3)+ym*zm-yp*zp
                qrtmat(4,4)=qrtmat(4,4)+xp*xp+yp*yp+zm*zm
             endif
          end do
       endif
    end do
    ! symmetrize
    do m=1,3
       do k=m+1,4
          qrtmat(m,k)=qrtmat(k,m)
       end do
    end do
    do m=1,4
       do k=1,4
          qrtmat(k,m)=qrtmat(k,m)/totm
       end do
    end do
    call DSYEV('V','U',4,qrtmat,4,dumv,work,16,i)
    do m=0,3
       qrtn(m)=qrtmat(m+1,1)   ! eigenvector for the minimum eigenvalue
    end do
    totm=0.0e0
    do m=0,3
       totm=totm+qrtn(m)*qrtn(m)
    end do
    totm=sqrt(totm)
    do m=0,3
       qrtn(m)=qrtn(m)/totm
    end do
    return
  end subroutine dispref
!
!
  subroutine refsdev(strchr,refmol,caltype)
    use engmain, only: nummol,maxsite,numatm,numsite,specatm,sitepos
    implicit none
    integer refmol,rfyn,i,m,k,ati,rfi,sid,stmax
    character*6 caltype
    real strchr,centg(3),cenrf(3),qrtn(0:3),rotmat(3,3)
    real totm,factor,xsm(3)
    select case(caltype)
    case('system')
       call dispref(centg,cenrf,qrtn,refsys)
    case('extend')
       call dispref(centg,cenrf,qrtn,reftot)
    case default
       call insrt_stop('bug')
    end select
    call getrot(qrtn,rotmat)
    strchr=0.0e0
    totm=0.0e0
    do i=1,nummol
       call getrfyn(i,rfi,rfyn,refmol)
       if(rfyn.eq.1) then
          stmax=numsite(i)
          do sid=1,stmax
             if(refsatm(sid,i).eq.rfi) then
                ati=specatm(sid,i)
                totm=totm+1.0e0
                do m=1,3
                   factor=0.0e0
                   do k=1,3
                      factor=factor+rotmat(m,k)*(sitepos(k,ati)-centg(k))
                   end do
                   xsm(m)=factor+cenrf(m)
                end do
                do m=1,3
                   factor=xsm(m)-refspos(m,ati)
                   strchr=strchr+factor*factor
                end do
             endif
          end do
       endif
    end do
    strchr=sqrt(strchr/totm)
    return
  end subroutine refsdev
!
!
  subroutine getrfyn(i,rfi,rfyn,refmol)
    use engmain, only: nummol,sluvid,refmlid
    implicit none
    integer i,rfi,rfyn,refmol
    rfi=refmlid(i)
    rfyn=0
    select case(refmol)
    case(1,2)
       if(rfi.eq.refmol) rfyn=1
    case(reftot)
       if(rfi.ne.0) rfyn=1
    case(refsys)
       if((rfi.ne.0).and.(sluvid(i).le.1)) rfyn=1
    case default
       call insrt_stop('bug')
    end select
    return
  end subroutine getrfyn

  ! fit to reference structure
  subroutine reffit(ligmol)
    use engmain, only: refmlid, numsite, nummol, &
         mol_begin_index, mol_end_index, &
         sitepos, sitemass, bfcoord
    use quaternion, only: rotate
    use bestfit, only: fit_a_rotate_b, fit
    implicit none
    integer, intent(in) :: ligmol

    real, allocatable, save :: ref_solv_pos(:, :)  ! solvent structure in reference
    real, allocatable, save :: ref_lig_pos(:, :)   ! ligand structure in reference
    real, allocatable, save :: works(:, :), workl(:, :)
    real, allocatable, save :: solv_mass(:) ! masked mass vector;
    real, allocatable, save :: lig_mass(:)  ! atoms not used for fitting is zeroed

    real :: crd(3)
    integer :: natom_solv, natom_lig
    integer :: solvmol
    integer :: solv_begin, solv_end
    integer :: lig_begin, lig_end
    real :: com_solv(3), com_lig(3)
    integer :: i
    
    solvmol = -1
    do i = 1, nummol
       if(refmlid(i) == 1) solvmol = i
       if(refmlid(i) == 2 .and. i /= ligmol) stop "insmol / refmlid is inconsistent"
    end do
    if(solvmol == -1) stop "reffit: failed to find reference solvent"

    natom_solv = numsite(solvmol)
    natom_lig = numsite(ligmol)
    
    solv_begin = mol_begin_index(solvmol)
    solv_end = mol_end_index(solvmol)

    lig_begin = mol_begin_index(ligmol)
    lig_end = mol_end_index(ligmol)

    if(.not. allocated(ref_solv_pos)) then
       ! first time call
       allocate(ref_solv_pos(3, natom_solv), solv_mass(natom_solv))
       allocate(ref_lig_pos(3, natom_lig), lig_mass(natom_lig))
       allocate(works(3, natom_solv), workl(3, natom_lig))
       
       ! initialize mass first
       solv_mass(:) = sitemass(solv_begin:solv_end)
       lig_mass(:) = sitemass(lig_begin:lig_end)

       open(unit=refio, file=reffile, status='old')
       call load_structure(natom_solv, ref_solv_pos, solv_mass)
       call load_structure(natom_lig, ref_lig_pos, lig_mass)
       close(refio)
    endif
    
    ! first fit reference solvent to current solvent, to get "reference ligand position"
    call fit_a_rotate_b(natom_solv, &
         sitepos(1:3, solv_begin:solv_end), ref_solv_pos, solv_mass, &
         natom_lig, &
         ref_lig_pos, workl)
    ! then fit ligand structure read from file
    call fit(natom_lig, workl, bfcoord, lig_mass, sitepos(1:3, lig_begin:lig_end))
  contains
    subroutine load_structure(n, position, mass)
      integer, intent(in) :: n
      real, intent(out) :: position(3, n)
      real, intent(inout) :: mass(n)
      real :: crd(3)
      character(len=6) :: header
      integer :: i, ix

      do ix = 1, n
         do
            ! skip until ATOM/HETATM lines
            read(refio, '(A6)', advance='no') header
            if(header == 'ATOM  ' .or. header == 'HETATM') exit
            read(refio, *)
         end do
         read(refio, '(24X, 3F8.3)') (crd(i), i = 1, 3)
         position(1:3, ix) = crd(1:3)
      end do

      ! set mass = 0 to mask fitting
      ! user can implement his own special selection rule
      ! (e.g. by using B-factor, etc.)
      ! default: hydrogen is masked
      do i = 1, n
         if(mass(i) < 1.1) mass(i) = 0.0
      end do
    end subroutine load_structure
  end subroutine reffit

  integer function refsatm(site, mol)
    use engmain, only: specatm
    implicit none
    integer, intent(in) :: site, mol
    
    refsatm = refsatm_impl(specatm(site, mol))
    
  end function refsatm

end module
