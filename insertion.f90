
module ptinsrt
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
  !   refsatm : specification of the reference site
  !             --- 0 : not reference
  !                 1 : reference solvent  2 : reference solute
  !   refspos : coordiantes of interaction site for reference structure
  !   sltcen : coordinate of the solute center
  !   sltqrn : quarternion for the solute orientation
  !   movmax : number of Monte Carlo moves
  !   trmax : maximum of translational Monte Carlo move
  !   agmax : maximum of orientational Monte Carlo move
  integer, dimension(:,:), allocatable, save :: refsatm
  real, dimension(:,:),    allocatable, save :: refspos
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
  !
  subroutine instslt(wgtslcf,caltype)
    use engmain, only: nummol,slttype,inscnd,inscfg,numslt,sltlist,iseed
    integer insml,m
    real wgtslcf,pcom(3),qrtn(0:3)
    character*4 caltype
    character*3 insyn

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

    if(inscnd.le.2) then
       if(inscfg.ne.2) then
          insyn='not'
          do while(insyn.eq.'not')
             call sltpstn(m,pcom,'insert',insml)          ! center of mass
             if(inscfg.ne.1) call rndmvec('q',qrtn,0.0e0) ! random orientation
             call coordinate(insml,pcom,qrtn)             ! site coordinate
             call insscheme(insml,insyn)                  ! user-defined scheme
          end do
       endif
       if(inscfg.eq.2) call coordinate(insml,pcom,qrtn) ! site coordinate
    endif
    !
    if(inscnd.eq.3) call refmc('inst')                 ! MC with reference
    !
    return
  end subroutine instslt
  !
  ! FIXME: cleanup
  subroutine sltpstn(sltstat,pcom,type,tagslt)

    use engmain, only: nummol,maxsite,numatm,&
         boxshp,inscnd,inscfg,hostspec,&
         moltype,numsite,specatm,sitepos,cell,invcl,&
         lwreg,upreg
    use mpiproc, only: halt_with_error
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
!
! user-defined scheme to specify inserted molecule
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
    use engmain, only: nummol,maxsite,inscfg,numsite,bfcoord
    use setconf, only: molcen
    use OUTname, only: iofmt,bxiso,toptp,skpio,OUTconfig,OUTskip, solute_trajectory
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
    call mpi_bcast(psite, 3 * stmax, mpi_double_precision, &
         0, mpi_comm_activeprocs, ierror)
#endif
    goto 3195
    !
    if(toptp.eq.'int') then
       if(iofmt.eq.'yes') read(slcnf,*,END=3199) m
       if(iofmt.eq.'not') read(slcnf,END=3199) m
    endif
    if(toptp.eq.'rlv') then
       if(iofmt.eq.'yes') read(slcnf,*,END=3199) factor
       if(iofmt.eq.'not') read(slcnf,END=3199) factor
    endif
    if(toptp.eq.'rsg') then
       if(iofmt.eq.'yes') read(slcnf,*,END=3199) sglfct
       if(iofmt.eq.'not') read(slcnf,END=3199) sglfct
    endif
    if(toptp.eq.'chr') then
       if(iofmt.eq.'yes') read(slcnf,*,END=3199) rddum
       if(iofmt.eq.'not') read(slcnf,END=3199) rddum
    endif
    !
    backspace(slcnf)

    goto 3195
3199 rewind(slcnf)
    call OUTskip(slcnf,iofmt,skpio)
3195 continue
!
    if(inscfg.ne.2) call molcen(i,psite,xst,'com')
    do sid=1,stmax
       do m=1,3
          if(inscfg.ne.2) bfcoord(m,sid)=psite(m,sid)-xst(m)
          if(inscfg.eq.2) bfcoord(m,sid)=psite(m,sid)
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
!
  subroutine refmc(caltype)
    use engmain, only: nummol,maxsite,numatm,slttype,numsite,refmlid,&
         lwreg,upreg,bfcoord,specatm,sitepos
    integer i,k,m,q,ati,rfi,sid,stmax
    real xst(3),centg(3),cenrf(3)
    real factor,bfqrn(0:3),rtmbf(3,3),sltsite(3,maxsite)
    character*4 caltype
    character*8 atmtype,dump
    character eletype
    !
    if(caltype.eq.'init') then
       allocate( refsatm(maxsite,nummol),refspos(3,numatm) )
       do i=1,nummol
          do sid=1,maxsite
             refsatm(sid,i)=0
          end do
       end do
       do ati=1,numatm
          do m=1,3
             refspos(m,ati)=0.0e0
          end do
       end do
       !
       !  read the reference structure
       open(unit=refio,file=reffile,status='old')
       do i=1,nummol
          rfi=refmlid(i)
          if(rfi.ne.0) then
             stmax=numsite(i)
             do sid=1,stmax
                read(refio,*) dump,k,atmtype,dump,q,&
                     (xst(m), m=1,3),factor,factor
                eletype=atmtype(1:1)
                if(eletype.eq.'H') refsatm(sid,i)=0      ! hydrogen atom
                if(eletype.ne.'H') refsatm(sid,i)=rfi    ! heavy atom
                ati=specatm(sid,i)
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
       call refcen(xst,refsatm,sitepos,2)      ! solute center in bfcoord
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
  subroutine refcen(cen,lstatm,posatm,refmol)
    use engmain, only: nummol,maxsite,numatm,numsite,specatm
    integer refmol,rfyn,i,m,ati,rfi,sid,stmax,lstatm(maxsite,nummol)
    real cen(3),posatm(3,numatm),totwgt
    totwgt=0.0e0

    cen(:)=0.0e0
    do i=1,nummol
       call getrfyn(i,rfi,rfyn,refmol)
       if(rfyn.eq.1) then
          stmax=numsite(i)
          do sid=1,stmax
             if(lstatm(sid,i).eq.rfi) then
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
       ax1=mod(int(abs(unrn)),4)
       ax2=ax1
       do while(ax1.eq.ax2)
          call URAND(rdum)
          ax2=mod(int(abs(unrn)),4)
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
    call refcen(centg,refsatm,sitepos,refmol)
    call refcen(cenrf,refsatm,refspos,refmol)
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
!
end module
