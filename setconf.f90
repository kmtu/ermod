module OUTglobal
  implicit none
  !  renaming global variables in package program
  integer GLBnrun,GLBntype,GLBens,GLBpbc,GLBclt,GLBspl
  real GLBtmp,GLBelc,GLBvdc,GLBsrc
  integer GLBew1,GLBew2,GLBew3
#ifndef trjctry
#ifdef DLPOLY
contains
  subroutine DLglobal(nstinit,nstrun,ntpmls,keyens,imcon,keyfce,&
       temp,rcut,rvdw,alpha,nospl,kmax1,kmax2,kmax3)
    integer nstinit,nstrun,ntpmls,keyens,imcon,keyfce,nospl
    real temp,rcut,rvdw,alpha
    integer kmax1,kmax2,kmax3
    GLBnrun=nstrun-nstinit ; GLBntype=ntpmls
    GLBens=keyens ; GLBpbc=imcon ; GLBclt=keyfce ; GLBspl=nospl
    GLBtmp=temp ; GLBelc=rcut ; GLBvdc=rvdw ; GLBsrc=alpha
    GLBew1=kmax1 ; GLBew2=kmax2 ; GLBew3=kmax3
    return
  end subroutine DLglobal
#endif
#endif
end module OUTglobal
!
!
module OUTname
!  renaming outside parameters and parameters to avoid conflict
  use OUTglobal
  use trajectory, only: handle
!
  implicit none
  integer, parameter :: iotrj=99                 ! trajectory file IO
  character(*), parameter :: trjfile='HISTORY'   ! trajectory filename
  integer, parameter :: mdinf=89                 ! MD info file IO
  character(*), parameter :: inffile='MDinfo'    ! MD info filename

  character*3, save :: iofmt      ! formatted or not for trajectory file
  character*3, save :: cltrd      ! seperate file for cell trajectory file
  character*3, save :: mdird      ! seperate file for MD info
  character*3, save :: bxiso      ! PBC adopted for single-molecule MD
  character*3, save :: toptp      ! type of first variable in each line
  real, save :: sgmcnv,chgcnv,engcnv,lencnv      ! unit conversion factor

  integer OUTens,OUTbxs,OUTcmb,OUTclt,OUTspo
  real OUTtemp,OUTelc,OUTlwl,OUTupl,OUTscr
  integer OUTew1,OUTew2,OUTew3,OUTms1,OUTms2,OUTms3

  integer OUTnrun,OUTntype
  integer, dimension(:), allocatable :: OUTnmol,OUTsite
  real, dimension(:), allocatable :: OUTstmass,OUTcharge
  real, dimension(:), allocatable :: OUTljene,OUTljlen
  integer TotAtm
  logical :: use_mdlib

  type(handle) :: history_trajectory
  type(handle) :: solute_trajectory

contains

  subroutine OUTinitial
    iofmt='yes' ; cltrd='not' ; mdird='yes'      ! default
    bxiso='not' ; toptp='int'                    ! default
    sgmcnv=2.0e0**(5.0e0/6.0e0)                  ! from Rmin/2 to sigma
    chgcnv=1.0e0/1.60217653e-19                  ! from C to elementary
    engcnv=6.0221415e23/4.184e3                  ! from J to kcal/mol
    lencnv=1.0e1                                 ! from nm to Angstrom
  end subroutine OUTinitial


  subroutine opentrj
    use trajectory, only: open_trajectory
    call open_trajectory(history_trajectory, trjfile)
  end subroutine opentrj

  subroutine initconf
    call OUTinitial
  end subroutine initconf

  subroutine closetrj
    use trajectory, only: close_trajectory
    call close_trajectory(history_trajectory)
  end subroutine closetrj

  subroutine finiconf
  end subroutine finiconf
!
  ! Initialization 2nd phase - read topologies from mother MD program, or read MDinfo
  subroutine OUTrename
    implicit none
    integer i

    open(unit=mdinf,file=inffile,status='old')
    read(mdinf,*) OUTnrun,OUTntype
    allocate( OUTnmol(OUTntype),OUTsite(OUTntype) )

    read(mdinf,*) (OUTnmol(i), i=1,OUTntype)
    read(mdinf,*) (OUTsite(i), i=1,OUTntype)
!
    TotAtm=0
    do i=1,OUTntype
       TotAtm=TotAtm+OUTnmol(i)*OUTsite(i)
    end do
    allocate( OUTstmass(TotAtm),OUTcharge(TotAtm) )
    allocate( OUTljene(TotAtm),OUTljlen(TotAtm) )

    close(mdinf)
    return
  end subroutine OUTrename

  ! read system setup (coulomb rule, LJ, etc..)
  subroutine OUTintprm
  end subroutine OUTintprm

  ! Get molecular configuration
  subroutine OUTconfig(OUTpos,OUTcell,OUTatm,OUTbox,OUTtrj,rdconf)
    use mpiproc, only: halt_with_error
    use trajectory, only: read_trajectory, open_trajectory, close_trajectory
    implicit none
    integer, intent(in) :: OUTatm, OUTbox, OUTtrj
    real, intent(out) :: OUTpos(3, OUTatm), OUTcell(3, 3)
    character(len=3), intent(in) :: rdconf
    integer :: status
!
    
    OUTcell(:, :) = 0.0e0
    if(OUTtrj == 0) then
       call read_trajectory(history_trajectory, OUTatm, (OUTbox == 1), OUTpos, OUTcell, status)
       if(status /= 0) call halt_with_error("trj")
    else
       call read_trajectory(solute_trajectory, OUTatm, .false., OUTpos, OUTcell, status)
       if(status /= 0) then
          call close_trajectory(solute_trajectory)
          call open_trajectory(solute_trajectory, "SltConf")
          call read_trajectory(solute_trajectory, OUTatm, .false., OUTpos, OUTcell, status)
          if(status /= 0) then
             stop "Failed to reload solute trajectory"
          endif
       endif
    end if
    return
  end subroutine
!
end module
!
!
module setconf
  implicit none

contains
!
!  setting molecular simulation parameters
!
  subroutine iniparam
    use engmain, only: init_params,&
         iseed,&
         skpcnf,corrcal,&
         slttype,sltpick,refpick,wgtslf,wgtins,&
         estype,boxshp,inscnd,inscfg,hostspec,ljformat,&
         inptemp,temp,&
         engdiv,maxins,&
         lwreg,upreg,&
         intprm,elecut,lwljcut,upljcut,&
         cmbrule,cltype,screen,ewtoler,splodr,plmode,&
         ew1max,ew2max,ew3max,ms1max,ms2max,ms3max,&
         block_threshold,&
         force_calculation
    use OUTname, only: OUTintprm,&                            ! from outside
         OUTens,OUTbxs,OUTtemp,&                ! from outside
         OUTelc,OUTlwl,OUTupl,&                 ! from outside
         OUTcmb,OUTclt,OUTscr,OUTspo,&          ! from outside
         OUTew1,OUTew2,OUTew3,&                 ! from outside
         OUTms1,OUTms2,OUTms3                   ! from outside
    use mpiproc                                                      ! MPI
    real, parameter :: tiny=1.0e-20
    real :: real_seed
    character*3 scrtype
    call mpi_info                                                    ! MPI

    intprm=1                    ! trajectory reading
    if(intprm.eq.1) then        ! default settings for trajectory reading
       estype=2                                     ! constant pressure
       boxshp=1                                     ! periodic boundary
       cltype=2 ; ms1max=64                         ! PME
       inptemp=300.0e0                              ! Kelvin
       engdiv=1                                     ! number of divisions
       screen=0.0e0 ; ewtoler=0.0e0
       ewtoler=1.0e-6 ; elecut=12.0e0
       splodr=4 ; scrtype='dis'
       upljcut=elecut ; lwljcut=upljcut-2.0e0
    endif
    block_threshold = 4.0 ! block-wise calculation, corresponds to 13 atoms / box
    force_calculation = .false.
    ! only part of constants set here
    call init_params()

    ! default settings
    skpcnf=1                    ! no skip for trajectory reading
    
    if(slttype.eq.1) corrcal=0
    if(slttype.ge.2) corrcal=1
    wgtslf=0 ; wgtins=0
    if((slttype.ge.2).and.(cltype.ne.0)) wgtslf=1  ! Ewald and PME
    sltpick=0 ; refpick=0 ; hostspec=1 ; ljformat=1
    maxins=1000 ; inscnd=0 ; inscfg=0 ; lwreg=0.0e0 ; upreg=5.0e0
    if(intprm.ne.0) then
       cmbrule=0
       ew2max=ew1max ; ew3max=ew1max
       ms2max=ms1max ; ms3max=ms1max
    endif
    ! default settings done

    ! read again for non-default constants
    call init_params()

    ! initialize random seed
    if(iseed == 0) then
       CALL RANDOM_SEED
       CALL RANDOM_NUMBER(real_seed)
       iseed = 100000 + int(899999 * real_seed)
    endif
    
    temp=inptemp*8.314510e-3/4.184e0               ! kcal/mol
    if((screen.le.tiny).and.(cltype.ne.0)) then    ! Ewald and PME
       if(ewtoler.le.tiny) call set_stop('ewa')
       screen=getscrn(ewtoler,elecut,scrtype)
    endif
    if(cltype.eq.1) then                           ! Ewald
       if(ew1max*ew2max*ew3max.eq.0) call set_stop('ewa')
    endif
    if(cltype.eq.2) then                           ! PME
       if(ms1max*ms2max*ms3max.eq.0) call set_stop('ewa')
    endif
    if((slttype.lt.1).or.(slttype.gt.3)) then      ! check slttype parameter
       call set_stop('slt')
    endif
    if(boxshp.eq.0) then        ! check the consistency in parameters
       if((estype.eq.2).or.(cltype.ne.0)) call set_stop('prs')
    endif
    if(wgtins.eq.1) then        ! check the consistency of insertion scheme
       if(slttype.ne.3) call set_stop('ins')
    endif
    if(slttype.eq.1) inscfg=0   ! inscfg is effective only for insertion
    if((inscfg.lt.0).or.(inscfg.gt.2)) then        ! check inscfg parameter
       call set_stop('ins')
    endif
    if(inscfg.ne.0) then        ! check the consistency against inscfg
       if(slttype.eq.1) call set_stop('prs')
       if(inscnd.eq.3) call set_stop('ins')
    endif
    plmode=2                    ! energy calculation parallelization mode

    return
  end subroutine iniparam

  real function getscrn(ewtoler,elecut,scrtype)
    character*3 scrtype
    real ewtoler,elecut,ewasml,ewalrg,scrfac,factor
    real, parameter :: error=1.0e-20
    factor=error+1.0e0 ; ewasml=0.0e0 ; ewalrg=1.0e3
    do while(factor.gt.error)
       scrfac=(ewasml+ewalrg)/2.0e0
       factor=erfc(scrfac*elecut)
       if(scrtype.eq.'dis') factor=factor/elecut
       if(factor.gt.ewtoler) ewasml=scrfac
       if(factor.le.ewtoler) ewalrg=scrfac
       factor=abs(factor-ewtoler)
    end do
    getscrn=scrfac
    return
  end function getscrn

  ! Calls OUTinitial / iniparam / OUTrename, and sets parameters
  subroutine setparam
    use engmain, only: numtype,nummol,maxsite,numatm,maxcnf,&
         slttype,sltpick,refpick,inscfg,ljformat,&
         moltype,numsite,sluvid,refmlid,&
         bfcoord,sitemass,charge,ljene,ljlen,&
         specatm,sitepos, mol_begin_index, belong_to, mol_charge,&
         CAL_SOLN, CAL_REFS_RIGID, CAL_REFS_FLEX
    use OUTname, only: OUTinitial,OUTrename,&                 ! from outside
         OUTntype,OUTnmol,OUTsite,OUTnrun,&     ! from outside
         OUTstmass,OUTcharge,OUTljene,OUTljlen  ! from outside
    use mpiproc, only: halt_with_error
    implicit none
    integer, parameter :: large=1000000
    ! only integer power is allowed as the initialization expression (7.1.6.1)
    real, parameter :: sgmcnv=1.7817974362806784e0 ! from Rmin/2 to sigma, 2.0**(5.0/6.0)
    real, parameter :: lencnv=1.0e1                ! from nm to Angstrom
    real, parameter :: engcnv=1.0e0/4.184e0        ! from kJ/mol to kcal/mol
    integer pti,stmax,uvtype,rftype,cmin,cmax,sid,i,ati,m
    integer :: solute_index
    real factor,xst(3)
    integer, dimension(:), allocatable :: pttype,ptcnt
    real, dimension(:,:), allocatable :: psite
    character*8 atmtype
    character*7 molfile
    character(*), parameter :: sltfile='SltInfo'
    character(*), parameter :: prmfile='MolPrm'
    character(*), parameter :: numbers='123456789'
    integer, parameter :: sltio=71                 ! IO for sltfile
    integer, parameter :: molio=72                 ! IO for molfile
    integer, parameter :: PT_PHYSICAL = 0, PT_TEST = 1

    call OUTinitial                ! initialization of OUTname module
    call iniparam                  ! initialization of parameters
    call OUTrename                 ! matching with outside variables

    maxcnf=OUTnrun                                     ! from outside
    if(slttype == CAL_SOLN) then
       numtype=OUTntype                  ! from outside
    else
       numtype=OUTntype+1                ! from outside
    end if

    ! pttype is particle type for each molecule group
    allocate( pttype(numtype),ptcnt(numtype) )
    pttype(:) = PT_PHYSICAL                  ! Default is physical particle, i.e. coordinate exists in (HISTORY) trajectory
    ptcnt(1:OUTntype) = OUTnmol(1:OUTntype)  ! Default: same as written in MDinfo
    if(slttype /= CAL_SOLN) then
       pttype(numtype) = PT_TEST             ! Last particle is test particle (for insertion)
       ptcnt(numtype) = 1                    ! Test particle can only be one molecule
    endif

    nummol = sum(ptcnt(1:numtype))

    allocate( moltype(nummol),numsite(nummol) )
    allocate( sluvid(nummol),refmlid(nummol) )

    ! make mapping from molecule no. [1..nummol] to particle type [1..numtype]
    cmin=1
    cmax=0
    do pti = 1, numtype
       cmax=cmax+ptcnt(pti)
       moltype(cmin:cmax) = pti ! sequential identification
       cmin=cmax+1
    end do
    if(cmax.ne.nummol) call set_stop('num')

    ! Determine which molecule is solute?
    if(slttype == CAL_SOLN) then
       solute_index = 1                   ! default
       if((1.le.sltpick).and.(sltpick.le.numtype)) solute_index = sltpick
    else
       solute_index = numtype             ! default
    endif
    if((slttype /= CAL_SOLN).and.(refpick == numtype)) call set_stop('ref')

    do i = 1, nummol
       pti = moltype(i)
       if(pttype(pti) == PT_PHYSICAL) stmax=OUTsite(pti)      ! from outside
       if(pttype(pti) == PT_TEST) then                        ! read from file
          open(unit=sltio,file=sltfile,status='old')
          stmax=0
          do sid=1,large
             if(slttype == CAL_REFS_RIGID) read(sltio,*,END=1219) m,atmtype,&
                  (xst(m), m=1,3),(xst(m), m=1,3)
             if(slttype == CAL_REFS_FLEX)  read(sltio,*,END=1219) m,atmtype,&
                  (xst(m), m=1,3)
             stmax=stmax+1
          end do
1219      continue
          close(sltio)
       endif

       if(pti.ne.solute_index) then                  ! solvent
          uvtype=0
          rftype=0                             ! default
          if(pti.eq.refpick) rftype=1          ! superposition reference
       endif
       if(pti.eq.solute_index) then                  ! solute
          uvtype=slttype
          rftype=2
       endif
       numsite(i)=stmax
       sluvid(i)=uvtype
       refmlid(i)=rftype
    end do
    deallocate( pttype,ptcnt )

    ! set max and total no. of atoms 
    maxsite = maxval(numsite(1:nummol))
    numatm = sum(numsite(1:nummol))

    allocate( bfcoord(3,maxsite),sitemass(numatm) )
    allocate( charge(numatm),ljene(numatm),ljlen(numatm) )
    allocate(sitepos(3,numatm))
    allocate(mol_begin_index(nummol + 1))
    allocate(mol_charge(nummol))
    allocate(belong_to(numatm))

    ! initial setting to zero
    bfcoord(1:3,1:maxsite) = 0.0e0
    sitemass(1:numatm) = 0.0e0
    charge(1:numatm) = 0.0e0
    ljene(1:numatm) = 0.0e0
    ljlen(1:numatm) = 0.0e0

    ! initialize mol_begin_index
    ! mol_begin_index(i) .. (mol_begin_index(i+1) - 1) will be the index range for i-th molecule
    mol_begin_index(1) = 1
    do i = 1, nummol
       mol_begin_index(i + 1) = mol_begin_index(i) + numsite(i)
    end do
    if(mol_begin_index(nummol + 1) /= numatm + 1) call halt_with_error("bug")

    ! initialize belong_to(map from atom number to molecule no)
    do i = 1, nummol
       belong_to(mol_begin_index(i):(mol_begin_index(i + 1) - 1)) = i
    end do

    allocate( psite(3,maxsite) )         ! temporary set of coordinates
    do i=1,nummol
       uvtype=sluvid(i)
       stmax=numsite(i)

       if(uvtype.eq.0) then                       ! solvent
          pti=moltype(i)
          m=pti
          do sid=1,nummol
             if((sluvid(sid).ge.1).and.(moltype(sid).lt.pti)) m=m-1
          end do
          molfile=prmfile//numbers(m:m)
       endif
       if(uvtype.ge.1) molfile=sltfile            ! solute
       open(unit=molio,file=molfile,status='old')
       do sid=1,stmax
          if(uvtype.eq.2) read(molio,*) m,atmtype,(xst(m), m=1,3),&
               (psite(m,sid), m=1,3)
          if(uvtype.ne.2) read(molio,*) m,atmtype,(xst(m), m=1,3)
          call getmass(factor,atmtype)
          ati=specatm(sid,i)
          sitemass(ati)=factor
          charge(ati)=xst(1)
          if(ljformat.eq.1) xst(3)=sgmcnv*xst(3)
          if((ljformat.eq.3).or.(ljformat.eq.4).and.(xst(3).ne.0.0)) then
             factor=(xst(2)/xst(3))**(1.0e0/6.0e0)
             xst(2)=xst(3)/(4.0e0*(factor**6.0e0))
             xst(3)=factor
          endif
          if((ljformat.eq.2).or.(ljformat.eq.4)) then
             xst(2)=engcnv*xst(2)
             xst(3)=lencnv*xst(3)
          endif
          ljene(ati)=xst(2)
          ljlen(ati)=xst(3)
       end do
       close(molio)

       if(uvtype.eq.2) then        
          if(inscfg.eq.2) bfcoord(1:3,1:stmax)=psite(1:3,1:stmax)
          if(inscfg.ne.2) then 
             ! setting the center of mass to zero
             call molcen(i,psite,xst,'com')
             do sid=1,stmax
                bfcoord(1:3,sid)=psite(1:3,sid)-xst(1:3)
             end do
          endif
       endif
    end do
    deallocate( psite )

    ! conversion to (kcal/mol angstrom)^(1/2)
    ! == sqrt(e^2 * coulomb const * avogadro / (kcal / mol angstrom))
    charge(1:numatm)=18.22261721e0 * charge(1:numatm)

    ! get molecule-wise charges
    do i = 1, nummol
       mol_charge(i) = sum(charge(mol_begin_index(i):(mol_begin_index(i+1)-1)))
    end do

    ! set L-J epsilon value to be square-rooted, because these values are only used in sqrt-form.
    ljene(:) = sqrt(ljene(:))

  end subroutine setparam


! Reads up to maxread coordinates serially and distributes frames
! returns number of frames read (EXCLUDING skipped frames)
  subroutine getconf_parallel(maxread, actual_read)
    use engmain, only: nummol,numatm,boxshp,&
         numsite,sluvid,sitepos,cell,skpcnf
    use OUTname, only: OUTconfig                     ! from outside
    use mpiproc
    implicit none
    integer, intent(in) :: maxread
    integer, intent(out) :: actual_read
    
    real, dimension(:,:), allocatable :: OUTpos,OUTcell
    integer i,m,k,OUTatm, iproc, nread
    character*3 rdconf
    
#ifdef trjctry
    rdconf='trj'                                     ! reading from file
#else
    rdconf='fly'                                     ! on-the-fly reading
#endif
    OUTatm=0
    do i=1,nummol
       if(sluvid(i).le.1) OUTatm=OUTatm+numsite(i) ! only solvent & solute; no test particles
    end do
    allocate( OUTpos(3,OUTatm),OUTcell(3,3) )
    
    nread = min(nprocs, maxread)

    if(myrank == 0) then
       do iproc = 1, nread
          ! get configuration and store in OUTpos / OUTcell
          do i = 1, skpcnf
             call OUTconfig(OUTpos,OUTcell,OUTatm,boxshp,0,rdconf)
          end do
          
          if(iproc == 1) then ! copy to self
             sitepos(:,1:OUTatm)=OUTpos(:,:)
             cell(:,:)=OUTcell(:,:)
          else
#ifndef noMPI
             ! FIXME: rewrite with grouping and scatter?
             call mpi_send(OUTpos, 3 * OUTatm, &
                  mpi_double_precision, iproc - 1, tag_coord, mpi_comm_world, ierror)
             call mpi_send(OUTcell, 3 * 3, &
                  mpi_double_precision, iproc - 1, tag_cell, mpi_comm_world, ierror)
#endif
          endif
       end do
    elseif(myrank < nread) then
#ifndef noMPI
       call mpi_recv(sitepos, 3 * OUTatm, &
            mpi_double_precision, 0, tag_coord, mpi_comm_world, mpistatus, ierror)
       call mpi_recv(cell, 3 * 3, &
            mpi_double_precision, 0, tag_cell, mpi_comm_world, mpistatus, ierror)
#endif
    endif

    deallocate( OUTpos,OUTcell )
    actual_read = nread
  end subroutine getconf_parallel

      subroutine set_stop(type)
      use engmain, only: io6
      use mpiproc                                                      ! MPI
      character*3 type
      if(type.eq.'slt') write(io6,991)
      if(type.eq.'num') write(io6,992)
      if(type.eq.'ref') write(io6,993)
      if(type.eq.'prs') write(io6,994)
      if(type.eq.'ins') write(io6,995)
      if(type.eq.'ewa') write(io6,996)
991   format(' The solute type is incorrectly set')
992   format(' The number of molecules or atoms is incorrectly set')
993   format(' The reference structure of solvent is incorrectly set')
994   format(' The system parameters are incorrectly set')
995   format(' The insertion parameters are incorrectly set')
996   format(' The Ewald parameters are incorrectly set')
      call mpi_setup('stop')                                           ! MPI
      stop
    end subroutine set_stop
!
!
    subroutine getmass(stmass,atmtype)
!
      real stmass
      real massH,massC,massO,massN,massS,massP
      real massNa,massCa,massZn,massFe,massCu,massF,massCl,massBr
      character*1 eltp1
      character*2 eltp2
      character*3 eltp3
      character*5 atmtype
!
      massH=1.00794e0             ! mass number (hydrogen)
      massC=12.0107e0             ! mass number (carbon)
      massO=15.9994e0             ! mass number (oxygen)
      massN=14.00674e0            ! mass number (nitrogen)
      massS=32.066e0              ! mass number (sulfur)
      massP=30.973761e0           ! mass number (phosphorus)
      massNa=22.989770e0          ! mass number (sodium)
      massCa=40.078e0             ! mass number (calcium)
      massZn=65.409e0             ! mass number (zinc)
      massFe=55.845e0             ! mass number (iron)
      massCu=63.546e0             ! mass number (copper)
      massF=18.9984032e0          ! mass number (fluorine)
      massCl=35.4527e0            ! mass number (chlorine)
      massBr=79.904e0             ! mass number (bromine)
!
      eltp1=atmtype(1:1)
      if(eltp1.eq.'H') stmass=massH
      if(eltp1.eq.'C') stmass=massC
      if(eltp1.eq.'O') stmass=massO
      if(eltp1.eq.'N') stmass=massN
      if(eltp1.eq.'S') stmass=massS
      if(eltp1.eq.'F') stmass=massF
      if(eltp1.eq.'P') stmass=massP
      eltp2=atmtype(1:2)
      if(eltp2.eq.'Na') stmass=massNa
      if(eltp2.eq.'Ca') stmass=massCa
      if(eltp2.eq.'Zn') stmass=massZn
      if(eltp2.eq.'Fe') stmass=massFe
      if(eltp2.eq.'Cu') stmass=massCu
      if(eltp2.eq.'Cl') stmass=massCl
      if(eltp2.eq.'Br') stmass=massBr
      if(eltp2.eq.'CH') stmass=massC+massH
      eltp3=atmtype(1:3)
      if(eltp3.eq.'CH2') stmass=massC+2.0e0*massH
      if(eltp3.eq.'CH3') stmass=massC+3.0e0*massH
!
      return
    end subroutine getmass
!
!
    subroutine molcen(i,psite,cen,caltype)       ! getting molecular center
      use engmain, only: nummol,maxsite,numatm,numsite,sitemass,specatm
      integer i,ati,sid,stmax,m
      real psite(3,maxsite),cen(3),wgt,sitm
      character*3 caltype
      wgt=0.0e0
      do 3331 m=1,3
        cen(m)=0.0e0
3331  continue
      stmax=numsite(i)
      do 3332 sid=1,stmax
        ati=specatm(sid,i)
        if(caltype.eq.'cnt') sitm=1.0e0            ! centroid
        if(caltype.eq.'com') sitm=sitemass(ati)    ! center of mass
        wgt=wgt+sitm
        do 3333 m=1,3
          cen(m)=cen(m)+sitm*psite(m,sid)
3333    continue
3332  continue
      do 3334 m=1,3
        cen(m)=cen(m)/wgt
3334  continue
      return
    end subroutine
!
end module
