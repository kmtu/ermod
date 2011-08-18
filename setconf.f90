module OUTglobal
  implicit none
  !  renaming global variables in package program
  integer GLBnrun,GLBntype,GLBens,GLBpbc,GLBclt,GLBspl
  real GLBtmp,GLBelc,GLBvdc,GLBsrc
  integer GLBew1,GLBew2,GLBew3
end module OUTglobal
!
!
module OUTname
!  renaming outside parameters and parameters to avoid conflict
  use OUTglobal
  implicit none
  integer, parameter :: iotrj=99                 ! trajectory file IO
  character(*), parameter :: trjfile='HISTORY'   ! trajectory filename
  integer, parameter :: cltrj=98                 ! cell trajectory file IO
  character(*), parameter :: celfile='HISTCELL'  ! cell trajectory filename
  integer, parameter :: mdinf=89                 ! MD info file IO
  character(*), parameter :: inffile='MDinfo'    ! MD info filename

  character*3, save :: iofmt      ! formatted or not for trajectory file
  character*3, save :: cltrd      ! seperate file for cell trajectory file
  character*3, save :: mdird      ! seperate file for MD info
  character*3, save :: bxiso      ! PBC adopted for single-molecule MD
  character*3, save :: toptp      ! type of first variable in each line
  integer, save :: skpio          ! number of first skipped lines
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

#ifdef GROMACS
  character(len=80) :: buffer
#ifdef MDLIB
  integer(8) :: gmxhandle
#endif
#endif
#ifdef VMDPLUGINS
  integer(8) :: vmdhisthandle
#endif

contains

  subroutine OUTinitial
    iofmt='yes' ; cltrd='not' ; mdird='yes'      ! default
    bxiso='not' ; toptp='int' ; skpio=0          ! default
    sgmcnv=2.0e0**(5.0e0/6.0e0)                  ! from Rmin/2 to sigma
    chgcnv=1.0e0/1.60217653e-19                  ! from C to elementary
    engcnv=6.0221415e23/4.184e3                  ! from J to kcal/mol
    lencnv=1.0e10                                ! from meter to Angstrom
#ifdef GROMACS
    bxiso='yes' ; toptp='chr'                                 ! GROMACS
    lencnv=1.0e1                                 ! from nm to Angstrom
#endif
#ifdef NAMD
    iofmt='not' ; cltrd='yes' ; toptp='rsg' ; skpio=3         ! NAMD
#endif
    return
  end subroutine OUTinitial


  subroutine opentrj
#ifdef VMDPLUGINS
    integer :: status
    external vmdfio_open_traj
#endif

#if defined(GROMACS) && defined(MDLIB)
!     GROMACS + MDLIB
    integer :: status
    external open_gmtraj
#endif
    use_mdlib = .false.
#if defined(GROMACS) && defined(MDLIB)
!     GROMACS + MDLIB
    call open_gmtraj(gmxhandle, 0, status) ! 0 == HISTORY
    if(status == 0) then
       use_mdlib = .true.
    endif
#endif

#ifdef VMDPLUGINS
    call vmdfio_open_traj(vmdhisthandle, trjfile, len_trim(trjfile), status)
    if(status /= 0) then
       stop "vmdfio_open_traj_: unable to open trajectory. HISTORY must be a symlink"
    endif
    use_mdlib = .true.
#endif
    if(.not.(use_mdlib)) then
       if(iofmt.eq.'yes') open(unit=iotrj,file=trjfile,status='old')
       if(iofmt.eq.'not') open(unit=iotrj,file=trjfile,status='old',&
            form='unformatted')
    endif
    if(cltrd.eq.'yes') open(unit=cltrj,file=celfile,status='old')
  end subroutine opentrj

  subroutine initconf
    call OUTinitial

    if(mdird.eq.'yes') open(unit=mdinf,file=inffile,status='old')
    return
  end subroutine initconf

  subroutine closetrj
#if defined(GROMACS) && defined(MDLIB)
!     GROMACS + MDLIB
    external close_gmtraj
    if(use_mdlib) then
       call close_gmtraj(gmxhandle)
    endif
#endif
#ifdef VMDPLUGINS
    external vmdfio_close_traj
    call vmdfio_close_traj(vmdhisthandle)
    use_mdlib = .true.
#endif
    if(.not.use_mdlib) then
       close(iotrj)                                 ! trajectory file
    endif
    if(cltrd.eq.'yes') close(cltrj)                ! cell trajectory file
    use_mdlib = .false.                            ! Clear for next
    return
  end subroutine closetrj

  subroutine finiconf
    if(mdird.eq.'yes') close(mdinf)                ! MD info
  end subroutine finiconf
!
  ! Initialization 2nd phase - read topologies from mother MD program, or read MDinfo
  subroutine OUTrename
    integer i,pti,sid,ctm
    real trjene,trjlen
!
#ifdef trjctry
    if(mdird.eq.'yes') read(mdinf,*) OUTnrun,OUTntype
    if(mdird.eq.'not') read(iotrj,*) OUTnrun,OUTntype         ! CHARMM
    if(.not.use_mdlib) then
       call OUTskip(iotrj,iofmt,skpio)
    endif
    if(cltrd.eq.'yes') call OUTskip(cltrj,'yes',3)            ! NAMD
#endif
!
    allocate( OUTnmol(OUTntype),OUTsite(OUTntype) )
!
#ifdef trjctry
    if(mdird.eq.'yes') then
       read(mdinf,*) (OUTnmol(i), i=1,OUTntype)
       read(mdinf,*) (OUTsite(i), i=1,OUTntype)
    endif
#endif
!
    TotAtm=0
    do i=1,OUTntype
       TotAtm=TotAtm+OUTnmol(i)*OUTsite(i)
    end do
    allocate( OUTstmass(TotAtm),OUTcharge(TotAtm) )
    allocate( OUTljene(TotAtm),OUTljlen(TotAtm) )
!
    return
  end subroutine OUTrename


  ! read system setup (coulomb rule, LJ, etc..)
  subroutine OUTintprm
    return
  end subroutine OUTintprm

  ! Get molecular configuration
  subroutine OUTconfig(OUTpos,OUTcell,OUTatm,OUTbox,OUTtrj,rdconf)
    implicit none
    integer OUTatm,OUTbox,OUTtrj,trjID,i,m,k
    real OUTpos(3,OUTatm),OUTcell(3,3),xst(9),factor
    character*3 rdconf
    character*8 dumchr
#ifdef NAMD
    real*4, dimension(:), allocatable :: snglcrd   ! used to read DCD file
#endif
#ifdef VMDPLUGINS
    integer :: vmdstatus
      
    external vmdfio_read_traj_step
#endif
!
    if(OUTtrj.eq.0) trjID=iotrj
    if(OUTtrj.ne.0) trjID=OUTtrj
    OUTcell(:, :) = 0.0e0

!
      if(rdconf.eq.'trj') then
#ifdef VMDPLUGINS
        if(trjID == iotrj) then
          call vmdfio_read_traj_step(vmdhisthandle, OUTpos, OUTcell, OUTatm, vmdstatus)
          return
        endif
! for sltconf keep traditional I/O
#endif
#ifdef GROMACS
        buffer = ""
        do while((buffer /= "POSITIONRED").and.(buffer /= "POSITION"))! GROMACS
          read(trjID,*) buffer                                  ! GROMACS
          buffer = trim(buffer)                                 ! GROMACS
        enddo                                                   ! GROMACS
        do 7711 i=1,OUTatm                                      ! GROMACS
          if(buffer == "POSITION") then                         ! GROMACS
            read(trjID,*) k,(dumchr, m=1,2),k,(xst(m), m=1,3)   ! GROMACS
          else                                                  ! GROMACS
            read(trjID,*) (xst(m), m=1,3)                       ! GROMACS
          endif                                                 ! GROMACS
          do 7712 m=1,3                                         ! GROMACS
            OUTpos(m,i)=lencnv*xst(m)                           ! GROMACS
7712      continue                                              ! GROMACS
7711    continue                                                ! GROMACS
        call OUTskip(trjID,iofmt,1)                             ! GROMACS
        if(OUTbox.ne.0) then                                    ! GROMACS
          do while(trim(buffer) /= "BOX")                       ! GROMACS
            read(trjID,*) buffer                                ! GROMACS
          enddo                                                 ! GROMACS
          read(trjID,*) (xst(m), m=1,3)                         ! GROMACS
          do 7713 m=1,3                                         ! GROMACS
            OUTcell(m,m)=lencnv*xst(m)                          ! GROMACS
7713      continue                                              ! GROMACS
          call OUTskip(trjID,iofmt,1)                           ! GROMACS
        endif                                                   ! GROMACS
#endif
#ifdef NAMD
        if(OUTbox.ne.0) then                                    ! NAMD
          read(cltrj,*) i,(xst(m), m=1,9)                       ! NAMD
          do 7711 k=1,3                                         ! NAMD
            do 7712 m=1,3                                       ! NAMD
              OUTcell(m,k)=xst(3*(k-1)+m)                       ! NAMD
7712        continue                                            ! NAMD
7711      continue                                              ! NAMD
          call OUTskip(trjID,iofmt,1)                           ! NAMD
        endif                                                   ! NAMD
        allocate( snglcrd(OUTatm) )                             ! NAMD
        do 7713 m=1,3                                           ! NAMD
          read(trjID) (snglcrd(i), i=1,OUTatm)                  ! NAMD
          do 7714 i=1,OUTatm                                    ! NAMD
            OUTpos(m,i)=dble(snglcrd(i))                        ! NAMD
7714      continue                                              ! NAMD
7713    continue                                                ! NAMD
        deallocate( snglcrd )                                   ! NAMD
#endif
      endif
!
      return
  end subroutine
!
! Utility function to skip trajectories (FIXME: what's the call from insertion?)
  subroutine OUTskip(fileio,ioform,skpmax)
    use mpiproc, only: myrank
    integer fileio,skpcnt,skpmax
    character*3 ioform
    if(skpmax.gt.0 .and. myrank == 0) then
       do skpcnt=1,skpmax
          if(ioform.eq.'yes') read(fileio,*)
          if(ioform.eq.'not') read(fileio)
       end do
    endif
    return
  end subroutine OUTskip
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
!
#ifndef trjctry
      intprm=0                    ! combined with parent MD program
#else
      intprm=1                    ! trajectory reading
#endif
      if(intprm.eq.0) then        ! from parent MD setup
        call OUTintprm
        estype=OUTens ; boxshp=OUTbxs ; inptemp=OUTtemp
        elecut=OUTelc ; lwljcut=OUTlwl ; upljcut=OUTupl
        cmbrule=OUTcmb ; cltype=OUTclt ; screen=OUTscr ; splodr=OUTspo
        ew1max=OUTew1 ; ew2max=OUTew2 ; ew3max=OUTew3
        ms1max=OUTms1 ; ms2max=OUTms2 ; ms3max=OUTms3
      endif
      if(intprm.eq.1) then        ! default settings for trajectory reading
        estype=2                                     ! constant pressure
        boxshp=1                                     ! periodic boundary
        cltype=2 ; ms1max=64                         ! PME
        inptemp=300.0e0                              ! Kelvin
        engdiv=1                                     ! number of divisions
        screen=0.0e0 ; ewtoler=0.0e0
#ifdef GROMACS
        ewtoler=1.0e-5 ; elecut=10.0e0                          ! GROMACS
        splodr=6 ; scrtype='rel'                                ! GROMACS
        upljcut=14.0e0 ; lwljcut=upljcut                        ! GROMACS
#endif
#ifdef NAMD
        ewtoler=1.0e-6 ; elecut=12.0e0                          ! NAMD
        splodr=4 ; scrtype='dis'                                ! NAMD
        upljcut=elecut ; lwljcut=upljcut-2.0e0                  ! NAMD
#endif
      endif
      block_threshold = 4.0 ! block-wise calculation
      force_calculation = .false.
!     only part of constants set here
      call init_params()
!
!  default settings
#ifndef trjctry
      if(slttype.eq.1) skpcnf=5
      if(slttype.ge.2) skpcnf=50
#else
      skpcnf=1                    ! no skip for trajectory reading
#endif
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
!  default settings done
!
!     read again for non-default constants
      call init_params()

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
!
      return
      end subroutine
!
!
      real function getscrn(ewtoler,elecut,scrtype)
      character*3 scrtype
      real ewtoler,elecut,ewasml,ewalrg,scrfac,factor
      real, parameter :: error=1.0e-20
      factor=error+1.0e0 ; ewasml=0.0e0 ; ewalrg=1.0e3
      do 5511 while(factor.gt.error)
        scrfac=(ewasml+ewalrg)/2.0e0
        factor=erfc(scrfac*elecut)
        if(scrtype.eq.'dis') factor=factor/elecut
        if(factor.gt.ewtoler) ewasml=scrfac
        if(factor.le.ewtoler) ewalrg=scrfac
        factor=abs(factor-ewtoler)
5511  continue
      getscrn=scrfac
      return
    end function
!
!
    ! Calls OUTinitial / iniparam / OUTrename, and sets parameters
    subroutine setparam
      use engmain, only: numtype,nummol,maxsite,numatm,maxcnf,&
                         slttype,sltpick,refpick,inscfg,ljformat,&
                         moltype,numsite,sluvid,refmlid,&
                         bfcoord,sitemass,charge,ljene,ljlen,&
                         specatm,sitepos, mol_begin_index, belong_to, mol_charge
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
      real factor,xst(3)
      integer, dimension(:), allocatable :: pttype,ptcnt
      real, dimension(:,:), allocatable :: psite
      character*8 prmread,atmtype
      character*7 molfile
      character(*), parameter :: sltfile='SltInfo'
      character(*), parameter :: prmfile='MolPrm'
      character(*), parameter :: numbers='123456789'
      integer, parameter :: sltio=71                 ! IO for sltfile
      integer, parameter :: molio=72                 ! IO for molfile
!
      call OUTinitial                ! initialization of OUTname module
      call iniparam                  ! initialization of parameters
      call OUTrename                 ! matching with outside variables
      maxcnf=OUTnrun                                     ! from outside
      if(slttype.eq.1) numtype=OUTntype                  ! from outside
      if(slttype.ge.2) numtype=OUTntype+1                ! from outside
!
      allocate( pttype(numtype),ptcnt(numtype) )
      do 1101 pti=1,numtype
        ati=9
        if((slttype.eq.1).or.(pti.lt.numtype))  ati=0  ! physical particle
        if((slttype.ge.2).and.(pti.eq.numtype)) ati=1  ! test particle
        if((ati.ne.0).and.(ati.ne.1)) call set_stop('slt')
        pttype(pti)=ati
        if(pti.le.OUTntype) ptcnt(pti)=OUTnmol(pti)      ! from outside
        if(pti.gt.OUTntype) ptcnt(pti)=1
1101  continue
!
      nummol=0
      do 1102 pti=1,numtype
        nummol=nummol+ptcnt(pti)
1102  continue
!
      allocate( moltype(nummol),numsite(nummol) )
      allocate( sluvid(nummol),refmlid(nummol) )
!
      cmin=1
      cmax=0
      do 1201 pti=1,numtype
        cmax=cmax+ptcnt(pti)
        do 1202 i=cmin,cmax                ! sequential identification
          moltype(i)=pti
1202    continue
        cmin=cmax+1
1201  continue
      if(cmax.ne.nummol) call set_stop('num')
!
      do 1251 i=1,nummol
        pti=moltype(i)
        if(pttype(pti).eq.0) stmax=OUTsite(pti)          ! from outside
        if(pttype(pti).eq.1) then                        ! read from file
          open(unit=sltio,file=sltfile,status='old')
          stmax=0
          do 1211 sid=1,large
            if(slttype.eq.2) read(sltio,*,END=1219) m,atmtype,&
                                       (xst(m), m=1,3),(xst(m), m=1,3)
            if(slttype.eq.3) read(sltio,*,END=1219) m,atmtype,&
                                                       (xst(m), m=1,3)
            stmax=stmax+1
1211      continue
1219      continue
          close(sltio)
        endif
        if(slttype.eq.1) then
          m=1                                  ! default
          if((1.le.sltpick).and.(sltpick.le.numtype)) m=sltpick
        endif
        if(slttype.ge.2) m=numtype             ! default
        if((slttype.ge.2).and.(refpick.eq.numtype)) call set_stop('ref')
        if(pti.ne.m) then                  ! solvent
          uvtype=0
          rftype=0                             ! default
          if(pti.eq.refpick) rftype=1          ! superposition reference
        endif
        if(pti.eq.m) then                  ! solute
          uvtype=slttype
          rftype=2
        endif
        numsite(i)=stmax
        sluvid(i)=uvtype
        refmlid(i)=rftype
1251  continue
      deallocate( pttype,ptcnt )
!
      maxsite=0
      numatm=0
      do 1281 i=1,nummol
        stmax=numsite(i)
        if(stmax.gt.maxsite) maxsite=stmax
        numatm=numatm+stmax
1281  continue
!
      allocate( bfcoord(3,maxsite),sitemass(numatm) )
      allocate( charge(numatm),ljene(numatm),ljlen(numatm) )
      allocate(sitepos(3,numatm))
      allocate(mol_begin_index(nummol + 1))
      allocate(mol_charge(nummol))
      allocate(belong_to(numatm))

      do 7301 sid=1,maxsite                ! initial setting to zero
        do 7302 m=1,3
          bfcoord(m,sid)=0.0e0
7302    continue
7301  continue
      do 7305 ati=1,numatm
        sitemass(ati)=0.0e0
        charge(ati)=0.0e0
        ljene(ati)=0.0e0
        ljlen(ati)=0.0e0
7305  continue

      ! initialize mol_begin_index
      ! mol_begin_index(i) .. (mol_begin_index(i+1) - 1) will be the index range for i-th molecule
      mol_begin_index(1) = 1
      do i = 1, nummol
         mol_begin_index(i + 1) = mol_begin_index(i) + numsite(i)
      end do
      if(mol_begin_index(nummol + 1) /= numatm + 1) call halt_with_error("bug")

      ! initialize belong_to
      do i = 1, nummol
         belong_to(mol_begin_index(i):(mol_begin_index(i + 1) - 1)) = i
      end do

      ati=0                                ! specifying the site in molecule
      do 1501 i=1,nummol
        stmax=numsite(i)
        do 1512 sid=1,stmax
          ati=ati+1
1512    continue
1501  continue
      if(ati.ne.numatm) call set_stop('num')
!
      allocate( psite(3,maxsite) )         ! temporary set of coordinates
      do 5001 i=1,nummol
        uvtype=sluvid(i)
#ifndef trjctry
        if(uvtype.le.1) prmread='internal'           ! physical particle
        if(uvtype.ge.2) prmread='external'           ! test particle
#else
        prmread='external'
#endif
        stmax=numsite(i)
        if(prmread.eq.'internal') then     ! read through translation module
          do 5101 sid=1,stmax
            ati=specatm(sid,i)                       ! from outside
            sitemass(ati)=OUTstmass(ati)             ! from outside
            charge(ati)=OUTcharge(ati)               ! from outside
            ljene(ati)=OUTljene(ati)                 ! from outside
            ljlen(ati)=OUTljlen(ati)                 ! from outside
5101      continue
        endif
        if(prmread.eq.'external') then     ! read from file
          if(uvtype.eq.0) then                       ! solvent
            pti=moltype(i)
            m=pti
            do 5087 sid=1,nummol
              if((sluvid(sid).ge.1).and.(moltype(sid).lt.pti)) m=m-1
5087        continue
            molfile=prmfile//numbers(m:m)
          endif
          if(uvtype.ge.1) molfile=sltfile            ! solute
          open(unit=molio,file=molfile,status='old')
          do 7111 sid=1,stmax
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
7111      continue
          close(molio)
        endif
!
        if(uvtype.eq.2) then        ! setting the center of mass to zero
          if(inscfg.ne.2) call molcen(i,psite,xst,'com')
          do 5251 sid=1,stmax
            do 5252 m=1,3
              if(inscfg.ne.2) bfcoord(m,sid)=psite(m,sid)-xst(m)
              if(inscfg.eq.2) bfcoord(m,sid)=psite(m,sid)
5252        continue
5251      continue
        endif
5001  continue
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

      return
    end subroutine


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
