c  program enganal.f
c
c   angstrom and kcal/mol are taken as the units of length and energy.
c
c  names of parameters
c   numtype : number of molcular types
c   nummol : total number of molecules
c   maxsite : maximum number of interaction sites in molecule
c   numatm : total number of atomic sites
c   maxcnf : maximum number of configurations in MD
c   engdiv : number of divisions of the total simulation length
c   skpcnf : interval to skip the configuration examined
c   corrcal : calculation of the correlation matrix
c               0 : no calculation  1 : calculation performed
c   slttype : type of solute treatment
c               1 : physical
c               2 : test particle (rigid)  3 : test particle (flexible)
c   sltpick : specifying the solute species
c               1 <= sltpick <= numtype (default = 1) if slttype = 1
c               sltpick = numtype if slttype >= 2
c   refpick : specifying the mixed solvent species for superposition
c               1 <= refpick <= numtype    when slttype = 1
c               1 <= refpick <= numtype-1  when slttype >= 2
c   wgtslf : weighting by the self-energy  --- 0 : no  1 : yes
c   wgtins : weight of the solute intramolecular configuration
c               0 : no  1 : yes (can be = 1 only when slttype = 3)
c   estype : type of system
c               1 : constant volume  2 : constant pressure
c   boxshp : shape of the unit cell box
c               0 : non-periodic  1 : periodic and parallelepiped
c   inscnd : geometrical condition of the solute configuration
c               0 : random  1 : spherical  2 : slab  3 : reference
c   inscfg : position and orientation for the inserted solute
c               0 : only the intramolecular configuration is from the file
c               1 : orientation is fixed from the file with random position
c               2 : position and orientation are also fixed from the file
c            The file for the solute configuration
c            is SltInfo when slttype = 2 and is SltConf when slttype = 3
c              default = 0
c   hostspec : type of molecule forming host (micelle, membrane, or protein)
c              active only when inscnd is not equal to 0
c               1 <= hostspec <= numtype    when slttype = 1
c               1 <= hostspec <= numtype-1  when slttype >= 2
c   ljformat : input-file format for the LJ energy and length parameters
c               0 : epsilon (kcal/mol) and sigma (A)
c               1 : epsilon (kcal/mol) and Rmin/2 (A)
c               2 : epsilon (kJ/mol) and sigma (nm)
c               3 : A (kcal/mol A^12) and C (kcal/mol A^6)
c               4 : C12 (kJ/mol nm^12) and C6 (kJ/mol nm^6)
c              default = 1
c   iseed : seed parameter for uniform random number
c   inptemp : temperature of the system in Kelvin
c   temp  : temperature of the system in kcal/mol
c   io6   : standard output
c
c  names of constants and variables for trajectory generation
c   moltype : type of the molecule numbering between 1 and numtype
c   numsite : number of sites in a molecule
c   sluvid : solvent or solute
c               0 : solvent  1 : solute
c               2 : test particle (rigid)  3 : test particle (flexible)
c   refmlid : superposition reference among solvent species
c               0 : not reference
c               1 : reference solvent  2 : reference solute
c             default = 0 unless 1 <= refpick <= numtype for solvent species
c   bfcoord : coordinate of a rigid molecule in body-fixed frame
c             used only for test particle and kept NULL for physical particle
c   sitemass : mass of an interaction site in a molecule
c   charge : partial charges on the sites in a molecule
c   ljene : energy parameter for Lennard-Jones potential in a molecule
c   ljlen : length parameter for Lennard-Jones potential in a molecule
c   intprm : whether the intereaction paramters given below
c                    (from elecut to ms1max,ms2max,ms3max)
c                    and estype, boxshp, and inptemp
c                    are read from the parent MD program
c      default = 0 in the case of on-the-fly calculation
c      default = 1 in the case of trajectory reading
c   elecut : cutoff of the real-space part of electrostatic interaction
c   lwljcut : lower limit of the LJ cutoff tapering function
c   upljcut : upper limit of the LJ cutoff tapering function
c   cmbrule : combination rule for LJ interaction
c        0 : arithmetic mean is used for LJ sigma as for AMBER and CHARMM
c        1 : geometric mean is used for LJ sigma as for OPLS
c      default = 0
c      geometric mean is always used for LJ epsilon
c   cltype : treatment of Coulomb interaction   0 : bare  1 : Ewald  2 : PME
c   screen : screening constant in Ewald summation
c   ewtoler : Ewald and PME tolerance to calculate the screen parameter
c      when screen is given, screen has the priority
c   splodr : order of spline function used in PME
c   ew1max,ew2max,ew3max : number of reciprocal vectors along one direction
c   ms1max,ms2max,ms3max : number of meshes in PME along one direction
c   plmode : parallelization mode for calculation of solute-solvent interaction
c        0 : parallel over solvent molecules and serial over insertion
c        1 : serial over solvent molecules and parallel over insertion
c      default = 0 if slttype = 1 or ( slttype >= 2 and cltype = 0 or 1)
c              = 1 if slttype >= 2 and cltype = 2
c   specatm : specification of the site
c   sitepos : coordiantes of interaction site
c   cell : unit cell vector
c   invcl : inversion of the cell matrix
c   volume : system volume
c
c  names of constants and variables for energy distribution
c   ermax : size of the energy-represented distribution functions
c   numslv : number of solvent species
c   uvmax : number of discretization for each solvent species
c   uvsoft : number of discretization in soft interaction region
c   esmax : number of discretization of the solute self-energy
c   maxins : maximum number of insertions for test solute particle
c   uvspec : assignment to the number representing the solvent species
c   numslt : number of solute molecules
c   sltlist : list of solute molecules
c   engnorm : normalization factor
c   engsmpl : number of samplings
c   voffset : offset value for the self-energy
c   lwreg : lower bound for the region of solute position
c   upreg : upper bound for the region of solute position
c
c  constants defining the discretized energy coordinate
c     these appear only in the enginit subroutine
c     and are used for each of the solute and solvent species
c   peread : determines whether the parameters are read from a separate file
c       0 : parameters are read from param_eng
c       1 : parameters are read separately from a file EcdInfo
c     default = 0
c   pemax : number of discretization of the solute-solvent energy
c   pesoft : number of discretization in the soft interaction region
c   pecore : number of discretization in the core interaction region
c            pemax and pesoft are constructed from other parameters
c            pemax = pesoft + pecore
c   ecdmin : minimum value of the solute-solvent energy
c   ecfmns : smaller side of the finely dicretized solute-solvent energy
c   ecmns0 : smaller side of the very finely dicretized energy near ecdcen
c   ecdcen : central value of the energy coordinate, typically zero
c   ecpls0 : larger side of the very finely dicretized energy near ecdcen
c   ecfpls : larger side of the finely dicretized solute-solvent energy
c   eccore : the solute-solvent energy at which the dicretization is changed
c   ecdmax : maximum value of the solute-solvent energy
c   eclbin : linear mesh for the solute-solvent energy
c   ecfbin : fine linear mesh for the solute-solvent energy
c   ec0bin : very fine linear mesh for the solute-solvent energy near 0
c
c
      module engmain
c
      implicit none
c
      integer :: numtype,nummol,maxsite,numatm,maxcnf,engdiv,skpcnf,corrcal
      integer :: slttype,sltpick,refpick,wgtslf,wgtins
      integer :: estype,boxshp,inscnd,inscfg,hostspec,ljformat,iseed
      real :: inptemp,temp
      integer, parameter :: io6=6             ! standard output
c
      integer, dimension(:),   allocatable, save :: moltype,numsite
      integer, dimension(:),   allocatable, save :: sluvid,refmlid
      real, dimension(:,:),    allocatable, save :: bfcoord
      real, dimension(:),      allocatable, save :: sitemass
      real, dimension(:),      allocatable, save :: charge,ljene,ljlen
      integer, dimension(:,:), allocatable, save :: specatm
      real, dimension(:,:),    allocatable, save :: sitepos
      real, dimension(3,3),                 save :: cell,invcl
      real,    save :: volume
c
      real,    save :: elecut,lwljcut,upljcut,screen,ewtoler
      integer, save :: intprm,cmbrule,cltype,splodr,plmode
      integer, save :: ew1max,ew2max,ew3max,ms1max,ms2max,ms3max
c
      integer, save :: ermax,numslv,esmax,maxins
      integer, dimension(:),   allocatable, save :: uvmax,uvsoft,uvspec
      real, dimension(:),      allocatable, save :: uvcrd,edens
      real, dimension(:,:),    allocatable, save :: ecorr
      real, dimension(:),      allocatable, save :: escrd,eself
      real, dimension(:,:),    allocatable, save :: aveuv
      real, dimension(:),      allocatable, save :: slnuv
      real, dimension(:,:),    allocatable, save :: avediv
      real,    save :: avslf
      real, dimension(:),      allocatable, save :: minuv,maxuv
      integer, save :: numslt
      integer, dimension(:),   allocatable, save :: sltlist
      real,    save :: engnorm,engsmpl,voffset
c
      real,    save :: lwreg,upreg
c
      end module
c
c
      module mpiproc                                                   ! MPI
#ifndef noMPI
c                                                                      ! MPI
      include "mpif.h"
#endif
      integer ierror,myrank,nprocs                                     ! MPI
      contains
      subroutine mpi_setup(type)                                       ! MPI
      character*4 type                                                 ! MPI
#ifndef noMPI
      if(type.eq.'init') call mpi_init(ierror)                         ! MPI
      if(type.eq.'stop') call mpi_finalize(ierror)                     ! MPI
#endif
      return
      end subroutine                                                   ! MPI
      subroutine mpi_info                                              ! MPI
      nprocs=1
      myrank=0
#ifndef noMPI
      call mpi_comm_size(mpi_comm_world,nprocs,ierror)                 ! MPI
      call mpi_comm_rank(mpi_comm_world,myrank,ierror)                 ! MPI
#endif
      return
      end subroutine                                                   ! MPI
      end module                                                       ! MPI
c
c
#     include "setconf.f"
#     include "insertion.f"
c
c
      module engproc
      contains
c
c  procedure for constructing energy distribution functions
c
      subroutine enginit
c
      use engmain, only: numtype,nummol,engdiv,corrcal,slttype,
     #                   moltype,sluvid,
     #                   ermax,numslv,uvmax,uvsoft,esmax,uvspec,
     #                   uvcrd,edens,ecorr,escrd,eself,
     #                   aveuv,slnuv,avediv,minuv,maxuv,numslt,sltlist
c
      real ecdmin,ecfmns,ecmns0,ecdcen,ecpls0,ecfpls,eccore,ecdmax
      real eclbin,ecfbin,ec0bin,finfac,ectmvl
      integer peread,pemax,pesoft,pecore,sltmltp
      character(*), parameter :: ecdfile='EcdInfo'
      integer, parameter :: ecdio=51       ! IO for ecdfile
      real, parameter :: infty=1.0e50      ! essentially equal to infinity
c
      integer, parameter :: rglmax=5, large=10000
      real, parameter :: tiny=1.0e-30
      integer iduv,i,k,q,pti,regn,prvmx,curmx,uprgcd(rglmax+1)
      real factor,incre,cdrgvl(0:rglmax+1),rgcnt(rglmax),ecpmrd(10)
      integer, dimension(:), allocatable :: tplst
      real, dimension(:,:), allocatable  :: ercrd
c
      allocate( tplst(nummol) )
      numslt=0
      do 3701 i=1,nummol
        if(sluvid(i).gt.0) then
          numslt=numslt+1
          tplst(numslt)=i
          sltmltp=moltype(i)
        endif
3701  continue
      do 3702 i=1,nummol
        if((sluvid(i).gt.0).and.(moltype(i).ne.sltmltp)) then
          call eng_stop('typ')
        endif
3702  continue
      iduv=0
      if(numslt.eq.0) iduv=9
      if((slttype.ge.2).and.(numslt.ne.1)) iduv=9
      if(iduv.ne.0) call eng_stop('num')
      allocate( sltlist(numslt) )
      do 3703 i=1,numslt
        sltlist(i)=tplst(i)
3703  continue
      if((slttype.ge.2).and.(sltlist(1).ne.nummol)) call eng_stop('ins')
      deallocate( tplst )
c
      if(numslt.eq.1) numslv=numtype-1
      if(numslt.gt.1) numslv=numtype
c
      allocate( uvspec(nummol) )
      do 3531 i=1,nummol
        pti=moltype(i)
        if(sluvid(i).eq.0) then            ! solvent
          if(pti.lt.sltmltp) uvspec(i)=pti
          if(pti.eq.sltmltp) call eng_stop('typ')
          if(pti.gt.sltmltp) uvspec(i)=pti-1
        endif
        if(sluvid(i).ne.0) then            ! solute
          if(pti.ne.sltmltp) call eng_stop('typ')
          if(numslt.eq.1) uvspec(i)=0
          if(numslt.gt.1) uvspec(i)=numtype
        endif
3531  continue
c
      allocate( uvmax(numslv),uvsoft(numslv),ercrd(large,0:numslv) )
c
      peread=0
      ermax=0
      do 3001 pti=0,numslv
#       include "param_eng"
        if(peread.eq.1) then   ! read coordinate parameters from separate file
          open(unit=ecdio,file=ecdfile,status='old')
          read(ecdio,*)        ! comment line
          do 3101 i=1,large
            read(ecdio,*,END=3109) q
            if(q.eq.pti) then
              backspace(ecdio)
              if(pti.eq.0) then
                read(ecdio,*) k,(ecpmrd(k),k=1,8)
                pecore=0       ! no core region for solute self-energy
              endif
              if(pti.gt.0) then
                read(ecdio,*) k,(ecpmrd(k),k=1,9),pecore
                ecdmax=ecpmrd(9)
              endif
              eclbin=ecpmrd(1) ; ecfbin=ecpmrd(2) ; ec0bin=ecpmrd(3)
              finfac=ecpmrd(4) ; ecdmin=ecpmrd(5) ; ecfmns=ecpmrd(6)
              ecdcen=ecpmrd(7) ; eccore=ecpmrd(8)
              if(eccore.lt.tiny) pecore=0
              if(pecore.eq.1) call eng_stop('ecd')
              goto 3109
            endif
3101      continue
3109      continue
          close(ecdio)
        endif
        ectmvl=finfac*ecfbin ; ecdmin=ecdmin-ectmvl
        ecfmns=ecfmns-ectmvl ; ecmns0=ecdcen-ectmvl
        ecpls0=2.0e0*ecdcen-ecmns0 ; ecfpls=2.0e0*ecdcen-ecfmns
        eccore=eccore+ecfpls-ecdcen
c
        rgcnt(1)=(ecfmns-ecdmin)/eclbin
        rgcnt(2)=(ecmns0-ecfmns)/ecfbin
        rgcnt(3)=(ecpls0-ecmns0)/ec0bin
        rgcnt(4)=(ecfpls-ecpls0)/ecfbin
        rgcnt(5)=(eccore-ecfpls)/eclbin
        pesoft=0
        do 3151 regn=1,rglmax
          factor=rgcnt(regn)
          if(int(factor).lt.1) call eng_stop('ecd')
          pesoft=pesoft+nint(factor)
3151    continue
        pemax=pesoft+pecore
        if(pemax.gt.large) call eng_stop('siz')
c
        cdrgvl(0)=ecdmin
        cdrgvl(1)=ecfmns
        cdrgvl(2)=ecmns0
        cdrgvl(3)=ecpls0
        cdrgvl(4)=ecfpls
        cdrgvl(5)=eccore
        cdrgvl(6)=ecdmax
        do 3351 regn=1,rglmax-1
          if((regn.eq.1).or.(regn.eq.5)) factor=eclbin
          if((regn.eq.2).or.(regn.eq.4)) factor=ecfbin
          if(regn.eq.3) factor=ec0bin
          iduv=nint((cdrgvl(regn)-cdrgvl(regn-1))/factor)
          if(regn.eq.1) uprgcd(regn)=iduv
          if(regn.gt.1) uprgcd(regn)=uprgcd(regn-1)+iduv
3351    continue
        uprgcd(rglmax)=pesoft
        uprgcd(rglmax+1)=pemax-1
c
        do 3352 regn=0,rglmax+1
          if(regn.eq.0) iduv=0
          if(regn.gt.0) iduv=uprgcd(regn)
          ercrd(iduv+1,pti)=cdrgvl(regn)
3352    continue
c
        if(pecore.eq.0) i=rglmax       ! no core region
        if(pecore.gt.0) i=rglmax+1     ! explicit treatment of core region
        curmx=0
        do 3331 regn=1,i
          prvmx=curmx
          curmx=uprgcd(regn)
          if((regn.eq.1).or.(regn.eq.5)) factor=eclbin
          if((regn.eq.2).or.(regn.eq.4)) factor=ecfbin
          if(regn.eq.3) factor=ec0bin
          if(regn.eq.(rglmax+1)) then
            incre=log(ercrd(pemax,pti)/ercrd(pesoft+1,pti))
            factor=incre/real(pecore-1)
          endif
          do 3332 iduv=prvmx+1,curmx
            incre=factor*real(iduv-prvmx-1)
            if(regn.le.rglmax) ercrd(iduv,pti)=ercrd(prvmx+1,pti)+incre
            if(regn.eq.(rglmax+1)) then
              ercrd(iduv,pti)=ercrd(prvmx+1,pti)*exp(incre)
            endif
3332      continue
3331    continue
c
        if(pti.eq.0) esmax=pesoft      ! solute self-energy
        if(pti.gt.0) then              ! solute-solvent interaction energy
          uvmax(pti)=pemax
          uvsoft(pti)=pesoft
          ermax=ermax+pemax
        endif
3001  continue
c
      allocate( uvcrd(ermax),edens(ermax) )
      if(corrcal.eq.1) allocate( ecorr(ermax,ermax) )
      allocate( escrd(esmax),eself(esmax) )
      if(slttype.eq.1) allocate( aveuv(engdiv,numslv),slnuv(numslv) )
      allocate( avediv(engdiv,2) )
      allocate( minuv(0:numslv),maxuv(0:numslv) )
c
      i=0
      do 3731 pti=1,numslv
        pemax=uvmax(pti)
        do 3732 iduv=1,pemax
          i=i+1
          uvcrd(i)=ercrd(iduv,pti)
3732    continue
3731  continue
      do 3751 iduv=1,esmax
        escrd(iduv)=ercrd(iduv,0)
3751  continue
      deallocate( ercrd )
c
      do 3771 pti=0,numslv
        minuv(pti)=infty
        maxuv(pti)=-infty
3771  continue
c
      call engclear
c
      return
      end subroutine
c
c
      subroutine engclear
      use engmain, only: corrcal,slttype,ermax,numslv,esmax,
     #                   edens,ecorr,eself,slnuv,avslf,engnorm,engsmpl
      integer iduv,iduvp,pti
      do 3501 iduv=1,ermax
        edens(iduv)=0.0e0
3501  continue
      if(corrcal.eq.1) then
        do 3511 iduv=1,ermax
         do 3512 iduvp=1,ermax
           ecorr(iduvp,iduv)=0.0e0
3512     continue
3511    continue
      endif
      do 3521 iduv=1,esmax
        eself(iduv)=0.0e0
3521  continue
      if(slttype.eq.1) then
        do 3531 pti=1,numslv
          slnuv(pti)=0.0e0
3531    continue
      endif
      avslf=0.0e0
      engnorm=0.0e0
      engsmpl=0.0e0
      return
      end subroutine
c
c
      subroutine engconst(stnum)
c
      use engmain, only: nummol,maxcnf,skpcnf,corrcal,slttype,wgtslf,
     #                   estype,sluvid,temp,volume,plmode,
     #                   maxins,ermax,numslv,esmax,uvspec,
     #                   edens,ecorr,eself,
     #                   slnuv,avslf,minuv,maxuv,numslt,sltlist,
     #                   engnorm,engsmpl,voffset
      use ptinsrt, only: instslt
      use mpiproc                                                      ! MPI
      integer, parameter :: flcio=91                    ! IO unit for flcuv
      integer stnum,cntdst,maxdst,tagslt,slvmax,i,pti,iduv,iduvp,k,q
      integer ptinit,ptskip,dsinit,dsskip
      real engnmfc,pairep,wgtslcf,factor
      integer, dimension(:), allocatable :: insdst,engdst,tagpt,tplst
      real, dimension(:),    allocatable :: uvengy,flceng,svfl
      call mpi_info                                                    ! MPI
c
      if((slttype.eq.1).and.(myrank.eq.0).and.(stnum.eq.skpcnf)) then
        open(unit=flcio,file='flcuv.tt',status='new')   ! open flcuv file
      endif
c
      call sltcnd(q,0,'sys')
      if((slttype.eq.1).and.(q.ne.1)) q=9
      if((slttype.ge.2).and.(q.ne.2)) q=9
      if(q.eq.9) call eng_stop('par')
      if(slttype.eq.1) maxdst=numslt
      if(slttype.ge.2) maxdst=maxins
c
      if(plmode.eq.0) then
        ptinit=myrank ; ptskip=nprocs
        dsinit=0 ; dsskip=1
      endif
      if(plmode.eq.1) then
        ptinit=0 ; ptskip=1
        dsinit=myrank ; dsskip=nprocs
      endif
c
      allocate( tplst(nummol) )
      slvmax=0
      do 3001 i=1+ptinit,nummol,ptskip
         if((slttype.eq.1) .or. 
     #      ((slttype.ge.2).and.(sluvid(i).eq.0))) then
            slvmax=slvmax+1
            tplst(slvmax)=i
        endif
3001  continue
      allocate( insdst(ermax),engdst(ermax),tagpt(slvmax) )
      allocate( uvengy(0:slvmax),flceng(numslv),svfl(numslv) )
      do 3002 k=1,slvmax
        tagpt(k)=tplst(k)
3002  continue
      deallocate( tplst )
c
      call cellinfo
c
      if(stnum.eq.skpcnf) then                ! Ewald and PME initialization
        call recpcal(0,0,factor,slvmax,tagpt,'alloct')
      endif
      call recpcal(0,0,factor,slvmax,tagpt,'preeng')
c
      do 3101 k=1,slvmax
        i=tagpt(k)
        call recpcal(i,i,factor,slvmax,tagpt,'charge')
3101  continue
c
      do 90000 cntdst=1,maxdst
        if(slttype.eq.1) then                 ! solution system
          tagslt=sltlist(cntdst)
          call sltcnd(q,tagslt,'pos')
          if(q.eq.0) go to 99999
          if((q.ne.0).and.(q.ne.1)) call eng_stop('slt')
        endif
        if(slttype.ge.2) then                 ! solute insertion
          tagslt=sltlist(1)
          if((stnum.eq.skpcnf).and.(cntdst.eq.1)) then
            call instslt(wgtslcf,'init')
          endif
          call instslt(wgtslcf,'proc')
          if((stnum.eq.maxcnf).and.(cntdst.eq.maxdst)) then
            call instslt(wgtslcf,'last')
          endif
          if(mod(cntdst-1,dsskip).ne.dsinit) go to 99999
        endif
c
        call recpcal(tagslt,tagslt,factor,slvmax,tagpt,'charge')
        do 1101 k=0,slvmax
          if(k.eq.0) i=tagslt                 ! solute self
          if(k.gt.0) then                     ! solute-solvent pair
            i=tagpt(k)
            if(i.eq.tagslt) goto 1199
          endif
          call realcal(tagslt,i,pairep)
          call recpcal(tagslt,i,factor,slvmax,tagpt,'energy')
          pairep=pairep+factor
          uvengy(k)=pairep
1199      continue
1101    continue
c
        if(wgtslf.eq.0) engnmfc=1.0e0
        if(wgtslf.eq.1) then
          factor=uvengy(0)
          if((stnum.eq.skpcnf).and.(cntdst.eq.(1+dsinit))) then
            if(dsinit.eq.0) voffset=factor
#ifndef noMPI
            if(plmode.eq.1) call mpi_bcast(voffset,1,                  ! MPI
     #                mpi_double_precision,0,mpi_comm_world,ierror)    ! MPI
#endif
          endif
          factor=factor-voffset               ! shifted by offset
          if(slttype.eq.1) engnmfc=exp(factor/temp)
          if(slttype.ge.2) engnmfc=exp(-factor/temp)
        endif
        if(estype.eq.2) call volcorrect(engnmfc)
        if(slttype.ge.2) engnmfc=engnmfc*wgtslcf
c
        engnorm=engnorm+engnmfc               ! normalization factor
        engsmpl=engsmpl+1.0e0                 ! number of sampling
        if(estype.le.1) avslf=avslf+1.0e0
        if(estype.eq.2) avslf=avslf+volume
c
        do 1011 iduv=1,ermax
          insdst(iduv)=0
1011    continue
        do 1012 pti=1,numslv
          flceng(pti)=0.0e0                   ! sum of solute-solvent energy
1012    continue
        do 1001 k=0,slvmax
          if(k.eq.0) pti=0                    ! solute self
          if(k.gt.0) then                     ! solute-solvent pair
            i=tagpt(k)
            if(i.eq.tagslt) goto 1099
            pti=uvspec(i)
            if(pti.le.0) call eng_stop('eng')
          endif
          pairep=uvengy(k)
          call getiduv(pti,pairep,iduv)
          if(k.eq.0) eself(iduv)=eself(iduv)+engnmfc
          if(k.gt.0) then
            insdst(iduv)=insdst(iduv)+1
            flceng(pti)=flceng(pti)+pairep    ! sum of solute-solvent energy
          endif
          if(pairep.lt.minuv(pti)) minuv(pti)=pairep        ! minimum energy
          if(pairep.gt.maxuv(pti)) maxuv(pti)=pairep        ! maximum energy
1099      continue
1001    continue
c
#ifndef noMPI
        if(plmode.eq.0) then
          call mpi_allreduce(insdst,engdst,ermax,                      ! MPI
     #                mpi_integer,mpi_sum,mpi_comm_world,ierror)       ! MPI
          call mpi_allreduce(flceng,svfl,numslv,mpi_double_precision,  ! MPI
     #                mpi_sum,mpi_comm_world,ierror)                   ! MPI
          do 1331 iduv=1,ermax                                         ! MPI
            insdst(iduv)=engdst(iduv)                                  ! MPI
1331      continue                                                     ! MPI
          do 1332 pti=1,numslv                                         ! MPI
            flceng(pti)=svfl(pti)                                      ! MPI
1332      continue                                                     ! MPI
        endif
#endif
        if(slttype.eq.1) then
          do 1151 pti=1,numslv
            slnuv(pti)=slnuv(pti)+flceng(pti)*engnmfc
1151      continue
          if(myrank.eq.0) then
            write(flcio,911) stnum,(flceng(pti), pti=1,numslv)
          endif
911       format(i9,999999999f15.5)
        endif
c
        do 1501 iduv=1+ptinit,ermax,ptskip
          k=insdst(iduv)
          if(k.gt.0) edens(iduv)=edens(iduv)+engnmfc*real(k)
1501    continue
        if(corrcal.eq.1) then
          do 1701 iduv=1+ptinit,ermax,ptskip
            k=insdst(iduv)
            if(k.gt.0) then
              do 1702 iduvp=1,ermax
                q=insdst(iduvp)
                if(q.gt.0) then
                  factor=engnmfc*real(k)*real(q)
                  ecorr(iduvp,iduv)=ecorr(iduvp,iduv)+factor
                endif
1702          continue
            endif
1701      continue
        endif
c
99999 continue
90000 continue
c
      if((slttype.eq.1).and.(myrank.eq.0).and.(stnum.eq.maxcnf)) then
        endfile(flcio) ; close(flcio)                   ! close flcuv file
      endif
c
      deallocate( insdst,engdst,tagpt,uvengy,flceng,svfl )
c
      return
      end subroutine
c
c
      subroutine engstore(stnum)
c
      use engmain, only: nummol,maxcnf,engdiv,corrcal,slttype,wgtslf,
     #                   plmode,ermax,numslv,esmax,temp,
     #                   edens,ecorr,eself,
     #                   aveuv,slnuv,avediv,avslf,minuv,maxuv,
     #                   engnorm,engsmpl,voffset
      use mpiproc                                                      ! MPI
      integer stnum,i,pti,j,iduv,iduvp,k,q,cntdst
      character*10, parameter :: numbers='0123456789'
      character*9 engfile
      character*3 suffeng
      real factor
      real, parameter :: tiny=1.0e-30
      real, dimension(:), allocatable :: sve1,sve2
      call mpi_info                                                    ! MPI
c
#ifndef noMPI
      if(plmode.eq.1) then                                             ! MPI
        call mpi_reduce(avslf,factor,1,                                ! MPI
     #       mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)     ! MPI
        avslf=factor                                                   ! MPI
        call mpi_reduce(engnorm,factor,1,                              ! MPI
     #       mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)     ! MPI
        engnorm=factor                                                 ! MPI
        call mpi_reduce(engsmpl,factor,1,                              ! MPI
     #       mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)     ! MPI
        engsmpl=factor                                                 ! MPI
        allocate( sve1(esmax) )                                        ! MPI
        do 7701 iduv=1,esmax                                           ! MPI
          sve1(iduv)=eself(iduv)                                       ! MPI
7701    continue                                                       ! MPI
        call mpi_reduce(sve1,eself,esmax,                              ! MPI
     #       mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)     ! MPI
        deallocate( sve1 )                                             ! MPI
      endif                                                            ! MPI
      allocate( sve1(0:numslv),sve2(0:numslv) )                        ! MPI
      do 7731 pti=0,numslv                                             ! MPI
        sve1(pti)=minuv(pti)                                           ! MPI
        sve2(pti)=maxuv(pti)                                           ! MPI
7731  continue                                                         ! MPI
      call mpi_reduce(sve1,minuv,numslv+1,                             ! MPI
     #     mpi_double_precision,mpi_min,0,mpi_comm_world,ierror)       ! MPI
      call mpi_reduce(sve2,maxuv,numslv+1,                             ! MPI
     #     mpi_double_precision,mpi_max,0,mpi_comm_world,ierror)       ! MPI
      deallocate( sve1,sve2 )                                          ! MPI
      allocate( sve1(ermax) )                                          ! MPI
      do 7502 iduv=1,ermax                                             ! MPI
        sve1(iduv)=edens(iduv)                                         ! MPI
7502  continue                                                         ! MPI
      call mpi_reduce(sve1,edens,ermax,                                ! MPI
     #     mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)       ! MPI
      deallocate( sve1 )                                               ! MPI
#endif
      do 7501 iduv=1,ermax
        edens(iduv)=edens(iduv)/engnorm
7501  continue
      if(corrcal.eq.1) then
        do 7511 iduv=1,ermax
#ifndef noMPI
          allocate( sve1(ermax),sve2(ermax) )                          ! MPI
          do 7513 iduvp=1,ermax                                        ! MPI
            sve1(iduvp)=ecorr(iduvp,iduv)                              ! MPI
7513      continue                                                     ! MPI
          call mpi_reduce(sve1,sve2,ermax,                             ! MPI
     #         mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror)   ! MPI
          do 7514 iduvp=1,ermax                                        ! MPI
            ecorr(iduvp,iduv)=sve2(iduvp)                              ! MPI
7514      continue                                                     ! MPI
          deallocate( sve1,sve2 )                                      ! MPI
#endif
          do 7512 iduvp=1,ermax
            ecorr(iduvp,iduv)=ecorr(iduvp,iduv)/engnorm
7512      continue
7511    continue
      endif
      do 7521 iduv=1,esmax
        eself(iduv)=eself(iduv)/engnorm
7521  continue
      avslf=avslf/engnorm
c
      if(myrank.ne.0) go to 7999                                       ! MPI
      i=stnum/(maxcnf/engdiv)
      if(slttype.eq.1) then
        do 7551 pti=1,numslv
          aveuv(i,pti)=slnuv(pti)/engnorm
7551    continue
      endif
      avediv(i,1)=engnorm/engsmpl
      if(slttype.eq.1) avediv(i,2)=voffset-temp*log(avslf)
      if(slttype.ge.2) avediv(i,2)=voffset+temp*log(avslf)
c
      if(stnum.eq.maxcnf) then
        if(slttype.eq.1) then
          open(unit=75,file='aveuv.tt',status='new')
          do 7553 k=1,engdiv
            write(75,751) k,(aveuv(k,pti), pti=1,numslv)
7553      continue
          endfile(75)
          close(75)
751       format(i5,99999f15.5)
        endif
        if(slttype.eq.1) open(unit=73,file='weight_soln',status='new')
        if(slttype.ge.2) open(unit=73,file='weight_refs',status='new')
        do 7311 k=1,engdiv
          if(wgtslf.eq.0) write(73,731) k,avediv(k,1)
          if(wgtslf.eq.1) write(73,732) k,avediv(k,1),avediv(k,2)
7311    continue
        endfile(73)
        close(73)
731     format(i5,f15.3)
732     format(i5,f15.3,g17.5)
        open(77,file='uvrange.tt',status='new')
        write(77,771)
        do 7711 pti=0,numslv
          factor=maxuv(pti)
          if(factor.lt.1.0e5) write(77,772) pti,minuv(pti),factor
          if(factor.ge.1.0e5) write(77,773) pti,minuv(pti),factor
7711    continue
        endfile(77)
        close(77)
771     format(' species     minimum        maximum')
772     format(i5,2f15.5)
773     format(i5,f15.5,g18.5)
      endif
c
      j=i/10
      k=i-10*j
      if(engdiv.eq.1) suffeng='.tt'
      if(engdiv.gt.1) suffeng='.'//numbers(j+1:j+1)//numbers(k+1:k+1)
      do 7011 cntdst=1,3
        if(cntdst.eq.1) then
          if(slttype.eq.1) engfile='engsln'//suffeng
          if(slttype.ge.2) engfile='engref'//suffeng
        endif
        if((cntdst.eq.2).and.(corrcal.ne.1)) goto 7199
        if((cntdst.eq.2).and.(corrcal.eq.1)) then
          if(slttype.eq.1) engfile='corsln'//suffeng
          if(slttype.ge.2) engfile='corref'//suffeng
        endif
        if(cntdst.eq.3) engfile='slfeng'//suffeng
        open(unit=71,file=engfile,status='new')
        if(cntdst.eq.1) then
          do 7101 iduv=1,ermax
            call repval(iduv,factor,pti,'intn')
            write(71,'(g15.7,i5,g25.15)') factor,pti,edens(iduv)
7101      continue
        endif
        if(cntdst.eq.2) then
          do 7111 iduv=1,ermax
           do 7112 iduvp=1,ermax
             factor=ecorr(iduvp,iduv)
             if(factor.gt.tiny) write(71,'(g25.15)') factor
             if(factor.le.tiny) write(71,'(i1)') 0
7112       continue
7111      continue
        endif
        if(cntdst.eq.3) then
          do 7151 iduv=1,esmax
            call repval(iduv,factor,pti,'self')
            write(71,'(g15.7,g25.15)') factor,eself(iduv)
7151      continue
        endif
        endfile(71)
        close(71)
7199    continue
7011  continue
7999  continue                                                         ! MPI
c
      return
      end subroutine
c
c
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
991   format(' The number of solute types is incorrectly set')
992   format(' The number of solute molecules is incorrectly set')
993   format(' The solute numbering is incorrect for insertion')
994   format(' The input parameter is incorrectly set')
995   format(' The input parameter is incorrect for solute')
996   format(' The coordinate system is incorrectly set')
997   format(' Inconsistency is present in the program')
998   format(' The number of energy-coordinate meshes is too large')
999   format(' The minimum of the energy coordinate is too large')
981   format(' The energy-coordinate system is inconsistent')
      call mpi_setup('stop')                                           ! MPI
      stop
      end subroutine
c
c
      subroutine realcal(i,j,pairep)
c
      use engmain, only:  nummol,maxsite,numatm,boxshp,numsite,
     #                    elecut,lwljcut,upljcut,cmbrule,cltype,screen,
     #                    charge,ljene,ljlen,specatm,sitepos,
     #                    cell,invcl,volume
      integer i,j,is,js,ismax,jsmax,ati,atj,m,k
      real pi,reelcut,pairep,ljeps,ljsgm,chr2,rst,dis2,rtp1,rtp2
      real eplj,epcl,xst(3),clm(3),swth
      real, parameter :: infty=1.0e50      ! essentially equal to infinity
c
      pi=real(4)*atan(real(1))
      if(i.eq.j) reelcut=infty             ! no cutoff for self-interaction
      if(i.ne.j) then
        if(cltype.eq.0) then                  ! bare Coulomb
          if(boxshp.eq.0) reelcut=infty
          if(boxshp.ne.0) reelcut=elecut
        endif
        if(cltype.ne.0) reelcut=elecut        ! Ewald and PME
      endif
c
      pairep=0.0e0
      ismax=numsite(i)
      jsmax=numsite(j)
c
      do 1501 is=1,ismax
       do 1502 js=1,jsmax
         if((i.eq.j).and.(is.gt.js)) go to 1599
         ati=specatm(is,i)
         atj=specatm(js,j)
         do 5501 m=1,3
           xst(m)=sitepos(m,ati)-sitepos(m,atj)
5501     continue
         if(boxshp.ne.0) then              ! when the system is periodic
           do 5551 k=1,3
             rst=0.0e0
             do 5552 m=1,3
               rst=rst+invcl(k,m)*xst(m)
5552         continue
             clm(k)=real(nint(rst))
5551       continue
           do 5553 m=1,3
             rst=0.0e0
             do 5554 k=1,3
               rst=rst+cell(m,k)*clm(k)
5554         continue
             xst(m)=xst(m)-rst             ! get the nearest distance between i,j
5553       continue
         endif
         dis2=xst(1)*xst(1)+xst(2)*xst(2)+xst(3)*xst(3)
         rst=sqrt(dis2)
         if((rst.gt.upljcut).or.(i.eq.j)) eplj=0.0e0
         if((rst.le.upljcut).and.(i.ne.j)) then
           ljeps=sqrt(ljene(ati)*ljene(atj))
           if(cmbrule.eq.0) ljsgm=(ljlen(ati)+ljlen(atj))/2.0e0
           if(cmbrule.eq.1) ljsgm=sqrt(ljlen(ati)*ljlen(atj))
           rtp1=ljsgm*ljsgm/dis2
           rtp2=rtp1*rtp1*rtp1
           eplj=4.0e0*ljeps*rtp2*(rtp2-1.0e0)
           if(rst.gt.lwljcut) then    ! CHARMM form of switching function
             rtp1=lwljcut*lwljcut
             rtp2=upljcut*upljcut
             swth=(2.0e0*dis2+rtp2-3.0e0*rtp1)*(dis2-rtp2)*(dis2-rtp2)
     #           /(rtp2-rtp1)/(rtp2-rtp1)/(rtp2-rtp1)
             eplj=swth*eplj
           endif
         endif
         if(rst.gt.reelcut) epcl=0.0e0
         if(rst.le.reelcut) then
           chr2=charge(ati)*charge(atj)
           if((i.eq.j).and.(is.eq.js)) then
             if(cltype.eq.0) rtp1=real(0)                  ! bare Coulomb
             if(cltype.ne.0) rtp1=screen                   ! Ewald and PME
             epcl=-chr2*rtp1/sqrt(pi)
           endif
           if((i.ne.j).or.(is.ne.js)) then
             if(cltype.eq.0) rtp1=real(0)                  ! bare Coulomb
             if(cltype.ne.0) rtp1=screen*rst               ! Ewald and PME
             if(i.eq.j) epcl=-chr2*derf(rtp1)/rst
             if(i.ne.j) epcl=chr2*(1.0e0-derf(rtp1))/rst
           endif
         endif
         pairep=pairep+(eplj+epcl)
1599     continue
1502   continue
1501  continue
c
      if(cltype.ne.0) then                                 ! Ewald and PME
        rtp1=0.0e0
        do 2501 is=1,ismax
          ati=specatm(is,i)
          rtp1=rtp1+charge(ati)
2501    continue
        rtp2=0.0e0
        do 2502 js=1,jsmax
          atj=specatm(js,j)
          rtp2=rtp2+charge(atj)
2502    continue
        epcl=pi*rtp1*rtp2/screen/screen/volume
        if(i.eq.j) epcl=epcl/2.0e0                         ! self-interaction
        pairep=pairep-epcl
      endif
c
      return
      end subroutine
c
c
      subroutine volcorrect(engnmfc)
      use engmain, only:  nummol,maxsite,numatm,temp,numsite,sluvid,
     #                    cltype,screen,charge,specatm,volume
      integer i,ati,sid,stmax
      real pi,factor,engnmfc
      engnmfc=volume*engnmfc
      if(cltype.ne.0) then                                 ! Ewald and PME
        pi=real(4)*atan(real(1))
        factor=0.0e0
        do 1011 i=1,nummol
          if(sluvid(i).le.1) then                          ! physical particle
            stmax=numsite(i)
            do 1012 sid=1,stmax
              ati=specatm(sid,i)
              factor=factor+charge(ati)
1012        continue
          endif
1011    continue
        factor=pi*factor*factor/screen/screen/volume/2.0e0
        engnmfc=engnmfc*exp(factor/temp)
      endif
      return
      end subroutine
c
c
      subroutine cellinfo
      use engmain, only:  cell,invcl,volume
      integer m,k
      volume=cell(1,1)*cell(2,2)*cell(3,3)
     #      +cell(1,2)*cell(2,3)*cell(3,1)+cell(1,3)*cell(2,1)*cell(3,2)
     #      -cell(1,3)*cell(2,2)*cell(3,1)
     #      -cell(1,2)*cell(2,1)*cell(3,3)-cell(1,1)*cell(2,3)*cell(3,2)
      invcl(1,1)=cell(2,2)*cell(3,3)-cell(2,3)*cell(3,2)
      invcl(1,2)=cell(1,3)*cell(3,2)-cell(1,2)*cell(3,3)
      invcl(1,3)=cell(1,2)*cell(2,3)-cell(1,3)*cell(2,2)
      invcl(2,1)=cell(2,3)*cell(3,1)-cell(2,1)*cell(3,3)
      invcl(2,2)=cell(1,1)*cell(3,3)-cell(1,3)*cell(3,1)
      invcl(2,3)=cell(1,3)*cell(2,1)-cell(1,1)*cell(2,3)
      invcl(3,1)=cell(2,1)*cell(3,2)-cell(2,2)*cell(3,1)
      invcl(3,2)=cell(1,2)*cell(3,1)-cell(1,1)*cell(3,2)
      invcl(3,3)=cell(1,1)*cell(2,2)-cell(1,2)*cell(2,1)
      do 3111 m=1,3
       do 3112 k=1,3
         invcl(k,m)=invcl(k,m)/volume
3112   continue
3111  continue
      return
      end subroutine
c
c
      subroutine recpcal(i,j,pairep,slvmax,tagpt,scheme)
c
      use engmain, only:  nummol,maxsite,numatm,numsite,sluvid,
     #                    cltype,screen,splodr,charge,
     #                    ew1max,ew2max,ew3max,ms1max,ms2max,ms3max,
     #                    specatm,sitepos,invcl,volume
      use spline, only: spline_init, spline_value
      use fft_iface, only: fft_init_ctc, fft_init_inplace, 
     #                     fft_ctc, fft_inplace,
     #                     fft_set_size
      integer i,j,ptrnk,slvmax,tagpt(slvmax)
      integer svi,svj,uvi,uvj,ati,sid,stmax,m,k
      integer rc1,rc2,rc3,rci,rcimax,spi,cg1,cg2,cg3
      real pi,pairep,chr,xst(3),inm(3),rtp2,cosk,sink,factor
      complex rcpi,rcpj,nmfact(3)
      character*6 scheme
c
      integer, save :: rc1min,rc1max,rc2min,rc2max,rc3min,rc3max
      integer, dimension(:),       allocatable, save :: slvtag
      real, dimension(:,:,:),      allocatable, save :: engfac
      complex, dimension(:,:,:,:), allocatable, save :: rcpslv
      complex, dimension(:,:,:),   allocatable, save :: rcpslt
      real, dimension(:,:,:),      allocatable, save :: splslv
      integer, dimension(:,:),     allocatable, save :: grdslv
      complex, dimension(:,:,:),   allocatable, save :: cnvslt
      real, dimension(:),          allocatable, save :: splint
c
      real, dimension(:,:,:),      allocatable :: splval
      integer, dimension(:,:),     allocatable :: grdval
      integer :: si, gridsize(3)
      complex, allocatable, save :: fft_buf(:, :, :)
c
      if(cltype.eq.0) then                                   ! bare Coulomb
        pairep=0.0e0
        return
      endif
c
      pi=real(4)*atan(real(1))
c
      if(scheme.eq.'alloct') then
        allocate( slvtag(nummol) )
        do 3109 m=1,nummol
          slvtag(m)=-1
3109    continue
        ptrnk=0
        do 3101 k=1,slvmax
          m=tagpt(k)
          slvtag(m)=ptrnk+1
          if(cltype.eq.1) ptrnk=ptrnk+1                      ! Ewald
          if(cltype.eq.2) ptrnk=ptrnk+numsite(m)             ! PME
3101    continue
        if(cltype.eq.1) then                                 ! Ewald
          rc1min=-ew1max ; rc1max=ew1max
          rc2min=-ew2max ; rc2max=ew2max
          rc3min=-ew3max ; rc3max=ew3max
          allocate( rcpslv(rc1min:rc1max,rc2min:rc2max,
     #                                   rc3min:rc3max,ptrnk) )
        endif
        if(cltype.eq.2) then                                 ! PME
          rc1min=0 ; rc1max=ms1max-1
          rc2min=0 ; rc2max=ms2max-1
          rc3min=0 ; rc3max=ms3max-1
          call spline_init(splodr)
          allocate( splslv(0:splodr-1,3,ptrnk),grdslv(3,ptrnk) )
          allocate( cnvslt(rc1min:rc1max,rc2min:rc2max,rc3min:rc3max) )
          allocate(splint(1:(splodr-1)))
          do si = 1, splodr - 1
             splint(si) = spline_value(dble(si))
          enddo
          ! allocate fft-buffers
          allocate( fft_buf(rc1min:rc1max,rc2min:rc2max,rc3min:rc3max) )
          gridsize(1) = ms1max
          gridsize(2) = ms2max
          gridsize(3) = ms3max
          call fft_set_size(gridsize)
        endif
        allocate( engfac(rc1min:rc1max,rc2min:rc2max,rc3min:rc3max) )
        allocate( rcpslt(rc1min:rc1max,rc2min:rc2max,rc3min:rc3max) )
        if(cltype .eq. 2) then
          ! init fft
          call fft_init_inplace(rcpslt)
          call fft_init_ctc(fft_buf, cnvslt)
        endif
      endif
c
      if(scheme.eq.'preeng') then
        do 3201 rc3=rc3min,rc3max
         do 3202 rc2=rc2min,rc2max
          do 3203 rc1=rc1min,rc1max
            factor=0.0e0
            m=rc1*rc1+rc2*rc2+rc3*rc3
            if(m.eq.0) go to 3219
            if(cltype.eq.1) then                             ! Ewald
              rcimax=min(ew1max,ew2max,ew3max)
              if(m.gt.rcimax*rcimax) go to 3219
            endif
            do 3211 m=1,3
              if(m.eq.1) rci=rc1
              if(m.eq.2) rci=rc2
              if(m.eq.3) rci=rc3
              if(cltype.eq.1) inm(m)=real(rci)               ! Ewald
              if(cltype.eq.2) then                           ! PME
                if(m.eq.1) rcimax=ms1max
                if(m.eq.2) rcimax=ms2max
                if(m.eq.3) rcimax=ms3max
                if((mod(splodr,2).eq.1).and.(2*abs(rci).eq.rcimax)) then
                  go to 3219
                endif
                if(rci.le.rcimax/2) inm(m)=real(rci)
                if(rci.gt.rcimax/2) inm(m)=real(rci-rcimax)
                rcpi=(0.0e0,0.0e0)
                do 3215 spi=0,splodr-2
                  chr=splint(spi+1)
                  rtp2=2.0e0*pi*real(spi*rci)/real(rcimax)
                  cosk=chr*cos(rtp2)
                  sink=chr*sin(rtp2)
                  rcpi=rcpi+cmplx(cosk,sink)
3215            continue
                nmfact(m)=rcpi
              endif
3211        continue
            do 3221 m=1,3
              factor=0.0e0
              do 3222 k=1,3
                factor=factor+invcl(k,m)*inm(k)
3222          continue
              xst(m)=factor
3221        continue
            rtp2=xst(1)*xst(1)+xst(2)*xst(2)+xst(3)*xst(3)
            chr=pi*pi*rtp2/screen/screen
            factor=exp(-chr)/rtp2/pi/volume
            if(cltype.eq.2) then                             ! PME
              rtp2=abs(nmfact(1)*nmfact(2)*nmfact(3))
              factor=factor/rtp2/rtp2
            endif
3219        continue
            engfac(rc1,rc2,rc3)=factor
3203      continue
3202     continue
3201    continue
      endif
c
      if(scheme.eq.'charge') then
        if(i.ne.j) return
        uvi=sluvid(i)
        if(uvi.eq.0) svi=slvtag(i)
        if((uvi.eq.0).and.(svi.le.0)) return
        stmax=numsite(i)
        if(cltype.eq.1) then                                 ! Ewald
          do 5101 rc3=rc3min,rc3max
           do 5102 rc2=rc2min,rc2max
            do 5103 rc1=rc1min,rc1max
              cosk=0.0e0
              sink=0.0e0
              do 5111 sid=1,stmax
                ati=specatm(sid,i)
                do 5112 m=1,3
                  xst(m)=sitepos(m,ati)
5112            continue
                do 5113 k=1,3
                  factor=0.0e0
                  do 5114 m=1,3
                    factor=factor+invcl(k,m)*xst(m)
5114              continue
                  inm(k)=factor
5113            continue
                rtp2=real(rc1)*inm(1)+real(rc2)*inm(2)+real(rc3)*inm(3)
                rtp2=2.0e0*pi*rtp2
                cosk=cosk+charge(ati)*cos(rtp2)
                sink=sink-charge(ati)*sin(rtp2)
5111          continue
              rcpi=cmplx(cosk,sink)
              if(uvi.eq.0) rcpslv(rc1,rc2,rc3,svi)=rcpi
              if(uvi.gt.0) rcpslt(rc1,rc2,rc3)=rcpi
5103        continue
5102       continue
5101      continue
        endif
        if(cltype.eq.2) then                                 ! PME
          allocate( splval(0:splodr-1,3,stmax),grdval(3,stmax) )
          do 5201 sid=1,stmax
            ati=specatm(sid,i)
            do 5202 m=1,3
              xst(m)=sitepos(m,ati)
5202        continue
            do 5203 k=1,3
              factor=0.0e0
              do 5204 m=1,3
                factor=factor+invcl(k,m)*xst(m)
5204          continue
              if(factor.lt.0.0e0) factor=factor+1.0e0
              if(factor.gt.1.0e0) factor=factor-1.0e0
              if((factor.lt.0.0e0).or.(factor.gt.1.0e0)) then
                call eng_stop('crd')
              endif
              inm(k)=factor
5203        continue
            do 5205 m=1,3
              if(m.eq.1) rcimax=ms1max
              if(m.eq.2) rcimax=ms2max
              if(m.eq.3) rcimax=ms3max
              factor=inm(m)*real(rcimax)
              rci=int(factor)
              do 5206 spi=0,splodr-1
                rtp2=factor-real(rci-spi)
                splval(spi,m,sid)=spline_value(rtp2)
5206          continue
              grdval(m,sid)=rci
5205        continue
5201      continue
          if(uvi.eq.0) then
            do 5211 sid=1,stmax
              ptrnk=svi+sid-1
              do 5212 m=1,3
                do 5213 spi=0,splodr-1
                  splslv(spi,m,ptrnk)=splval(spi,m,sid)
5213            continue
                grdslv(m,ptrnk)=grdval(m,sid)
5212          continue
5211        continue
          endif
          if(uvi.gt.0) then
            rcpslt(:, :, :)=(0.0e0,0.0e0)
            do 5251 sid=1,stmax
              ati=specatm(sid,i)
              chr=charge(ati)
              do 5231 cg3=0,splodr-1
               do 5232 cg2=0,splodr-1
                do 5233 cg1=0,splodr-1
                  rc1=modulo(grdval(1,sid)-cg1,ms1max)
                  rc2=modulo(grdval(2,sid)-cg2,ms2max)
                  rc3=modulo(grdval(3,sid)-cg3,ms3max)
                  factor=chr*splval(cg1,1,sid)*splval(cg2,2,sid)
     #                                        *splval(cg3,3,sid)
                  rcpi=cmplx(factor,0.0e0)
                  rcpslt(rc1,rc2,rc3)=rcpslt(rc1,rc2,rc3)+rcpi
5233            continue
5232           continue
5231          continue
5251        continue
            ! FIXME: rewrite to real-to-complex transform
            call fft_inplace(rcpslt)                         ! 3D-FFT
            do 5241 rc3=rc3min,rc3max
             do 5242 rc2=rc2min,rc2max
              do 5243 rc1=rc1min,rc1max
                rcpi=cmplx(engfac(rc1,rc2,rc3),0.0e0)
                fft_buf(rc1,rc2,rc3)=rcpi*conjg(rcpslt(rc1,rc2,rc3))
5243          continue
5242         continue
5241        continue
            call fft_ctc(fft_buf, cnvslt)                    ! 3D-FFT
          endif
          deallocate( splval,grdval )
        endif
      endif
c
      if(scheme.eq.'energy') then
        pairep=0.0e0
        uvi=sluvid(i)
        uvj=sluvid(j)
        if(uvi.eq.0) svi=slvtag(i)
        if(uvj.eq.0) svj=slvtag(j)
        if((uvi.eq.0).and.(svi.le.0)) return
        if((uvj.eq.0).and.(svj.le.0)) return
        if((i.eq.j).and.(uvi.eq.0).and.(uvj.eq.0)) call eng_stop('eng')
        if(cltype.eq.1) then                                 ! Ewald
          do 7101 rc3=rc3min,rc3max
           do 7102 rc2=rc2min,rc2max
            do 7103 rc1=rc1min,rc1max
              if(uvi.eq.0) rcpi=rcpslv(rc1,rc2,rc3,svi)
              if(uvi.gt.0) rcpi=rcpslt(rc1,rc2,rc3)
              if(uvj.eq.0) rcpj=rcpslv(rc1,rc2,rc3,svj)
              if(uvj.gt.0) rcpj=rcpslt(rc1,rc2,rc3)
              pairep=pairep+engfac(rc1,rc2,rc3)*real(rcpi*conjg(rcpj))
7103        continue
7102       continue
7101      continue
          if(i.eq.j) pairep=pairep/2.0e0
        endif
        if(cltype.eq.2) then                                 ! PME
          if(uvi.eq.0) call eng_stop('eng')   ! particle i needs to be solute
          if(i.eq.j) then
            do 7201 rc3=rc3min,rc3max
             do 7202 rc2=rc2min,rc2max
              do 7203 rc1=rc1min,rc1max
                rcpi=rcpslt(rc1,rc2,rc3)
                pairep=pairep+engfac(rc1,rc2,rc3)*real(rcpi*conjg(rcpi))
7203          continue
7202         continue
7201        continue
            pairep=pairep/2.0e0
          endif
          if(i.ne.j) then
            stmax=numsite(j)
            do 7251 sid=1,stmax
              ptrnk=svj+sid-1
              ati=specatm(sid,j)
              chr=charge(ati)
              do 7271 cg3=0,splodr-1
               do 7272 cg2=0,splodr-1
                do 7273 cg1=0,splodr-1
                  factor=chr*splslv(cg1,1,ptrnk)*splslv(cg2,2,ptrnk)
     #                                          *splslv(cg3,3,ptrnk)
                  rc1=modulo(grdslv(1,ptrnk)-cg1,ms1max)
                  rc2=modulo(grdslv(2,ptrnk)-cg2,ms2max)
                  rc3=modulo(grdslv(3,ptrnk)-cg3,ms3max)
                  pairep=pairep+factor*real(cnvslt(rc1,rc2,rc3))
7273            continue
7272           continue
7271          continue
7251        continue
          endif
        endif
      endif
c
      return
      end subroutine
c
c
      subroutine getiduv(pti,factor,iduv)
      use engmain, only: ermax,numslv,uvmax,uvcrd,esmax,escrd,io6
      integer pti,iduv,k,idpick,idmax
      real factor,egcrd
      if(pti.eq.0) idmax=esmax               ! solute self-energy
      if(pti.gt.0) idmax=uvmax(pti)          ! solute-solvent interaction
      idpick=0
      if(pti.gt.1) then
        do 1051 k=1,pti-1
          idpick=idpick+uvmax(k)
1051    continue
      endif
      iduv=idpick
      do 1011 k=1,idmax
        if(pti.eq.0) egcrd=escrd(k)          ! solute self-energy
        if(pti.gt.0) egcrd=uvcrd(k+idpick)   ! solute-solvent interaction
        if(factor.ge.egcrd) iduv=iduv+1
        if(factor.lt.egcrd) go to 1099
1011  continue
1099  continue
      if(iduv.le.idpick) then
        iduv=idpick+1                                 ! smallest energy mesh
        write(io6,199) factor,pti
199     format('  energy of ',g12.4,' for ',i3,'-th species')
        call eng_stop('min')
      endif
      if(iduv.gt.(idpick+idmax)) iduv=idpick+idmax    ! largest energy mesh
      return
      end subroutine
c
c
      subroutine sltcnd(systype,tagslt,type)
      use engmain, only: nummol,sluvid
      use ptinsrt, only: sltpstn
      integer systype,tagslt,i,uvi,cntuv(2)
      real xst(3)
      character*3 type
      if(type.eq.'sys') then
        systype=9
        do 5553 uvi=1,2
          cntuv(uvi)=0
5553    continue
        do 5551 i=1,nummol
          uvi=sluvid(i)
          if(uvi.eq.3) uvi=2
          if(uvi.gt.0) cntuv(uvi)=cntuv(uvi)+1
5551    continue
        if((cntuv(1).ne.0).and.(cntuv(2).eq.0)) systype=1
        if((cntuv(1).eq.0).and.(cntuv(2).ne.0)) systype=2
      endif
      if(type.eq.'pos') call sltpstn(systype,xst,'solutn',tagslt)
      return
      end subroutine
c
c
      subroutine repval(iduv,factor,pti,caltype)
      use engmain, only: ermax,numslv,uvmax,uvsoft,uvcrd,esmax,escrd
      integer iduv,idpt,pti,cnt,idpick,idmax,idsoft
      real factor
      character*4 caltype
      if(caltype.eq.'self') then
        if(iduv.lt.esmax) factor=(escrd(iduv)+escrd(iduv+1))/2.0e0
        if(iduv.eq.esmax) factor=escrd(esmax)
      endif
      if(caltype.eq.'intn') then
        idpick=0
        do 3331 cnt=1,numslv
          idpick=idpick+uvmax(cnt)
          if(iduv.le.idpick) goto 3339
3331    continue
3339    continue
        pti=cnt
        idpick=idpick-uvmax(pti)
        idsoft=uvsoft(pti)
        idmax=uvmax(pti)
        idpt=iduv-idpick
        if((idpt.lt.0).or.(idpt.gt.idmax)) call eng_stop('ecd')
        if(idpt.le.idsoft) factor=(uvcrd(iduv)+uvcrd(iduv+1))/2.0e0
        if((idpt.gt.idsoft).and.(idpt.lt.idmax)) then
          factor=sqrt(uvcrd(iduv)*uvcrd(iduv+1))
        endif
        if(idpt.eq.idmax) factor=uvcrd(iduv)
      endif
      return
      end subroutine
c
      end module
c
c
c  structural analysis performed separately from energetic analysis
#     include "stranal.f"
c
c
c  connection to the main routine of trajectory generation is done in
c   setparam for parameter setting and getconf for configuration reading
      subroutine enganal(stnum)
      use engmain, only: maxcnf,engdiv,skpcnf,inscnd
      use engproc, only: enginit,engclear,engconst,engstore
      use setconf, only: setparam,getconf
      use ptinsrt, only: refmc
      use stranal, only: spacedst
      integer stnum
      if(stnum.eq.1) call setparam
#ifndef trjctry
      if(mod(stnum,skpcnf).ne.0) return
#endif
      if(stnum.eq.skpcnf) call enginit
      if(mod(stnum,(maxcnf/engdiv)).eq.skpcnf) call engclear
      call getconf
#ifdef trjctry
      if(mod(stnum,skpcnf).ne.0) return
#endif
      if((stnum.eq.skpcnf).and.(inscnd.eq.3)) call refmc('init')
      call engconst(stnum)
      call spacedst(stnum)
      if(mod(stnum,(maxcnf/engdiv)).eq.0) call engstore(stnum)
      return
      end subroutine
c
c
#ifdef trjctry
      program trjmain
      use engmain, only: maxcnf
      use OUTname, only: opentrj,closetrj
      use mpiproc                                                      ! MPI
      integer stnum
      integer, parameter :: large=100000000
      call mpi_setup('init')                                           ! MPI
      call opentrj
      do 99999 stnum=1,large
        call enganal(stnum)
        if(stnum.eq.maxcnf) go to 1111
99999 continue
1111  call closetrj
      call mpi_setup('stop')                                           ! MPI
      stop
      end
#endif
