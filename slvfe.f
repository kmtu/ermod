      module sysvars
c
      implicit none
c
      character*5 :: clcond='merge'
c
      integer :: pecore=200,numprm=15
      integer :: numsln=10,numref=5,numdiv=-1
c
      character*3 :: peread='not',uvread='yes',slfslt='yes'
      character*4 :: zerosft='mxco',wgtfnform='harm'
      character*3 :: refmerge='yes',extsln='lin'
      character*3 :: wgtf2smpl='yes',slncor='not'
      character*3 :: normalize='yes',showdst='not'
      character*3 :: wrtzrsft='not',readwgtfl='yes'
c
      real :: inptemp=300.0e0                          ! Kelvin
      integer :: pickgr=3
      integer :: maxmesh=30000,large=500000,itrmax=100
      real :: error=1.0e-8,tiny=1.0e-10
c
      character(len=1024) :: wgtslnfl='soln/weight_soln'
      character(len=1024) :: wgtreffl='refs/weight_refs'
      character(len=1024) :: slndnspf='soln/engsln.'
      character(len=1024) :: slncorpf='soln/corsln.'
      character(len=1024) :: refdnspf='refs/engref.'
      character(len=1024) :: refcorpf='refs/corref.'
      character(*), parameter :: numbers='0123456789'
c
      integer prmmax,maxsln,maxref,numrun
      integer numslv,ermax
      real temp,kT,slfeng
c
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

      namelist /fevars/ clcond, pecore, numprm, numsln, numref, numdiv, 
     #                  peread, uvread, slfslt, zerosft, wgtfnform, refmerge, extsln,
     #                  wgtf2smpl, slncor, normalize, showdst, wrtzrsft, readwgtfl,
     #                  inptemp, pickgr, maxmesh, large, itrmax, error, tiny,
     #                  wgtslnfl, wgtreffl, slndnspf, slncorpf, refdnspf, refcorpf

      contains

      subroutine init_sysvars
      character(len=*), parameter :: parmfname = 'parameters_fe'
      integer, parameter :: iounit = 191
      integer :: iostat
      
      open(unit=iounit, file= parmfname, action='read', iostat=iostat)
      if(iostat /= 0) goto 99
      read(iounit, nml=fevars)
      close(iounit)

 99   if(numdiv == -1) numdiv = numsln
      
      end subroutine init_sysvars
c
      end module
c
c
c
      module sysread
      use sysvars, only: clcond,uvread,slfslt,slncor,tiny,
     #                   maxsln,maxref,numrun,
     #                   numslv,ermax,slfeng,nummol,
     #                   rdcrd,rddst,rddns,rdslc,rdcor,rdspec,
     #                   aveuv,uvene,blkuv,wgtsln,wgtref
      character*85 engfile(5)
      contains
c
      subroutine defcond
c
      use sysvars, only: peread,readwgtfl,wgtslnfl,wgtreffl,
     #                   numprm,prmmax,numsln,numref,numdiv,
     #                   inptemp,temp,kT,
     #                   pecore,maxmesh,large,
     #                   rduvmax,rduvcore,
     #                   chmpt,svgrp,svinf
c
      integer group,inft,prmcnt,iduv,i,j,k,m,pti,cnt
      real factor
      character*85 opnfile
c
      if((clcond.ne.'basic').and.(clcond.ne.'range')
     #                      .and.(clcond.ne.'merge')) then
        write(6,*) ' The clcond parameter is incrrect ' ; stop
      endif
c
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
781     format(' What is the energy correlation in solution?')
782     format(' What is the energy correlation in solution?')
783     format(' What is the energy density for insertion?')
784     format(' What is the energy density for insertion?')
785     format(' Which file describes energy-coordinate parameters?')
        maxsln=1 ; maxref=1 ; numrun=1
      endif
      if(clcond.eq.'merge') then
        maxsln=numsln ; maxref=numref ; numrun=numdiv
      endif
      if(clcond.eq.'basic') prmmax=1
      if(clcond.ne.'basic') prmmax=numprm
c
      if(clcond.ne.'merge') opnfile=engfile(1)
      if(clcond.eq.'merge') opnfile='soln/engsln.01'
      open(unit=71,file=opnfile,status='old')
      ermax=0 ; numslv=0 ; k=0
      do 7801 iduv=1,large
        read(71,*,END=7809) factor,pti,factor
        if(pti.ne.k) then
          numslv=numslv+1 ; k=pti
        endif
        ermax=ermax+1
7801  continue
7809  continue
      close(71)
c
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
c
      if(peread.ne.'yes') then
        cnt=ermax/numslv
        rduvmax(:)=cnt ; rduvcore(:)=pecore
      endif
      if(peread.eq.'yes') then
        open(unit=71,file=opnfile,status='old')
        k=0
        do 7811 iduv=1,ermax
          read(71,*) factor,pti,factor
          if(pti.ne.k) then
            rduvmax(pti)=1 ; k=pti
          else
            rduvmax(pti)=rduvmax(pti)+1
          endif
7811    continue
        close(71)
        do 7821 pti=1,numslv
          rduvcore(pti)=pecore
          if(clcond.ne.'merge') opnfile=engfile(5)
          if(clcond.eq.'merge') opnfile='soln/EcdInfo'
          open(unit=71,file=opnfile,status='old')
          read(71,*)        ! comment line
          read(71,*)        ! line for solute-self energy
          do 7822 i=1,large
            read(71,*,END=7829) k
            if(k.eq.pti) then
              backspace(71)
              read(71,*) m,(factor,j=1,9),rduvcore(pti)
              goto 7829
            endif
7822      continue
7829      continue
          close(71)
7821    continue
        k=0
        do 7825 pti=1,numslv
          k=k+rduvmax(pti)
7825    continue
        if(k.ne.ermax) then
          write(6,*) ' The file format is incorrect' ; stop
        endif
      endif
c
      if(ermax.gt.maxmesh) then
        write(6,*) ' The number of meshes is too large' ; stop
      endif
c
      if(clcond.eq.'basic') then
        write(6,*) ' How many data are grouped into one?'
        read(5,*) group
        write(6,*) ' How many large-energy meshes are merged ? (in %)'
        read(5,*) inft
        svgrp(1)=group ; svinf(1)=inft
        write(6,*) ' What is the temperature in Kelvin?'
        read(5,*) temp
      endif
      if(clcond.ne.'basic') then
        do 1001 prmcnt=1,prmmax                 ! inft in % (percentage)
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
          svgrp(prmcnt)=group ; svinf(prmcnt)=inft
1001    continue
        temp=inptemp
      endif
      kT=temp*8.314510e-3/4.184e0               ! kcal/mol
c
      if(uvread.eq.'yes') then
        if(clcond.ne.'merge') then
          write(6,891)
891       format(' What is average solute-solvent energy in solution?')
          read(5,*) (aveuv(pti), pti=1,numslv)
        endif
        if(clcond.eq.'merge') then
          opnfile='soln/aveuv.tt'
          open(unit=82,file=opnfile,status='old')
          do 8951 k=1,maxsln
            read(82,*) m,(uvene(pti,k), pti=1,numslv)
8951      continue
          close(82)
        endif
      endif
c
      do 8891 cnt=1,2
        if(cnt.eq.1) j=maxsln
        if(cnt.eq.2) j=maxref
        do 8892 i=1,j
          if(cnt.eq.1) wgtsln(i)=1.0e0
          if(cnt.eq.2) wgtref(i)=1.0e0
8892    continue
        if((clcond.eq.'merge').and.(readwgtfl.eq.'yes')) then
          if(cnt.eq.1) opnfile=wgtslnfl
          if(cnt.eq.2) opnfile=wgtreffl
          open(unit=81,file=opnfile,status='old')
          do 8893 i=1,j
            read(81,*) k,factor
            if(cnt.eq.1) wgtsln(i)=factor
            if(cnt.eq.2) wgtref(i)=factor
8893      continue
          close(81)
          factor=0.0e0
          do 8894 i=1,j
            if(cnt.eq.1) factor=factor+wgtsln(i)
            if(cnt.eq.2) factor=factor+wgtref(i)
8894      continue
          do 8895 i=1,j
            if(cnt.eq.1) wgtsln(i)=wgtsln(i)/factor
            if(cnt.eq.2) wgtref(i)=wgtref(i)/factor
8895      continue
        endif
8891  continue
c
      if(slfslt.eq.'yes') then
        if(clcond.ne.'merge') then
          write(6,881) ; read(5,*) slfeng
881       format(' What is the solute self-energy?')
        endif
        if(clcond.eq.'merge') then
          if(readwgtfl.eq.'not') then
            write(6,882) ; stop
882         format(' readwgtfl needs to be yes when slfslt is yes')
          endif
          slfeng=0.0e0
          open(unit=81,file=wgtreffl,status='old')
          do 8871 i=1,maxref
            read(81,*) k,factor,factor
            slfeng=slfeng+wgtref(i)*factor
8871      continue
          close(81)
        endif
      endif
c
      return
      end subroutine
c
c
      subroutine datread(cntrun)
c
      use sysvars, only: refmerge,
     #                   slndnspf,slncorpf,refdnspf,refcorpf,numbers
      integer cntrun,slnini,slnfin,refini,reffin,ecmin,ecmax
      integer iduv,iduvp,i,k,m,pti,cnt
      real factor,ampl
      character*85 opnfile
      character*2 suffnum
c
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
c
      rddst(:)=0.0e0
      if(slncor.eq.'yes') rdslc(:,:)=0.0e0
      if((cntrun.eq.1).or.(refmerge.eq.'not')) then
        rddns(:)=0.0e0 ; rdcor(:,:)=0.0e0
      endif
c
      do 7100 cnt=1,4
        if((cnt.eq.2).and.(slncor.ne.'yes')) goto 7100
        if((cnt.ge.3).and.(cntrun.gt.1)
     #               .and.(refmerge.eq.'yes')) goto 7899
c
        if(cnt.le.2) ecmin=slnini
        if(cnt.le.2) ecmax=slnfin
        if(cnt.ge.3) ecmin=refini
        if(cnt.ge.3) ecmax=reffin
c
        factor=0.0e0
        do 7108 i=ecmin,ecmax
          if(cnt.le.2) factor=factor+wgtsln(i)
          if(cnt.ge.3) factor=factor+wgtref(i)
7108    continue
        do 7109 i=ecmin,ecmax
          if(cnt.le.2) wgtsln(i)=wgtsln(i)/factor
          if(cnt.ge.3) wgtref(i)=wgtref(i)/factor
7109    continue
c
        do 7105 i=ecmin,ecmax
          if(clcond.ne.'merge') opnfile=engfile(cnt)
          if(clcond.eq.'merge') then
            m=i/10 ; k=i-10*m
            suffnum=numbers(m+1:m+1)//numbers(k+1:k+1)
            if(cnt.eq.1) opnfile=trim(slndnspf)//suffnum
            if(cnt.eq.2) opnfile=trim(slncorpf)//suffnum
            if(cnt.eq.3) opnfile=trim(refdnspf)//suffnum
            if(cnt.eq.4) opnfile=trim(refcorpf)//suffnum
          endif
          open(unit=71,file=opnfile,status='old')
          if((cnt.eq.1).or.(cnt.eq.3)) then
            k=0 ; m=0
          endif
          do 7101 iduv=1,ermax
            if((cnt.eq.1).or.(cnt.eq.3)) then
              read(71,*) rdcrd(iduv),pti,factor
              if(pti.ne.k) then
                if(factor.gt.tiny) then
                  write(6,*) ' Incorrect energy range with species ',pti
                  stop
                endif
                k=pti ; m=m+1
              endif
              if(cnt.eq.1) rddst(iduv)=rddst(iduv)+wgtsln(i)*factor
              if(cnt.eq.3) rddns(iduv)=rddns(iduv)+wgtref(i)*factor
              rdspec(iduv)=m
            endif
            if((cnt.eq.2).or.(cnt.eq.4)) then
              do 7102 iduvp=1,ermax
                read(71,*) factor
                if(cnt.eq.2) then
                  rdslc(iduvp,iduv)=rdslc(iduvp,iduv)+wgtsln(i)*factor
                endif
                if(cnt.eq.4) then
                  rdcor(iduvp,iduv)=rdcor(iduvp,iduv)+wgtref(i)*factor
                endif
7102          continue
            endif
7101      continue
          close(71)
7105    continue
7899    continue
7100  continue
c
      if(cntrun.eq.1) write(6,*)
      do 7120 pti=1,numslv
        factor=0.0e0 ; ampl=0.0e0
        do 7121 iduv=1,ermax
          if(rdspec(iduv).eq.pti) then
            factor=factor+rddst(iduv) ; ampl=ampl+rddns(iduv)
          endif
7121    continue
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
718     format('  Incorrect normalization at ',i4)
719     format('  Number of the ',i3,'-th solvent  = ',i12)
7120  continue
      if(cntrun.eq.1) write(6,*)
c
      if((uvread.eq.'yes').and.(clcond.eq.'merge')) then
        do 7231 pti=1,numslv
          factor=0.0e0
          do 7232 i=slnini,slnfin
            factor=factor+wgtsln(i)*uvene(pti,i)
7232      continue
          aveuv(pti)=factor
7231    continue
        factor=0.0e0
        do 7233 pti=1,numslv
          blkuv(pti,cntrun)=aveuv(pti)
          factor=factor+aveuv(pti)
7233    continue
        if(slfslt.eq.'yes') factor=factor+slfeng
        blkuv(0,cntrun)=factor
      endif
c
      return
      end subroutine
c
      end module
c
c
c
      module sfecalc
      use sysvars, only: zerosft,wgtfnform,extsln,slncor,
     #                   numslv,ermax,nummol,kT,itrmax,error,tiny
      integer, dimension(:), allocatable :: idrduv,uvmax
      real, dimension(:),    allocatable :: uvcrd,edist,edens
      real, dimension(:,:),  allocatable :: edscr,ecorr
      integer, dimension(:), allocatable :: uvspec
      real, dimension(:),    allocatable :: slncv,inscv,sdrcv
      real, dimension(:),    allocatable :: zrsln,zrref,zrsdr
      integer gemax
      contains
c
      subroutine chmpot(prmcnt,cntrun)
c
      use sysvars, only: uvread,slfslt,
     #                   normalize,showdst,wrtzrsft,slfeng,
     #                   rduvmax,rduvcore,
     #                   rdcrd,rddst,rddns,rdslc,rdcor,rdspec,
     #                   chmpt,aveuv,svgrp,svinf
c
      integer prmcnt,cntrun,group,inft
      integer iduv,iduvp,pti,cnt,j,k,m
      real factor,ampl,slvfe,uvpot,lcsln,lcref
      integer, dimension(:), allocatable :: gpnum
c
      group=svgrp(prmcnt) ; inft=svinf(prmcnt)
c
      allocate( idrduv(ermax),uvmax(numslv) )
c
      gemax=0
      do 1701 pti=1,numslv
        k=rduvcore(pti)*inft/100
        m=(rduvmax(pti)-k)/group
        uvmax(pti)=m ; gemax=gemax+m
        cnt=1 ; j=0
        if(pti.gt.1) then
          do 1702 k=1,pti-1
            cnt=cnt+rduvmax(k) ; j=j+uvmax(k)
1702      continue
        endif
        do 1703 iduv=cnt,cnt+rduvmax(pti)-1
          k=(iduv-cnt)/group+1
          if(k.ge.m) k=m
          idrduv(iduv)=k+j
1703    continue
1701  continue
c
      allocate( uvcrd(gemax),edist(gemax),edens(gemax),uvspec(gemax) )
      uvcrd(:)=0.0e0 ; edist(:)=0.0e0 ; edens(:)=0.0e0
c
      if(slncor.eq.'yes') then
        allocate( edscr(gemax,gemax) )
        edscr(:,:)=0.0e0
      endif
      allocate( ecorr(gemax,gemax) )
      ecorr(:,:)=0.0e0
c
      allocate( slncv(gemax),inscv(gemax),zrsln(numslv),zrref(numslv) )
      slncv(:)=0.0e0 ; inscv(:)=0.0e0 ; zrsln(:)=0.0e0 ; zrref(:)=0.0e0
      if(slncor.eq.'yes') then
        allocate( sdrcv(gemax),zrsdr(numslv) )
        sdrcv(:)=0.0e0 ; zrsdr(:)=0.0e0
      endif
c
      allocate( gpnum(gemax) )
      gpnum(:)=0
      do 1111 iduv=1,ermax
        k=idrduv(iduv)
        uvcrd(k)=uvcrd(k)+rdcrd(iduv) ; gpnum(k)=gpnum(k)+1
        uvspec(k)=rdspec(iduv)
1111  continue
      do 1112 k=1,gemax
        if(gpnum(k).gt.0) uvcrd(k)=uvcrd(k)/real(gpnum(k))
1112  continue
      cnt=rduvmax(1) ; k=uvmax(1)
      do 1113 pti=1,numslv
        if(pti.gt.1) cnt=cnt+rduvmax(pti)
        if(pti.gt.1) k=k+uvmax(pti)
        uvcrd(k)=rdcrd(cnt)
1113  continue
      deallocate( gpnum )
c
      do 1115 cnt=1,2
        do 1116 iduv=1,ermax
          k=idrduv(iduv)
          if(cnt.eq.1) edist(k)=edist(k)+rddst(iduv)
          if(cnt.eq.2) edens(k)=edens(k)+rddns(iduv)
1116    continue
        if((cnt.eq.1).and.(slncor.ne.'yes')) goto 1115
        do 1117 iduv=1,ermax
         do 1118 iduvp=1,ermax
           k=idrduv(iduv) ; m=idrduv(iduvp)
           if(cnt.eq.1) edscr(m,k)=edscr(m,k)+rdslc(iduvp,iduv)
           if(cnt.eq.2) ecorr(m,k)=ecorr(m,k)+rdcor(iduvp,iduv)
1118     continue
1117    continue
1115  continue
c
      if(normalize.eq.'yes') call distnorm
      if(showdst.eq.'yes')   call distshow
c
      call getslncv
      call getinscv
c
      do 5001 pti=1,numslv
        uvpot=0.0e0 ; slvfe=0.0e0
        do 5002 iduv=1,gemax
          if(uvspec(iduv).eq.pti) then
            if((edist(iduv).le.tiny).and.(edens(iduv).le.tiny)) goto 5009
            uvpot=uvpot+uvcrd(iduv)*edist(iduv)
            slvfe=slvfe-kT*(edist(iduv)-edens(iduv))
            factor=-(slncv(iduv)+zrsln(pti)+uvcrd(iduv))
            if((slncor.eq.'yes').and.
     #         (edist(iduv).gt.tiny).and.(edens(iduv).le.tiny)) then
              ampl=factor*edens(iduv)/edist(iduv)
              factor=ampl-(zrsln(pti)+uvcrd(iduv))
     #                   *(1.0e0-edens(iduv)/edist(iduv))
            endif
            slvfe=slvfe+factor*edist(iduv)
            do 5005 cnt=1,2
              if(cnt.eq.1) lcsln=pyhnc(slncv(iduv),cnt)
              if(cnt.eq.2) lcref=pyhnc(inscv(iduv),cnt)
              if((slncor.eq.'yes').and.(cnt.eq.1).and.
     #           (edist(iduv).gt.tiny).and.(edens(iduv).le.tiny)) then
                lcsln=pyhnc(sdrcv(iduv)+zrsdr(pti),3)
              endif
5005        continue
            ampl=sfewgt(edist(iduv),edens(iduv))
            factor=ampl*lcsln+(1.0e0-ampl)*lcref
            slvfe=slvfe+kT*factor*(edist(iduv)-edens(iduv))
5009        continue
          endif
5002    continue
        if((uvread.ne.'yes').and.(prmcnt.eq.1)) aveuv(pti)=uvpot
        chmpt(pti,prmcnt,cntrun)=slvfe+aveuv(pti)
5001  continue
c
      slvfe=0.0e0
      do 5199 pti=1,numslv
        slvfe=slvfe+chmpt(pti,prmcnt,cntrun)
5199  continue
      if(slfslt.eq.'yes') slvfe=slvfe+slfeng
      chmpt(0,prmcnt,cntrun)=slvfe
c
      if(wrtzrsft.eq.'yes') then
        do 5181 cnt=1,3
          if(cnt.eq.1) write(6,381) (zrsln(pti),pti=1,numslv)
          if(cnt.eq.2) write(6,382) (zrref(pti),pti=1,numslv)
          if((cnt.eq.3).and.(slncor.eq.'yes')) then
            write(6,383) (zrsdr(pti),pti=1,numslv)
          endif
5181    continue
      endif
381   format('  Zero shift for solution             = ', 999999f12.4)
382   format('  Zero shift for reference solvent    = ', 999999f12.4)
383   format('  Zero shift for solution correlation = ', 999999f12.4)
c
      deallocate( slncv,inscv,zrsln,zrref )
      if(slncor.eq.'yes') deallocate( sdrcv,zrsdr,edscr )
      deallocate( uvcrd,edist,edens,ecorr,uvspec,idrduv,uvmax )
c
      return
      end subroutine
c
c
      subroutine getslncv
c
      integer iduv,iduvp,pti,j,k,m
      real factor,ampl,lcsln,lcref
      real, dimension(:), allocatable :: work
c
      do 3111 iduv=1,gemax
        if((edist(iduv).gt.tiny).and.(edens(iduv).gt.tiny)) then
          factor=edist(iduv)/edens(iduv)
          slncv(iduv)=-kT*log(factor)-uvcrd(iduv)
        endif
3111  continue
c
      do 3211 iduv=1,gemax
        if((edist(iduv).le.tiny).or.(edens(iduv).le.tiny)) then
          if(edist(iduv).le.tiny) then
            slncv(iduv)=0.0e0 ; goto 3211
          endif
          pti=uvspec(iduv)
          m=1 ; k=gemax
          do 3212 iduvp=1,iduv-1
            if((uvspec(iduvp).eq.pti).and.(m.lt.iduvp).and.
     #         (edist(iduvp).gt.tiny).and.(edens(iduvp).gt.tiny)) then
              m=iduvp
            endif
3212      continue
          do 3213 iduvp=gemax,iduv+1,-1
            if((uvspec(iduvp).eq.pti).and.(k.gt.iduvp).and.
     #         (edist(iduvp).gt.tiny).and.(edens(iduvp).gt.tiny)) then
              k=iduvp
            endif
3213      continue
c
          if(extsln.eq.'sim') then
            if(abs(m-iduv).lt.abs(k-iduv)) factor=slncv(m)
            if(abs(m-iduv).gt.abs(k-iduv)) factor=slncv(k)
            if(abs(m-iduv).eq.abs(k-iduv)) then
              factor=(slncv(m)+slncv(k))/2.0e0
            endif
          else
            j=k ; if(abs(m-iduv).lt.abs(k-iduv)) j=m
            allocate( work(gemax) ) ; work(:)=0.0e0
            do 3221 iduvp=1,gemax
              if((uvspec(iduvp).eq.pti).and.
     #           (edist(iduvp).gt.tiny).and.(edens(iduvp).gt.tiny)) then
                if(iduvp.eq.iduv) then
                  write(6,*) ' A bug in program or data' ; stop
                endif
                factor=uvcrd(iduvp)-uvcrd(j)
                if(iduvp.lt.iduv) then
                  factor=-factor-2.0e0*(uvcrd(j)-uvcrd(iduv))
                endif
                ampl=wgtdst(iduvp,1,'extsl',wgtfnform)
                work(iduvp)=exp(-factor/kT)*ampl
              endif
3221        continue
            factor=0.0e0
            do 3222 iduvp=1,gemax
              if(work(iduvp).gt.tiny) factor=factor+work(iduvp)
3222        continue
            do 3223 iduvp=1,gemax
              work(iduvp)=work(iduvp)/factor
3223        continue
            factor=0.0e0 ; ampl=0.0e0 ; lcsln=0.0e0 ; lcref=0.0e0
            do 3224 iduvp=1,gemax
              if(work(iduvp).gt.tiny) then
                factor=factor+work(iduvp)*uvcrd(iduvp)
                ampl=ampl+work(iduvp)*uvcrd(iduvp)*uvcrd(iduvp)
                lcsln=lcsln+work(iduvp)*slncv(iduvp)
                lcref=lcref+work(iduvp)*uvcrd(iduvp)*slncv(iduvp)
              endif
3224        continue
            work(1)=(ampl*lcsln-factor*lcref)/(ampl-factor*factor)
            work(2)=(lcref-factor*lcsln)/(ampl-factor*factor)
            factor=work(1)+work(2)*uvcrd(iduv)
            deallocate( work )
          endif
          slncv(iduv)=factor
        endif
3211  continue
c
      do 3231 pti=1,numslv
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
            write(6,*) ' zerosft not properly set ' ; stop
        end select
        do 3232 iduv=1,gemax
          if(uvspec(iduv).eq.pti) slncv(iduv)=slncv(iduv)-factor
3232    continue
        zrsln(pti)=factor
3231  continue
c
      return
      end subroutine
c
c
      subroutine getinscv
c
      integer iduv,iduvp,pti,cnt,wrksz,k
      real factor,ampl,lcsln,lcref
      real, dimension(:),   allocatable :: work,egnvl,zerouv
      real, dimension(:,:), allocatable :: edmcr
c
      do 1151 cnt=1,2
        if((cnt.eq.1).and.(slncor.ne.'yes')) goto 1151
        wrksz=gemax*gemax
        allocate( work(wrksz),egnvl(gemax),edmcr(gemax,gemax) )
        if(cnt.eq.1) edmcr(:,:)=edscr(:,:)
        if(cnt.eq.2) edmcr(:,:)=ecorr(:,:)
        do 1152 iduv=1,gemax
         do 1153 iduvp=1,gemax
           if(cnt.eq.1) then
             factor=edist(iduvp) ; ampl=edist(iduv)
           endif
           if(cnt.eq.2) then
             factor=edens(iduvp) ; ampl=edens(iduv)
           endif
           lcref=edmcr(iduvp,iduv)-factor*ampl
           if((factor.le.tiny).or.(ampl.le.tiny)) then
             if(iduv.eq.iduvp) lcref=1.0e0
             if(iduv.ne.iduvp) lcref=0.0e0
           endif
           edmcr(iduvp,iduv)=lcref
1153     continue
1152    continue
c
        call DSYEV('V','U',gemax,edmcr,gemax,egnvl,work,wrksz,k)
        pti=numslv+1
        do 1161 iduv=pti,gemax
          factor=0.0e0
          do 1162 iduvp=1,gemax
            if(cnt.eq.1) ampl=edist(iduvp)
            if(cnt.eq.2) ampl=edens(iduvp)
            if(ampl.gt.tiny) factor=factor+(edist(iduvp)-edens(iduvp))
     #                                    *edmcr(iduvp,iduv)
1162      continue
          work(iduv)=factor/egnvl(iduv)
1161    continue
        do 1163 iduv=1,gemax
          factor=0.0e0
          do 1164 iduvp=pti,gemax
            factor=factor+edmcr(iduv,iduvp)*work(iduvp)
1164      continue
          if(cnt.eq.1) sdrcv(iduv)=-kT*factor
          if(cnt.eq.2) inscv(iduv)=-kT*factor
1163    continue
        deallocate( work,egnvl,edmcr )
c
        allocate( zerouv(numslv) )
        do 1182 pti=1,numslv
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
1182    continue
c
        do 1181 iduv=1,gemax
          if(cnt.eq.1) ampl=edist(iduv)
          if(cnt.eq.2) ampl=edens(iduv)
          if(ampl.gt.tiny) then
            factor=-kT*(edist(iduv)-edens(iduv))/ampl
            if(cnt.eq.1) lcref=factor-sdrcv(iduv)
            if(cnt.eq.2) lcref=factor-inscv(iduv)
          endif
          if(ampl.le.tiny) lcref=0.0e0
          if(cnt.eq.1) sdrcv(iduv)=lcref
          if(cnt.eq.2) inscv(iduv)=lcref
1181    continue
c
        do 1171 pti=1,numslv
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
              write(6,*) ' zerosft not properly set ' ; stop
          end select
          do 1172 iduv=1,gemax
            if(uvspec(iduv).eq.pti) then
              if(cnt.eq.1) sdrcv(iduv)=sdrcv(iduv)-factor
              if(cnt.eq.2) inscv(iduv)=inscv(iduv)-factor
            endif
1172      continue
          if(cnt.eq.1) zrsdr(pti)=factor
          if(cnt.eq.2) zrref(pti)=factor
1171    continue
c
        deallocate( zerouv )
1151  continue
c
      return
      end subroutine
c
c
      real function wgtmxco(pti)
      integer pti
      real numpt
      numpt=nummol(pti)
      wgtmxco=1.0e0/numpt
      return
      end function
c
c
      real function cvfcen(pti,cnt,systype,wgttype,engtype)
      integer pti,cnt,iduv,errtag
      real factor,ampl,minuv,cvfnc
      real, dimension(:), allocatable :: weight
      character*5 systype
      character*4 wgttype
      character*3 engtype
      allocate( weight(gemax) )
      weight(:)=0.0e0
      do 3151 iduv=1,gemax
        if(uvspec(iduv).eq.pti) then
          weight(iduv)=wgtdst(iduv,cnt,systype,wgttype)
        endif
3151  continue
      if(engtype.eq.'yes') then
        minuv=abs(uvcrd(1))
        do 3152 iduv=1,gemax
          if((uvspec(iduv).eq.pti).and.(weight(iduv).gt.tiny)) then
            if(abs(uvcrd(iduv)).lt.minuv) minuv=abs(uvcrd(iduv))
          endif
3152    continue
        do 3153 iduv=1,gemax
        if(uvspec(iduv).eq.pti) then
          ampl=abs(uvcrd(iduv))-minuv
c         ampl=nummol(pti)*(abs(uvcrd(iduv))-minuv)
          weight(iduv)=exp(-ampl/kT)*weight(iduv)
        endif
3153    continue
      endif
      factor=0.0e0 ; ampl=0.0e0 ; errtag=0
      do 3251 iduv=1,gemax
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
          if(errtag.ne.0) then
            write(6,*) ' Bug in the program' ; stop
          endif
          factor=factor+cvfnc*weight(iduv)
          ampl=ampl+weight(iduv)
        endif
3251  continue
      if(ampl.gt.tiny) factor=factor/ampl
      if(ampl.le.tiny) factor=0.0e0
      cvfcen=factor
      deallocate( weight )
      return
      end function
c
c
      integer function zeroec(pti,cnt)
      integer iduv,cnt,pti,k
      real factor,ampl,lcsln,lcref
      do 1170 iduv=1,gemax-1
        if(uvspec(iduv).eq.pti) then
          if((uvcrd(iduv).le.0).and.(uvcrd(iduv+1).ge.0)) then
            factor=abs(uvcrd(iduv)) ; ampl=uvcrd(iduv+1)
            if(ampl.gt.factor+tiny) k=iduv
            if(factor.gt.ampl+tiny) k=iduv+1
            if(abs(factor-ampl).le.tiny) then
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
1170  continue
      zeroec=k
      return
      end function
c
c
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
            if(factor.gt.tiny) wght=sqrt(factor)
          case default              ! corresponding to wgttype = 'harm'
            factor=fsln+fref
            if(factor.gt.tiny) wght=fsln*fref/factor
        end select
      endif
      if(errtag.ne.0) then
        write(6,*) ' Bug in the program' ; stop
      endif
      wgtdst=wght
      return
      end function
c
c
      real function sfewgt(fsln,fref)
      real fsln,fref,wght,factor
      if(fsln.ge.fref) wght=1.0e0
      if(fsln.lt.fref) then
        factor=(fsln-fref)/(fsln+fref)
        wght=1.0e0-factor*factor
      endif
      sfewgt=wght
      return
      end function
c
c
      real function pyhnc(indpmf,cnt)
      real intg,indpmf,factor
      integer cnt
      factor=indpmf/kT
      if(factor.lt.-tiny) then
        if(cnt.eq.1) intg=factor+factor/(exp(-factor)-1.0e0)
        if(cnt.eq.2) intg=(log(1.0e0-factor))*(1.0e0/factor-1.0e0)
        intg=intg+1.0e0
      endif
      if(factor.ge.-tiny) intg=factor/2.0e0
      if(cnt.eq.3) then
        if(factor.ge.tiny) then
          intg=1.0e0-(log(1.0e0+factor))*(1.0e0/factor+1.0e0)
        endif
        if(factor.lt.tiny) intg=-factor/2.0e0
      endif
      pyhnc=intg
      return
      end function
c
c
      subroutine distnorm
      real, dimension(:), allocatable :: correc
      integer iduv,iduvp,pti,cnt,itrcnt
      real factor,ampl,lcsln,lcref,errtmp
      allocate( correc(gemax) )
      do 2100 cnt=1,2
        do 2101 pti=1,numslv
          factor=0.0e0
          do 2102 iduv=1,gemax
            if(uvspec(iduv).eq.pti) then
              if(cnt.eq.1) factor=factor+edist(iduv)
              if(cnt.eq.2) factor=factor+edens(iduv)
            endif
2102      continue
          if(factor.gt.tiny) factor=nummol(pti)/factor
          if(factor.le.tiny) factor=0.0e0
          do 2103 iduv=1,gemax
            if(uvspec(iduv).eq.pti) then
              if(cnt.eq.1) edist(iduv)=factor*edist(iduv)
              if(cnt.eq.2) edens(iduv)=factor*edens(iduv)
            endif
2103      continue
2101    continue
c
        if((cnt.eq.1).and.(slncor.ne.'yes')) goto 2100
        errtmp=error+1.0e0 ; itrcnt=0
        do 2111 iduv=1,gemax
          correc(iduv)=1.0e0
2111    continue
        do 2121 while((errtmp.gt.error).and.(itrcnt.le.itrmax))
          do 2125 iduv=1,gemax
            lcsln=0.0e0
            do 2126 pti=1,numslv
              ampl=0.0e0
              do 2127 iduvp=1,gemax
                if(uvspec(iduvp).eq.pti) then
                  if(cnt.eq.1) lcref=edscr(iduvp,iduv)
                  if(cnt.eq.2) lcref=ecorr(iduvp,iduv)
                  ampl=ampl+correc(iduvp)*lcref
                endif
2127          continue
              if(ampl.gt.tiny) lcsln=lcsln+nummol(pti)/ampl
2126        continue
            lcsln=lcsln/real(numslv)
            if(cnt.eq.1) correc(iduv)=lcsln*edist(iduv)
            if(cnt.eq.2) correc(iduv)=lcsln*edens(iduv)
2125      continue
          do 2131 iduv=1,gemax
           do 2132 iduvp=1,gemax
             ampl=correc(iduv)*correc(iduvp)
             if(cnt.eq.1) edscr(iduvp,iduv)=ampl*edscr(iduvp,iduv)
             if(cnt.eq.2) ecorr(iduvp,iduv)=ampl*ecorr(iduvp,iduv)
2132       continue
2131      continue
          errtmp=0.0e0 ; itrcnt=itrcnt+1
          do 2135 iduv=1,gemax
            if(cnt.eq.1) ampl=edist(iduv)
            if(cnt.eq.2) ampl=edens(iduv)
            if(ampl.gt.tiny) then
              factor=abs(correc(iduv)-1.0e0)
              if(factor.gt.errtmp) errtmp=factor
            endif
2135      continue
          if(itrcnt.eq.itrmax) then
            write(6,*) ' The optimzation of the correlation matrix'
            write(6,*) '  did not converge with an error of ',errtmp
            stop
          endif
2121    continue
2100  continue
      deallocate( correc )
      return
      end subroutine
c
c
      subroutine distshow
      integer iduv,iduvp,pti,cnt,ecmin,ecmax,i,k
      real factor
      write(6,*)
      do 2150 pti=1,numslv
        do 2151 cnt=1,2
          if((cnt.eq.1).and.(numslv.gt.1)) write(6,211) pti
          if((cnt.eq.2).and.(numslv.gt.1)) write(6,212) pti
          if((cnt.eq.1).and.(numslv.eq.1)) write(6,*) 'SOLUTION'
          if((cnt.eq.2).and.(numslv.eq.1)) write(6,*) 'INSERTION'
211       format('SOLUTION for',i4,'-th species')
212       format('INSERTION for',i4,'-th species')
          ecmin=gemax ; ecmax=1
          do 2152 iduv=1,gemax
            if(uvspec(iduv).eq.pti) then
              i=0
              if((cnt.eq.1).and.(edist(iduv).gt.tiny)) i=1
              if((cnt.eq.2).and.(edens(iduv).gt.tiny)) i=1
              if(i.eq.1) then
                if(ecmin.gt.iduv) ecmin=iduv
                if(ecmax.lt.iduv) ecmax=iduv
              endif
            endif
2152      continue
          k=0 ; factor=0.0e0
          do 2153 iduv=ecmin,ecmax
            i=0
            if((cnt.eq.1).and.(edist(iduv).gt.tiny)) i=1
            if((cnt.eq.2).and.(edens(iduv).gt.tiny)) i=1
            if(i.eq.1) then
              k=k+1
              if(cnt.eq.1) factor=factor+edist(iduv)
              if(cnt.eq.2) factor=factor+edens(iduv)
            endif
            if((cnt.eq.2).and.(edens(iduv).le.tiny)) then
              write(6,215) iduv,uvcrd(iduv)
215           format('     No sampling at ',i5,' with energy ',g14.6)
            endif
2153      continue
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
2151    continue
2150  continue
      return
      end subroutine
c
      end module
c
c
c
      module opwrite
      use sysvars, only: clcond,uvread,slfslt,
     #                   prmmax,numrun,numslv,pickgr,slfeng,
     #                   chmpt,aveuv,blkuv,svgrp,svinf
      integer grref
      contains
c
      subroutine wrtresl
c
      integer prmcnt,pti,k
      real factor,valcp
c
      if(slfslt.eq.'yes') write(6,321) slfeng
321   format('  Self-energy of the solute   =   ',f12.4, '  kcal/mol')
c
      if(clcond.ne.'merge') then
        write(6,*)
        if(numslv.gt.1) then
          if(uvread.ne.'yes') write(6,331) (aveuv(pti), pti=1,numslv)
        endif
        factor=0.0e0
        do 3311 pti=1,numslv
          factor=factor+aveuv(pti)
3311    continue
        if(slfslt.eq.'yes') factor=factor+slfeng
        write(6,332) factor
331     format('  Solute-solvent energy       =   ', 999999f12.4)
332     format('  Total solvation energy      =   ',f12.4, '  kcal/mol')
      endif
c
      if(clcond.eq.'basic') then
        if(numslv.gt.1) write(6,351) (chmpt(pti,1,1), pti=1,numslv)
        write(6,352) chmpt(0,1,1)
351     format('  Solvation free energy       =   ', 999999f12.4)
352     format('  Total solvation free energy =   ',f12.4, '  kcal/mol')
      endif
c
      if(clcond.ne.'basic') then
        do 9975 prmcnt=1,prmmax
          if(svgrp(prmcnt).eq.pickgr) then
            grref=prmcnt ; goto 9976
          endif
9975    continue
9976    continue
      endif
c
      if(clcond.eq.'range') then
        if(numslv.eq.1) k=0
        if(numslv.gt.1) k=numslv
        do 5195 pti=0,k
          write(6,*)
          if(pti.eq.0) write(6,591)
          if(pti.ne.0) write(6,592) pti
591       format('               chemical potential     difference')
592       format('               ',i3,'-th component       difference')
          factor=chmpt(pti,grref,1)
          do 5196 prmcnt=1,prmmax
            valcp=chmpt(pti,prmcnt,1)
            write(6,661) svgrp(prmcnt),svinf(prmcnt),valcp,valcp-factor
5196      continue
5195    continue
661     format(i4,i7,f17.5,f18.5)
      endif
c
      if(clcond.eq.'merge') call wrtmerge
c
      return
      end subroutine
c
c
      subroutine wrtmerge
c
      use sysvars, only: large
      integer prmcnt,cntrun,group,inft,pti,i,j,k,m
      real avecp,stdcp,avcp0,factor,slvfe,shcp(large)
      real, dimension(:,:), allocatable :: wrtdata
c
      allocate( wrtdata(0:numslv,numrun) )
      if(uvread.eq.'yes') then
        do 8811 pti=0,numslv
          do 8812 cntrun=1,numrun
            wrtdata(pti,cntrun)=blkuv(pti,cntrun)
8812      continue
8811    continue
        write(6,*) ; write(6,*) ; write(6,773)
        call wrtcumu(wrtdata) ; write(6,*)
      endif
773   format(' cumulative average & 95% error for solvation energy')
c
      do 9983 pti=0,numslv
        if((numslv.eq.1).and.(pti.ne.0)) goto 9984
        avcp0=0.0e0
        do 9978 cntrun=1,numrun
          avcp0=avcp0+chmpt(pti,grref,cntrun)
9978    continue
        avcp0=avcp0/real(numrun)
        do 9985 prmcnt=1,prmmax
          group=svgrp(prmcnt)
          inft=svinf(prmcnt)
          avecp=0.0e0
          stdcp=0.0e0
          do 9986 cntrun=1,numrun
            slvfe=chmpt(pti,prmcnt,cntrun)
            avecp=avecp+slvfe
            stdcp=stdcp+slvfe*slvfe
9986      continue
          factor=real(numrun)
          avecp=avecp/factor
          stdcp=sqrt(factor/(factor-1.0e0))
     #         *sqrt(stdcp/factor-avecp*avecp)
          stdcp=2.0e0*stdcp/sqrt(factor)
          if(prmcnt.eq.1) then
            write(6,*)
            if(pti.eq.0) write(6,670)
670         format(' group  inft  solvation free energy     error',
     #             '          difference')
            if(numslv.gt.1) then
              if(pti.eq.0) write(6,661)
              if(pti.ge.1) write(6,662) pti
661           format('  total solvation free energy')
662           format('  contribution from ',i2,'-th solvent component')
            endif
          endif
          write(6,671) group,inft,avecp,stdcp,avecp-avcp0
671       format(i4,i7,f17.5,2f18.5)
9985    continue
9984    continue
9983  continue
c
      write(6,*) ; write(6,*)
      do 9981 pti=0,numslv
        if((numslv.eq.1).and.(pti.ne.0)) goto 9982
        do 9988 prmcnt=1,prmmax
          group=svgrp(prmcnt)
          inft=svinf(prmcnt)
          do 9987 cntrun=1,numrun
            shcp(cntrun)=chmpt(pti,prmcnt,cntrun)
9987      continue
          if(prmcnt.eq.1) then
            if(numslv.eq.1) write(6,667)
            if(numslv.gt.1) then
              if(pti.eq.0) write(6,681)
              if(pti.ge.1) then
                write(6,*)
                write(6,682) pti
              endif
            endif
667         format(' group  inft   Estimated free energy (kcal/mol)')
681         format(' group  inft   Estimated free energy:',
     #             ' total (kcal/mol)')
682         format(' group  inft   Estimated free energy:',
     #               i2,'-th solvent contribution (kcal/mol)')
          endif
          k=(numrun-1)/5
          if(k.eq.0) write(6,121) group,inft,(shcp(m), m=1,numrun)
          if(k.ge.1) then
            write(6,125) group,inft,(shcp(m), m=1,5)
            if(k.gt.1) then
              do i=1,k-1
                write(6,126) (shcp(5*i+m), m=1,5)
              enddo
            endif
            j=numrun-5*k
            if(j.gt.0) write(6,127) (shcp(m), m=5*k+1,numrun)
          endif
121       format(i4,i7,999999f13.4)
125       format(i4,i7,5f13.4)
126       format('           ',5f13.4)
127       format('           ',999999f13.4)
9988    continue
9982    continue
9981  continue
c
      do 9911 pti=0,numslv
        do 9912 cntrun=1,numrun
          wrtdata(pti,cntrun)=chmpt(pti,grref,cntrun)
9912    continue
9911  continue
      write(6,*) ; write(6,*) ; write(6,770)
770   format(' cumulative average & 95% error ',
     #       'for solvation free energy')
      call wrtcumu(wrtdata)
      deallocate( wrtdata )
c
      return
      end subroutine
c
c
      subroutine wrtcumu(wrtdata)
      use sysvars, only: large,tiny
      integer cntrun,pti
      real avecp,factor,slvfe,shcp(large),wrtdata(0:numslv,numrun)
      real, dimension(:),   allocatable :: runcp,runer
      allocate( runcp(0:numslv),runer(0:numslv) )
      if(numslv.eq.2) write(6,771)
771   format('              total             1st component',
     #       '         2nd component')
      do 7773 pti=0,numslv
        runcp(pti)=0.0e0 ; runer(pti)=0.0e0
7773  continue
      do 7771 cntrun=1,numrun
        factor=real(cntrun)
        do 7774 pti=0,numslv
          slvfe=wrtdata(pti,cntrun)
          runcp(pti)=runcp(pti)+slvfe
          runer(pti)=runer(pti)+slvfe*slvfe
          avecp=runcp(pti)/factor ; shcp(2*pti+1)=avecp
          if(cntrun.gt.1) then
            slvfe=runer(pti)/factor-avecp*avecp
            if(slvfe.lt.tiny) shcp(2*pti+2)=0.0e0
            if(slvfe.ge.tiny) shcp(2*pti+2)=(2.0e0/sqrt(factor))
     #                          *sqrt(factor/(factor-1.0e0))*sqrt(slvfe)
          endif
7774    continue
        if(cntrun.eq.1) then
          do 7777 pti=0,numslv
            shcp(pti+1)=shcp(2*pti+1)
7777      continue
          if(numslv.eq.1) write(6,772) cntrun,shcp(1)
          if(numslv.gt.1) write(6,773) cntrun,shcp(1),
     #                                 (shcp(pti), pti=2,numslv+1)
        endif
        if(cntrun.gt.1) then
          if(numslv.eq.1) write(6,774) cntrun,(shcp(pti), pti=1,2)
          if(numslv.gt.1) write(6,775) cntrun,
     #                                 (shcp(pti), pti=1,2*numslv+2)
        endif
7771  continue
772   format(i3,f11.4)
773   format(i3,f11.4,999999f22.4)
774   format(i3,2f11.4)
775   format(i3,999999f11.4)
      deallocate( runcp,runer )
      return
      end subroutine
c
      end module
c
c
c
      program sfemain
      use sysvars, only: numrun,prmmax,init_sysvars
      use sysread, only: defcond,datread
      use sfecalc, only: chmpot
      use opwrite, only: wrtresl
      integer cntrun,prmcnt
      call init_sysvars
      call defcond
      do 99999 cntrun=1,numrun
        call datread(cntrun)
        do 88888 prmcnt=1,prmmax
          call chmpot(prmcnt,cntrun)
88888   continue
99999 continue
      call wrtresl
      stop
      end
