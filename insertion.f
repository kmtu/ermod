      module ptinsrt
c
c  test particle insertion of the solute
c
      real, save :: unrn
c
c  single-solute trajectrory file           used only when slttype = 3
      character(*), parameter :: slttrj='SltConf'    ! solute filename
      integer, parameter :: slcnf=31                 ! solute file ID
      character(*), parameter :: sltwgt='SltWght'    ! solute weight filename
      integer, parameter :: swinf=32                 ! solute weight ID
c
c  insertion against reference structure    used only when inscnd = 3
c   refmlid : superposition reference among solvent species
c             --- 0 : not reference
c                 1 : reference solvent  2 : reference solute
c             value set in subroutine setparam
c   refsatm : specification of the reference site
c             --- 0 : not reference
c                 1 : reference solvent  2 : reference solute
c   refspos : coordiantes of interaction site for reference structure
c   sltcen : coordinate of the solute center
c   sltqrn : quarternion for the solute orientation
c   movmax : number of Monte Carlo moves
c   trmax : maximum of translational Monte Carlo move
c   agmax : maximum of orientational Monte Carlo move
      integer, dimension(:,:), allocatable, save :: refsatm
      real, dimension(:,:),    allocatable, save :: refspos
      real, save :: sltcen(3),sltqrn(0:3)
      integer, parameter :: movmax=10
      real, parameter :: trmax=0.20e0
      real, parameter :: agmax=0.10e0
c   file for reference structure
      character(*), parameter :: reffile='RefInfo'  ! reference structure
      integer, parameter :: refio=71                ! reference structure IO
c   specifier to treat the reference as the total or as the system part only
      integer, parameter :: reftot=99, refsys=98
c
c
      contains
c
      subroutine instslt(wgtslcf,caltype)
c
      use engmain, only: nummol,slttype,inscnd,inscfg,numslt,sltlist,
     #                   iseed
      integer insml,m
      real wgtslcf,pcom(3),qrtn(0:3)
      character*4 caltype
      character*3 insyn
c
      call instwgt(wgtslcf,caltype)
c
      if(caltype.eq.'init') unrn=real(iseed)
      if((caltype.eq.'init').or.(caltype.eq.'last')) then
        if(slttype.eq.3) call getsolute(0,caltype)       ! flexible solute
        return
      endif
c
      m=0
      if(numslt.ne.1) m=9
      if((numslt.eq.1).and.(sltlist(1).ne.nummol)) m=9
      if(m.ne.0) call insrt_stop('set')    ! incorrect solute specification
c
      insml=nummol                                       ! inserted solute
      if(slttype.eq.3) call getsolute(insml,'conf')      ! flexible solute
c
      if(inscnd.le.2) then
        if(inscfg.ne.2) then
          insyn='not'
          do 2101 while(insyn.eq.'not')
            call sltpstn(m,pcom,'insert',insml)          ! center of mass
            if(inscfg.ne.1) call rndmvec('q',qrtn,0.0e0) ! random orientation
            call coordinate(insml,pcom,qrtn)             ! site coordinate
            call insscheme(insml,insyn)                  ! user-defined scheme
2101      continue
        endif
        if(inscfg.eq.2) call coordinate(insml,pcom,qrtn) ! site coordinate
      endif
c
      if(inscnd.eq.3) call refmc('inst')                 ! MC with reference
c
      return
      end subroutine
c
c
      subroutine sltpstn(sltstat,pcom,type,tagslt)
c
      use engmain, only: nummol,maxsite,numatm,
     #                   boxshp,inscnd,inscfg,hostspec,
     #                   moltype,numsite,specatm,sitepos,cell,invcl,
     #                   lwreg,upreg
      integer sltstat,tagslt,stmax,sid,ati,pti,i,m,k,centag(numatm)
      real rdum,clm(3),pcom(3),qrtn(0:3),rst,dis,syscen(3),elen
      character*6 type
c
      if(inscfg.eq.2) then
        sltstat=1
        return
      endif
c
      if(type.eq.'solutn') then
        do 1701 ati=1,numatm
          centag(ati)=0
1701    continue
        stmax=numsite(tagslt)
        do 1702 sid=1,stmax
          ati=specatm(sid,tagslt)
          centag(ati)=1
1702    continue
        call getcen(centag,pcom)
      endif
c
      if(inscnd.eq.0) then   ! solute with random position
        if(type.eq.'insert') then
          if(boxshp.eq.0) call insrt_stop('geo')  ! system has to be periodic
          do 1501 k=1,3
            call URAND(rdum)
            clm(k)=rdum-0.50e0
1501      continue
          do 1502 m=1,3
            rst=0.0e0
            do 1503 k=1,3
              rst=rst+cell(m,k)*clm(k)
1503        continue
            pcom(m)=rst
1502      continue
        endif
      endif
c
      if((inscnd.eq.1).or.(inscnd.eq.2)) then     ! system center
        do 1751 i=1,nummol
          pti=moltype(i)
          if(pti.eq.hostspec) m=1
          if(pti.ne.hostspec) m=0
          stmax=numsite(i)
          do 1752 sid=1,stmax
            ati=specatm(sid,i)
            centag(ati)=m
1752      continue
1751    continue
        call getcen(centag,syscen)
      endif
c
      if(inscnd.eq.1) then   ! solute in spherical object or isolated droplet
        if(type.eq.'insert') then
          call rndmvec('p',qrtn,(lwreg/upreg))
          do 1511 m=1,3
            pcom(m)=upreg*qrtn(m)+syscen(m)
1511      continue
        endif
        dis=0.0e0
        do 1512 m=1,3
          rst=pcom(m)-syscen(m)
          dis=dis+rst*rst
1512    continue
        dis=sqrt(dis)
      endif
c
      if(inscnd.eq.2) then   ! solute in slab geometry
        if(boxshp.eq.0) call insrt_stop('geo')    ! system has to be periodic
        elen=0.0e0
        do 1521 m=1,3
          elen=elen+cell(m,3)*cell(m,3)
1521    continue
        elen=sqrt(elen)
        if(type.eq.'insert') then
          do 1522 k=1,2
            call URAND(rdum)
            clm(k)=rdum-0.50e0
1522      continue
          do 1523 m=1,3
            rst=0.0e0
            do 1524 k=1,2
              rst=rst+cell(m,k)*clm(k)
1524        continue
            pcom(m)=rst
1523      continue
          call URAND(rdum)
          rst=lwreg+rdum*(upreg-lwreg)
          call URAND(rdum)
          if(rdum.le.0.50e0) rst=-rst
          dis=0.0e0
          do 1525 m=1,3
            dis=dis+invcl(3,m)*syscen(m)
1525      continue
          dis=dis+rst/elen
          do 1526 m=1,3
            pcom(m)=pcom(m)+dis*cell(m,3)
1526      continue
        endif
        rst=0.0e0
        do 1527 m=1,3
          rst=rst+invcl(3,m)*(pcom(m)-syscen(m))
1527    continue
        dis=abs(rst)*elen
      endif
c
      if(inscnd.eq.3) then   ! solute against reference structure
        if(type.eq.'insert') call insrt_stop('bug')
        call refsdev(dis,2,'system')
      endif
c
      if(inscnd.eq.0) sltstat=1
      if(inscnd.gt.0) then
        if((lwreg.le.dis).and.(dis.le.upreg)) then
          sltstat=1
        else
          if(type.eq.'solutn') sltstat=0
          if(type.eq.'insert') call insrt_stop('bug')
        endif
      endif
c
      return
      end subroutine
c
c
      subroutine getcen(centag,cen)             ! getting the center of mass
      use engmain, only: numatm,sitemass,sitepos
      integer ati,m,centag(numatm)
      real wgt,cen(3),sitm
      wgt=0.0e0
      do 3331 m=1,3
        cen(m)=0.0e0
3331  continue
      do 3332 ati=1,numatm
        if(centag(ati).eq.1) then
          sitm=sitemass(ati)
          wgt=wgt+sitm
          do 3333 m=1,3
            cen(m)=cen(m)+sitm*sitepos(m,ati)
3333      continue
        endif
3332  continue
      do 3334 m=1,3
        cen(m)=cen(m)/wgt
3334  continue
      return
      end subroutine
c
c
      subroutine rndmvec(vectp,qrtn,lwbnd)
      character vectp
      integer m,inim
      real qrtn(0:3),lwbnd,rdum,factor
      if(vectp.eq.'p') inim=1
      if(vectp.eq.'q') inim=0
      factor=2.0e0
      do 1515 m=0,3
        qrtn(m)=0.0e0
1515  continue
      do 1511 while((factor.gt.1.0e0).or.(factor.lt.lwbnd))
        factor=0.0e0
        do 1512 m=inim,3
          call URAND(rdum)
          qrtn(m)=2.0e0*rdum-1.0e0
          factor=factor+qrtn(m)*qrtn(m)
1512    continue
        factor=sqrt(factor)
1511  continue
      if(vectp.eq.'q') then
        do 1513 m=0,3
          qrtn(m)=qrtn(m)/factor
1513    continue
      endif
      return
      end subroutine
c
c
      subroutine coordinate(i,pcom,qrtn)
      use engmain, only: nummol,maxsite,numatm,inscfg,
     #                   numsite,bfcoord,specatm,sitepos
      integer stmax,sid,ati,m,k,i
      real pcom(3),qrtn(0:3),rotmat(3,3),rst
      stmax=numsite(i)
      if(inscfg.eq.0) call getrot(qrtn,rotmat)        ! rotation matrix
      do 3501 sid=1,stmax
        ati=specatm(sid,i)
        do 3511 m=1,3
          if(inscfg.eq.0) then                        ! translation + rotation
            rst=0.0e0
            do 3512 k=1,3
              rst=rst+rotmat(k,m)*bfcoord(k,sid)
3512        continue
            rst=rst+pcom(m)
          endif
          if(inscfg.eq.1) rst=pcom(m)+bfcoord(m,sid)  ! translation only
          if(inscfg.eq.2) rst=bfcoord(m,sid)          ! no change
          sitepos(m,ati)=rst
3511    continue
3501  continue
      return
      end subroutine
c
c
      subroutine getrot(qrtn,rotmat)
      real qrtn(0:3),rotmat(3,3)
      rotmat(1,1)=qrtn(0)*qrtn(0)+qrtn(1)*qrtn(1)
     #                           -qrtn(2)*qrtn(2)-qrtn(3)*qrtn(3)
      rotmat(2,2)=qrtn(0)*qrtn(0)-qrtn(1)*qrtn(1)
     #                           +qrtn(2)*qrtn(2)-qrtn(3)*qrtn(3)
      rotmat(3,3)=qrtn(0)*qrtn(0)-qrtn(1)*qrtn(1)
     #                           -qrtn(2)*qrtn(2)+qrtn(3)*qrtn(3)
      rotmat(1,2)=2.0e0*(qrtn(1)*qrtn(2)+qrtn(0)*qrtn(3))
      rotmat(2,1)=2.0e0*(qrtn(1)*qrtn(2)-qrtn(0)*qrtn(3))
      rotmat(1,3)=2.0e0*(qrtn(1)*qrtn(3)-qrtn(0)*qrtn(2))
      rotmat(3,1)=2.0e0*(qrtn(1)*qrtn(3)+qrtn(0)*qrtn(2))
      rotmat(2,3)=2.0e0*(qrtn(2)*qrtn(3)+qrtn(0)*qrtn(1))
      rotmat(3,2)=2.0e0*(qrtn(2)*qrtn(3)-qrtn(0)*qrtn(1))
      return
      end subroutine
c
c
      subroutine insscheme(insml,insyn)                ! user-defined scheme
      use engmain, only: nummol,maxsite,numatm,numsite,specatm,sitepos
      integer insml,stmax,sid,ati
      character*3 insyn
      insyn='yes'
      stmax=numsite(insml)
      return
      end subroutine
c
c
      subroutine insrt_stop(type)
      use engmain, only: io6
      use mpiproc                                                      ! MPI
      character*3 type
      if(type.eq.'set') write(io6,791)
      if(type.eq.'geo') write(io6,792)
      if(type.eq.'bug') write(io6,793)
791   format(' The solute specification is incorrectly set')
792   format(' The system geometry is incorrectly set')
793   format(' There is a bug in the insertion program')
      call mpi_setup('stop')                                           ! MPI
      stop
      end subroutine
c
c
      subroutine getsolute(i,caltype)
c
      use engmain, only: nummol,maxsite,inscfg,numsite,bfcoord
      use setconf, only: molcen
      use OUTname, only: iofmt,bxiso,toptp,skpio,OUTconfig,OUTskip
      character*4 caltype
      integer i,sid,stmax,m
      real xst(3),dumcl(3,3),factor
      real*4 sglfct
      character rddum
      real, dimension(:,:), allocatable :: psite
c
      if(caltype.eq.'init') then
        if(iofmt.eq.'yes') open(unit=slcnf,file=slttrj,status='old')
        if(iofmt.eq.'not') open(unit=slcnf,file=slttrj,status='old',
     #                                           form='unformatted')
        call OUTskip(slcnf,iofmt,skpio)
        return
      endif
      if(caltype.eq.'last') then
        close(slcnf)
        return
      endif
c
      stmax=numsite(i)
      allocate( psite(3,stmax) )
      if(bxiso.eq.'not') sid=0
      if(bxiso.eq.'yes') sid=1
      call OUTconfig(psite,dumcl,stmax,sid,slcnf,'trj')
c
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
c
      backspace(slcnf)
      goto 3195
3199  rewind(slcnf)
      call OUTskip(slcnf,iofmt,skpio)
3195  continue
c
      if(inscfg.ne.2) call molcen(i,psite,xst,'com')
      do 3131 sid=1,stmax
        do 3132 m=1,3
          if(inscfg.ne.2) bfcoord(m,sid)=psite(m,sid)-xst(m)
          if(inscfg.eq.2) bfcoord(m,sid)=psite(m,sid)
3132    continue
3131  continue
      deallocate( psite )
c
      return
      end subroutine
c
c
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
2199    rewind(swinf)
2195    continue
      endif
      return
      end subroutine
c
c
      subroutine URAND(rndm)          ! uniform random number generator
      real rndm,mtpl,mdls
      integer gen1,gen2,rat
      parameter(gen1=16807,gen2=2147483647)
      mdls=real(gen2)
      mtpl=unrn*real(gen1)
      rat=int(mtpl/mdls)
      mtpl=mtpl-mdls*real(rat)
      unrn=mtpl
      rndm=unrn/mdls
      return
      end subroutine
c
c
      subroutine refmc(caltype)
c
      use engmain, only: nummol,maxsite,numatm,slttype,numsite,refmlid,
     #                   lwreg,upreg,bfcoord,specatm,sitepos
      integer i,k,m,q,ati,rfi,sid,stmax
      real xst(3),centg(3),cenrf(3)
      real factor,bfqrn(0:3),rtmbf(3,3),sltsite(3,maxsite)
      character*4 caltype
      character*8 atmtype,dump
      character eletype
c
      if(caltype.eq.'init') then
        allocate( refsatm(maxsite,nummol),refspos(3,numatm) )
        do 1311 i=1,nummol
          do 1312 sid=1,maxsite
            refsatm(sid,i)=0
1312      continue
1311    continue
        do 1313 ati=1,numatm
          do 1314 m=1,3
            refspos(m,ati)=0.0e0
1314      continue
1313    continue
c
c  read the reference structure
        open(unit=refio,file=reffile,status='old')
        do 1001 i=1,nummol
          rfi=refmlid(i)
          if(rfi.ne.0) then
            stmax=numsite(i)
            do 1011 sid=1,stmax
              read(refio,*) dump,k,atmtype,dump,q,
     #                      (xst(m), m=1,3),factor,factor
              eletype=atmtype(1:1)
              if(eletype.eq.'H') refsatm(sid,i)=0      ! hydrogen atom
              if(eletype.ne.'H') refsatm(sid,i)=rfi    ! heavy atom
              ati=specatm(sid,i)
              do 1012 m=1,3
                refspos(m,ati)=xst(m)
1012          continue
1011        continue
          endif
1001    continue
        close(refio)
c
c  build the initial configuration of the solute
        if(slttype.ge.2) then
          call dispref(centg,cenrf,bfqrn,1)     ! matching reference
          call getrot(bfqrn,rtmbf)
          stmax=numsite(nummol)                 ! inserted solute molecule
          do 1511 sid=1,stmax
            ati=specatm(sid,nummol)
            do 1521 m=1,3
              factor=0.0e0
              do 1522 k=1,3
                factor=factor+rtmbf(k,m)*(refspos(k,ati)-cenrf(k))
1522          continue
              sitepos(m,ati)=factor+centg(m)
1521        continue
1511      continue
          call dispref(sltcen,cenrf,sltqrn,2)   ! solute center & orientation
        endif
      endif
c
      if(caltype.eq.'inst') then                ! inserted solute molecule
        stmax=numsite(nummol)
        do 3101 sid=1,stmax
          ati=specatm(sid,nummol)
          do 3102 m=1,3
            sltsite(m,sid)=sitepos(m,ati)
3102      continue
          do 3103 m=1,3
            sitepos(m,ati)=bfcoord(m,sid)
3103      continue
3101    continue
        call refcen(xst,refsatm,sitepos,2)      ! solute center in bfcoord
        do 3111 sid=1,stmax
          ati=specatm(sid,nummol)
          do 3112 m=1,3
            bfcoord(m,sid)=bfcoord(m,sid)-xst(m)
3112      continue
3111    continue
        call dispref(centg,cenrf,bfqrn,2)       ! matching solute
        call getrot(bfqrn,rtmbf)
        k=0
        do 3001 q=1,movmax
          call sltmove(sltcen,sltqrn,rtmbf,'frwd')
          call refsdev(factor,2,'extend')
          if((factor.lt.lwreg).or.(factor.gt.upreg)) then
            k=k+1
            call sltmove(sltcen,sltqrn,rtmbf,'back')
          endif
3001    continue
        if(k.eq.movmax) then                    ! when all MC are rejected
          do 3211 sid=1,stmax
            ati=specatm(sid,nummol)
            do 3212 m=1,3
              sitepos(m,ati)=sltsite(m,sid)
3212        continue
3211      continue
        endif
      endif
c
      return
      end subroutine
c
c
      subroutine refcen(cen,lstatm,posatm,refmol)
      use engmain, only: nummol,maxsite,numatm,numsite,specatm
      integer refmol,rfyn,i,m,ati,rfi,sid,stmax,lstatm(maxsite,nummol)
      real cen(3),posatm(3,numatm),totwgt
      totwgt=0.0e0
      do 1101 m=1,3
        cen(m)=0.0e0
1101  continue
      do 1001 i=1,nummol
        call getrfyn(i,rfi,rfyn,refmol)
        if(rfyn.eq.1) then
          stmax=numsite(i)
          do 1011 sid=1,stmax
            if(lstatm(sid,i).eq.rfi) then
              ati=specatm(sid,i)
              totwgt=totwgt+1.0e0
              do 1021 m=1,3
                cen(m)=cen(m)+posatm(m,ati)
1021          continue
            endif
1011      continue
        endif
1001  continue
      do 1201 m=1,3
        cen(m)=cen(m)/totwgt
1201  continue
      return
      end subroutine
c
c
      subroutine sltmove(pcen,qrtn,rtmbf,mvtype)
      use engmain, only: nummol,maxsite,numatm,
     #                   numsite,bfcoord,specatm,sitepos
      integer m,k,q,ati,sid,stmax,ax1,ax2
      real pcen(3),qrtn(0:3),rtmbf(3,3)
      real rdum,factor,rotmat(3,3),rtmmov(3,3),axis(3)
      character*4 mvtype
      real, save :: odcn(3),oqrn(0:3)
      stmax=numsite(nummol)                       ! inserted solute molecule
      do 1001 m=1,3
        if(mvtype.eq.'frwd') odcn(m)=pcen(m)
        if(mvtype.eq.'back') pcen(m)=odcn(m)
1001  continue
      do 1002 m=0,3
        if(mvtype.eq.'frwd') oqrn(m)=qrtn(m)
        if(mvtype.eq.'back') qrtn(m)=oqrn(m)
1002  continue
      if(mvtype.eq.'frwd') then
        do 3001 m=1,3                             ! translation
          call URAND(rdum)
          pcen(m)=pcen(m)+trmax*(rdum+rdum-1.0e0)
3001    continue
        factor=2.0e0                              ! rotation
        do 3101 while(factor.gt.1.0e0)
          do 3102 m=1,3
            call URAND(rdum)
            axis(m)=rdum+rdum-1.0e0
3102      continue
          factor=axis(1)*axis(1)+axis(2)*axis(2)+axis(3)*axis(3)
3101    continue
        factor=sqrt(factor)
        do 3103 m=1,3
          axis(m)=axis(m)/factor
3103    continue
        call URAND(rdum)
        ax1=mod(int(abs(unrn)),4)
        ax2=ax1
        do 3104 while(ax1.eq.ax2)
          call URAND(rdum)
          ax2=mod(int(abs(unrn)),4)
3104    continue
        call URAND(rdum)
        factor=agmax*(rdum+rdum-1.0e0)
        qrtn(ax1)=cos(factor)*oqrn(ax1)-sin(factor)*oqrn(ax2)
        qrtn(ax2)=sin(factor)*oqrn(ax1)+cos(factor)*oqrn(ax2)
        factor=0.0e0
        do 3105 m=0,3
          factor=factor+qrtn(m)*qrtn(m)
3105    continue
        factor=sqrt(factor)
        do 3106 m=0,3
          qrtn(m)=qrtn(m)/factor
3106    continue
      endif
      call getrot(qrtn,rtmmov)                    ! site coordinate
      do 3201 m=1,3
       do 3202 k=1,3
         factor=0.0e0
         do 3203 q=1,3
           factor=factor+rtmmov(q,m)*rtmbf(q,k)
3203     continue
         rotmat(m,k)=factor
3202   continue
3201  continue
      do 3211 sid=1,stmax
        ati=specatm(sid,nummol)
        do 3212 m=1,3
          factor=0.0e0
          do 3213 k=1,3
            factor=factor+rotmat(m,k)*bfcoord(k,sid)
3213      continue
          sitepos(m,ati)=factor+pcen(m)
3212    continue
3211  continue
      return
      end subroutine
c
c
      subroutine dispref(centg,cenrf,qrtn,refmol)
      use engmain, only: nummol,maxsite,numatm,numsite,specatm,sitepos
      integer refmol,rfyn,i,m,k,ati,rfi,sid,stmax
      real centg(3),cenrf(3),qrtn(0:3),qrtmat(4,4),xsm(3),xrf(3)
      real totm,xp,yp,zp,xm,ym,zm,dumv(4),work(16)
      call refcen(centg,refsatm,sitepos,refmol)
      call refcen(cenrf,refsatm,refspos,refmol)
      totm=0.0e0
      do 3211 m=1,4
       do 3212 k=1,4
         qrtmat(k,m)=0.0e0
3212   continue
3211  continue
      do 3201 i=1,nummol
        call getrfyn(i,rfi,rfyn,refmol)
        if(rfyn.eq.1) then
          stmax=numsite(i)
          do 3202 sid=1,stmax
            if(refsatm(sid,i).eq.rfi) then
              ati=specatm(sid,i)
              totm=totm+1.0e0
              do 3203 m=1,3
                xsm(m)=sitepos(m,ati)-centg(m)
3203          continue
              do 3204 m=1,3
                xrf(m)=refspos(m,ati)-cenrf(m)
3204          continue
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
3202      continue
        endif
3201  continue
      do 3213 m=1,3
       do 3214 k=m+1,4
         qrtmat(m,k)=qrtmat(k,m)
3214   continue
3213  continue
      do 3215 m=1,4
       do 3216 k=1,4
         qrtmat(k,m)=qrtmat(k,m)/totm
3216   continue
3215  continue
      call DSYEV('V','U',4,qrtmat,4,dumv,work,16,i)
      do 3221 m=0,3
        qrtn(m)=qrtmat(m+1,1)   ! eigenvector for the minimum eigenvalue
3221  continue
      totm=0.0e0
      do 3222 m=0,3
        totm=totm+qrtn(m)*qrtn(m)
3222  continue
      totm=sqrt(totm)
      do 3223 m=0,3
        qrtn(m)=qrtn(m)/totm
3223  continue
      return
      end subroutine
c
c
      subroutine refsdev(strchr,refmol,caltype)
      use engmain, only: nummol,maxsite,numatm,numsite,specatm,sitepos
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
      do 5001 i=1,nummol
        call getrfyn(i,rfi,rfyn,refmol)
        if(rfyn.eq.1) then
          stmax=numsite(i)
          do 5002 sid=1,stmax
            if(refsatm(sid,i).eq.rfi) then
              ati=specatm(sid,i)
              totm=totm+1.0e0
              do 5011 m=1,3
                factor=0.0e0
                do 5012 k=1,3
                  factor=factor+rotmat(m,k)*(sitepos(k,ati)-centg(k))
5012            continue
                xsm(m)=factor+cenrf(m)
5011          continue
              do 5013 m=1,3
                factor=xsm(m)-refspos(m,ati)
                strchr=strchr+factor*factor
5013          continue
            endif
5002      continue
        endif
5001  continue
      strchr=sqrt(strchr/totm)
      return
      end subroutine
c
c
      subroutine getrfyn(i,rfi,rfyn,refmol)
      use engmain, only: nummol,sluvid,refmlid
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
      end subroutine
c
      end module
