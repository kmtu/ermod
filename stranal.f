      module stranal
      contains
      subroutine spacedst(stnum)
c
      use engmain, only: nummol,maxsite,numatm,
     #                   maxcnf,slttype,skpcnf,
     #                   estype,boxshp,inscnd,hostspec,
     #                   moltype,numsite,sluvid,refmlid,
     #                   specatm,sitepos,invcl,volume
      use engproc, only: volcorrect
      use ptinsrt, only: getcen,refsdev
      use mpiproc                                                      ! MPI
c
      real, parameter    :: dstbin=1.0e-2,maxdst=50.0e0       ! angstrom
      integer, save      :: dstmax,dstmin
      integer, parameter :: hostsite=7     ! number of sites treated for host
      integer, save      :: hostmin,hostmax,numlgnd,numrmsd
      integer, dimension(:), allocatable, save :: lgndspec
      real, dimension(:,:), allocatable, save  :: hostdist
      real, dimension(:,:), allocatable, save  :: lgnddist
      real, dimension(:,:), allocatable, save  :: rmsddist
      real, save :: dstnorm,avearea
c
      character(*), parameter :: sltevl='trjslt.tt'
      integer, parameter :: utrio=59
      real, parameter :: tiny=1.0e-20
      character, save :: sltrank
c
      integer stnum,i,uvi,rfi,lgi,ati,sid,stmax,m
      integer iddst,idtyp,centag(numatm)
      real pcom(3),syscen(3),dstnmfc,dis,rst
      real, dimension(:), allocatable :: svg,rcg
c
      if((boxshp.eq.0).or.(inscnd.eq.0)) return
      call mpi_info                                                    ! MPI
c
      if(stnum.eq.skpcnf) then
        hostmin=nummol
        hostmax=0
        sid=0
        do 2101 i=1,nummol
          if(inscnd.le.2) then
            lgi=moltype(i)
            m=hostspec
          endif
          if(inscnd.eq.3) then
            lgi=refmlid(i)
            m=1
          endif
          if(lgi.eq.m) then
            if(i.lt.hostmin) hostmin=i
            if(i.gt.hostmax) hostmax=i
            sid=sid+1
          endif
2101    continue
        if(sid.ne.(hostmax-hostmin+1)) then
          if(inscnd.le.2) call spd_stop('hst')
          if(inscnd.eq.3) call spd_stop('ref')
        endif
c
        if(inscnd.eq.3) then
          if(slttype.eq.1) numrmsd=2
          if(slttype.ge.2) numrmsd=1
          do 2131 i=1,nummol
            if(sluvid(i).le.1) then
              rfi=refmlid(i)
              if((rfi.lt.0).and.(rfi.gt.numrmsd)) call spd_stop('rfs')
            endif
2131      continue
        endif
c
        allocate( lgndspec(nummol) )
        numlgnd=0
        do 2501 i=1,nummol
          uvi=sluvid(i)
          if(i.lt.hostmin) then
            if(uvi.le.1) lgi=moltype(i)
            if(uvi.ge.2) lgi=-1
          endif
          if((hostmin.le.i).and.(i.le.hostmax)) lgi=0
          if(i.gt.hostmax) then
            if(uvi.le.1) then
              idtyp=0
              if(hostmin.gt.1) then
                do 2506 m=1,hostmin-1
                  if(idtyp.lt.moltype(m)) idtyp=moltype(m)
2506            continue
              endif
              ati=moltype(hostmin)
              do 2507 m=hostmin,hostmax
                if(ati.lt.moltype(m)) ati=moltype(m)
2507          continue
              lgi=idtyp+moltype(i)-ati
              if(lgi.le.0) call spd_stop('mol')
            endif
            if(uvi.ge.2) lgi=-1
          endif
          lgndspec(i)=lgi
          if(lgi.gt.numlgnd) numlgnd=lgi
2501    continue
c
        dstmax=nint(maxdst/dstbin)
        if(inscnd.eq.1) dstmin=1
        if(inscnd.eq.2) dstmin=-dstmax
        if(inscnd.eq.3) dstmin=1
        if(inscnd.le.2) allocate( hostdist(dstmin:dstmax,hostsite) )
        allocate( lgnddist(dstmin:dstmax,numlgnd) )
        if(inscnd.eq.3) allocate( rmsddist(dstmin:dstmax,numrmsd) )
        do 2511 m=1,3
          if(m.eq.1) stmax=hostsite
          if(m.eq.2) stmax=numlgnd
          if(m.eq.3) stmax=numrmsd
          if((m.eq.1).and.(inscnd.eq.3)) goto 2599
          if((m.eq.3).and.(inscnd.le.2)) goto 2599
          do 2512 idtyp=1,stmax
            do 2513 iddst=dstmin,dstmax
              if(m.eq.1) hostdist(iddst,idtyp)=0.0e0
              if(m.eq.2) lgnddist(iddst,idtyp)=0.0e0
              if(m.eq.3) rmsddist(iddst,idtyp)=0.0e0
2513        continue
2512      continue
2599      continue
2511    continue
        dstnorm=0.0e0
        avearea=0.0e0
c
        if(slttype.eq.1) then
          sltrank='n'
          do 2531 i=myrank+1,nummol,nprocs                             ! MPI
            if(sluvid(i).eq.1) sltrank='y'
2531      continue
          if(sltrank.eq.'y') open(unit=utrio,file=sltevl,status='new')
        endif
      endif
c
c
      dstnmfc=1.0e0
      if(estype.eq.2) call volcorrect(dstnmfc)
      dstnorm=dstnorm+dstnmfc
      rst=0.0e0
      do 3811 m=1,3
        rst=rst+invcl(3,m)*invcl(3,m)
3811  continue
      avearea=avearea+dstnmfc*volume*sqrt(rst)
      do 3101 i=1,nummol                         ! setting the system center
        lgi=lgndspec(i)
        if(lgi.eq.0) m=1
        if(lgi.ne.0) m=0
        stmax=numsite(i)
        do 3102 sid=1,stmax
          ati=specatm(sid,i)
          centag(ati)=m
3102    continue
3101  continue
      call getcen(centag,syscen)
c
      do 2001 i=1+myrank,nummol,nprocs                                 ! MPI
        uvi=sluvid(i)
        lgi=lgndspec(i)
        stmax=numsite(i)
        if((inscnd.le.2).and.(lgi.eq.0).and.(uvi.eq.0)) then
          do 2011 sid=1,stmax
            call hostid(sid,idtyp)
            if((idtyp.gt.0).and.(idtyp.le.hostsite)) then
              ati=specatm(sid,i)
              if(inscnd.eq.1) then
                dis=0.0e0
                do 2012 m=1,3
                  rst=sitepos(m,ati)-syscen(m)
                  dis=dis+rst*rst
2012            continue
                rst=sqrt(dis)
                iddst=int(rst/dstbin)+1
              endif
              if(inscnd.eq.2) then
                rst=sitepos(3,ati)-syscen(3)
                iddst=nint(rst/dstbin)
              endif
              if((dstmin.le.iddst).and.(iddst.le.dstmax)) then
                hostdist(iddst,idtyp)=hostdist(iddst,idtyp)+dstnmfc
              endif
            endif
2011      continue
        endif
c
        if(lgi.gt.0) then
          do 2021 ati=1,numatm
            centag(ati)=0
2021      continue
          stmax=numsite(i)
          do 2022 sid=1,stmax
            ati=specatm(sid,i)
            centag(ati)=1
2022      continue
          call getcen(centag,pcom)
          if((sltrank.eq.'y').and.(uvi.eq.1)) then
            write(utrio,'(2i10,3f12.3)') stnum,i,(pcom(m), m=1,3)
          endif
          if((inscnd.eq.1).or.(inscnd.eq.3)) then
            dis=0.0e0
            do 2023 m=1,3
              rst=pcom(m)-syscen(m)
              dis=dis+rst*rst
2023        continue
            rst=sqrt(dis)
            iddst=int(rst/dstbin)+1
          endif
          if(inscnd.eq.2) iddst=nint((pcom(3)-syscen(3))/dstbin)
          idtyp=lgndspec(i)
          if((dstmin.le.iddst).and.(iddst.le.dstmax)) then
            lgnddist(iddst,idtyp)=lgnddist(iddst,idtyp)+dstnmfc
          endif
        endif
c
        if((inscnd.eq.3).and.(uvi.le.1)) then
          rfi=refmlid(i)
          if((rfi.gt.0).and.(rfi.le.numrmsd)) then
            call refsdev(rst,rfi,'system')
            iddst=int(rst/dstbin)+1
            if((dstmin.le.iddst).and.(iddst.le.dstmax)) then
              if(rfi.eq.1) rst=1.0e0/real(hostmax-hostmin+1)
              if(rfi.eq.2) rst=1.0e0
              rmsddist(iddst,rfi)=rmsddist(iddst,rfi)+rst*dstnmfc
            endif
          endif
        endif
2001  continue
c
      if(stnum.eq.maxcnf) then
        if((slttype.eq.1).and.(sltrank.eq.'y')) then
          endfile(utrio)
          close(utrio)
        endif
c
        avearea=avearea/dstnorm
        allocate( svg(dstmin:dstmax),rcg(dstmin:dstmax) )
        do 2701 m=1,3
          if(m.eq.1) stmax=hostsite
          if(m.eq.2) stmax=numlgnd
          if(m.eq.3) stmax=numrmsd
          if((m.eq.1).and.(inscnd.eq.3)) goto 2799
          if((m.eq.3).and.(inscnd.le.2)) goto 2799
          do 2702 idtyp=1,stmax
            do 2703 iddst=dstmin,dstmax
              if(m.eq.1) svg(iddst)=hostdist(iddst,idtyp)
              if(m.eq.2) svg(iddst)=lgnddist(iddst,idtyp)
              if(m.eq.3) svg(iddst)=rmsddist(iddst,idtyp)
2703        continue
#ifndef noMPI
            call mpi_reduce(svg,rcg,(dstmax-dstmin+1),                 ! MPI
     #           mpi_double_precision,mpi_sum,0,mpi_comm_world,ierror) ! MPI
            do 2705 iddst=dstmin,dstmax                                ! MPI
              svg(iddst)=rcg(iddst)                                    ! MPI
2705        continue                                                   ! MPI
#endif
            do 2704 iddst=dstmin,dstmax
              if(m.eq.1) hostdist(iddst,idtyp)=svg(iddst)/dstnorm
              if(m.eq.2) lgnddist(iddst,idtyp)=svg(iddst)/dstnorm
              if(m.eq.3) rmsddist(iddst,idtyp)=svg(iddst)/dstnorm
2704        continue
2702      continue
c
          if(myrank.ne.0) go to 2799                                   ! MPI
          if(m.eq.1) open(unit=58,file='hostdst.tt',status='new')
          if(m.eq.2) open(unit=58,file='lgnddst.tt',status='new')
          if(m.eq.3) open(unit=58,file='rmsddst.tt',status='new')
          do 2711 iddst=dstmin,dstmax
            if(inscnd.eq.1) then
              rst=(real(iddst)-0.50e0)*dstbin
              dis=4.1887902e0*dstbin*dstbin*dstbin
     #                       *real(3*iddst*(iddst-1)+1)
            endif
            if(inscnd.eq.2) then
              rst=real(iddst)*dstbin
              dis=dstbin*avearea
            endif
            if(inscnd.eq.3) then
              rst=(real(iddst)-0.50e0)*dstbin
              if(m.eq.2) dis=4.1887902e0*dstbin*dstbin*dstbin
     #                                  *real(3*iddst*(iddst-1)+1)
              if(m.eq.3) dis=dstbin
            endif
            do 2712 idtyp=1,stmax
              if(m.eq.1) svg(idtyp)=hostdist(iddst,idtyp)/dis
              if(m.eq.2) svg(idtyp)=lgnddist(iddst,idtyp)/dis
              if(m.eq.3) svg(idtyp)=rmsddist(iddst,idtyp)/dis
              if(svg(idtyp).le.tiny) svg(idtyp)=0.0e0
2712        continue
            write(58,581) rst,(svg(idtyp), idtyp=1,stmax)
581         format(f8.3,99999999g10.3)
2711      continue
          endfile(58)
          close(58)
2799      continue                                                     ! MPI
2701    continue
        deallocate( svg,rcg )
      endif
c
      return
      end subroutine
c
c
      subroutine hostid(sid,idtyp)
      integer sid,idtyp
      idtyp=0
c
c  DMPC
      if(sid.eq.1) idtyp=1                                    ! N
      if((sid.eq.2).or.(sid.eq.6).or.(sid.eq.10)) idtyp=2     ! C
      if((sid.eq.14).or.(sid.eq.17)) idtyp=2                  ! C
      if(sid.eq.20) idtyp=3                                   ! P
      if((sid.eq.21).or.(sid.eq.22).or.(sid.eq.24)) idtyp=4   ! O (PO3-)
      if(sid.eq.23) idtyp=5                                   ! O (glycerol)
      if((sid.eq.25).or.(sid.eq.28).or.(sid.eq.36)) idtyp=5   ! C (glycerol)
      if((sid.eq.30).or.(sid.eq.39)) idtyp=5                  ! O (glycerol)
      if((sid.ge.31).and.(sid.le.33)) idtyp=6                 ! CO, CH2 (chain)
      if((sid.ge.40).and.(sid.le.42)) idtyp=6                 ! CO, CH2 (chain)
      if((sid.ge.45).and.(sid.le.77)) then                    ! chain
        if(mod(sid-45,3).eq.0) idtyp=6                          ! C (methylene)
      endif
      if((sid.ge.82).and.(sid.le.114)) then                   ! chain
        if(mod(sid-82,3).eq.0) idtyp=6                          ! C (methylene)
      endif
      if((sid.eq.78).or.(sid.eq.115)) idtyp=7                 ! chain (methyl)
c
      return
      end subroutine
c
c
      subroutine spd_stop(type)
      use engmain, only: io6
      use mpiproc                                                      ! MPI
      character*3 type
      if(type.eq.'hst') write(io6,591)
      if(type.eq.'ref') write(io6,592)
      if(type.eq.'rfs') write(io6,593)
      if(type.eq.'mol') write(io6,594)
591   format(' The numbering of host molecules is not consecutive')
592   format(' The numbering of reference molecules is not consecutive')
593   format(' The refmlid variable is incorrectly set')
594   format(' The numbering of molecules is not consecutive')
      call mpi_setup('stop')                                           ! MPI
      stop
      end subroutine
c
      end module
