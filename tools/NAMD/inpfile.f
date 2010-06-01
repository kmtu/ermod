      program inpfile
c
      integer, parameter :: large=1000000,ljform=1
      real, parameter :: sgmcnv=2.0e0**(5.0e0/6.0e0) ! from Rmin/2 to sigma
      character*80 nfile,oufile
      character*7 rdcr
      character*4 dumc,atmtp(large),psatp(large)
      character*1 topc
      character*9, parameter :: numbers='123456789'
      real eps(large),rminh(large),chrg(large),dumf
      integer sid,i,j,k,m,pt,nummol(large),numsite(large)
      integer numtype,sltpick,slvcnt
c
      write(6,*) ' What is the name of the parameter file?'
      read(5,*) nfile
      open(2,file=nfile,status='old')
      do i=1,large
        read(2,*,ERR=7777) rdcr
        if(rdcr(1:7).eq.'NONBOND') exit
7777    continue
      enddo
c
      read(2,*)
      sid=0
      do i=1,large
        k=0
        read(2,*,ERR=9999) topc
        if(topc.ne.'!') then
          backspace(2)
          read(2,*,ERR=9999) rdcr
          if(rdcr(1:5).eq.'HBOND') then
            exit
          else
            k=1
          endif
        endif
        if(k.eq.1) then
          backspace(2)
          sid=sid+1
          read(2,*) atmtp(sid),dumf,eps(sid),rminh(sid)
          eps(sid)=-eps(sid)
          if(ljform.eq.0) rminh(sid)=sgmcnv*rminh(sid)
        endif
      enddo
9999  continue
      close(2)
c
      write(6,*) ' Where is the MDinfo file?'
      read(5,*) nfile
      open(3,file=nfile,status='old')
      read(3,*) i,numtype
      read(3,*) (nummol(i), i=1,numtype)
      read(3,*) (numsite(i), i=1,numtype)
      close(3)
      write(6,*) ' Which is the solute?'
      read(5,*) sltpick
c
      write(6,*) ' What is the name of the PSF file?'
      read(5,*) nfile
      open(9,file=nfile,status='old')
      do k=1,2 ; read(9,*) ; enddo
      read(9,*) i
      do k=1,i+2 ; read(9,*) ; enddo
      slvcnt=0
      do i=1,numtype ; do j=1,nummol(i)
        do k=1,numsite(i)
          read(9,*) pt,topc,m,dumc,dumc,psatp(pt),chrg(pt)
        enddo
        if(j.gt.2) then
          k=0
          if(psatp(pt).ne.psatp(pt-numsite(i))) k=1
          if(chrg(pt).ne.chrg(pt-numsite(i))) k=1
          if(k.eq.1) then
            write(6,*) ' Inconsistency in PSF file' ; stop
          endif
        endif
        if(j.eq.nummol(i)) then
          if(i.eq.sltpick) oufile='SltInfo'
          if(i.ne.sltpick) then
            slvcnt=slvcnt+1
            oufile='MolPrm'//numbers(slvcnt:slvcnt)
          endif
          open(31,file=oufile,status='new')
          do k=1,numsite(i) ; backspace(9) ; enddo
          do k=1,numsite(i)
            read(9,*) pt,topc,m,dumc,dumc,psatp(pt),chrg(pt)
            do m=1,sid
              if(psatp(pt).eq.atmtp(m)) then
                dumc=atmtp(m)
                write(31,351) k,atmtp(m),chrg(pt),eps(m),rminh(m)
              endif
351           format(i4,'   ',a4,3f12.6)
            enddo
          enddo
          endfile(31) ; close(31)
        endif
      enddo ; enddo
c
      stop ; end
