      program sltref
      character*80 nfile1,nfile2,nfile3,nfile4
      integer i,im,j,k,np,sltpick
      integer, dimension(:), allocatable :: mnum,anum
      real, dimension(:,:), allocatable :: cord,prm
      character*4, dimension(:), allocatable :: atm
      character*30 dmc
c
      write(6,*) ' What is the name of the PDB file?'
      read(5,*) nfile1
      write(6,*) ' What is the name of the SltInfo file for soln?'
      read(5,*) nfile2
      write(6,*) ' What is the name of the SltInfo file for refs?'
      read(5,*) nfile3
      write(6,*) ' What is the name of the MDinfo file for soln?'
      read(5,*) nfile4
      write(6,*) ' Which is the solute?'
      read(5,*) sltpick
c
      open(2,file=nfile4,status='old')
      read(2,*) i,im ; allocate( mnum(im),anum(im) )
      read(2,*) (mnum(i), i=1,im)
      read(2,*) (anum(i), i=1,im) ; close(2)
c
      np=anum(sltpick)
      allocate( cord(3,np),prm(3,np),atm(np) )
c
      open(2,file=nfile1,status='old')
      read(2,*)
      do i=1,im
        if(i.ne.sltpick) then
          do j=1,mnum(i) ; do k=1,anum(i) ; read(2,*) ; enddo ; enddo
        endif
        if(i.eq.sltpick) then
          if(mnum(i).ne.1) then ; write(6,*) 'Error' ; stop ; endif
          do k=1,np
            read(2,'(a30,3f8.3)') dmc,cord(1:3,k)
          enddo
        endif
      enddo
      close(2)
c
      open(3,file=nfile2,status='old')
      do k=1,np
        read(3,*) j,atm(k),prm(1:3,k)
      enddo
      close(3)
c
      open(31,file=nfile3,status='new')
      do k=1,np
        write(31,311) k,atm(k),prm(1:3,k),cord(1:3,k)
      enddo
      endfile(31) ; close(31)
311   format(i4,'   ',a4,3f12.6,3f9.3)
      stop ; end
