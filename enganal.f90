! -*- F90 -*-

subroutine enganal_init()
  use setconf, only: setparam
  use engproc, only: enginit
  implicit none

  call setparam
  call enginit
end subroutine enganal_init

! FIXME: recover routine which runs as "combined with MD program"
!  connection to the main routine of trajectory generation is done in
!  setparam for parameter setting and getconf for configuration reading
subroutine enganal(stnum, nactiveproc)
  use engmain, only: maxcnf,engdiv,skpcnf,inscnd
  use engproc, only: engclear,engconst,engstore
  use ptinsrt, only: refmc
  implicit none
  integer, intent(in) :: stnum, nactiveproc

  if((inscnd.eq.3)) call refmc('init')
  call engconst(stnum, nactiveproc)
end subroutine enganal
!     
program trjmain
  use engmain, only: maxcnf, skpcnf, engdiv
  use OUTname, only: opentrj,closetrj, OUTinitial
  use setconf, only: getconf_parallel
  use vmdfio_interface, only: init_vmdplugins, finish_vmdplugins
  use engproc, only: engclear,engstore, engproc_cleanup
  use mpiproc               ! MPI
  implicit none
  integer :: stnum, iproc, iskip, idiv, frames_per_div, nread, iframe

  call mpi_setup('init')    ! MPI

#ifdef VMDPLUGINS
  call init_vmdplugins()
#endif
  call opentrj()

  ! initialize
  call enganal_init()

  stnum = 0
  frames_per_div = maxcnf / skpcnf / engdiv
  if(frames_per_div <= 0) call halt_with_error("par")

  do idiv = 1, engdiv
     call engclear

     do iframe = 1, frames_per_div, nprocs
        call getconf_parallel(frames_per_div - iframe + 1, nread)
        call enganal(stnum + myrank + 1, nread)
        stnum = stnum + nread
     end do
     call engstore(stnum)
  end do
  call engproc_cleanup

  call closetrj
#ifdef VMDPLUGINS
  call finish_vmdplugins()
#endif

  call mpi_setup('stop')    ! MPI
  stop
end program trjmain
