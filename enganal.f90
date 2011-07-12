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
subroutine enganal(stnum)
  use engmain, only: maxcnf,engdiv,skpcnf,inscnd
  use engproc, only: engclear,engconst,engstore
  use ptinsrt, only: refmc
  implicit none
  integer, intent(in) :: stnum

  if(mod(stnum,(maxcnf/engdiv)).eq.skpcnf) call engclear

  if((stnum.eq.skpcnf).and.(inscnd.eq.3)) call refmc('init')
  call engconst(stnum)
  if(mod(stnum,(maxcnf/engdiv)).eq.0) call engstore(stnum)
  return
end subroutine enganal
!     
program trjmain
  use engmain, only: maxcnf, skpcnf
  use OUTname, only: opentrj,closetrj
  use setconf, only: getconf
  use mpiproc               ! MPI
  implicit none
  integer stnum
  integer, parameter :: large=100000000
#ifdef VMDPLUGINS
  external vmdfio_init_traj, vmdfio_fini_traj
  call vmdfio_init_traj
#endif
  call mpi_setup('init')    ! MPI
  call opentrj()

  ! initialize
  call enganal_init()

  do stnum=1,maxcnf
     call getconf
     if(mod(stnum,skpcnf) == 0) call enganal(stnum)
  end do
  call closetrj
  call mpi_setup('stop')    ! MPI
#ifdef VMDPLUGINS
  call vmdfio_fini_traj
#endif
  stop
end program trjmain
