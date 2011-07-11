! -*- F90 -*-
!  connection to the main routine of trajectory generation is done in
!  setparam for parameter setting and getconf for configuration reading
subroutine enganal(stnum)
  use engmain, only: maxcnf,engdiv,skpcnf,inscnd
  use engproc, only: enginit,engclear,engconst,engstore
  use setconf, only: setparam,getconf
  use ptinsrt, only: refmc
  implicit none
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
  if(mod(stnum,(maxcnf/engdiv)).eq.0) call engstore(stnum)
  return
end subroutine enganal
!     
#ifdef trjctry
program trjmain
  use engmain, only: maxcnf
  use OUTname, only: opentrj,closetrj
  use mpiproc               ! MPI
  implicit none
  integer stnum
#ifdef VMDPLUGINS
  external vmdfio_init_traj, vmdfio_fini_traj
  call vmdfio_init_traj
#endif
  call mpi_setup('init')    ! MPI
  call opentrj
  do stnum=1,maxconf
     call enganal(stnum)
  end do
  call closetrj
  call mpi_setup('stop')    ! MPI
#ifdef VMDPLUGINS
  call vmdfio_fini_traj
#endif
  stop
end program trjmain
#endif
