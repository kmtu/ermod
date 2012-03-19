! -*- F90 -*-
! ERMOD - Eneregy Representation Module
! Copyright (C) 2000-2012 Nobuyuki Matubayasi
! Copyright (C) 2010-2012 Shun Sakuraba
! 
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION ""
#endif

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
  use OUTname, only: opentrj, closetrj, OUTinitial, initconf, finiconf
  use setconf, only: getconf_parallel
  use vmdfio_interface, only: init_vmdplugins, finish_vmdplugins
  use engproc, only: engclear,engstore, engproc_cleanup
  use mpiproc               ! MPI
  implicit none
  integer :: stnum, iproc, iskip, idiv, frames_per_div, nread, iframe

  call mpi_setup('init')    ! MPI
  if(myrank == 0) then
     print *, "ERMOD " // PACKAGE_VERSION // ", Copyright (C) 2000-2012 Nobuyuki Matubayasi"
     print *, "                           2010-2012 Shun Sakuraba"
     print *, "ERMOD comes with ABSOLUTELY NO WARRANTY."
     print *, "This is free software, and you can redistribute it"
     print *, "and/or modify it under certain conditions."
     print *, "See LICENSE file for details."
  end if
  call initconf()

  if(myrank == 0) then
#ifdef VMDPLUGINS
     call init_vmdplugins()
#endif
     call opentrj()
  end if

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

  if(myrank == 0) then
     call closetrj
#ifdef VMDPLUGINS
     call finish_vmdplugins()
#endif
  end if

  call finiconf
  call mpi_setup('stop')    ! MPI
  stop
end program trjmain
