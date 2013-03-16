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

module ptinsrt
  !  test particle insertion of the solute
  implicit none
  real, save :: unrn
  !  single-solute trajectrory file           used only when slttype = 3
  character(*), parameter :: slttrj = 'SltConf'  ! solute filename
  character(*), parameter :: sltwgt = 'SltWght'  ! solute weight filename
  integer, parameter :: sltwgt_io = 31           ! solute weight ID
  !
  !  insertion against reference structure
  !   refmlid : superposition reference among solvent species
  !             --- 0 : not reference
  !                 1 : reference solvent  2 : reference solute
  !             value set in subroutine setparam
  !  file for reference structure
  character(*), parameter :: reffile='RefInfo'  ! reference structure
  integer, parameter :: refio=71                ! reference structure IO
  !
  !
  ! parameters and variable used only in the legacy part  starting here
  !   refsatm_impl : specification of the reference site
  !             --- 0 : not reference
  !                 1 : reference solvent  2 : reference solute
  !   refspos : coordiantes of interaction site for reference structure
  !   sltcen : coordinate of the solute center
  !   sltqrn : quarternion for the solute orientation
  !   movmax : number of Monte Carlo moves
  !   trmax : maximum of translational Monte Carlo move
  !   agmax : maximum of orientational Monte Carlo move
  !   inscnd and inscfg are deprecated and used only in the legacy parts
  integer, dimension(:),   allocatable :: refsatm_impl
  real, dimension(:,:),    allocatable :: refspos
  real, save :: sltcen(3),sltqrn(0:3)
  integer, parameter :: movmax=10
  real, parameter :: trmax=0.20e0
  real, parameter :: agmax=0.10e0
  !   specifier to treat the reference as the total or as the system part only
  integer, parameter :: reftot=99, refsys=98
  ! parameters and variable used only in the legacy part  ending here
  !
  !
contains
  subroutine instslt(caltype, stat_weight_solute)
    use engmain, only: nummol, slttype, numslt, sltlist, iseed, CAL_REFS_FLEX
    use mpiproc, only: halt_with_error
    implicit none
    character(len=4), intent(in) :: caltype
    real, optional, intent(out) :: stat_weight_solute
    integer, save :: insml
    logical :: reject
    
    if(.not. present(stat_weight_solute)) then
       select case(caltype)
       case('init')
          ! sanity check of solute specification
          reject = .false.
          if(numslt /= 1) reject = .true.
          if((numslt == 1).and.(sltlist(1) /= nummol)) reject = .true.
          if(reject) call halt_with_error('ins_set')
          ! inserted solute is set to the last molecule in the system
          insml = sltlist(1)
          ! initialize the random number used for solute insertion
          call urand_init(iseed)
          ! opening the file for coordinate of flexible solute
          if(slttype == CAL_REFS_FLEX) call getsolute(caltype)
          return
       case('last')
          ! closing the file for coordinate of flexible solute
          if(slttype == CAL_REFS_FLEX) call getsolute(caltype)
          return
       case('proc')
          call halt_with_error('ins_bug')
       case default
          stop "Incorrect caltype in instslt"
       end select
    endif

    if(.not. present(stat_weight_solute)) call halt_with_error('ins_bug')

    if(slttype == CAL_REFS_FLEX) then
       ! get the coordinate and structure-specific weight of flexible solute
       call getsolute(caltype, insml, stat_weight_solute)
    else
       ! no structure-specific weight when the solute is rigid
       stat_weight_solute = 1.0e0
    endif

    reject = .true.
    do while(reject)
       call set_solute_origin(insml)
       call set_shift_com(insml, stat_weight_solute)
       call apply_orientation(insml)
       ! user-defined scheme to apply change / reject the solute configuration
       call insscheme(insml,reject)
    end do
    
    return
  end subroutine instslt

  subroutine set_solute_origin(insml)
    use engmain, only: insorigin, numsite, mol_begin_index, mol_end_index, &
                   bfcoord, sitepos, &
                   INSORG_ORIGIN, INSORG_NOCHANGE, INSORG_AGGCEN, INSORG_REFSTR
    use bestfit, only: com_aggregate
    implicit none
    integer, intent(in) :: insml
    integer :: molb, mole, nsite
    real :: syscen(3)

    nsite = numsite(insml)
    molb = mol_begin_index(insml)
    mole = mol_end_index(insml) 
    sitepos(1:3, molb:mole) = bfcoord(1:3, 1:nsite)

    select case(insorigin)
    case(INSORG_ORIGIN)
       syscen(:) = (/ 0., 0., 0. /)
       call set_solute_com(insml, syscen)     ! set solute COM to (0,0,0)
    case(INSORG_NOCHANGE)
       ! do nothing and use the coordinate as read from the file
    case(INSORG_AGGCEN)
       call com_aggregate(syscen)             ! get aggregate center (syscen)
       call set_solute_com(insml, syscen)     ! set solute COM to syscen
    case(INSORG_REFSTR)
       call reffit(insml)
    case default
       stop "Unknown insorigin"
    end select
  end subroutine set_solute_origin

  subroutine set_shift_com(insml, weight)
    use engmain, only: insposition, lwreg, upreg, &
                       cell, celllen, boxshp, SYS_NONPERIODIC, &
                       INSPOS_RANDOM, INSPOS_NOCHANGE, &
                       INSPOS_SPHERE, INSPOS_SLAB, INSPOS_GAUSS
    use mpiproc, only: halt_with_error
    implicit none
    integer, intent(in) :: insml
    real, intent(inout) :: weight
    
    integer :: i
    real :: com(3), syscen(3), r(3), norm, dir, t, s

    select case(insposition)
    case(INSPOS_RANDOM)
       ! fully random position within periodic box
       if(boxshp == SYS_NONPERIODIC) then    ! system has to be periodic
          call halt_with_error('ins_geo')
       endif
       do i = 1, 3
          call urand(r(i))
       end do
       com(:) = matmul(cell(:, :), r(:))
       call set_solute_com(insml, com)
       return
    case(INSPOS_NOCHANGE)
       ! fixed position as read from the file
       ! do nothing
       return
    case(INSPOS_SPHERE)
       ! spherical random position
       ! rejection method. Probability may not be good esp. lwreg ~ upreg.
       do
          do i = 1, 3
             call urand(r(i))
          end do
          r(:) = r(:) * 2.0 - 1.0
          norm = sum(r ** 2)
          if(norm < 1 .and. norm > (lwreg/upreg) ** 2) exit
       end do
       com(:) = r(:) * upreg
       call shift_solute_com(insml, com)
       return
    case(INSPOS_SLAB)
       ! slab random position
       if(boxshp == SYS_NONPERIODIC) then    ! system has to be periodic
          call halt_with_error('ins_geo')
       endif
       call urand(r(1))
       call urand(r(2))
       
       call urand(t)
       call urand(dir)
       
       t = t * (upreg - lwreg) + lwreg
       if(dir > 0.5) t = -t
       r(3) = t / celllen(3)
       
       com(:) = matmul(cell(:, :), r(:))
       call shift_solute_com(insml, com)
       return
    case(INSPOS_GAUSS)
       ! 50% mixture of uniform distribution and weighted distribution
       call uniform_gauss_mixture(com, weight)
       call shift_solute_com(insml, com)
       return
    case default
       stop "Unknown insposition"
    end select

  contains
    subroutine shift_solute_com(insml, com)
      use engmain, only: numsite, mol_begin_index, mol_end_index, sitepos
      implicit none
      integer, intent(in) :: insml
      real, intent(in) :: com(3)
      integer :: insb, inse, i, n
      n = numsite(insml)
      insb = mol_begin_index(insml)
      inse = mol_end_index(insml)
      do i = 1, 3
         sitepos(i, insb:inse) = sitepos(i, insb:inse) + com(i)
      end do
    end subroutine shift_solute_com

    ! FIXME: this routine does not work for skewed periodic box
    subroutine uniform_gauss_mixture(com, weight)
      use mpiproc, only: myrank
      use engmain, only: PI, celllen, upreg
      implicit none
      real, intent(out) :: com(3), weight

      real, parameter :: uniform_ratio = 0.5

      integer :: i
      real :: scaled_coord(3), sqsum
      real :: r

      real, save :: l_of_sigma(3), z0, z1
      logical, save :: use_uniform = .false.
      logical, save :: first_time = .true.

      if(first_time) then
         l_of_sigma(:) = (celllen(:) / 2) / upreg
         if(myrank == 0) print *, "Lx/2 / sigma = ", l_of_sigma(:)
         z0 = 8 * l_of_sigma(1) * l_of_sigma(2) * l_of_sigma(3)
         z1 = (sqrt(PI) ** 3) * erf(l_of_sigma(1)) * erf(l_of_sigma(2)) * erf(l_of_sigma(3))
         first_time = .false.
      endif

      use_uniform = .not. use_uniform
      if(use_uniform) then
         ! use random position
         do i = 1, 3
            call urand(r)
            scaled_coord(i) = (2 * r - 1) * l_of_sigma(i)
         end do

      else
         ! for weighted insertion
         ! W = x^2
         ! get three N(0, 1) values
         do i = 1, 3
            do 
               r = nrand() / sqrt(2.0)
               if(abs(r) < l_of_sigma(i)) exit
            end do
            scaled_coord(i) = r
         end do

      endif
      sqsum = sum(scaled_coord(:) ** 2)
      ! from WHAM
      weight = 1.0 / (uniform_ratio * 1 / z0 + (1 - uniform_ratio) * exp(-sqsum) / z1)
      com(:) = matmul(cell(:, :) , scaled_coord(:) / (l_of_sigma(:) * 2))
    end subroutine uniform_gauss_mixture
  end subroutine set_shift_com

  subroutine set_solute_com(insml, com)
    use engmain, only: numsite, mol_begin_index, mol_end_index, &
                       sitepos, sitemass
    use bestfit, only: com_shift, com_unshift
    implicit none
    integer, intent(in) :: insml
    real, intent(in) :: com(3)
    integer :: insb, inse, i, n
    real :: tempcom(3)
    
    n = numsite(insml)
    insb = mol_begin_index(insml)
    inse = mol_end_index(insml)
    
    call com_shift(n, sitepos(1:3, insb:inse), sitemass(insb:inse), tempcom)
    do i = 1, 3
       sitepos(i, insb:inse) = sitepos(i, insb:inse) + com(i)
    end do
  end subroutine set_solute_com

  subroutine apply_orientation(insml)
    use engmain, only: insorient, INSROT_RANDOM, INSROT_NOCHANGE, &
                       numsite, mol_begin_index, mol_end_index, &
                       sitepos, sitemass
    use quaternion, only: rotate_inplace
    use bestfit, only: com_shift, com_unshift
    implicit none
    integer, intent(in) :: insml
    integer :: insb, inse, i, n
    real :: com(3), randq(0:3)
    
    n = numsite(insml)
    insb = mol_begin_index(insml)
    inse = mol_end_index(insml)

    select case(insorient)
    case(INSROT_RANDOM)     ! random orientation
       call com_shift(n, sitepos(1:3, insb:inse), sitemass(insb:inse), com)

       ! get random quaternion (prob. success ~ 0.3)
       ! if this is really bad implement:
       ! Marsaglia, G. "Choosing a Point from the Surface of a Sphere."
       !   Ann. Math. Stat. 43, 645-646, 1972.
       do
          do i = 0, 3
             call urand(randq(i))
          end do
          randq(:) = randq(:) * 2.0 - 1.0
          if(sum(randq ** 2) < 1) exit
       end do

       randq(:) = randq(:) / sqrt(sum(randq**2)) ! set on unit hyper-sphere surface
       call rotate_inplace(numsite(insml), sitepos(1:3, insb:inse), randq)

       call com_unshift(n, sitepos(1:3, insb:inse), sitemass(insb:inse), com)
       return
    case(INSROT_NOCHANGE)   ! no orientational change
       return
    case default
       stop "Unknown insorient"
    end select
  end subroutine apply_orientation
  !
  ! user-defined scheme to specify inserted molecule
  ! user may reject snapshot by specifying out_of_range to .true.
  ! or set coordinate in sitepos manually
  subroutine insscheme(insml, out_of_range)
    use engmain, only: nummol, numsite, specatm, sitepos, cell
    use bestfit, only: center_of_mass, rmsd
    implicit none
    integer, intent(in) :: insml
    logical, intent(out) :: out_of_range
    integer :: stmax, sid, ati
    out_of_range = .false.
    stmax=numsite(insml)
    return
  end subroutine insscheme


  subroutine getsolute(caltype, insml, stat_weight)
    use trajectory, only: open_trajectory, close_trajectory
    use engmain, only: slttype, wgtins, numsite, bfcoord, stdout, CAL_REFS_FLEX
    use OUTname, only: OUTconfig, solute_trajectory
    use mpiproc
    implicit none
    character(len=4), intent(in) :: caltype
    integer, optional, intent(in) :: insml
    real, optional, intent(out) :: stat_weight
    logical, save :: read_weight
    integer :: stmax, dumint, ioerr
    real :: dumcl(3,3)
    real, dimension(:,:), allocatable :: psite
!
    if(slttype /= CAL_REFS_FLEX) call halt_with_error('ins_bug')
!
    if((.not. present(insml)) .and. &
       (.not. present(stat_weight))) then
       select case(caltype)
       case('init')
          if(wgtins == 1) then
             read_weight = .true.
          else
             read_weight = .false.
          endif
          if(myrank /= 0) return
          call open_trajectory(solute_trajectory, slttrj)
          if(read_weight) open(unit=sltwgt_io, file=sltwgt, status='old')
          return
       case('last')
          if(myrank /= 0) return
          call close_trajectory(solute_trajectory)
          if(read_weight) close(unit=sltwgt_io)
          return
       case('proc')
          call halt_with_error('ins_bug')
       case default
          stop "Incorrect caltype in getsolute"
       end select
    endif
!
    if((.not. present(insml)) .or. &
       (.not. present(stat_weight))) call halt_with_error('ins_bug')
!
    stmax=numsite(insml)
    allocate( psite(3,stmax) )
    if(myrank == 0) then
       call OUTconfig(psite,dumcl,stmax,0,'solute','trjfl_read')
       if(read_weight) then
          read(sltwgt_io, *, iostat = ioerr) dumint, stat_weight
          if(ioerr /= 0) then      ! wrap around
             rewind(sltwgt_io)
             read(sltwgt_io, *, iostat = ioerr) dumint, stat_weight
             if(ioerr /= 0) then
                write(stdout,*) " The weight file (", sltwgt, ") is ill-formed"
                call mpi_setup('stop')
                stop
             endif
          endif
       else
          stat_weight = 1.0e0
       endif
    endif

#ifndef noMPI
    ! distribute to non rank-0 nodes
    call mpi_bcast(psite, 3 * stmax, mpi_double_precision, &
                   0, mpi_comm_activeprocs, ierror)
    call mpi_bcast(stat_weight, 1, mpi_double_precision, &
                   0, mpi_comm_activeprocs, ierror)
#endif
    bfcoord(1:3,1:stmax)=psite(1:3,1:stmax)
    deallocate( psite )

  end subroutine getsolute


  ! returns random value from [0,1)
  ! Any sane compiler implements random_number
  ! (which is included in fortran 95 standards)
  subroutine urand(rndm)          ! uniform random number generator
    implicit none
    real, intent(out) :: rndm
    call random_number(rndm)
  end subroutine urand

  ! Normal random variable N(0,1)
  ! uses Box-Muller method
  real function nrand()
    use engmain, only: PI
    implicit none
    real :: r1, r2
    call urand(r1)
    call urand(r2)
    ! get (0,1] instead of [0, 1)
    r1 = 1 - r1
    nrand = sqrt(-2.0 * log(r1)) * cos(2 * PI * r2)
  end function nrand

  subroutine urand_init(seed)
    use mpiproc, only: myrank
    implicit none
    integer, intent(in) :: seed
    integer :: seedsize
    integer, allocatable :: seedarray(:)

    call random_seed(size = seedsize)
    allocate(seedarray(seedsize))
    seedarray(:) = 1

    seedarray(1) = myrank + seed
    if(seed == 0) call system_clock(count = seedarray(1))

    call random_seed(put = seedarray)
    deallocate(seedarray)
  end subroutine urand_init
!
!
! fit to reference structure
  subroutine reffit(ligmol)
    use engmain, only: refmlid, numsite, nummol, &
         mol_begin_index, mol_end_index, &
         sitepos, sitemass, bfcoord
    use quaternion, only: rotate
    use bestfit, only: fit_a_rotate_b, fit
    implicit none
    integer, intent(in) :: ligmol

    real, allocatable, save :: ref_solv_pos(:, :)  ! solvent structure in reference
    real, allocatable, save :: ref_lig_pos(:, :)   ! ligand structure in reference
    real, allocatable, save :: works(:, :), workl(:, :)
    real, allocatable, save :: solv_mass(:) ! masked mass vector;
    real, allocatable, save :: lig_mass(:)  ! atoms not used for fitting is zeroed

    integer :: natom_solv, natom_lig
    integer :: solvmol
    integer :: solv_begin, solv_end
    integer :: lig_begin, lig_end
    integer :: i
    
    solvmol = -1
    do i = 1, nummol
       if(refmlid(i) == 1) solvmol = i
       if(refmlid(i) == 2 .and. i /= ligmol) stop "insmol / refmlid is inconsistent"
    end do
    if(solvmol == -1) stop "reffit: failed to find reference solvent"

    natom_solv = numsite(solvmol)
    natom_lig = numsite(ligmol)
    
    solv_begin = mol_begin_index(solvmol)
    solv_end = mol_end_index(solvmol)

    lig_begin = mol_begin_index(ligmol)
    lig_end = mol_end_index(ligmol)

    if(.not. allocated(ref_solv_pos)) then
       ! first time call
       allocate(ref_solv_pos(3, natom_solv), solv_mass(natom_solv))
       allocate(ref_lig_pos(3, natom_lig), lig_mass(natom_lig))
       allocate(works(3, natom_solv), workl(3, natom_lig))
       
       ! initialize mass first
       solv_mass(:) = sitemass(solv_begin:solv_end)
       lig_mass(:) = sitemass(lig_begin:lig_end)

       open(unit=refio, file=reffile, status='old')
       call load_structure(natom_solv, ref_solv_pos, solv_mass)
       call load_structure(natom_lig, ref_lig_pos, lig_mass)
       close(refio)
    endif
    
    ! first fit reference solvent to current solvent, to get "reference ligand position"
    call fit_a_rotate_b(natom_solv, &
         sitepos(1:3, solv_begin:solv_end), ref_solv_pos, solv_mass, &
         natom_lig, &
         ref_lig_pos, workl)
    ! then fit ligand structure read from file
    call fit(natom_lig, workl, bfcoord, lig_mass, sitepos(1:3, lig_begin:lig_end))
  contains
    subroutine load_structure(n, position, mass)
      implicit none
      integer, intent(in) :: n
      real, intent(out) :: position(3, n)
      real, intent(inout) :: mass(n)
      real :: crd(3)
      character(len=6) :: header
      integer :: i, ix

      do ix = 1, n
         do
            ! skip until ATOM/HETATM lines
            read(refio, '(A6)', advance='no') header
            if(header == 'ATOM  ' .or. header == 'HETATM') exit
            read(refio, *)
         end do
         read(refio, '(24X, 3F8.3)') (crd(i), i = 1, 3)
         position(1:3, ix) = crd(1:3)
      end do

      ! set mass = 0 to mask fitting
      ! user can implement his own special selection rule
      ! (e.g. by using B-factor, etc.)
      ! default: hydrogen is masked
      do i = 1, n
         if(mass(i) < 1.1) mass(i) = 0.0
      end do
    end subroutine load_structure
  end subroutine reffit





! legacy part follows
! coordinate: make the solute coordinate through translation and orientation
  subroutine coordinate(i,pcom,qrtn)
    use engmain, only: nummol,numatm,inscfg,numsite,bfcoord,specatm,sitepos
    implicit none
    integer stmax,sid,ati,m,k,i
    real pcom(3),qrtn(0:3),rotmat(3,3),rst
    stmax=numsite(i)
    if(inscfg.eq.0) call getrot(qrtn,rotmat)        ! rotation matrix
    do  sid=1,stmax
       ati=specatm(sid,i)
       do  m=1,3
          if(inscfg.eq.0) then                        ! translation + rotation
             rst=dot_product(rotmat(:, m), bfcoord(:, sid))
             rst=rst+pcom(m)
          endif
          if(inscfg.eq.1) rst=pcom(m)+bfcoord(m,sid)  ! translation only
          if(inscfg.eq.2) rst=bfcoord(m,sid)          ! no change
          sitepos(m,ati)=rst
       end do
    end do
    return
  end subroutine coordinate
! getrot: convert quaternion to rotation matrix,
!         used only in coordinate, refmc, sltmove, and refsdev
  subroutine getrot(qrtn,rotmat)
    implicit none
    real, intent(in) :: qrtn(0:3)
    real, intent(out) :: rotmat(3,3)
    rotmat(1,1)=qrtn(0)*qrtn(0)+qrtn(1)*qrtn(1)&
         -qrtn(2)*qrtn(2)-qrtn(3)*qrtn(3)
    rotmat(2,2)=qrtn(0)*qrtn(0)-qrtn(1)*qrtn(1)&
         +qrtn(2)*qrtn(2)-qrtn(3)*qrtn(3)
    rotmat(3,3)=qrtn(0)*qrtn(0)-qrtn(1)*qrtn(1)&
         -qrtn(2)*qrtn(2)+qrtn(3)*qrtn(3)
    rotmat(1,2)=2.0e0*(qrtn(1)*qrtn(2)+qrtn(0)*qrtn(3))
    rotmat(2,1)=2.0e0*(qrtn(1)*qrtn(2)-qrtn(0)*qrtn(3))
    rotmat(1,3)=2.0e0*(qrtn(1)*qrtn(3)-qrtn(0)*qrtn(2))
    rotmat(3,1)=2.0e0*(qrtn(1)*qrtn(3)+qrtn(0)*qrtn(2))
    rotmat(2,3)=2.0e0*(qrtn(2)*qrtn(3)+qrtn(0)*qrtn(1))
    rotmat(3,2)=2.0e0*(qrtn(2)*qrtn(3)-qrtn(0)*qrtn(1))
  end subroutine getrot
!
!
! sltpstn: identifying the solute position according to the inscnd value
  ! FIXME: cleanup
  subroutine sltpstn(sltstat,pcom,type,tagslt)
    use engmain, only: nummol,numatm,boxshp,inscnd,inscfg,hostspec,&
         moltype,numsite,specatm,sitepos,cell,invcl,lwreg,upreg
    use mpiproc, only: halt_with_error
    use bestfit, only: com_aggregate
    implicit none
    integer sltstat,tagslt,stmax,sid,ati,pti,i,m,k,centag(numatm)
    real rdum,clm(3),pcom(3),qrtn(0:3),rst,dis,syscen(3),elen
    character*6 type

    if(inscfg.eq.2) then
       sltstat=1
       return
    endif

    if(type.ne.'insert') stop "BUG"
    !
    if(inscnd.eq.0) then   ! solute with random position
       if(boxshp.eq.0) then    ! system has to be periodic
          call halt_with_error('ins_geo')
       endif
       do k=1,3
          call URAND(rdum)
          clm(k)=rdum-0.50e0
       end do
       do m=1,3
          pcom(m)=dot_product(cell(m, :), clm(:))
       end do
    endif
    !
    if((inscnd.eq.1).or.(inscnd.eq.2)) then     ! system center
       call com_aggregate(syscen)
    endif
    !
    if(inscnd.eq.1) then   ! solute in spherical object or isolated droplet
       call rndmvec('p',qrtn,(lwreg/upreg))

       pcom(:)=upreg*qrtn(1:3)+syscen(:)
       dis=0.0e0
       do m=1,3
          rst=pcom(m)-syscen(m)
          dis=dis+rst*rst
       end do
       dis=sqrt(dis)
    endif
    !
    if(inscnd.eq.2) then   ! solute in slab geometry
       if(boxshp.eq.0) then    ! system has to be periodic
          call halt_with_error('ins_geo')
       endif
       elen=0.0e0
       do m=1,3
          elen=elen+cell(m,3)*cell(m,3)
       end do
       elen=sqrt(elen)

       do k=1,2
          call URAND(rdum)
          clm(k)=rdum-0.50e0
       end do
       do m=1,3
          rst=0.0e0
          do  k=1,2
             rst=rst+cell(m,k)*clm(k)
          end do
          pcom(m)=rst
       end do
       call URAND(rdum)
       rst=lwreg+rdum*(upreg-lwreg)
       call URAND(rdum)
       if(rdum.le.0.50e0) rst=-rst
       dis=0.0e0
       do m=1,3
          dis=dis+invcl(3,m)*syscen(m)
       end do
       dis=dis+rst/elen
       do m=1,3
          pcom(m)=pcom(m)+dis*cell(m,3)
       end do

       rst=0.0e0
       do m=1,3
          rst=rst+invcl(3,m)*(pcom(m)-syscen(m))
       end do
       dis=abs(rst)*elen
    endif
    !
    if(inscnd.eq.3) then   ! solute against reference structure
       if(type.eq.'insert') call halt_with_error('ins_bug')
       call refsdev(dis,2,'system')
    endif
    !
    if(inscnd.eq.0) sltstat=1
    if(inscnd.gt.0) then
       if((lwreg.le.dis).and.(dis.le.upreg)) then
          sltstat=1
       else
          call halt_with_error('ins_bug')
       endif
    endif
    !
    return
  end subroutine sltpstn
!
! rndmvec: generator of random vector on sphere, used only in sltpstn
  subroutine rndmvec(vectp,qrtn,lwbnd)
    implicit none
    character vectp
    integer m,inim
    real qrtn(0:3),lwbnd,rdum,factor
    if(vectp.eq.'p') inim=1
    if(vectp.eq.'q') inim=0
    factor=2.0e0
    qrtn(:) = 0.0e0
    do while((factor.gt.1.0e0).or.(factor.lt.lwbnd))
       do m=inim,3
          call URAND(rdum)
          qrtn(m)=2.0e0*rdum-1.0e0
       end do
       factor=sqrt(sum(qrtn(inim:3) ** 2))
    end do
    if(vectp.eq.'q') then
       do m=0,3
          qrtn(m)=qrtn(m)/factor
       end do
    endif
  end subroutine rndmvec
!
!
! refmc: subroutine to generate an insertion configuration at inscnd = 3
  subroutine refmc
    use engmain, only: nummol,numatm,slttype,numsite,refmlid,&
         lwreg,upreg,bfcoord,specatm,sitepos
    use bestfit
    integer i,k,m,q,ati,rfi,sid,stmax
    real xst(3),centg(3),cenrf(3)
    real factor,bfqrn(0:3),rtmbf(3,3)
    real, dimension(:,:), allocatable :: sltsite
    character*4 caltype
    character*8 atmtype,dump
    character eletype

    if(caltype.eq.'init') then
       allocate( refsatm_impl(numatm),refspos(3,numatm) )
       refsatm_impl(:) = 0
       refspos(:, :) = 0.0e0

       !  read the reference structure
       open(unit=refio,file=reffile,status='old')
       do i=1,nummol
          rfi=refmlid(i)
          if(rfi.ne.0) then
             stmax=numsite(i)
             do sid=1,stmax
                read(refio, '(12X, A4, 14X, 3F8.3)') atmtype, (xst(m), m=1,3)
                ! FIXME: this element selection looks very ugly
                atmtype = trim(atmtype)
                eletype = atmtype(1:1)
                ! atom name can be like "1HG"
                if('0' <= eletype .and. eletype <= '9') eletype = atmtype(2:2)
                ati=specatm(sid,i)
                if(eletype.eq.'H') refsatm_impl(ati) = 0      ! hydrogen atom
                if(eletype.ne.'H') refsatm_impl(ati) = rfi    ! heavy atom
                do m=1,3
                   refspos(m,ati)=xst(m)
                end do
             end do
          endif
       end do
       close(refio)
       !
       !  build the initial configuration of the solute
       if(slttype.ge.2) then
          call dispref(centg,cenrf,bfqrn,1)     ! matching reference
          call getrot(bfqrn,rtmbf)
          stmax=numsite(nummol)                 ! inserted solute molecule
          do  sid=1,stmax
             ati=specatm(sid,nummol)
             do  m=1,3
                factor=0.0e0
                do  k=1,3
                   factor=factor+rtmbf(k,m)*(refspos(k,ati)-cenrf(k))
                end do
                sitepos(m,ati)=factor+centg(m)
             end do
          end do
          call dispref(sltcen,cenrf,sltqrn,2)   ! solute center & orientation
       endif
    endif
    !
    if(caltype.eq.'inst') then                ! inserted solute molecule
       stmax=numsite(nummol)
       allocate( sltsite(3,stmax) )
       do sid=1,stmax
          ati=specatm(sid,nummol)
          do m=1,3
             sltsite(m,sid)=sitepos(m,ati)
          end do
          do m=1,3
             sitepos(m,ati)=bfcoord(m,sid)
          end do
       end do
       call refcen(xst,sitepos,2)      ! solute center in bfcoord
       do sid=1,stmax
          ati=specatm(sid,nummol)
          do m=1,3
             bfcoord(m,sid)=bfcoord(m,sid)-xst(m)
          end do
       end do
       call dispref(centg,cenrf,bfqrn,2)       ! matching solute
       call getrot(bfqrn,rtmbf)
       k=0
       do q=1,movmax
          call sltmove(sltcen,sltqrn,rtmbf,'frwd')
          call refsdev(factor,2,'extend')
          if((factor.lt.lwreg).or.(factor.gt.upreg)) then
             k=k+1
             call sltmove(sltcen,sltqrn,rtmbf,'back')
          endif
       end do
       if(k.eq.movmax) then                    ! when all MC are rejected
          do sid=1,stmax
             ati=specatm(sid,nummol)
             sitepos(:,ati)=sltsite(:,sid)
          end do
       endif
       deallocate( sltsite )
    endif
    !
    return
  end subroutine refmc
!
! refcen: center of reference, used only in refmc and dispref
  subroutine refcen(cen,posatm,refmol)
    use engmain, only: nummol,numatm,numsite,specatm
    implicit none
    integer refmol,rfyn,i,m,ati,rfi,sid,stmax
    real cen(3),posatm(3,numatm),totwgt
    totwgt=0.0e0
    cen(:)=0.0e0
    do i=1,nummol
       call getrfyn(i,rfi,rfyn,refmol)
       if(rfyn.eq.1) then
          stmax=numsite(i)
          do sid=1,stmax
             if(refsatm(sid,i).eq.rfi) then
                ati=specatm(sid,i)
                totwgt=totwgt+1.0e0
                cen(:) = cen(:) + posatm(:,ati)
             endif
          end do
       endif
    end do
    cen(:) = cen(:) / totwgt
    return
  end subroutine refcen
!
!
  integer function refsatm(site, mol)
    use engmain, only: specatm
    implicit none
    integer, intent(in) :: site, mol
    refsatm = refsatm_impl(specatm(site, mol))
  end function refsatm
!
! sltmove: moving the solute at inscnd = 3, used only in refmc
  subroutine sltmove(pcen,qrtn,rtmbf,mvtype)
    use engmain, only: nummol,numatm,numsite,bfcoord,specatm,sitepos
    implicit none
    integer m,k,q,ati,sid,stmax,ax1,ax2
    real pcen(3),qrtn(0:3),rtmbf(3,3)
    real rdum,factor,rotmat(3,3),rtmmov(3,3),axis(3)
    character*4 mvtype
    real, save :: odcn(3),oqrn(0:3)
    stmax=numsite(nummol)                       ! inserted solute molecule
    do m=1,3
       if(mvtype.eq.'frwd') odcn(m)=pcen(m)
       if(mvtype.eq.'back') pcen(m)=odcn(m)
    end do
    do m=0,3
       if(mvtype.eq.'frwd') oqrn(m)=qrtn(m)
       if(mvtype.eq.'back') qrtn(m)=oqrn(m)
    end do
    if(mvtype.eq.'frwd') then
       do m=1,3                             ! translation
          call URAND(rdum)
          pcen(m)=pcen(m)+trmax*(rdum+rdum-1.0e0)
       end do
       factor=2.0e0                              ! rotation
       do while(factor.gt.1.0e0)
          do m=1,3
             call URAND(rdum)
             axis(m)=rdum+rdum-1.0e0
          end do
          factor=axis(1)*axis(1)+axis(2)*axis(2)+axis(3)*axis(3)
       end do
       factor=sqrt(factor)
       do m=1,3
          axis(m)=axis(m)/factor
       end do
       call URAND(rdum)
       ax1=floor(rdum * 4)
       ax2=ax1
       do while(ax1.eq.ax2)
          call URAND(rdum)
          ax2=floor(rdum * 4)
       end do
       call URAND(rdum)
       factor=agmax*(rdum+rdum-1.0e0)
       qrtn(ax1)=cos(factor)*oqrn(ax1)-sin(factor)*oqrn(ax2)
       qrtn(ax2)=sin(factor)*oqrn(ax1)+cos(factor)*oqrn(ax2)
       factor=0.0e0
       do m=0,3
          factor=factor+qrtn(m)*qrtn(m)
       end do
       factor=sqrt(factor)
       do m=0,3
          qrtn(m)=qrtn(m)/factor
       end do
    endif
    call getrot(qrtn,rtmmov)                    ! site coordinate
    do m=1,3
       do k=1,3
          factor=0.0e0
          do q=1,3
             factor=factor+rtmmov(q,m)*rtmbf(q,k)
          end do
          rotmat(m,k)=factor
       end do
    end do
    do sid=1,stmax
       ati=specatm(sid,nummol)
       do m=1,3
          factor=0.0e0
          do k=1,3
             factor=factor+rotmat(m,k)*bfcoord(k,sid)
          end do
          sitepos(m,ati)=factor+pcen(m)
       end do
    end do
    return
  end subroutine sltmove
!
! dispref: best fit to the reference structure, used only in refmc and refsdev
  subroutine dispref(centg,cenrf,qrtn,refmol)
    use engmain, only: nummol,numatm,numsite,specatm,sitepos
    implicit none
    integer refmol,rfyn,i,m,k,ati,rfi,sid,stmax
    real centg(3),cenrf(3),qrtn(0:3),qrtmat(4,4),xsm(3),xrf(3)
    real totm,xp,yp,zp,xm,ym,zm,dumv(4),work(16)
    call refcen(centg,sitepos,refmol)
    call refcen(cenrf,refspos,refmol)
    totm=0.0e0
    qrtmat(:, :) = 0.0e0
    do i=1,nummol
       call getrfyn(i,rfi,rfyn,refmol)
       if(rfyn.eq.1) then
          stmax=numsite(i)
          do sid=1,stmax
             if(refsatm(sid,i).eq.rfi) then
                ati=specatm(sid,i)
                totm=totm+1.0e0
                do m=1,3
                   xsm(m)=sitepos(m,ati)-centg(m)
                end do
                do m=1,3
                   xrf(m)=refspos(m,ati)-cenrf(m)
                end do
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
          end do
       endif
    end do
    ! symmetrize
    do m=1,3
       do k=m+1,4
          qrtmat(m,k)=qrtmat(k,m)
       end do
    end do
    do m=1,4
       do k=1,4
          qrtmat(k,m)=qrtmat(k,m)/totm
       end do
    end do
    call DSYEV('V','U',4,qrtmat,4,dumv,work,16,i)
    do m=0,3
       qrtn(m)=qrtmat(m+1,1)   ! eigenvector for the minimum eigenvalue
    end do
    totm=0.0e0
    do m=0,3
       totm=totm+qrtn(m)*qrtn(m)
    end do
    totm=sqrt(totm)
    do m=0,3
       qrtn(m)=qrtn(m)/totm
    end do
    return
  end subroutine dispref
!
! refsdev: RMSD calculation, used only in refmc and sltpstn
  subroutine refsdev(strchr,refmol,caltype)
    use engmain, only: nummol,numatm,numsite,specatm,sitepos
    use mpiproc, only: halt_with_error
    implicit none
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
       call halt_with_error('ins_bug')
    end select
    call getrot(qrtn,rotmat)
    strchr=0.0e0
    totm=0.0e0
    do i=1,nummol
       call getrfyn(i,rfi,rfyn,refmol)
       if(rfyn.eq.1) then
          stmax=numsite(i)
          do sid=1,stmax
             if(refsatm(sid,i).eq.rfi) then
                ati=specatm(sid,i)
                totm=totm+1.0e0
                do m=1,3
                   factor=0.0e0
                   do k=1,3
                      factor=factor+rotmat(m,k)*(sitepos(k,ati)-centg(k))
                   end do
                   xsm(m)=factor+cenrf(m)
                end do
                do m=1,3
                   factor=xsm(m)-refspos(m,ati)
                   strchr=strchr+factor*factor
                end do
             endif
          end do
       endif
    end do
    strchr=sqrt(strchr/totm)
    return
  end subroutine refsdev
!
! getrfyn: judging whether the molecules is employed for reference calculation,
!          used only in refmc, dispref, and refsdev
  subroutine getrfyn(i,rfi,rfyn,refmol)
    use engmain, only: nummol,sluvid,refmlid
    use mpiproc, only: halt_with_error
    implicit none
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
       call halt_with_error('ins_bug')
    end select
    return
  end subroutine getrfyn
!
end module
