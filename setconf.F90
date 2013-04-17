! -*- F90 -*-
! ERmod - Eneregy Representation Module
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
!
!
! renaming outside parameters and parameters to avoid conflict
module OUTname
  use trajectory, only: handle
!
  implicit none
  integer, parameter :: iotrj=99                 ! trajectory file IO
  character(*), parameter :: trjfile='HISTORY'   ! trajectory filename
  integer, parameter :: mdinf=89                 ! MD info file IO
  character(*), parameter :: inffile='MDinfo'    ! MD info filename

  integer OUTens,OUTbxs,OUTcmb,OUTclt,OUTspo
  real OUTtemp,OUTelc,OUTlwl,OUTupl,OUTscr
  integer OUTew1,OUTew2,OUTew3,OUTms1,OUTms2,OUTms3

  integer OUTnrun,OUTntype
  integer, dimension(:), allocatable :: OUTnmol,OUTsite
  real, dimension(:), allocatable :: OUTstmass,OUTcharge
  real, dimension(:), allocatable :: OUTljene,OUTljlen
  logical :: use_mdlib

  type(handle) :: history_trajectory
  type(handle) :: solute_trajectory

contains

  subroutine OUTinitial
    real, save :: sgmcnv,chgcnv,engcnv,lencnv    ! unit conversion factor
    sgmcnv=2.0e0**(5.0e0/6.0e0)                  ! from Rmin/2 to sigma
    chgcnv=1.0e0/1.60217653e-19                  ! from C to elementary
    engcnv=6.0221415e23/4.184e3                  ! from J to kcal/mol
    lencnv=1.0e1                                 ! from nm to Angstrom
  end subroutine OUTinitial

  subroutine opentrj
    use trajectory, only: open_trajectory
    call open_trajectory(history_trajectory, trjfile)
  end subroutine opentrj

  subroutine initconf
    call OUTinitial
  end subroutine initconf

  subroutine closetrj
    use trajectory, only: close_trajectory
    call close_trajectory(history_trajectory)
  end subroutine closetrj

  subroutine finiconf
  end subroutine finiconf
!
  ! Initialization 2nd phase - read MDinfo
  ! when ermod is built into MD program, read topologies from mother MD program
  subroutine OUTrename
    implicit none
    integer TotAtm
    open(unit=mdinf,file=inffile,status='old')
    read(mdinf,*) OUTnrun,OUTntype
    allocate( OUTnmol(OUTntype),OUTsite(OUTntype) )
    read(mdinf,*) OUTnmol(1:OUTntype)
    read(mdinf,*) OUTsite(1:OUTntype)
    TotAtm = sum(OUTnmol(1:OUTntype) * OUTsite(1:OUTntype))
    allocate( OUTstmass(TotAtm),OUTcharge(TotAtm) )
    allocate( OUTljene(TotAtm),OUTljlen(TotAtm) )
    close(mdinf)
    return
  end subroutine OUTrename

  ! read system setup (coulomb rule, LJ, etc..)
  ! used when the ermod program is built into MD program and runs on-the-fly
  subroutine OUTintprm
  end subroutine OUTintprm

  ! Get molecular configuration
  subroutine OUTconfig(OUTpos,OUTcell,OUTatm,OUTbox,particle_type,calc_type)
    use mpiproc, only: halt_with_error
    use trajectory, only: read_trajectory, open_trajectory, close_trajectory
    implicit none
    integer, intent(in) :: OUTatm, OUTbox
    character*6, intent(in) :: particle_type
    character*10, intent(in) :: calc_type
    real, intent(out) :: OUTpos(3, OUTatm), OUTcell(3, 3)
    integer :: status
    !
    select case(calc_type)
    case('trjfl_read')              ! reading from a trajectory file
       ! do nothing
    case('on_the_fly')              ! on-the-fly calculation during MD
       stop "Unavailable in the present version"
    case default
       stop "Unknown calc_type in OUTconfig"
    end select
    !
    OUTcell(:, :) = 0.0e0
    !
    select case(particle_type)
    case('system')                  ! reading the HISTORY file
       call read_trajectory(history_trajectory, OUTatm, (OUTbox == 1), OUTpos, OUTcell, status)
       if(status /= 0) call halt_with_error("set_trj")
    case('solute')                  ! reading the SltConf file
       call read_trajectory(solute_trajectory, OUTatm, .false., OUTpos, OUTcell, status)
       if(status /= 0) then
          ! wrap around
          call close_trajectory(solute_trajectory)
          call open_trajectory(solute_trajectory, "SltConf")
          call read_trajectory(solute_trajectory, OUTatm, .false., OUTpos, OUTcell, status)
          if(status /= 0) then
             stop "Failed to reload solute trajectory"
          endif
       endif
    case default
       stop "Unknown particle_type"
    end select
    !
    return
  end subroutine
!
end module
!
!
module setconf
  implicit none
  character(len=*), parameter :: weight_file = "SysWght"
  integer, parameter :: weight_io = 33
  character(len=*), parameter :: perm_file = "PermIndex"
  integer, parameter :: perm_io = 75
contains
!
!  setting molecular simulation parameters
!
  subroutine iniparam
    use engmain, only: init_params, &
         iseed, &
         skpcnf, corrcal, selfcal, &
         slttype, wgtslf, wgtins, wgtsys, estype, boxshp, &
         sltspec, hostspec, refspec, ljformat, ljswitch, &
         insorigin, insposition, insorient, insstructure, &
         sltpick, refpick, inscnd, inscfg, &     ! deprecated   
         lwreg, upreg, lwstr, upstr, &
         inptemp, temp, &
         engdiv, maxins, &
         intprm, elecut, lwljcut, upljcut, &
         cmbrule, cltype, screen, ewtoler, splodr, plmode, &
         ew1max, ew2max, ew3max, ms1max, ms2max, ms3max, &
         block_threshold, &
         force_calculation, &
         SYS_NONPERIODIC, SYS_PERIODIC, EL_PME, ES_NPT, &
         CAL_SOLN, CAL_REFS_RIGID, CAL_REFS_FLEX, &
         INSORG_ORIGIN, INSORG_NOCHANGE, INSORG_AGGCEN, INSORG_REFSTR, &
         INSPOS_RANDOM, INSPOS_NOCHANGE, &
         INSPOS_SPHERE, INSPOS_SLAB_GENERIC, INSPOS_SLAB_SYMMETRIC, &
         INSPOS_RMSD, INSPOS_GAUSS, &
         INSROT_RANDOM, INSROT_NOCHANGE, &
         INSSTR_NOREJECT, INSSTR_RMSD
    use OUTname, only: OUTintprm, &              ! from outside
         OUTens, OUTbxs, OUTtemp, &              ! from outside
         OUTelc, OUTlwl, OUTupl, &               ! from outside
         OUTcmb, OUTclt, OUTscr, OUTspo, &       ! from outside
         OUTew1, OUTew2, OUTew3, &               ! from outside
         OUTms1, OUTms2, OUTms3                  ! from outside
    use mpiproc, only: mpi_info, halt_with_error                     ! MPI
    implicit none
    real, parameter :: tiny = 1.0e-20
    real :: real_seed
    character(len=3) :: scrtype
    call mpi_info                                                    ! MPI

    intprm=1                                     ! trajectory reading
    select case(intprm)
    case(0)  ! on-the-fly reading from parent MD setup (currently ineffective)
      call OUTintprm
      estype = OUTens  ; boxshp = OUTbxs  ; inptemp = OUTtemp
      elecut = OUTelc  ; lwljcut = OUTlwl ; upljcut = OUTupl
      cmbrule = OUTcmb ; cltype = OUTclt  ; screen = OUTscr  ; splodr = OUTspo
      ew1max = OUTew1  ; ew2max = OUTew2  ; ew3max = OUTew3
      ms1max = OUTms1  ; ms2max = OUTms2  ; ms3max = OUTms3
    case(1)  ! default settings for trajectory reading
      estype = ES_NPT                            ! constant pressure
      boxshp = SYS_PERIODIC                      ! periodic boundary
      cltype = EL_PME  ; ms1max = 64             ! PME
      inptemp = 300.0                            ! Kelvin
      engdiv = 1                                 ! number of divisions
      screen = 0.0     ; ewtoler = 0.0  
      ewtoler = 1.0e-6 ; elecut = 12.0
      splodr = 4       ; scrtype = 'dis'
      upljcut = elecut ; lwljcut = upljcut - 2.0
    case default
       stop "Unknown intprm"
    end select
    
    sltpick = 0 ; refpick = 0 ; inscnd = 0 ; inscfg = 0    ! deprecated
    insorigin = INSORG_ORIGIN
    insposition = INSPOS_RANDOM
    insorient = INSROT_RANDOM
    insstructure = INSSTR_NOREJECT

    ! block-wise calculation, corresponds to 13 atoms / box
    block_threshold = 4.0

    force_calculation = .false.

    ! only part of constants set here
    call init_params()

    ! default settings
    skpcnf = 1                     ! no skip for trajectory reading
    
    selfcal = 0                    ! no construction of self-energy distribution
    wgtslf = 0 ; wgtins = 0 ; wgtsys = 0
    select case(slttype)
    case(CAL_SOLN)
       corrcal = 0                 ! no calculation of correlation matrix
    case(CAL_REFS_RIGID, CAL_REFS_FLEX)
       corrcal = 1                 ! calculation of correlation matrix
       if(cltype /= 0) wgtslf = 1  ! Ewald and PME
    case default
       stop "Unknown slttype"
    end select
    sltspec = 1 ; hostspec = 1 ; refspec = 0 ; ljformat = 1 ; ljswitch = 0
    maxins = 1000 ; lwreg = 0.0 ; upreg = 5.0 ; lwstr = 0.0 ; upstr = 2.0
    if(intprm /= 0) then
      cmbrule = 0                  ! arithmetic mean of LJ sigma
      ew2max = ew1max ; ew3max = ew1max
      ms2max = ms1max ; ms3max = ms1max
    endif

    if(sltpick > 0) sltspec = sltpick                      ! deprecated
    if(refpick > 0) refspec = refpick                      ! deprecated

    select case(inscnd)
    case(0)    ! random
       insorigin = INSORG_ORIGIN ; insposition = INSPOS_RANDOM
    case(1)    ! spherical
       insorigin = INSORG_AGGCEN ; insposition = INSPOS_SPHERE
    case(2)    ! slab (symmetric bilayer)
       insorigin = INSORG_AGGCEN ; insposition = INSPOS_SLAB_SYMMETRIC
    case(3)    ! reference
       insorigin = INSORG_REFSTR ; insposition = INSPOS_RMSD
    case default
       stop "Unknown inscnd"
    end select

    select case(inscfg)
    case(0)    ! only the intramolecular configuration is from the file
       insorient = INSROT_RANDOM
    case(1)    ! orientation is fixed from the file with random position
       insorient = INSROT_NOCHANGE
    case(2)    ! position and orientation are both fixed from the file
       insorient = INSROT_NOCHANGE ; insposition = INSPOS_NOCHANGE
    case default
       stop "Unknown inscfg"
    end select

    ! default settings done

    ! read again for non-default constants
    call init_params()

    ! initialize random seed
    if(iseed == 0) then
       CALL RANDOM_SEED
       CALL RANDOM_NUMBER(real_seed)
       iseed = 100000 + int(899999 * real_seed)
    endif
    
    ! temperature converted into the unit of kcal/mol
    temp=inptemp*8.314510e-3/4.184e0
    ! get the screening parameter in Ewald and PME
    if((screen <= tiny).and.(cltype /= 0)) then
       if(ewtoler <= tiny) call halt_with_error('set_ewa')
       screen = getscrn(ewtoler, elecut, scrtype)
    endif
    ! check Ewald parameters, not effective in the current version
    if(cltype == 1) then
       if(ew1max*ew2max*ew3max == 0) call halt_with_error('set_ewa')
    endif
    ! check PME parameters
    if(cltype == EL_PME) then
       if(ms1max*ms2max*ms3max == 0) call halt_with_error('set_ewa')
    endif
    ! ljswitch parameter
    if((ljswitch < 0).or.(ljswitch > 2)) call halt_with_error('set_prs')
    ! check slttype parameter
    if((slttype < CAL_SOLN).or.(slttype > CAL_REFS_FLEX)) then
       call halt_with_error('set_slt')
    endif
    ! check the consistency in parameters for non-periodic system
    if(boxshp == SYS_NONPERIODIC) then
       if((estype == 2).or.(cltype /= 0)) call halt_with_error('set_prs')
    endif
    ! check the consistency of insertion with structure-dependent weight
    if(wgtins == 1) then
       if(slttype /= 3) call halt_with_error('set_ins')
    endif

    plmode=2                    ! energy calculation parallelization mode

    return
  end subroutine iniparam

  subroutine check_insparam
    use engmain, only: numtype, maxins, slttype, hostspec, refspec, &
         insorigin, insposition, insorient, &
         CAL_SOLN, CAL_REFS_RIGID, CAL_REFS_FLEX, &
         INSORG_ORIGIN, INSORG_NOCHANGE, INSORG_AGGCEN, INSORG_REFSTR, &
         INSPOS_RANDOM, INSPOS_NOCHANGE, &
         INSPOS_SPHERE, INSPOS_SLAB_GENERIC, INSPOS_SLAB_SYMMETRIC, &
         INSPOS_RMSD, INSPOS_GAUSS, &
         INSROT_RANDOM, INSROT_NOCHANGE
    use mpiproc, only: halt_with_error, warning
    implicit none
    integer :: valmin, valmax
    logical :: check_ins

    ! insorigin and insorient is effective only for insertion
    if(slttype == CAL_SOLN) then
       insorigin = INSORG_ORIGIN ; insorient = INSROT_RANDOM
    endif

    check_ins = .true.
    ! when restrained relative to aggregate
    if(insorigin == INSORG_AGGCEN) then
       if(slttype == CAL_SOLN) then
          if((hostspec < 1) .or. (hostspec > numtype)) check_ins = .false.
       endif
       if((slttype == CAL_REFS_RIGID) .or. (slttype == CAL_REFS_FLEX)) then
          if((hostspec < 1) .or. (hostspec > numtype - 1)) check_ins = .false.
       endif
    endif
    ! when restrained against reference
    if(insorigin == INSORG_REFSTR) then
       if(slttype == CAL_SOLN) then
          if((refspec < 1) .or. (refspec > numtype)) check_ins = .false.
       endif
       if((slttype == CAL_REFS_RIGID) .or. (slttype == CAL_REFS_FLEX)) then
          if((refspec < 1) .or. (refspec > numtype - 1)) check_ins = .false.
       endif
    endif
    ! when fit to reference
    if((insposition == INSPOS_RMSD) .or. (insposition == INSPOS_GAUSS)) then
       if(insorigin /= INSORG_REFSTR) check_ins = .false.
    else  ! insposition /= INSPOS_RMSD .and. insposition /= INSPOS_GAUSS
       if(insorigin == INSORG_REFSTR) check_ins = .false.
    endif

    ! check the consistency for insorigin, insposition, and insorient
    valmin = min(INSORG_ORIGIN, INSORG_NOCHANGE, INSORG_AGGCEN, INSORG_REFSTR)
    valmax = max(INSORG_ORIGIN, INSORG_NOCHANGE, INSORG_AGGCEN, INSORG_REFSTR)
    if((insorigin < valmin) .or. (insorigin > valmax)) check_ins = .false.

    valmin = min(INSPOS_RANDOM, INSPOS_NOCHANGE, &
                 INSPOS_SPHERE, INSPOS_SLAB_GENERIC, INSPOS_SLAB_SYMMETRIC, &
                 INSPOS_RMSD, INSPOS_GAUSS)
    valmax = max(INSPOS_RANDOM, INSPOS_NOCHANGE, &
                 INSPOS_SPHERE, INSPOS_SLAB_GENERIC, INSPOS_SLAB_SYMMETRIC, &
                 INSPOS_RMSD, INSPOS_GAUSS)
    if((insposition < valmin) .or. (insposition > valmax)) check_ins = .false.

    valmin = min(INSROT_RANDOM, INSROT_NOCHANGE)
    valmax = max(INSROT_RANDOM, INSROT_NOCHANGE)
    if((insorient < valmin) .or. (insorient > valmax)) check_ins = .false.

    if((insorigin == INSORG_NOCHANGE .and. insposition /= INSPOS_NOCHANGE).or.&
       (insorigin /= INSORG_NOCHANGE .and. insposition == INSPOS_NOCHANGE)) then
       check_ins = .false.
    endif
    if(.not. check_ins) call halt_with_error('set_ins')

    ! maxins > 1 in insertion makes no sense if coordinate is used as is read
    if((slttype == CAL_REFS_RIGID .or. slttype == CAL_REFS_FLEX) .and. &
       insposition == INSPOS_NOCHANGE .and. insorient == INSROT_NOCHANGE .and. &
       maxins /= 1) then
       call warning('insu')
    endif

    return
  end subroutine check_insparam


  real function getscrn(ewtoler, elecut, scrtype)
    implicit none
    character(len=3), intent(in) :: scrtype
    real, intent(in) :: ewtoler, elecut
    real :: ewasml, ewalrg, scrfac, factor
    real, parameter :: error=1.0e-20
    factor = error + 1.0 ; ewasml = 0.0 ; ewalrg = 1.0e3
    do while(factor > error)
       scrfac = (ewasml + ewalrg) / 2.0
       factor = erfc(scrfac * elecut)
       if(scrtype == 'dis') factor = factor / elecut
       if(factor > ewtoler) then
          ewasml=scrfac
       else
          ewalrg=scrfac
       endif
       factor = abs(factor - ewtoler)
    end do
    getscrn = scrfac
    return
  end function getscrn

  ! Calls OUTinitial / iniparam / OUTrename, and sets parameters
  subroutine setparam
    use engmain, only: numtype, nummol, numatm, maxcnf, &
         slttype, sltspec, ljformat, &
         moltype, numsite, sluvid, &
         bfcoord, sitemass, charge, &
         ljene_mat, ljlensq_mat, ljtype, ljtype_max, cmbrule, &
         specatm, sitepos, mol_begin_index, belong_to, mol_charge, &
         CAL_SOLN, CAL_REFS_RIGID, CAL_REFS_FLEX, &
         PT_SOLVENT, PT_SOLUTE, PT_TEST_RIGID, PT_TEST_FLEX
    use OUTname, only: OUTinitial, OUTrename, &    ! from outside
         OUTntype, OUTnmol, OUTsite, OUTnrun, &    ! from outside
         OUTstmass, OUTcharge, OUTljene, OUTljlen  ! from outside
    use mpiproc, only: halt_with_error
    use utility, only: itoa
    implicit none
    ! only integer power is allowed as the initialization expression (7.1.6.1)
    real, parameter :: sgmcnv = 1.7817974362806784e0 ! from Rmin/2 to sigma, 2.0**(5.0/6.0)
    real, parameter :: lencnv = 1.0e1                ! from nm to Angstrom
    real, parameter :: engcnv = 1.0e0/4.184e0        ! from kJ/mol to kcal/mol
    integer :: pti, stmax, maxsite, uvtype, cmin, cmax, sid, i, ati, m
    integer :: solute_index, cur_solvent, prev_solvent_type, cur_atom
    real :: factor, xst(3)
    real, dimension(:), allocatable :: sitemass_temp, charge_temp
    integer, allocatable :: ljtype_temp(:)
    real, dimension(:), allocatable :: ljlen_temp, ljene_temp
    real, dimension(:), allocatable :: ljlen_temp_table, ljene_temp_table
    integer :: ljtype_found
    logical :: lj_is_new
    integer, dimension(:), allocatable :: pttype, ptcnt, ptsite
    real, dimension(:,:), allocatable :: psite
    character(len=5) :: atmtype
    character(len=80) :: molfile
    character(len=*), parameter :: sltfile = 'SltInfo'
    character(len=*), parameter :: prmfile = 'MolPrm'
    character(len=*), parameter :: ljtablefile = 'LJTable'
    integer, parameter :: molio = 71                 ! IO for molfile
    integer, parameter :: ljtableio = 70             ! IO for LJ table

    call OUTinitial                ! initialization of OUTname module
    call iniparam                  ! initialization of parameters
    call OUTrename                 ! matching with outside variables

    maxcnf=OUTnrun                                   ! from outside
    if(slttype == CAL_SOLN) then
       numtype=OUTntype                              ! from outside
    else
       numtype=OUTntype+1                            ! from outside
    end if

    ! pttype is particle type for each molecule group
    allocate(pttype(numtype), ptcnt(numtype), ptsite(numtype))
    ! Default is physical particle, coordinate existing in (HISTORY) trajectory
    pttype(:) = PT_SOLVENT
    ! Default: same as written in MDinfo
    ptcnt(1:OUTntype) = OUTnmol(1:OUTntype)
    ptsite(1:OUTntype) = OUTsite(1:OUTntype)

    if(slttype == CAL_SOLN) then
       ! Determine which molecule is solute?
       solute_index = 1            ! default for soln
       if((1 <= sltspec).and.(sltspec <= numtype)) solute_index = sltspec
       pttype(solute_index) = PT_SOLUTE
    else
       ! Test particle information will be defined outside
       ptcnt(numtype) = 1          ! Test particle can only be one molecule

       open(unit = molio, file = sltfile, status='old')
       ! here only counts no. of lines
       stmax = 0
       do
          read(molio, *, end = 99) m
          stmax = stmax + 1
       end do
99     close(molio)
       ptsite(numtype) = stmax
       pttype(numtype) = slttype   ! Test particle is the last (for insertion)       
       solute_index = numtype
    endif

    ! count up number of mols, 
    ! set max and total no. of atoms 
    nummol = sum(ptcnt(1:numtype))
    numatm = sum(ptcnt(1:numtype) * ptsite(1:numtype))

    allocate( moltype(nummol), numsite(nummol), sluvid(nummol) )

    ! make mapping from molecule no. [1..nummol] to particle type [1..numtype]
    cmin=1
    cmax=0
    do pti = 1, numtype
       cmax=cmax+ptcnt(pti)
       moltype(cmin:cmax) = pti ! sequential identification
       cmin=cmax+1
    end do
    if(cmax /= nummol) call halt_with_error("set_num")

    ! Read solute specification
    do i = 1, nummol
       pti = moltype(i)
       numsite(i) = ptsite(pti)
       sluvid(i) = pttype(pti)
    end do

    ! check if all the solute molecules have the same number of atoms
    stmax = -1
    do i = 1, nummol
       if(moltype(i) == solute_index) then        ! solute
          if(stmax == -1) then                 ! initialize
             stmax = numsite(i)
          else
             if(stmax /= numsite(i)) call halt_with_error("set_slt")
          endif
       endif
    end do

    if(numatm /= sum(numsite(1:nummol))) stop "something is wrong in setconf::setparam, numatm"


    ! check the insertion parameters
    call check_insparam

    allocate(bfcoord(3,stmax), sitemass(numatm))
    allocate(charge(numatm), ljtype(numatm))
    allocate(sitepos(3, numatm))
    allocate(mol_begin_index(nummol + 1))
    allocate(mol_charge(nummol))
    allocate(belong_to(numatm))

    ! initial setting to zero
    bfcoord(:,:) = 0.0e0
    sitemass(:) = 0.0e0
    charge(:) = 0.0e0

    ! initialize mol_begin_index
    ! mol_begin_index(i) .. (mol_begin_index(i+1) - 1) will be the index range for i-th molecule
    mol_begin_index(1) = 1
    do i = 1, nummol
       mol_begin_index(i + 1) = mol_begin_index(i) + numsite(i)
    end do
    if(mol_begin_index(nummol + 1) /= numatm + 1) call halt_with_error("set_bug")

    ! initialize belong_to(map from atom number to molecule no)
    do i = 1, nummol
       belong_to(mol_begin_index(i):(mol_begin_index(i + 1) - 1)) = i
    end do

    ! large enough LJ table size
    allocate( ljlen_temp_table(1:sum(ptsite(:))), &
              ljene_temp_table(1:sum(ptsite(:))) )

    ! temporary set of LJ & coordinates
    maxsite = maxval(ptsite(1:numtype))
    allocate( psite(3,maxsite), sitemass_temp(maxsite), charge_temp(maxsite), &
              ljtype_temp(maxsite), ljlen_temp(maxsite), ljene_temp(maxsite) )

    cur_solvent = 0
    cur_atom = 1
    ljtype_max = 0

    ! read molecules specification
    do pti = 1, numtype
       uvtype = pttype(pti)
       if(uvtype == PT_SOLVENT) then            ! solvent
          cur_solvent = cur_solvent + 1
          molfile = prmfile//trim(adjustl(itoa(cur_solvent)))
       else
          if(ptcnt(pti) > 1) cur_solvent = cur_solvent + 1
          molfile = sltfile            ! solute / test particle
       endif
       stmax = ptsite(pti)

       ! This part is a bit complicated due to backward compatibility.
       ! for ljtype /= 5, read the table and make table by program
       open(unit = molio, file = molfile, status='old')
       do sid = 1, stmax
          if(uvtype == CAL_REFS_RIGID) then
             read(molio,*) m, atmtype, xst(1:3), psite(1:3,sid)
          else
             read(molio,*) m, atmtype, xst(1:3)
          endif
          call getmass(sitemass_temp(sid), atmtype)

          charge_temp(sid) = xst(1)
          if(ljformat == 1) xst(3) = sgmcnv * xst(3)
          if((ljformat == 3).or.(ljformat == 4)) then
             if(xst(3) /= 0.0) then
                factor = (xst(2)/xst(3))**(1.0/6.0)
                xst(2) = xst(3)/(4.0*(factor**6.0))
                xst(3) = factor
             else
                xst(2) = 0.0
             endif
          endif
          if((ljformat == 2).or.(ljformat == 4)) then
             xst(2) = engcnv * xst(2)
             xst(3) = lencnv * xst(3)
          endif
          ljene_temp(sid) = xst(2)
          ljlen_temp(sid) = xst(3)
       end do
       close(molio)

       if(ljformat == 5) then
          ! use numbers directly
          ! No sane system will have the problem with 
          ! string -> double -> int conversion ...
          ljtype_temp(1:stmax) = ljene_temp(1:stmax)
       else
          ! allocate LJ types
          do sid = 1, stmax
             lj_is_new = .true.
             do i = 1, ljtype_max
                ! linear search LJ table
                if((ljlen_temp_table(i) == ljlen_temp(sid)) .and. &
                   (ljene_temp_table(i) == ljene_temp(sid))) then
                   ljtype_found = i
                   lj_is_new = .false.
                   exit
                endif
             end do
             if(lj_is_new) then
                ! new LJ type
                ljtype_max = ljtype_max + 1
                ljlen_temp_table(ljtype_max) = ljlen_temp(sid)
                ljene_temp_table(ljtype_max) = ljene_temp(sid)
                ljtype_found = ljtype_max
             endif
             ljtype_temp(sid) = ljtype_found
          end do
       endif

       do i = 1, ptcnt(pti)
          ljtype(cur_atom:(cur_atom + stmax - 1)) = ljtype_temp(1:stmax)
          charge(cur_atom:(cur_atom + stmax - 1)) = charge_temp(1:stmax)
          sitemass(cur_atom:(cur_atom + stmax - 1)) = sitemass_temp(1:stmax)
          cur_atom = cur_atom + stmax
       end do

       if(uvtype == CAL_REFS_RIGID) bfcoord(1:3, 1:stmax) = psite(1:3, 1:stmax)
    end do

    deallocate( psite, sitemass_temp, charge_temp, &
                ljlen_temp, ljene_temp, ljtype_temp )

    ! Fill LJ table
    if(ljformat == 5) then
       ! From table (directly)
       open(unit = ljtableio, file = ljtablefile, status = 'old', action = 'read')
       read(ljtableio, *) ljtype_max
       allocate( ljlensq_mat(ljtype_max, ljtype_max), &
                 ljene_mat(ljtype_max, ljtype_max) )
       do i = 1, ljtype_max
          read (ljtableio, *) ljlensq_mat(i, 1:ljtype_max)
          ljlensq_mat(i, 1:ljtype_max) = ljlensq_mat(i, 1:ljtype_max) ** 2
       end do
       do i = 1, ljtype_max
          read (ljtableio, *) ljene_mat(i, 1:ljtype_max)
       end do
       close(ljtableio)
    else
       ! From LJ data
       allocate( ljlensq_mat(ljtype_max, ljtype_max), &
                 ljene_mat(ljtype_max, ljtype_max) )
       do i = 1, ljtype_max
          select case(cmbrule)
          case(0) ! arithmetic mean
             ljlensq_mat(1:ljtype_max, i) = (( ljlen_temp_table(1:ljtype_max) &
                                             + ljlen_temp_table(i) ) / 2.0) ** 2
          case(1) ! geometric mean
             ljlensq_mat(1:ljtype_max, i) = ljlen_temp_table(1:ljtype_max) &
                                          * ljlen_temp_table(i)
          case default
             stop "Incorrect cmbrule"
          end select
          ljene_mat(1:ljtype_max, i) = sqrt( ljene_temp_table(1:ljtype_max)  &
                                           * ljene_temp_table(i) )
       end do
    endif
    deallocate(ljlen_temp_table, ljene_temp_table)


    ! conversion to (kcal/mol angstrom)^(1/2)
    ! == sqrt(e^2 * coulomb const * avogadro / (kcal / mol angstrom))
    charge(1:numatm) = 18.22261721e0 * charge(1:numatm)

    ! get molecule-wise charges
    do i = 1, nummol
       mol_charge(i) = sum(charge(mol_begin_index(i):(mol_begin_index(i+1)-1)))
    end do

    deallocate( pttype, ptcnt, ptsite )
  end subroutine setparam


! Reads up to maxread coordinates serially and distributes frames
! returns number of frames read (EXCLUDING skipped frames)
  subroutine getconf_parallel(maxread, actual_read)
    use engmain, only: skpcnf, boxshp, numsite, sluvid, stdout, &
                       sitepos, cell, stat_weight_system, PT_SOLVENT, PT_SOLUTE
    use OUTname, only: OUTconfig                     ! from outside
    use mpiproc
    implicit none
    integer, intent(in) :: maxread
    integer, intent(out) :: actual_read
    
    real, dimension(:,:), allocatable :: OUTpos, OUTcell, readpos
    real :: readcell(3, 3)
    real :: weight, readweight
    integer i, OUTatm, iproc, nread

    logical, save :: first_time = .true.
    integer, allocatable, save :: permutation(:)
    integer :: stat

    ! sum over solvent & solute in trajectory file (HISTORY); no test particle
    OUTatm = sum( numsite, &
                  mask = ((sluvid == PT_SOLVENT) .or. (sluvid == PT_SOLUTE)) )
    if(myrank == 0) allocate(readpos(3, OUTatm))
    allocate( OUTpos(3,OUTatm), OUTcell(3,3) )

    ! first time setup: read index permutation
    if(first_time) then
       first_time = .false.
       open(file = perm_file, unit = perm_io, status = 'old', action = 'read', iostat = stat)
       if(stat == 0) then ! file exists and is successfully opened
          if(myrank == 0) write(stdout, *) "Reading permutation information"
          allocate(permutation(OUTatm))
          do i = 1, OUTatm
             permutation(i) = i
          end do
          do
             ! the i-th atom in the trajectory (HISTORY) file is set to
             ! the permutation(i)-th atom in the ermod program
             ! All the other variables such as mass and interaction parameters
             ! and the input files have the order of particles AFTER permutation
             read(perm_io, *, end = 99) i, permutation(i)
             if (i < 1 .or. i > OUTatm) then
                print *, "Incorrect range of particle number at ",i
                call halt_with_error('set_pmt')
             endif
             if (permutation(i) < 1 .or. permutation(i) > OUTatm) then
                print *, i, "-th permutation points to ", permutation(i)
                call halt_with_error('set_pmt')
             endif
          end do
99        close(perm_io)
          ! each particle appears once and only once in permutation
          do i = 1, OUTatm
             if( count(permutation == i) /= 1 ) call halt_with_error('set_pmt')
          end do
       endif
    end if
    
    nread = min(nprocs, maxread)

    if(myrank < nread) then
       if(myrank == 0) then              ! rank-0 to read from file
          do iproc = 1, nread
             ! get configuration and store in OUTpos / OUTcell
             do i = 1, skpcnf
                call OUTconfig(readpos,readcell,OUTatm,boxshp,'system','trjfl_read')
                call read_weight(readweight)
             end do
             
             if(iproc /= 1) then         ! send the data to other rank
#ifndef noMPI
                ! FIXME: rewrite with grouping and scatter?
                call mpi_send(readpos, 3 * OUTatm, mpi_double_precision, &
                              iproc - 1, tag_coord, mpi_comm_world, ierror)
                call mpi_send(readcell, 3 * 3, mpi_double_precision, &
                              iproc - 1, tag_cell, mpi_comm_world, ierror)
                call mpi_send(readweight, 1, mpi_double_precision, &
                              iproc - 1, tag_weight, mpi_comm_world, ierror)
#endif
             else                        ! rank-0 to use the data as read
                ! if this causes memory bottleneck, rewrite with in-place permutation
                OUTpos(:, :) = readpos(:, :)
                OUTcell(:, :) = readcell(:, :)
                weight = readweight
             endif
          end do
       else                              ! non-0 rank to receive the data
#ifndef noMPI
          call mpi_recv(OUTpos, 3 * OUTatm, mpi_double_precision, &
                        0, tag_coord, mpi_comm_world, mpistatus, ierror)
          call mpi_recv(OUTcell, 3 * 3, mpi_double_precision, &
                        0, tag_cell, mpi_comm_world, mpistatus, ierror)
          call mpi_recv(weight, 1, mpi_double_precision, &
                        0, tag_weight, mpi_comm_world, mpistatus, ierror)
#endif
       endif

       if(allocated(permutation)) then
          sitepos(:, permutation(1:OUTatm)) = OUTpos(:, 1:OUTatm)
       else
          sitepos(:, 1:OUTatm) = OUTpos(:, 1:OUTatm)
       endif
       cell(:, :) = OUTcell(:, :)
       stat_weight_system = weight
    endif

    if(myrank == 0) deallocate(readpos)
    deallocate( OUTpos,OUTcell )
    actual_read = nread

  end subroutine getconf_parallel

  subroutine read_weight(weight)
    use engmain, only: wgtsys, stdout
    use mpiproc, only: mpi_setup
    implicit none
    real, intent(out) :: weight
    integer :: dummy, ioerr    
    logical, save :: file_opened = .false.

    if(wgtsys /= 1) then
       weight = 1.0
       return
    endif

    if(.not. file_opened) then
       open(unit = weight_io, file = weight_file, action = 'read')
       file_opened = .true.
    endif
    
    read(weight_io, *, iostat = ioerr) dummy, weight
    if(ioerr /= 0) then 
       rewind(weight_io)
       read(weight_io, *, iostat = ioerr) dummy, weight
       if(ioerr /= 0) then
          write(stdout,*) " The weight file (", weight_file, ") is ill-formed"
          call mpi_setup('stop')
          stop
       endif
    endif
    
  end subroutine read_weight


  subroutine getmass(stmass,atmtype)
    implicit none
    real, parameter :: massH=1.00794e0          ! mass number (hydrogen)
    real, parameter :: massC=12.0107e0          ! mass number (carbon)
    real, parameter :: massO=15.9994e0          ! mass number (oxygen)
    real, parameter :: massN=14.00674e0         ! mass number (nitrogen)
    real, parameter :: massS=32.066e0           ! mass number (sulfur)
    real, parameter :: massP=30.973761e0        ! mass number (phosphorus)
    real, parameter :: massLi=6.941e0           ! mass number (lithium)
    real, parameter :: massNa=22.989770e0       ! mass number (sodium)
    real, parameter :: massK=39.0983e0          ! mass number (potassium)
    real, parameter :: massF=18.9984032e0       ! mass number (fluorine)
    real, parameter :: massCl=35.4527e0         ! mass number (chlorine)
    real, parameter :: massBr=79.904e0          ! mass number (bromine)
    real, parameter :: massCa=40.078e0          ! mass number (calcium)
    real, parameter :: massZn=65.409e0          ! mass number (zinc)
    real, parameter :: massFe=55.845e0          ! mass number (iron)
    real, parameter :: massCu=63.546e0          ! mass number (copper)

    real, intent(out) :: stmass
    character(len=5), intent(in) :: atmtype
    character(len=1) :: eltp1
    character(len=2) :: eltp2
    character(len=3) :: eltp3

    eltp1=atmtype(1:1)
    if(eltp1.eq.'H') stmass=massH
    if(eltp1.eq.'C') stmass=massC
    if(eltp1.eq.'O') stmass=massO
    if(eltp1.eq.'N') stmass=massN
    if(eltp1.eq.'S') stmass=massS
    if(eltp1.eq.'P') stmass=massP
    if(eltp1.eq.'K') stmass=massK
    if(eltp1.eq.'F') stmass=massF
    eltp2=atmtype(1:2)
    if(eltp2.eq.'Li') stmass=massLi
    if(eltp2.eq.'Na') stmass=massNa
    if(eltp2.eq.'Cl') stmass=massCl
    if(eltp2.eq.'Br') stmass=massBr
    if(eltp2.eq.'Ca') stmass=massCa
    if(eltp2.eq.'Zn') stmass=massZn
    if(eltp2.eq.'Fe') stmass=massFe
    if(eltp2.eq.'Cu') stmass=massCu
    if(eltp2.eq.'CH') stmass=massC+massH
    eltp3=atmtype(1:3)
    if(eltp3.eq.'CH2') stmass=massC+2.0e0*massH
    if(eltp3.eq.'CH3') stmass=massC+3.0e0*massH
    return
  end subroutine getmass
end module
