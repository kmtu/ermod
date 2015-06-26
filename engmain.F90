! -*- F90 -*-
! ERmod - Eneregy Representation Module
! Copyright (C) 2000-2015 Nobuyuki Matubayasi
! Copyright (C) 2010-2015 Shun Sakuraba
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


! engmain.f90: various parameters about energy calculation
! 
! Angstrom and kcal/mol are taken as the units of length and energy.
!
!  names of parameters
!   numtype : number of molcular types
!   nummol : total number of molecules
!   numatm : total number of atomic sites
!   maxcnf : maximum number of configurations in MD
!   engdiv : number of divisions of the total simulation length
!   skpcnf : interval to skip the configuration examined
!   corrcal : calculation of the correlation matrix
!               0 : no calculation  1 : calculation performed
!   selfcal : construction of the self-energy distribution
!               0 (default) : no construction  1 : constructed
!   slttype : type of solute treatment
!               1 : physical
!               2 : test particle (rigid)  3 : test particle (flexible)
!            The file for the solute configuration
!            is SltInfo when slttype = 2 and is SltConf when slttype = 3
!   wgtslf : weighting by the self-energy  --- 0 : no  1 : yes
!   wgtins : weight of the solute intramolecular configuration
!               0 : no  1 : yes (can be = 1 only when slttype = 3)
!   wgtsys : weight of the solution / solvent configuration
!               0 : no  1 : yes
!   boxshp : shape of the unit cell box
!               0 : non-periodic  1 : periodic and parallelepiped
!   estype : type of system
!               1 : constant volume  2 : constant pressure
!
!   sltspec : specifying the solute species
!               1 <= sltspec <= numtype (default = 1) when slttype = 1
!               sltspec = numtype when slttype >= 2
!             This parameter is effective as an input only in soln calculation.
!   hostspec : solvent spcies to act as a host and bind the guest solute
!              (micelle, membrane or protein)
!               1 <= hostspec <= numtype      when slttype = 1
!               1 <= hostspec <= numtype - 1  when slttype >= 2
!   refspec : specifying the mixed solvent species for superposition reference
!               1 <= refspec <= numtype       when slttype = 1
!               1 <= refspec <= numtype - 1   when slttype >= 2
!             This parameter is effective as an input only when insorigin = 3
!
!   insorigin : translational origin of the solute position
!               0 (default) : mass weighted center is moved to (0, 0, 0)
!               1 : no COM change from the value in the file to be read
!                   (error unless insposition = 1)
!               2 : mass weighted center is moved to aggregate center
!                   (species forming the aggregate is defined by hostspec)
!               3 : fit to reference structure.
!                   reference structure needs be given as RefInfo in PDB format.
!                   RefInfo contains the structure of the host species
!                   (species forming the reference host is defined by refspec)
!                   and the solute structure, in order.
!                   (error unless insposition = 5 or 6)
!   insposition : position for the solute
!               0 (default) : fully random position (within perodic bondary)
!               1 : no position change from the value in the file to be read
!                   (error unless insorigin = 1)
!               2 : spherically random position,
!                   with radius specified from lwreg to upreg.
!               3 : slab random position (generic case)
!                   slab geometry specified as z = com(aggregate) + dz with
!                   lwreg < dz < upreg for rectangular box periodic condition.
!                   Positioning is more complicated in parallelpiped cell.
!                   (see insertion.F90)
!               4 : slab random position (symmetric bilayer)
!                   slab geometry specified as z = com(aggregate) + dz with
!                   -upreg < dz < -lwreg or lwreg < dz < upreg
!                   for rectangular box periodic condition.
!                   Positioning is more complicated in parallelpiped cell.
!                   (see insertion.F90)
!               5 : random position relative to a reference structure
!                   solvent species identified with the refspec parameter
!                   is set to the reference structure accompanying
!                   the reference position of solute insertion
!                   and the solute is placed relative to that reference
!                   with condition of lwreg < RMSD < upreg
!                   (error unless insorigin = 3)
!               6 : (experimental) Gaussian random position.
!                   Position is given by displacing the reference coordinate,
!                   or coordinate fit to reference (insorigin = 3), with upreg.
!                   Solute weight is automatically adjusted
!                   (error unless insorigin = 3)
!   insorient : orientation for the solute
!               0 (default) : random orientation
!               1 : no orientation change from the value in the file to be read
!   insstructure : intramolecular structure of the solute
!               0 (default) : no restriction, used as is from trajectory or file
!               1 : only the structures with lwstr < RMSD < upstr is counted
!                     RefInfo needs to be prepared to determine RMSD
!
!   inscnd : (deprecated) geometrical condition of the solute configuration
!               0 : random            (insorigin = 0, insposition = 0)  default
!               1 : spherical         (insorigin = 2, insposition = 2)
!               2 : symmetric bilayer (insorigin = 2, insposition = 4)
!               3 : reference         (insorigin = 3, insposition = 5)
!   inscfg : (deprecated) position and orientation for the inserted solute
!               0 : only the intramolecular configuration is from the file.
!                   (insorient = 0)  default
!               1 : orientation is fixed from the file with random position
!                   (insorient = 1)
!               2 : position and orientation are also fixed from the file
!                   (insorient = 1, insposition = 1)
!
!   lwreg : lower bound of the region of solute position
!   upreg : upper bound of the region of solute position
!            effective only when insorigin = 2 or 3 and insposition >= 2
!   lwstr : lower bound of order parameter of solute intramolecular structure
!   upstr : upper bound of order parameter of solute intramolecular structure
!            effective only when insstructure = 1
!
!   ljformat : input-file format for the LJ energy and length parameters
!               0 : epsilon (kcal/mol) and sigma (A)
!               1 (default) : epsilon (kcal/mol) and Rmin/2 (A)
!               2 : epsilon (kJ/mol) and sigma (nm)
!               3 : A (kcal/mol A^12) and C (kcal/mol A^6)
!               4 : C12 (kJ/mol nm^12) and C6 (kJ/mol nm^6)
!               5 : Read from table, LJTable file
!                   (epsilon in kcal/mol and sigma in A)
!   ljswitch : switching function for smooth LJ truncation
!               0 (default) : energy switch in CHARMM form
!               1 : energy switch in GROMACS form
!               2 : force switch
!               tapering function is defined by lwljcut and upljcut variables
!   iseed : seed parameter for uniform random number
!   inptemp : temperature of the system in Kelvin
!   temp  : temperature of the system in kcal/mol
!   block_threshold : box size for cell-link list based method in realcal.F90
!   force_calculation: if set to .true.,
!                      the program continues to run even if there is a warning
!   stdout : standard output
!   
!
!  names of constants and variables for trajectory generation
!   moltype : type of the molecule numbering between 1 and numtype
!   numsite : number of sites in a molecule
!   sluvid : solvent or solute
!               0 : solvent  1 : solute
!               2 : test particle (rigid)  3 : test particle (flexible)
!   bfcoord : coordinate of a rigid molecule in body-fixed frame
!             used only for test particle and kept NULL for physical particle
!   sitemass : mass of an interaction site in a molecule
!   charge : partial charges on the sites in a molecule
!   ljene : energy parameter for Lennard-Jones potential in a molecule
!   ljlen : length parameter for Lennard-Jones potential in a molecule
!   intprm : whether the intereaction paramters given below
!                    (from elecut to ms1max,ms2max,ms3max)
!                    and boxshp, estype, and inptemp
!                    are read from the parent MD program
!      default = 0 in the case of on-the-fly calculation
!      default = 1 in the case of trajectory reading
!            on-the-fly calculation is not effective in the current version
!            and some modification of the programs is necessary
!   elecut : cutoff of the real-space part of electrostatic interaction
!   lwljcut : lower limit of the LJ cutoff tapering function
!   upljcut : upper limit of the LJ cutoff tapering function
!   cmbrule : combination rule for LJ interaction
!        0 : arithmetic mean is used for LJ sigma as for AMBER and CHARMM
!        1 : geometric mean is used for LJ sigma as for OPLS
!      default = 0
!      geometric mean is always used for LJ epsilon
!   cltype : treatment of Coulomb interaction   0 : bare  1 : Ewald  2 : PME
!   screen : screening constant in Ewald summation
!   ewtoler : Ewald and PME tolerance to calculate the screen parameter
!      when screen is given, screen has the priority
!   splodr : order of spline function used in PME
!   ew1max,ew2max,ew3max : number of reciprocal vectors along one direction
!   ms1max,ms2max,ms3max : number of meshes in PME along one direction
!   plmode : parallelization mode for calculation of solute-solvent interaction
!        1 : parallel over solvent molecules in each trajectory snapshot
!        2 : each trajectory snapshot is assigned to
!                   each processor and calculated in parallel
!        default = 2
!        plmode = 1 is not used any more
!   specatm : specification of the site, defined as an integer function
!   sitepos : coordiantes of interaction site
!   cell : unit cell vector
!   invcl : inversion of the cell matrix
!   volume : system volume
!
!  names of constants and variables for energy distribution
!   ermax : size of the energy-represented distribution functions
!   numslv : number of solvent species
!   uvmax : number of discretization for each solvent species
!   uvsoft : number of discretization in soft interaction region
!   esmax : number of discretization of the solute self-energy
!   maxins : maximum number of insertions for test solute particle
!             This parameter is effective as an input only in refs calculation.
!   uvspec : assignment to the number representing the solvent species
!   numslt : number of solute molecules
!   sltlist : list of solute molecules
!   engnorm : normalization factor
!   engsmpl : number of samplings
!   voffset : offset value for the self-energy
!
!  constants defining the discretized energy coordinate
!     these appear only in the enginit subroutine of engproc.F90
!     and are used for each of the solute and solvent species
!     ecmns0, ecpls0, and ecfpls are internally set within enginit
!     and cannot be changed in the parameter files
!   peread : determines whether the parameters are read from a separate file
!       0 : parameters are read from parameters_er (default)
!       1 : parameters are read separately from a file EcdInfo
!   pecore : number of discretization in the core interaction region
!   ecdmin : minimum value of the solute-solvent energy
!   ecfmns : smaller side of the finely dicretized solute-solvent energy
!   ecmns0 : smaller side of the very finely dicretized energy near ecdcen
!   ecdcen : central value of the energy coordinate, typically zero
!   ecpls0 : larger side of the very finely dicretized energy near ecdcen
!   ecfpls : larger side of the finely dicretized solute-solvent energy
!   eccore : the solute-solvent energy at which the dicretization is changed
!   ecdmax : maximum value of the solute-solvent energy
!   eclbin : linear mesh for the solute-solvent energy
!   ecfbin : fine linear mesh for the solute-solvent energy
!   ec0bin : very fine linear mesh for the solute-solvent energy near 0
!   finfac : additional "margin" is set in the low-energy domain
!                       by shifting ecdmin and ecfmns by finfac * ecfbin 
!
!
module engmain
!
  implicit none
  ! Note for optimization: any major compilers shall inline expand "parameter"s
  ! mathematical & physical constants
  real, parameter :: PI = 3.1415926535897932
  real, parameter :: cal_per_joule = 4.1840   ! thermochemical cal / J
!
  integer :: numtype, nummol, numatm, maxcnf, engdiv, skpcnf, corrcal, selfcal
  integer :: slttype, wgtslf, wgtins, wgtsys, boxshp, estype
  integer :: sltspec, hostspec, refspec
  integer :: insorigin, insposition, insorient, insstructure
  integer :: sltpick, refpick, inscnd, inscfg           ! deprecated
  real :: lwreg, upreg, lwstr, upstr
  integer :: ljformat, ljswitch, iseed
  real :: inptemp, temp
  real :: block_threshold
  logical :: force_calculation

  ! IO units
  integer, parameter :: stdout = 6                      ! standard output
  integer, parameter :: io_flcuv = 91                   ! IO unit for flcuv

  integer, dimension(:), allocatable :: moltype, numsite, sluvid
  real, dimension(:,:),  allocatable :: bfcoord
  real, dimension(:),    allocatable :: sitemass, charge, ljene, ljlen

  integer                            :: ljtype_max
  integer, dimension(:), allocatable :: ljtype
  real, dimension(:,:),  allocatable :: ljlensq_mat, ljene_mat
  
  real, dimension(:,:),  allocatable :: sitepos
  real, dimension(:),    allocatable :: mol_charge
  integer, dimension(:), allocatable :: mol_begin_index, belong_to
  real, dimension(3,3)               :: cell, invcl
  real, dimension(3)                 :: celllen
  real                               :: volume

  real    :: elecut, lwljcut, upljcut, screen, ewtoler
  integer :: intprm, cmbrule, cltype, splodr, plmode
  integer :: ew1max, ew2max, ew3max, ms1max, ms2max, ms3max
  
  integer :: ermax, numslv, esmax, maxins
  integer, dimension(:), allocatable :: uvmax, uvsoft, uvspec
  real, dimension(:),    allocatable :: uvcrd, edens
  real, dimension(:,:),  allocatable :: ecorr
  real, dimension(:),    allocatable :: escrd, eself
  real, dimension(:,:),  allocatable :: aveuv
  real, dimension(:),    allocatable :: slnuv
  real, dimension(:,:),  allocatable :: avediv
  real                               :: avslf
  real, dimension(:),    allocatable :: minuv, maxuv
  integer                            :: numslt
  integer, dimension(:), allocatable :: sltlist
  real :: stat_weight_system
  real :: engnorm, engsmpl, voffset
  logical :: voffset_initialized = .false.


  ! numeric constants reference
  integer, parameter :: NO = 0, YES = 1
  integer, parameter :: SYS_NONPERIODIC = 0, SYS_PERIODIC = 1
  integer, parameter :: ES_NVT = 1, ES_NPT = 2
  integer, parameter :: LJFMT_EPS_cal_SGM_nm = 0, LJFMT_EPS_Rminh = 1, &
                        LJFMT_EPS_J_SGM_A = 2, LJFMT_A_C = 3, &
                        LJFMT_C12_C6 = 4, LJFMT_TABLE = 5
  integer, parameter :: LJSWT_POT_CHM = 0, LJSWT_POT_GMX = 1, &
                        LJSWT_FRC_CHM = 2, LJSWT_FRC_GMX = 3
  integer, parameter :: LJCMB_ARITH = 0, LJCMB_GEOM = 1
  integer, parameter :: EL_COULOMB = 0, EL_EWALD = 1, EL_PME = 2
  integer, parameter :: SLT_SOLN = 1, SLT_REFS_RIGID = 2, SLT_REFS_FLEX = 3
  integer, parameter :: PT_SOLVENT = 0, &
                        PT_SOLUTE = SLT_SOLN, PT_TEST_RIGID = SLT_REFS_RIGID, &
                                              PT_TEST_FLEX = SLT_REFS_FLEX
  ! PT_SOLUTE to PT_TEST_FLEX should correspond to SLT_SOLN to SLT_REFS_FLEX
  integer, parameter :: INSORG_ORIGIN = 0, INSORG_NOCHANGE= 1, &
                        INSORG_AGGCEN = 2, INSORG_REFSTR = 3
  integer, parameter :: INSPOS_RANDOM = 0, INSPOS_NOCHANGE= 1, &
                        INSPOS_SPHERE = 2, &
                        INSPOS_SLAB_GENERIC = 3, INSPOS_SLAB_SYMMETRIC = 4, &
                        INSPOS_RMSD = 5, INSPOS_GAUSS = 6
  integer, parameter :: INSROT_RANDOM = 0, INSROT_NOCHANGE= 1
  integer, parameter :: INSSTR_NOREJECT = 0, INSSTR_RMSD = 1

  character(len=*), parameter :: ene_confname = "parameters_er"

  namelist /ene_param/ iseed, &
       skpcnf, corrcal, selfcal, &
       slttype, wgtslf, wgtins, wgtsys, boxshp, estype, &
       sltspec, hostspec, refspec, lwreg, upreg, lwstr, upstr, &
       insorigin, insposition, insorient, insstructure, &
       sltpick, refpick, inscnd, inscfg, &                  ! deprecated
       ljformat, ljswitch, &
       inptemp, temp, &
       engdiv, maxins, &
       intprm, elecut, lwljcut, upljcut, &
       cmbrule, cltype, screen, ewtoler, splodr, plmode, &
       ew1max, ew2max, ew3max, ms1max, ms2max, ms3max, &
       block_threshold, force_calculation

contains 
  subroutine init_params()
    implicit none
    integer, parameter :: unit = 191
    integer :: err
    
    err = 0
    open(unit = unit, file = ene_confname, action = "read", iostat = err)
    
    if(err == 0) then
       read(unit, nml = ene_param)
       close(unit)
    else
       stop "parameter file does not exist"
    end if

  end subroutine init_params

  ! returns atom no. for [i]th atom in [mol]th molecule
  integer function specatm(i, mol)
    implicit none
    integer, intent(in) :: i, mol
    specatm = mol_begin_index(mol) + (i - 1)
  end function specatm

  ! helper function that corresponds mol_begin_index
  integer function mol_end_index(mol)
    implicit none
    integer, intent(in) :: mol
    mol_end_index = mol_begin_index(mol + 1) - 1
  end function mol_end_index
end module engmain
