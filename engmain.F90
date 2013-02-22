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
!   sltpick : specifying the solute species
!               1 <= sltpick <= numtype (default = 1) if slttype = 1
!               sltpick = numtype if slttype >= 2
!   refpick : specifying the mixed solvent species for superposition
!               1 <= refpick <= numtype    when slttype = 1
!               1 <= refpick <= numtype-1  when slttype >= 2
!   wgtslf : weighting by the self-energy  --- 0 : no  1 : yes
!   wgtins : weight of the solute intramolecular configuration
!               0 : no  1 : yes (can be = 1 only when slttype = 3)
!   wgtsys : weight of the solution / solvent configuration
!               0 : no  1 : yes
!   estype : type of system
!               1 : constant volume  2 : constant pressure
!   boxshp : shape of the unit cell box
!               0 : non-periodic  1 : periodic and parallelepiped
!
!   insorigin : the origin of the insertion position
!               0 : (default) mass weighted center is moved to (0, 0, 0)
!               1 : mass weighted center is moved to aggregate center
!                   (species forming the aggregate is defined by hostspec)
!               2 : fit to reference structure.
!                   reference structure needs be given as RefInfo in PDB format.
!                   RefInfo should conatin the structure specified in refpick
!                   and the solute structure, in order.
!               3 : no COM change from the value in the file to be read
!                   (error unless insposition = 3)
!   insposition : position for the inserted solute
!               0 : (default) fully random position (within perodic bondary)
!               1 : spherically random position,
!                   with radius specified from lwreg to upreg.
!               2 : slab random position
!                   slab geometry specified as z = com(aggregate) + dz with
!                   (-upreg < dz < -lwreg AND lwreg < dz < upreg)
!                   for rectangular box periodic condition.
!                   Position is much more complex for parallelpiped structure.
!                   (see insertion.F90)
!               3 : no position change from the value in the file to be read
!                   (error unless insorigin = 3)
!               4 : (experimental) Gaussian random position.
!                   Position is given by displacing the reference coordinate,
!                   or coordinate fit to reference (insorigin = 2), with upreg.
!                   Solute weight is automatically adjusted
!   insorient : orientation for the inserted solute
!               0 : (default) random orientation
!               1 : no orientation change from the value in the file to be read
!
!   inscnd : (deprecated) geometrical condition of the solute configuration
!               0 : random    (insorigin = 0, insposition = 0) 
!               1 : spherical (insorigin = 1, insposition = 1)
!               2 : slab      (insorigin = 1, insposition = 2)
!               3 : reference (insorigin = 2, insposition = 4)
!   inscfg : (deprecated) position and orientation for the inserted solute
!               0 : only the intramolecular configuration is from the file.
!                   (insorient = 0)
!               1 : orientation is fixed from the file with random position
!                   (insorient = 1)
!               2 : position and orientation are also fixed from the file
!                   (insorient = 1, insposition = 3)
!              default = 0
!   hostspec : type of molecule forming host (micelle, membrane, or protein)
!              active only when insorigin = 1
!               1 <= hostspec <= numtype    when slttype = 1
!               1 <= hostspec <= numtype-1  when slttype >= 2
!
!   ljformat : input-file format for the LJ energy and length parameters
!               0 : epsilon (kcal/mol) and sigma (A)
!               1 (default) : epsilon (kcal/mol) and Rmin/2 (A)
!               2 : epsilon (kJ/mol) and sigma (nm)
!               3 : A (kcal/mol A^12) and C (kcal/mol A^6)
!               4 : C12 (kJ/mol nm^12) and C6 (kJ/mol nm^6)
!               5 : Read from table, LJTable file
!                   (epsilon as kcal/mol and sigma as A)
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
!   refmlid : superposition reference among solvent species
!               0 : not reference
!               1 : reference solvent  2 : reference solute
!             default = 0 unless 1 <= refpick <= numtype for solvent species
!   bfcoord : coordinate of a rigid molecule in body-fixed frame
!             used only for test particle and kept NULL for physical particle
!   sitemass : mass of an interaction site in a molecule
!   charge : partial charges on the sites in a molecule
!   ljene : energy parameter for Lennard-Jones potential in a molecule
!   ljlen : length parameter for Lennard-Jones potential in a molecule
!   intprm : whether the intereaction paramters given below
!                    (from elecut to ms1max,ms2max,ms3max)
!                    and estype, boxshp, and inptemp
!                    are read from the parent MD program
!      default = 0 in the case of on-the-fly calculation
!      default = 1 in the case of trajectory reading
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
!        2 : each trajectory snapshot is assigned to each processor and calculated in parallel
!        default = 2
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
!   uvspec : assignment to the number representing the solvent species
!   numslt : number of solute molecules
!   sltlist : list of solute molecules
!   engnorm : normalization factor
!   engsmpl : number of samplings
!   voffset : offset value for the self-energy
!   lwreg : lower bound for the region of solute position
!   upreg : upper bound for the region of solute position
!
!  constants defining the discretized energy coordinate
!     these appear only in the enginit subroutine
!     and are used for each of the solute and solvent species
!   peread : determines whether the parameters are read from a separate file
!       0 : parameters are read from param_eng
!       1 : parameters are read separately from a file EcdInfo
!     default = 0
!   pemax : number of discretization of the solute-solvent energy
!   pesoft : number of discretization in the soft interaction region
!   pecore : number of discretization in the core interaction region
!            pemax and pesoft are constructed from other parameters
!            pemax = pesoft + pecore
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
!
!
module engmain
!
  implicit none
  ! Note for optimization: any major compilers shall inline expand "parameter"s
  ! mathematical & physical constants
  real, parameter :: PI = 3.1415926535897932
  real, parameter :: cal_per_joule = 4.1840e0 ! thermochemical cal / J
!
  integer :: numtype,nummol,numatm,maxcnf,engdiv,skpcnf,corrcal,selfcal
  integer :: slttype, sltpick, refpick, wgtslf, wgtins, wgtsys
  integer :: estype,boxshp
  integer :: insorigin, insposition, insorient, inscnd, inscfg, hostspec
  integer :: ljformat, ljswitch, iseed
  real :: inptemp,temp
  real :: block_threshold
  logical :: force_calculation

  ! IO units
  integer, parameter :: stdout = 6                      ! standard output
  integer, parameter :: io_flcuv = 91                   ! IO unit for flcuv

  integer, dimension(:),   allocatable :: moltype,numsite
  integer, dimension(:),   allocatable :: sluvid,refmlid
  real, dimension(:,:),    allocatable :: bfcoord
  real, dimension(:),      allocatable :: sitemass
  real, dimension(:),      allocatable :: charge,ljene,ljlen

  integer :: ljtype_max
  integer, allocatable :: ljtype(:)
  real, allocatable :: ljlensq_mat(:, :), ljene_mat(:, :)
  
  real, dimension(:,:),    allocatable :: sitepos
  real, allocatable :: mol_charge(:)
  integer, allocatable :: mol_begin_index(:), belong_to(:)
  real, dimension(3,3)                 :: cell,invcl
  real :: celllen(3)
  real :: volume

  real :: elecut,lwljcut,upljcut,screen,ewtoler
  integer :: intprm,cmbrule,cltype,splodr,plmode
  integer :: ew1max,ew2max,ew3max,ms1max,ms2max,ms3max
  
  integer :: ermax,numslv,esmax,maxins
  integer, dimension(:),   allocatable :: uvmax,uvsoft,uvspec
  real, dimension(:),      allocatable :: uvcrd,edens
  real, dimension(:,:),    allocatable :: ecorr
  real, dimension(:),      allocatable :: escrd,eself
  real, dimension(:,:),    allocatable :: aveuv
  real, dimension(:),      allocatable :: slnuv
  real, dimension(:,:),    allocatable :: avediv
  real :: avslf
  real, allocatable :: minuv(:), maxuv(:)
  integer :: numslt
  integer, allocatable :: sltlist(:)
  real :: stat_weight_system
  real :: engnorm,engsmpl,voffset
  logical :: voffset_initialized = .false.
  real :: lwreg,upreg


  ! numeric constants reference
  integer, parameter :: SYS_NONPERIODIC = 0, SYS_PERIODIC = 1
  integer, parameter :: EL_COULOMB = 0, EL_PME = 2
  integer, parameter :: CAL_SOLN = 1, CAL_REFS_RIGID = 2, CAL_REFS_FLEX = 3
  integer, parameter :: ES_NVT = 1, ES_NPT = 2
  integer, parameter :: PT_SOLVENT = 0, PT_SOLUTE = 1, PT_TEST_RIGID = 2, PT_TEST_FLEX = 3
  ! note: PT_SOLUTE to PT_TEST_FLEX should correspond to CAL_SOLN .. CAL_REFS_FLEX

  character(len=*), parameter :: ene_confname = "parameters_er"

  namelist /ene_param/ iseed, &
       skpcnf,corrcal,selfcal, &
       slttype,sltpick,refpick,wgtslf,wgtins,wgtsys, &
       estype,boxshp, &
       insorigin, insposition, insorient, inscnd, inscfg, &
       hostspec, lwreg, upreg, &
       ljformat, ljswitch, &
       inptemp,temp, &
       engdiv,maxins, &
       intprm,elecut,lwljcut,upljcut, &
       cmbrule,cltype,screen,ewtoler,splodr,plmode, &
       ew1max,ew2max,ew3max,ms1max,ms2max,ms3max, &
       block_threshold, force_calculation

contains 
  subroutine init_params()
    implicit none
    integer, parameter :: unit = 191
    integer :: err
    
    err = 0
    open(unit = unit, file = ene_confname, action = "read", iostat = err)
    
    if(err == 0) then
       read(unit, nml=ene_param)
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
