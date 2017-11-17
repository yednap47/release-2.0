module Reaction_Sandbox_Dithionite_class

! # Id: reaction_sandbox_dithionite3, Tue 08 Mar 2016 04:43:33 PM MST pandeys #
!       - update s2o4, so3, hcro4, hco3, h, cr(oh)3(s), fe(oh)3(s), siderite
!       - full reactant inhibition, no product inhibition

! 1. Change all references to "Dithionite" as desired to rename the module and
!    and subroutines within the module.

  use Reaction_Sandbox_Base_class

  use Global_Aux_module
  use Reactive_Transport_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

! 2. Add module variables here.  Note that one must use the PETSc data types
!    PetscInt, PetscReal, PetscBool to declare variables of type integer
!    float/real*8, and logical respectively.  E.g.
!
! PetscReal, parameter :: formula_weight_of_water = 18.01534d0

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_dithionite_type
! 3. Add variables/arrays associated with new reaction.
    character(len=MAXWORDLENGTH) :: name_spec_h ! H+
    character(len=MAXWORDLENGTH) :: name_spec_o2 ! O2(aq)
    character(len=MAXWORDLENGTH) :: name_spec_cr6 ! CrO4-
    character(len=MAXWORDLENGTH) :: name_spec_cr3 ! Cr+++
    character(len=MAXWORDLENGTH) :: name_spec_s2o4 ! S2O4--
    character(len=MAXWORDLENGTH) :: name_spec_s2o3  ! S2O3--
    character(len=MAXWORDLENGTH) :: name_spec_so3 ! SO3--
    character(len=MAXWORDLENGTH) :: name_spec_so4 ! SO4--
    character(len=MAXWORDLENGTH) :: name_spec_fe3 ! Fe+++
!    character(len=MAXWORDLENGTH) :: name_spec_fe2 ! Fe++
    character(len=MAXWORDLENGTH) :: name_bound_fe2_fast ! bound_Fe++ FAST
    character(len=MAXWORDLENGTH) :: name_bound_fe2_slow ! bound_Fe++ SLOW
    character(len=MAXWORDLENGTH) :: name_mnrl_fe3 ! Fe(OH)3(s)
    PetscInt :: id_spec_h, id_spec_o2, id_spec_cr6, id_spec_cr3, id_spec_s2o4 
!    PetscInt :: id_spec_s2o3, id_spec_so3, id_spec_so4, id_spec_fe3, id_spec_fe2
    PetscInt :: id_spec_s2o3, id_spec_so3, id_spec_so4, id_spec_fe3
    PetscInt :: id_mnrl_fe3, id_bound_fe2_fast, id_bound_fe2_slow
    PetscReal :: k_s2o4_disp, k_s2o4_o2, k_s2o4_fe3
    PetscReal :: k_fe2_o2_fast, k_fe2_o2_slow, k_fe2_cr6_fast, k_fe2_cr6_slow
    PetscReal :: ssa_feoh3, rock_density, sitefrac, eps
  contains
    procedure, public :: ReadInput => DithioniteRead
    procedure, public :: Setup => DithioniteSetup
    procedure, public :: Evaluate => DithioniteReact
    procedure, public :: UpdateKineticState => DithioniteKineticState
    procedure, public :: Destroy => DithioniteDestroy
  end type reaction_sandbox_dithionite_type

  public :: DithioniteCreate

contains

! ************************************************************************** !

function DithioniteCreate()
  !
  ! Allocates dithionite reaction object.
  !
  ! Author: Sachin Pandey 
  ! Date: 09/08/16 
  !

  implicit none

  class(reaction_sandbox_dithionite_type), pointer :: DithioniteCreate

! 4. Add code to allocate object and initialized all variables to zero and
!    nullify all pointers. E.g.
  allocate(DithioniteCreate)
  ! Names
  DithioniteCreate%name_spec_h = '' ! H+
  DithioniteCreate%name_spec_o2 = '' ! O2(aq)
  DithioniteCreate%name_spec_cr6 = '' ! CrO4-
  DithioniteCreate%name_spec_cr3 = '' ! Cr+++
  DithioniteCreate%name_spec_s2o4 = '' ! S2O4--
  DithioniteCreate%name_spec_s2o3  = '' ! S2O3--
  DithioniteCreate%name_spec_so3 = '' ! SO3--
  DithioniteCreate%name_spec_so4 = '' ! SO4--
  DithioniteCreate%name_spec_fe3 = '' ! Fe+++
!  DithioniteCreate%name_spec_fe2 = '' ! Fe++
  DithioniteCreate%name_bound_fe2_slow = '' ! bound_Fe++ SLOW
  DithioniteCreate%name_bound_fe2_fast = '' ! bound_Fe++ FAST
  DithioniteCreate%name_mnrl_fe3 = '' ! Fe(OH)3(s)
  
  ! IDs
  DithioniteCreate%id_spec_h = 0 ! H+
  DithioniteCreate%id_spec_o2 = 0 ! O2(aq)
  DithioniteCreate%id_spec_cr6 = 0 ! CrO4-
  DithioniteCreate%id_spec_cr3 = 0 ! Cr+++
  DithioniteCreate%id_spec_s2o4 = 0 ! S2O4--
  DithioniteCreate%id_spec_s2o3  = 0 ! S2O3--
  DithioniteCreate%id_spec_so3 = 0 ! SO3--
  DithioniteCreate%id_spec_so4 = 0 ! SO4--
  DithioniteCreate%id_spec_fe3 = 0 ! Fe+++
  !DithioniteCreate%id_spec_fe2 = 0 ! Fe++
  DithioniteCreate%id_bound_fe2_slow = 0 ! bound_Fe++ SLOW
  DithioniteCreate%id_bound_fe2_fast = 0 ! bound_Fe++ FAST
  DithioniteCreate%id_mnrl_fe3 = 0 ! Fe(OH)3(s)

  ! Rate constants
  DithioniteCreate%k_s2o4_disp = 0.d0
  DithioniteCreate%k_s2o4_o2 = 0.d0
  DithioniteCreate%k_s2o4_fe3 = 0.d0
  DithioniteCreate%k_fe2_o2_fast = 0.d0
  DithioniteCreate%k_fe2_o2_slow = 0.d0
  DithioniteCreate%k_fe2_cr6_fast = 0.d0
  DithioniteCreate%k_fe2_cr6_slow = 0.d0

  ! Other constants
  DithioniteCreate%ssa_feoh3 = 0.d0
  DithioniteCreate%rock_density = 0.d0
  DithioniteCreate%sitefrac = 0.d0
  DithioniteCreate%eps = 0.d0

  nullify(DithioniteCreate%next)

end function DithioniteCreate

! ************************************************************************** !

subroutine DithioniteRead(this,input,option)
  !
  ! Reads input deck for dithionite reaction parameters (if any)
  !
  ! Author: Sachin Pandey
  ! Date: 09/08/16
  !

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_dithionite_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word, internal_units

  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,DITHIONITE_PARAMETERS')
    call StringToUpper(word)

    select case(trim(word))
      case('K_S2O4_DISP')
        call InputReadDouble(input,option,this%k_s2o4_disp)
        call InputErrorMsg(input,option,'K_S2O4_DISP', &
                           'CHEMISTRY,REACTION_SANDBOX,DITHIONITE_PARAMETERS')
      case('K_S2O4_O2')
        call InputReadDouble(input,option,this%k_s2o4_o2)
        call InputErrorMsg(input,option,'K_S2O4_O2', &
                           'CHEMISTRY,REACTION_SANDBOX,DITHIONITE_PARAMETERS')
      case('K_S2O4_FE3')
        call InputReadDouble(input,option,this%k_s2o4_fe3)
        call InputErrorMsg(input,option,'K_S2O4_FE3', &
                           'CHEMISTRY,REACTION_SANDBOX,DITHIONITE_PARAMETERS')
      case('K_FE2_O2_FAST')
        call InputReadDouble(input,option,this%k_fe2_o2_fast)
        call InputErrorMsg(input,option,'K_FE2_O2_FAST', &
                           'CHEMISTRY,REACTION_SANDBOX,DITHIONITE_PARAMETERS')
      case('K_FE2_O2_SLOW')
        call InputReadDouble(input,option,this%k_fe2_o2_slow)
        call InputErrorMsg(input,option,'K_FE2_O2_SLOW', &
                           'CHEMISTRY,REACTION_SANDBOX,DITHIONITE_PARAMETERS')
      case('K_FE2_CR6_FAST')
        call InputReadDouble(input,option,this%k_fe2_cr6_fast)
        call InputErrorMsg(input,option,'K_FE2_CR6_FAST', &
                           'CHEMISTRY,REACTION_SANDBOX,DITHIONITE_PARAMETERS')
      case('K_FE2_CR6_SLOW')
        call InputReadDouble(input,option,this%k_fe2_cr6_slow)
        call InputErrorMsg(input,option,'K_FE2_CR6_SLOW', &
                           'CHEMISTRY,REACTION_SANDBOX,DITHIONITE_PARAMETERS')
      case('SSA_FEOH3')
        call InputReadDouble(input,option,this%ssa_feoh3)
        call InputErrorMsg(input,option,'SSA_FEOH3', &
                           'CHEMISTRY,REACTION_SANDBOX,DITHIONITE_PARAMETERS')
      case('ROCK_DENSITY')
        call InputReadDouble(input,option,this%rock_density)
        call InputErrorMsg(input,option,'ROCK_DENSITY', &
                           'CHEMISTRY,REACTION_SANDBOX,DITHIONITE_PARAMETERS')
      case('ALPHA')
        call InputReadDouble(input,option,this%sitefrac)
        call InputErrorMsg(input,option,'ALPHA', &
                           'CHEMISTRY,REACTION_SANDBOX,DITHIONITE_PARAMETERS')
      case('EPS')
        call InputReadDouble(input,option,this%eps)
        call InputErrorMsg(input,option,'EPS', &
                           'CHEMISTRY,REACTION_SANDBOX,DITHIONITE_PARAMETERS')

      case default
        call InputKeywordUnrecognized(word, &
                     'CHEMISTRY,REACTION_SANDBOX,DITHIONITE_PARAMETERS',option)
    end select
  enddo

end subroutine DithioniteRead

! ************************************************************************** !

subroutine DithioniteSetup(this,reaction,option)
  !
  ! Sets up the dithionite reaction either with parameters either
  ! read from the input deck or hardwired.
  !
  ! Author: Sachin Pandey
  ! Date: 09/08/16
  !

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Reaction_Immobile_Aux_module, only : GetImmobileSpeciesIDFromName
  use Reaction_Mineral_Aux_module, only : GetMineralIDFromName
  use Option_module

  implicit none

  class(reaction_sandbox_dithionite_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

! 9. Add code to initialize
  this%name_spec_h = 'H+'
  this%name_spec_o2 = 'O2(aq)'
  this%name_spec_cr6 = 'CrO4--'
  this%name_spec_cr3 = 'Cr+++'
  this%name_spec_s2o4 = 'S2O4--'
  this%name_spec_s2o3  = 'S2O3--'
  this%name_spec_so3 = 'SO3--'
  this%name_spec_so4 = 'SO4--'
  this%name_spec_fe3 = 'Fe+++'
!  this%name_spec_fe2 = 'Fe++'
  this%name_bound_fe2_fast  = 'fast_Fe++'
  this%name_bound_fe2_slow  = 'slow_Fe++'
  this%name_mnrl_fe3  = 'Fe(OH)3(s)'

  ! Aqueous species
  this%id_spec_h = &
      GetPrimarySpeciesIDFromName(this%name_spec_h,reaction,option)
  this%id_spec_o2 = &
      GetPrimarySpeciesIDFromName(this%name_spec_o2,reaction,option)
  this%id_spec_cr6 = &
      GetPrimarySpeciesIDFromName(this%name_spec_cr6,reaction,option)
  this%id_spec_cr3 = &
      GetPrimarySpeciesIDFromName(this%name_spec_cr3,reaction,option)
  this%id_spec_s2o4 = &
    GetPrimarySpeciesIDFromName(this%name_spec_s2o4,reaction,option)
  this%id_spec_s2o3 = &
      GetPrimarySpeciesIDFromName(this%name_spec_s2o3,reaction,option)
  this%id_spec_so3 = &
      GetPrimarySpeciesIDFromName(this%name_spec_so3,reaction,option)
  this%id_spec_so4 = &
      GetPrimarySpeciesIDFromName(this%name_spec_so4,reaction,option)
  this%id_spec_fe3 = &
      GetPrimarySpeciesIDFromName(this%name_spec_fe3,reaction,option)
!  this%id_spec_fe2 = &
!      GetPrimarySpeciesIDFromName(this%name_spec_fe2,reaction,option)

  ! Immobile species
  this%id_bound_fe2_fast = &
    GetImmobileSpeciesIDFromName(this%name_bound_fe2_fast,reaction%immobile,option)
  this%id_bound_fe2_slow = &
    GetImmobileSpeciesIDFromName(this%name_bound_fe2_slow,reaction%immobile,option)

  ! Mineral species
  this%id_mnrl_fe3 = &
    GetMineralIDFromName(this%name_mnrl_fe3,reaction%mineral,option)

end subroutine DithioniteSetup

! ************************************************************************** !

subroutine DithioniteReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)
  !
  ! Evaluates reaction storing residual and/or Jacobian
  !
  ! Author: Sachin Pandey
  ! Date: 09/08/16
  !

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class

  implicit none

  class(reaction_sandbox_dithionite_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: L_water, FeII_fast, FeII_slow, FeIII
  PetscReal :: vf_feoh3, mv_feoh3, mw_feoh3
  PetscReal :: r_s2o4_disp, r_s2o4_o2, r_s2o4_fe3
  PetscReal :: r_fe2_o2_fast, r_fe2_o2_slow, r_fe2_cr6_fast, r_fe2_cr6_slow
  PetscInt :: id_bound_fe2_slow_offset, id_bound_fe2_fast_offset

  ! Info for Fe(oh)3(s) and bound_Fe++
  vf_feoh3 = rt_auxvar%mnrl_volfrac(this%id_mnrl_fe3) ! m^3/m^3_bulk
  mv_feoh3 = reaction%mineral%kinmnrl_molar_vol(this%id_mnrl_fe3) ! m^3/mol
  mw_feoh3 = reaction%mineral%kinmnrl_molar_wt(this%id_mnrl_fe3) ! m^3/mol
  id_bound_fe2_slow_offset = reaction%offset_immobile + this%id_bound_fe2_slow
  id_bound_fe2_fast_offset = reaction%offset_immobile + this%id_bound_fe2_fast

  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3 ! L_water from m^3_water

  ! Surface bound Fe(II): mol_reactant/g_sed from from mol_reactant/m^3_bulk
  FeII_fast = rt_auxvar%immobile(this%id_bound_fe2_fast)/(this%rock_density*1000)
  FeII_slow = rt_auxvar%immobile(this%id_bound_fe2_slow)/(this%rock_density*1000)

  ! Mineral species: g_fe(oh)3/g_sed from from m^3_mnrl/m^3_bulk
  FeIII = vf_feoh3/mv_feoh3/this%rock_density/1.d3*mw_feoh3

  ! mole/s
  r_s2o4_disp = this%k_s2o4_disp * rt_auxvar%pri_molal(this%id_spec_s2o4) * L_water
  r_s2o4_o2 = this%k_s2o4_o2 * rt_auxvar%pri_molal(this%id_spec_s2o4) * rt_auxvar%pri_molal(this%id_spec_o2) * L_water
  r_fe2_o2_fast = this%k_fe2_o2_fast * FeII_fast * rt_auxvar%pri_molal(this%id_spec_o2) * L_water
  r_fe2_o2_slow = this%k_fe2_o2_slow * FeII_slow * rt_auxvar%pri_molal(this%id_spec_o2) * L_water
  r_fe2_cr6_fast = this%k_fe2_cr6_fast * FeII_fast * rt_auxvar%pri_molal(this%id_spec_cr6) * L_water
  r_fe2_cr6_slow = this%k_fe2_cr6_slow * FeII_slow * rt_auxvar%pri_molal(this%id_spec_cr6) * L_water
  r_s2o4_fe3 = this%k_s2o4_fe3 * this%ssa_feoh3 * rt_auxvar%pri_molal(this%id_spec_s2o4) * FeIII * L_water

  if (r_s2o4_disp < this%eps) then
    r_s2o4_disp = 0.d0
  endif

  if (r_s2o4_o2 < this%eps) then
    r_s2o4_o2 = 0.d0
  endif
  
  if (r_fe2_o2_fast < this%eps) then
    r_fe2_o2_fast = 0.d0
  endif
  
  if (r_fe2_o2_slow < this%eps) then
    r_fe2_o2_slow = 0.d0
  endif
  
  if (r_fe2_cr6_fast < this%eps) then
    r_fe2_cr6_fast = 0.d0
  endif
  
  if (r_fe2_cr6_slow < this%eps) then
    r_fe2_cr6_slow = 0.d0
  endif
  
  if (r_s2o4_fe3 < this%eps) then
    r_s2o4_fe3 = 0.d0
  endif

  ! for debug
#if 0
   print *, 'mv_feoh3', mv_feoh3
   print *, 'mw_feoh3', mw_feoh3
   print *, 'rock_density', this%rock_density
   print *, 'k_s2o4_disp', this%k_s2o4_disp
   print *, 'k_s2o4_o2', this%k_s2o4_o2
   print *, 'k_s2o4_fe3', this%k_s2o4_fe3
   print *, 'k_fe2_o2_fast', this%k_fe2_o2_fast
   print *, 'k_fe2_o2_slow', this%k_fe2_o2_slow
   print *, 'k_fe2_cr6_fast', this%k_fe2_cr6_fast
   print *, 'k_fe2_cr6_slow', this%k_fe2_cr6_slow
#endif

  ! NOTES
  ! 1. Always subtract contribution from residual
  ! 2. Units of residual are moles/second  
  Residual(this%id_spec_h) = Residual(this%id_spec_h) - &
    (1.0*r_s2o4_disp +2.0*r_s2o4_o2 -1.0*r_fe2_o2_fast -1.0*r_fe2_o2_slow &
    -2.66*r_fe2_cr6_fast -2.66*r_fe2_cr6_slow -2.0*r_s2o4_fe3)

  Residual(this%id_spec_o2) = Residual(this%id_spec_o2) - &
    (-1.0*r_s2o4_o2 -0.25*r_fe2_o2_fast -0.25*r_fe2_o2_slow)

  Residual(this%id_spec_cr6) = Residual(this%id_spec_cr6) - &
    (-0.33*r_fe2_cr6_fast -0.33*r_fe2_cr6_slow)

  Residual(this%id_spec_cr3) = Residual(this%id_spec_cr3) - &
    (0.33*r_fe2_cr6_fast +0.33*r_fe2_cr6_slow)

  Residual(this%id_spec_s2o4) = Residual(this%id_spec_s2o4) - &
    (-1.0*r_s2o4_disp -1.0*r_s2o4_o2 -1.0*r_s2o4_fe3)

  Residual(this%id_spec_s2o3) = Residual(this%id_spec_s2o3) - &
    (0.5*r_s2o4_disp)

  Residual(this%id_spec_so3) = Residual(this%id_spec_so3) - &
    (1.0*r_s2o4_disp +1.0*r_s2o4_o2 +2.0*r_s2o4_fe3)

  Residual(this%id_spec_so4) = Residual(this%id_spec_so4) - &
    (1.0*r_s2o4_o2)

  Residual(this%id_spec_fe3) = Residual(this%id_spec_fe3) - & 
    (1.0*r_fe2_o2_fast +1.0*r_fe2_o2_slow +1.0*r_fe2_cr6_fast +1.0*r_fe2_cr6_slow)

!  Residual(this%id_spec_fe2) = Residual(this%id_spec_fe2) - 0.0

  Residual(id_bound_fe2_fast_offset) = Residual(id_bound_fe2_fast_offset) - &
    (-1.0*r_fe2_o2_fast -1.0*r_fe2_cr6_fast +2.0*r_s2o4_fe3*this%sitefrac)

  Residual(id_bound_fe2_slow_offset) = Residual(id_bound_fe2_slow_offset) - &
    (-1.0*r_fe2_o2_slow -1.0*r_fe2_cr6_slow +2.0*r_s2o4_fe3*(1-this%sitefrac))

  if (compute_derivative) then

! 11. If using an analytical Jacobian, add code for Jacobian evaluation

     option%io_buffer = 'NUMERICAL_JACOBIAN_RXN must always be used ' // &
                        'due to assumptions in Dithionite'
     call printErrMsg(option)

  endif

end subroutine DithioniteReact

! ************************************************************************** !

subroutine DithioniteKineticState(this,rt_auxvar,global_auxvar, &
                                  material_auxvar,reaction,option)
  !
  ! For update mineral volume fractions
  !
  ! Author: Sachin Pandey
  ! Date: 09/08/16
  !

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class

  implicit none

  class(reaction_sandbox_dithionite_type) :: this
  type(option_type) :: option
  type(reaction_type) :: reaction
  ! the following arrays must be declared after reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: FeII, FeIII, rate, cnv_L2bulk
  PetscReal :: vf_feoh3, mv_feoh3, mw_feoh3
  PetscReal :: delta_volfrac

  ! Unit of the residual must be in moles/second
  ! Unit of mnr_rate should be mol/sec/m^3 bulk
  ! 1.d3 converts m^3 water -> L water
  vf_feoh3 = rt_auxvar%mnrl_volfrac(this%id_mnrl_fe3) ! m^3/m^3_bulk
  mv_feoh3 = reaction%mineral%kinmnrl_molar_vol(this%id_mnrl_fe3) ! m^3/mol
  mw_feoh3 = reaction%mineral%kinmnrl_molar_wt(this%id_mnrl_fe3) ! m^3/mol

  ! L_water/m^3_bulk
  cnv_L2bulk = material_auxvar%porosity*global_auxvar%sat(iphase)*1.d3

  ! UNITS: m^2_reactant/g_sed from from m^3_mnrl/m^3_bulk
  FeIII = vf_feoh3/mv_feoh3/this%rock_density/1.d3*mw_feoh3

  ! mole/l-s
  rate = this%k_s2o4_fe3*this%ssa_feoh3* &
             rt_auxvar%pri_molal(this%id_spec_s2o4)* &
             FeIII

  if (rate < this%eps) then
    rate = 0.d0
  endif

  ! rate = mol/m^3/sec
  ! dvolfrac = m^3 mnrl/m^3 bulk = rate (mol mnrl/m^3 bulk/sec) * mol_vol (m^3 mnrl/mol mnrl)
  ! Update Fe(OH)3
  delta_volfrac = (-2.0*rate)*cnv_L2bulk* &
                  reaction%mineral%kinmnrl_molar_vol(this%id_mnrl_fe3)* &
                  option%tran_dt
  rt_auxvar%mnrl_volfrac(this%id_mnrl_fe3) = &
    rt_auxvar%mnrl_volfrac(this%id_mnrl_fe3) + delta_volfrac

end subroutine DithioniteKineticState

! ************************************************************************** !

subroutine DithioniteDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this
  ! module
  !
  ! Author: Sachin Pandey
  ! Date: 09/08/16
  !

  implicit none

  class(reaction_sandbox_dithionite_type) :: this

! 12. Add code to deallocate contents of the dithionite object

end subroutine DithioniteDestroy

end module Reaction_Sandbox_Dithionite_class
