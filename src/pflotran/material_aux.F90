module Material_Aux_class
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private


  PetscInt, parameter, public :: perm_xx_index = 1
  PetscInt, parameter, public :: perm_yy_index = 2
  PetscInt, parameter, public :: perm_zz_index = 3
  PetscInt, parameter, public :: perm_xy_index = 4
  PetscInt, parameter, public :: perm_yz_index = 5
  PetscInt, parameter, public :: perm_xz_index = 6
  
  PetscInt, parameter, public :: POROSITY_CURRENT = 0
  PetscInt, parameter, public :: POROSITY_MINERAL = 1
  
!  PetscInt, public :: soil_thermal_conductivity_index
!  PetscInt, public :: soil_heat_capacity_index
  PetscInt, public :: soil_compressibility_index
  PetscInt, public :: soil_reference_pressure_index
  PetscInt, public :: max_material_index
  
  type, public :: material_auxvar_type
    PetscInt :: id
    PetscReal :: volume
    PetscReal :: porosity_base
    PetscReal :: porosity
    PetscReal :: dporosity_dp
    PetscReal :: tortuosity
    PetscReal :: soil_particle_density
    PetscReal, pointer :: permeability(:)
    PetscReal, pointer :: sat_func_prop(:)
    PetscReal, pointer :: soil_properties(:) ! den, therm. cond., heat cap.
    type(fracture_auxvar_type), pointer :: fracture
    PetscInt :: creep_closure_id

!    procedure(SaturationFunction), nopass, pointer :: SaturationFunction
  contains
    procedure, public :: PermeabilityTensorToScalar => &
                           MaterialDiagPermTensorToScalar
  end type material_auxvar_type
  
  type, public :: fracture_auxvar_type
    PetscBool :: fracture_is_on
    PetscReal :: initial_pressure
    PetscReal :: properties(4)
    PetscReal :: vector(3) ! < 0. 0. 0. >
  end type fracture_auxvar_type
  
  type, public :: material_parameter_type
    PetscReal, pointer :: soil_residual_saturation(:,:)
    PetscReal, pointer :: soil_heat_capacity(:) ! MJ/kg rock-K
    PetscReal, pointer :: soil_thermal_conductivity(:,:) ! W/m-K
  end type material_parameter_type  
  
  type, public :: material_type
    PetscReal :: time_t, time_tpdt  
    PetscInt :: num_aux
    type(material_parameter_type), pointer :: material_parameter
    class(material_auxvar_type), pointer :: auxvars(:)
  end type material_type
  
  ! procedure pointer declarations
  procedure(MaterialCompressSoilDummy), pointer :: &
    MaterialCompressSoilPtr => null()
 
  ! interface blocks
  interface
    subroutine MaterialCompressSoilDummy(auxvar,pressure,compressed_porosity, &
                                         dcompressed_porosity_dp)
    import material_auxvar_type
    implicit none
    class(material_auxvar_type), intent(in) :: auxvar
    PetscReal, intent(in) :: pressure
    PetscReal, intent(out) :: compressed_porosity
    PetscReal, intent(out) :: dcompressed_porosity_dp
    end subroutine MaterialCompressSoilDummy
  end interface 
  
  interface MaterialCompressSoil
    procedure MaterialCompressSoilPtr
  end interface
  
  public :: MaterialCompressSoilDummy, &
            MaterialCompressSoilPtr, &
            MaterialCompressSoil, &
            MaterialCompressSoilBragflo, &
            MaterialCompressSoilLeijnse, &
            MaterialCompressSoilQuadratic
  
  public :: MaterialAuxCreate, &
            MaterialAuxVarInit, &
            MaterialAuxVarCopy, &
            MaterialAuxVarStrip, &
            MaterialAuxVarGetValue, &
            MaterialAuxVarSetValue, &
            MaterialAuxIndexToPropertyName, &
            MaterialAuxDestroy, &
            MaterialAuxVarFractureStrip
  
contains

! ************************************************************************** !

function MaterialAuxCreate()
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Option_module

  implicit none
  
  type(material_type), pointer :: MaterialAuxCreate
  
  type(material_type), pointer :: aux

  allocate(aux)
  nullify(aux%auxvars)
  allocate(aux%material_parameter)
  nullify(aux%material_parameter%soil_residual_saturation)
  nullify(aux%material_parameter%soil_heat_capacity)
  nullify(aux%material_parameter%soil_thermal_conductivity)
  aux%num_aux = 0
  aux%time_t = 0.d0
  aux%time_tpdt = 0.d0

  MaterialAuxCreate => aux
  
end function MaterialAuxCreate

! ************************************************************************** !

subroutine MaterialAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Option_module

  implicit none
  
  class(material_auxvar_type) :: auxvar
  type(option_type) :: option
  
  auxvar%id = UNINITIALIZED_INTEGER
  auxvar%volume = UNINITIALIZED_DOUBLE
  auxvar%porosity = UNINITIALIZED_DOUBLE
  auxvar%dporosity_dp = 0.d0
  auxvar%porosity_base = UNINITIALIZED_DOUBLE
  auxvar%tortuosity = UNINITIALIZED_DOUBLE
  auxvar%soil_particle_density = UNINITIALIZED_DOUBLE
  if (option%iflowmode /= NULL_MODE) then
    allocate(auxvar%permeability(3))
    auxvar%permeability = UNINITIALIZED_DOUBLE
  else
    nullify(auxvar%permeability)
  endif
  nullify(auxvar%sat_func_prop)
  nullify(auxvar%fracture)
  auxvar%creep_closure_id = 1
  
  if (max_material_index > 0) then
    allocate(auxvar%soil_properties(max_material_index))
    ! initialize these to zero for now
    auxvar%soil_properties = UNINITIALIZED_DOUBLE
  else
    nullify(auxvar%soil_properties)
  endif
  
end subroutine MaterialAuxVarInit

! ************************************************************************** !

subroutine MaterialAuxVarCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Option_module

  implicit none
  
  class(material_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option
  
  auxvar2%volume = auxvar%volume
  auxvar2%porosity = auxvar%porosity
  auxvar2%porosity_base = auxvar%porosity_base
  auxvar2%tortuosity = auxvar%tortuosity
  auxvar2%soil_particle_density = auxvar%soil_particle_density
  if (associated(auxvar%permeability)) then
    auxvar2%permeability = auxvar%permeability
  endif
  if (associated(auxvar%sat_func_prop)) then
    auxvar2%sat_func_prop = auxvar%sat_func_prop
  endif
  if (associated(auxvar%soil_properties)) then
    auxvar2%soil_properties = auxvar%soil_properties
  endif
  auxvar2%creep_closure_id = auxvar%creep_closure_id
end subroutine MaterialAuxVarCopy

! ************************************************************************** !

subroutine MaterialDiagPermTensorToScalar(material_auxvar,dist, &
                                      scalar_permeability)
  ! 
  ! Transforms a diagonal permeability tensor to a scalar through a dot 
  ! product.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 

  use Option_module
  use Utility_module, only : Equal

  implicit none
  
  class(material_auxvar_type) :: material_auxvar
  ! -1 = fraction upwind
  ! 0 = magnitude
  ! 1 = unit x-dir
  ! 2 = unit y-dir
  ! 3 = unit z-dir
  PetscReal, intent(in) :: dist(-1:3)
  PetscReal, intent(out) :: scalar_permeability

  PetscReal :: kx, ky, kz

  kx = material_auxvar%permeability(perm_xx_index)
  ky = material_auxvar%permeability(perm_yy_index)
  kz = material_auxvar%permeability(perm_zz_index)
#if 0
  if (Equal(kx,ky) .and. Equal(ky,kz)) then
    scalar_permeability = kx
    return
  endif
#endif
  scalar_permeability = kx*dabs(dist(1))+ky*dabs(dist(2))+kz*dabs(dist(3))

end subroutine MaterialDiagPermTensorToScalar

! ************************************************************************** !

function MaterialAuxVarGetValue(material_auxvar,ivar)
  ! 
  ! Returns the value of an entry in material_auxvar_type based on ivar.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/28/14
  ! 

  use Variables_module
  
  implicit none

  class(material_auxvar_type) :: material_auxvar 
  PetscInt :: ivar

  PetscReal :: MaterialAuxVarGetValue

  MaterialAuxVarGetValue = UNINITIALIZED_DOUBLE
  select case(ivar)
    case(VOLUME)
      MaterialAuxVarGetValue = material_auxvar%volume
    case(POROSITY)
      MaterialAuxVarGetValue = material_auxvar%porosity
    case(MINERAL_POROSITY)
      MaterialAuxVarGetValue = material_auxvar%porosity_base
    case(TORTUOSITY)
      MaterialAuxVarGetValue = material_auxvar%tortuosity
    case(PERMEABILITY_X)
      MaterialAuxVarGetValue = material_auxvar%permeability(perm_xx_index)
    case(PERMEABILITY_Y)
      MaterialAuxVarGetValue = material_auxvar%permeability(perm_yy_index)
    case(PERMEABILITY_Z)
      MaterialAuxVarGetValue = material_auxvar%permeability(perm_zz_index)
    case(PERMEABILITY_XY)
      MaterialAuxVarGetValue = material_auxvar%permeability(perm_xy_index)
    case(PERMEABILITY_YZ)
      MaterialAuxVarGetValue = material_auxvar%permeability(perm_yz_index)
    case(PERMEABILITY_XZ)
      MaterialAuxVarGetValue = material_auxvar%permeability(perm_xz_index)
    case(SOIL_COMPRESSIBILITY)
      MaterialAuxVarGetValue = material_auxvar% &
                                 soil_properties(soil_compressibility_index)
    case(SOIL_REFERENCE_PRESSURE)
      MaterialAuxVarGetValue = material_auxvar% &
                                 soil_properties(soil_reference_pressure_index)
  end select
  
end function MaterialAuxVarGetValue

! ************************************************************************** !

subroutine MaterialAuxVarSetValue(material_auxvar,ivar,value)
  ! 
  ! Sets the value of an entry in material_auxvar_type based on ivar.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/28/14
  ! 

  use Variables_module
  
  implicit none

  class(material_auxvar_type) :: material_auxvar 
  PetscInt :: ivar
  PetscReal :: value

  select case(ivar)
    case(VOLUME)
      material_auxvar%volume = value
    case(POROSITY)
      material_auxvar%porosity = value
    case(MINERAL_POROSITY)
      material_auxvar%porosity_base = value
    case(TORTUOSITY)
      material_auxvar%tortuosity = value
    case(PERMEABILITY_X)
      material_auxvar%permeability(perm_xx_index) = value
    case(PERMEABILITY_Y)
      material_auxvar%permeability(perm_yy_index) = value
    case(PERMEABILITY_Z)
      material_auxvar%permeability(perm_zz_index) = value
    case(PERMEABILITY_XY)
      material_auxvar%permeability(perm_xy_index) = value
    case(PERMEABILITY_YZ)
      material_auxvar%permeability(perm_yz_index) = value
    case(PERMEABILITY_XZ)
      material_auxvar%permeability(perm_xz_index) = value
    case(SOIL_COMPRESSIBILITY)
      material_auxvar%soil_properties(soil_compressibility_index) = value
    case(SOIL_REFERENCE_PRESSURE)
      material_auxvar%soil_properties(soil_reference_pressure_index) = value
  end select
  
end subroutine MaterialAuxVarSetValue

! ************************************************************************** !

subroutine MaterialCompressSoilLeijnse(auxvar,pressure, &
                                       compressed_porosity, &
                                       dcompressed_porosity_dp)
  ! 
  ! Calculates soil matrix compression based on Leijnse, 1992.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/14/14
  ! 

  implicit none

  class(material_auxvar_type), intent(in) :: auxvar
  PetscReal, intent(in) :: pressure
  PetscReal, intent(out) :: compressed_porosity
  PetscReal, intent(out) :: dcompressed_porosity_dp
  
  PetscReal :: compressibility
  PetscReal :: compression
  PetscReal :: tempreal
  
  compressibility = auxvar%soil_properties(soil_compressibility_index)
  compression = &
    exp(-1.d0 * compressibility * &
        (pressure - auxvar%soil_properties(soil_reference_pressure_index)))
  tempreal = (1.d0 - auxvar%porosity_base) * compression
  compressed_porosity = 1.d0 - tempreal
  dcompressed_porosity_dp = tempreal * compressibility
  
end subroutine MaterialCompressSoilLeijnse

! ************************************************************************** !

subroutine MaterialCompressSoilBRAGFLO(auxvar,pressure, &
                                       compressed_porosity, &
                                       dcompressed_porosity_dp)
  ! 
  ! Calculates soil matrix compression based on Eq. 9.6.9 of BRAGFLO
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/14/14
  ! 

  implicit none

  class(material_auxvar_type), intent(in) :: auxvar
  PetscReal, intent(in) :: pressure
  PetscReal, intent(out) :: compressed_porosity
  PetscReal, intent(out) :: dcompressed_porosity_dp
  
  PetscReal :: compressibility
  
  ! convert to pore compressiblity by dividing by base porosity
  compressibility = auxvar%soil_properties(soil_compressibility_index) / &
                    auxvar%porosity_base
  compressed_porosity = auxvar%porosity_base * &
    exp(compressibility * &
        (pressure - auxvar%soil_properties(soil_reference_pressure_index)))
  dcompressed_porosity_dp = compressibility * compressed_porosity
  
end subroutine MaterialCompressSoilBRAGFLO

! ************************************************************************** !

subroutine MaterialCompressSoilQuadratic(auxvar,pressure, &
                                         compressed_porosity, &
                                         dcompressed_porosity_dp)
  ! 
  ! Calculates soil matrix compression based on a quadratic model
  ! This is thedefaul model adopted in ECLIPSE 
  !
  ! Author: Paolo Orsini
  ! Date: 02/27/17
  ! 

  implicit none

  class(material_auxvar_type), intent(in) :: auxvar
  PetscReal, intent(in) :: pressure
  PetscReal, intent(out) :: compressed_porosity
  PetscReal, intent(out) :: dcompressed_porosity_dp
  
  PetscReal :: compressibility
  PetscReal :: compress_factor

  compressibility = auxvar%soil_properties(soil_compressibility_index)

  compress_factor = compressibility * &
          (pressure - auxvar%soil_properties(soil_reference_pressure_index))

  compressed_porosity = auxvar%porosity_base * &
          ( 1.0 + compress_factor + (compress_factor**2)/2.0 )
  
  dcompressed_porosity_dp = auxvar%porosity_base * &
          ( 1.0 + compress_factor) * compressibility  
  
end subroutine MaterialCompressSoilQuadratic

! ************************************************************************** !

function MaterialAuxIndexToPropertyName(i)
  ! 
  ! Returns the name of the soil property associated with an index
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/06/16
  ! 
  
  implicit none

  PetscInt :: i

  character(len=MAXWORDLENGTH) :: MaterialAuxIndexToPropertyName

  if (i == soil_compressibility_index) then
    MaterialAuxIndexToPropertyName = 'soil compressibility'
  else if (i == soil_reference_pressure_index) then
    MaterialAuxIndexToPropertyName = 'soil reference pressure'
!  else if (i == soil_thermal_conductivity_index) then
!    MaterialAuxIndexToPropertyName = 'soil thermal conductivity'
!  else if (i == soil_heat_capacity_index) then
!    MaterialAuxIndexToPropertyName = 'soil heat capacity'
  else
    MaterialAuxIndexToPropertyName = 'unknown property'
  end if

end function MaterialAuxIndexToPropertyName

! ************************************************************************** !

subroutine MaterialAuxVarFractureStrip(fracture)
  ! 
  ! Deallocates a fracture auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/14/17
  ! 
  use Utility_module, only : DeallocateArray
  
  implicit none

  type(fracture_auxvar_type), pointer :: fracture

  if (.not.associated(fracture)) return

  ! properties and vector are now static arrays.
  deallocate(fracture)
  nullify(fracture)
  
end subroutine MaterialAuxVarFractureStrip

! ************************************************************************** !

subroutine MaterialAuxVarStrip(auxvar)
  ! 
  ! Deallocates a material auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/09/14
  ! 
  use Utility_module, only : DeallocateArray
  
  implicit none

  class(material_auxvar_type) :: auxvar
  
  call DeallocateArray(auxvar%permeability)
  call DeallocateArray(auxvar%sat_func_prop)
  call DeallocateArray(auxvar%soil_properties)
  call MaterialAuxVarFractureStrip(auxvar%fracture)
  
end subroutine MaterialAuxVarStrip

! ************************************************************************** !

subroutine MaterialAuxDestroy(aux)
  ! 
  ! Deallocates a material auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/02/11
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  type(material_type), pointer :: aux
  
  PetscInt :: iaux
  
  if (.not.associated(aux)) return
  
  if (associated(aux%auxvars)) then
    do iaux = 1, aux%num_aux
      call MaterialAuxVarStrip(aux%auxvars(iaux))
    enddo  
    deallocate(aux%auxvars)
  endif
  nullify(aux%auxvars)
    
  if (associated(aux%material_parameter)) then
    call DeallocateArray(aux%material_parameter%soil_residual_saturation)
    call DeallocateArray(aux%material_parameter%soil_heat_capacity)
    call DeallocateArray(aux%material_parameter%soil_thermal_conductivity)
  endif
  nullify(aux%material_parameter)
  
  deallocate(aux)
  nullify(aux)

end subroutine MaterialAuxDestroy

end module Material_Aux_class
