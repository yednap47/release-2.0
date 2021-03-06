module Saturation_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  public :: SaturationUpdateCoupler
 
contains

! ************************************************************************** !

subroutine SaturationUpdateCoupler(coupler,option,grid,saturation_functions, &
                                   sat_func_id)
  ! 
  ! Computes the pressures for a saturation
  ! initial/boundary condition
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/07/11
  ! 

  use Option_module
  use Grid_module
  use Coupler_module
  use Condition_module
  use Connection_module
  use Region_module
  use Saturation_Function_module

  implicit none

  type(coupler_type) :: coupler
  type(option_type) :: option
  type(grid_type) :: grid
  type(saturation_function_ptr_type) :: saturation_functions(:)
  PetscInt :: sat_func_id(:)

  PetscInt :: local_id, ghosted_id, iconn
  PetscReal :: saturation
  PetscReal :: capillary_pressure
  PetscReal :: liquid_pressure
  
  type(flow_condition_type), pointer :: condition
  
  type(connection_set_type), pointer :: cur_connection_set
  
  condition => coupler%flow_condition

  if (option%iflowmode /= RICHARDS_MODE) then
    option%io_buffer = 'SaturationUpdateCoupler is not set up for this flow mode.'
    call printErrMsg(option)
  endif
  
  ! in this case, the saturation is stored within concentration dataset
  saturation = condition%saturation%dataset%rarray(1)

  do iconn = 1, coupler%connection_set%num_connections
    local_id = coupler%connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)
    call SatFuncGetCapillaryPressure(capillary_pressure,saturation, &
                     option%reference_temperature, &
                     saturation_functions(sat_func_id(ghosted_id))%ptr,option)
    liquid_pressure = option%reference_pressure - capillary_pressure
    coupler%flow_aux_real_var(1,iconn) = liquid_pressure
  enddo

end subroutine SaturationUpdateCoupler

end module Saturation_module

