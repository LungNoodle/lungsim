module pressure_resistance_flow_c
implicit none
private

contains

!!!###################################################################################

subroutine evaluate_prq_c(mesh_type,mesh_type_len,grav_dirn,grav_factor,bc_type,bc_type_len,inlet_bc, &
               outlet_bc) bind(C, name="evaluate_prq_c")

use iso_c_binding, only: c_ptr
use utils_c, only: strncpy
use other_consts, only: MAX_STRING_LEN
use arrays, only: dp
use pressure_resistance_flow, only: evaluate_prq
implicit none

type(c_ptr), value, intent(in) :: mesh_type,bc_type
integer,intent(in) :: mesh_type_len,bc_type_len,grav_dirn
character(len=MAX_STRING_LEN) :: mesh_type_f,bc_type_f
real(dp),intent(in) :: grav_factor,inlet_bc,outlet_bc

call strncpy(mesh_type_f, mesh_type, mesh_type_len)
call strncpy(bc_type_f, bc_type, bc_type_len)

#if defined _WIN32 && defined __INTEL_COMPILER
call so_evaluate_prq(mesh_type_f,grav_dirn,grav_factor,bc_type_f,inlet_bc,outlet_bc)
#else
call evaluate_prq(mesh_type_f,grav_dirn,grav_factor,bc_type_f,inlet_bc,outlet_bc)
#endif

end subroutine evaluate_prq_c

!###################################################################################
end module pressure_resistance_flow_c
