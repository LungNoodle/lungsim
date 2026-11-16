module pressure_resistance_flow_c
implicit none
private

contains

  subroutine set_perfusion_flows_c(element_ids_len,element_ids,flow_values_len,flow_values) &
      bind(C,name="set_perfusion_flows_c")
    use arrays,only: dp
    use pressure_resistance_flow,only: set_perfusion_flows
    implicit none
    integer,intent(in) :: element_ids_len,flow_values_len
    integer,intent(in) :: element_ids(element_ids_len)
    real(dp),intent(in) :: flow_values(flow_values_len)

    if(element_ids_len.ne.flow_values_len)then
      print *, 'ERROR: prescribed perfusion element and flow lists must have the same length'
      call exit(1)
    endif
    call set_perfusion_flows(element_ids,flow_values)
  end subroutine set_perfusion_flows_c

  subroutine clear_perfusion_flows_c() bind(C,name="clear_perfusion_flows_c")
    use pressure_resistance_flow,only: clear_perfusion_flows
    implicit none

    call clear_perfusion_flows()
  end subroutine clear_perfusion_flows_c



  !
  !###################################################################################
  !
  ! the main growing subroutine. Generates a volume-filling tree into a closed surface.
  subroutine occlusion_list_c(surface_elems_len, surface_elems) bind(C, name="occlusion_list_c")

    !use arrays,only: dp
    !use iso_c_binding, only: c_ptr
    !use utils_c, only: strncpy
    !use other_consts, only: MAX_FILENAME_LEN
    use pressure_resistance_flow,only: occlusion_list
    implicit none

    integer,intent(in) :: surface_elems_len
    integer,intent(in) :: surface_elems(surface_elems_len)

    call occlusion_list(surface_elems)

  end subroutine occlusion_list_c
  !
  !###################################################################################
  !
  ! makes the vessels given rigid.
  subroutine compliance_list_c(surface_elems_len, surface_elems) bind(C, name="compliance_list_c")

    !use arrays,only: dp
    !use iso_c_binding, only: c_ptr
    !use utils_c, only: strncpy
    !use other_consts, only: MAX_FILENAME_LEN
    use pressure_resistance_flow,only: compliance_list
    implicit none

    integer,intent(in) :: surface_elems_len
    integer,intent(in) :: surface_elems(surface_elems_len)

    call compliance_list(surface_elems)

  end subroutine compliance_list_c

!!!###################################################################################

subroutine evaluate_prq_c(mesh_type,mesh_type_len,vessel_type,vessel_type_len,grav_dirn,grav_factor,bc_type,bc_type_len,inlet_bc, &
               outlet_bc,remodeling_grade) bind(C, name="evaluate_prq_c")

  use iso_c_binding, only: c_ptr
  use utils_c, only: strncpy
  use other_consts, only: MAX_STRING_LEN
  use arrays, only: dp
  use pressure_resistance_flow, only: evaluate_prq
  implicit none

  type(c_ptr), value, intent(in) :: mesh_type,bc_type,vessel_type
  integer,intent(in) :: mesh_type_len,bc_type_len, vessel_type_len,grav_dirn
  character(len=MAX_STRING_LEN) :: mesh_type_f,bc_type_f,vessel_type_f
  real(dp),intent(in) :: grav_factor,inlet_bc,outlet_bc,remodeling_grade

  call strncpy(mesh_type_f, mesh_type, mesh_type_len)
  call strncpy(bc_type_f, bc_type, bc_type_len)
  call strncpy(vessel_type_f, vessel_type, vessel_type_len)

  call evaluate_prq(mesh_type_f,vessel_type_f,grav_dirn,grav_factor,bc_type_f,inlet_bc,outlet_bc,remodeling_grade)

end subroutine evaluate_prq_c

!###################################################################################
end module pressure_resistance_flow_c
