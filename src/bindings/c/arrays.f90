module arrays_c
!*Brief Description:* This module wraps part of the arrays module that require a c interface
!
!*LICENSE:*
!
!
!*Contributor(s):* Merryn Tawhai, Alys Clark
!
!*Full Description:*
!
!This module wraps part of the arrays module that require a c interface
  use arrays,only: dp

  implicit none
  public set_node_field_value_c,update_parameter_c

contains
  subroutine set_node_field_value_c(row, col, value) bind(C, name="set_node_field_value_c")
    use arrays, only: set_node_field_value
    implicit none

    integer, intent(in) :: row, col
    real(dp), intent(in) :: value

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_set_node_field_value(row, col, value)
#else
    call set_node_field_value(row, col, value)
#endif

  end subroutine set_node_field_value_c

  subroutine update_parameter_c(parameter_name, parameter_name_len, parameter_value) bind(C, name="update_parameter_c")
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: max_filename_len
    use arrays, only: update_parameter
    implicit none

    integer, intent(in) :: parameter_name_len
    real(dp), intent(in) :: parameter_value
    type(c_ptr),value, intent(in) :: parameter_name
    character(len=max_filename_len) :: parameter_name_f

    call strncpy(parameter_name_f, parameter_name, parameter_name_len)
#if defined _WIN32 && defined __INTEL_COMPILER
    call so_update_parameter(parameter_name_f, parameter_value)
#else
    call update_parameter(parameter_name_f, parameter_value)
#endif
  end subroutine update_parameter_c

    

end module arrays_c
