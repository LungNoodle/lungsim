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
  public set_node_field_value_c

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


end module arrays_c
