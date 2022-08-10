module field_utilities_c
  implicit none
  private

contains
  
!
!###################################################################################
!
  subroutine scale_flow_to_inlet_c(INLET_FLOW,VORQ,vq_len) bind(C, name="scale_flow_to_inlet_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use arrays, only: dp
    use field_utilities, only: scale_flow_to_inlet
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: vq_len
    real(dp),intent(in) :: INLET_FLOW
    type(c_ptr), value, intent(in) :: VORQ
    character(len=MAX_FILENAME_LEN) :: vq_f

    call strncpy(vq_f, VORQ, vq_len)

    call scale_flow_to_inlet(INLET_FLOW,vq_f)

  end subroutine scale_flow_to_inlet_c

end module field_utilities_c

