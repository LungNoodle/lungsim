module ventilation_c
  implicit none
  private

contains

!!!###################################################################################

  subroutine evaluate_vent_c(flow_filename,flow_filename_len, &
       pressure_filename,pressure_filename_len) bind(C,name="evaluate_vent_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use ventilation, only: evaluate_vent
    use other_consts, only: MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: flow_filename_len,pressure_filename_len
    type(c_ptr),value,intent(in) :: flow_filename,pressure_filename
    character(len=MAX_FILENAME_LEN) :: flow_filename_f,pressure_filename_f

    call strncpy(flow_filename_f,flow_filename,flow_filename_len)
    call strncpy(pressure_filename_f,pressure_filename,pressure_filename_len)
    call evaluate_vent(flow_filename_f,pressure_filename_f)

  end subroutine evaluate_vent_c


  !###################################################################################

  subroutine evaluate_uniform_flow_c() bind(C, name="evaluate_uniform_flow_c")

    use ventilation, only: evaluate_uniform_flow
    implicit none

    call evaluate_uniform_flow

  end subroutine evaluate_uniform_flow_c


!###################################################################################
end module ventilation_c
