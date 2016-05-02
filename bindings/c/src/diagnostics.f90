module diagnostics_c

  implicit none

  private

contains

  !!!######################################################################
  subroutine enter_exit_c(sub_name, sub_name_len, state) bind(C, name="enter_exit_c")
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_ENTER_EXIT_C" :: ENTER_EXIT_C
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use diagnostics, only: enter_exit
    use other_consts, only: MAX_STRING_LEN
    implicit none
    integer,intent(in) :: state, sub_name_len
    type(c_ptr), value, intent(in) :: sub_name
    character(len=MAX_STRING_LEN) :: sub_name_f

    call strncpy(sub_name_f, sub_name, sub_name_len)
    call enter_exit(sub_name_f, state)

  end subroutine enter_exit_c

  !!!######################################################################
  subroutine set_diagnostics_on_c(state) bind(C, name="set_diagnostics_on_c")
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_SET_DIAGNOSTICS_ON_C" :: SET_DIAGNOSTICS_ON_C
    use diagnostics, only: set_diagnostics_on
    implicit none

    logical, intent(in) :: state

    call set_diagnostics_on(state)

  end subroutine set_diagnostics_on_c

end module diagnostics_c
