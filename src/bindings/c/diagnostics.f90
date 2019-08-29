module diagnostics_c

  implicit none

  private

contains

  !!!######################################################################
  subroutine enter_exit_c(sub_name, sub_name_len, state) bind(C, name="enter_exit_c")
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use diagnostics, only: enter_exit
    use other_consts, only: MAX_STRING_LEN
    implicit none
    integer,intent(in) :: state, sub_name_len
    type(c_ptr), value, intent(in) :: sub_name
    character(len=MAX_STRING_LEN) :: sub_name_f

    call strncpy(sub_name_f, sub_name, sub_name_len)
#if defined _WIN32 && defined __INTEL_COMPILER
    call so_enter_exit(sub_name_f, state)
#else
    call enter_exit(sub_name_f, state)
#endif

  end subroutine enter_exit_c

  !!!######################################################################
  subroutine set_diagnostics_on_c(state) bind(C, name="set_diagnostics_on_c")
    use diagnostics, only: set_diagnostics_on
    implicit none

    logical, intent(in) :: state

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_set_diagnostics_on(state)
#else
    call set_diagnostics_on(state)
#endif

  end subroutine set_diagnostics_on_c

  !!!######################################################################
  subroutine get_diagnostics_on_c(state) bind(C, name="get_diagnostics_on_c")
    use diagnostics, only: get_diagnostics_on
    implicit none

    logical :: state

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_get_diagnostics_on(state)
#else
    call get_diagnostics_on(state)
#endif

  end subroutine get_diagnostics_on_c

end module diagnostics_c
