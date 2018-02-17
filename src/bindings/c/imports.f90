module imports_c
  implicit none
  private

contains
!
!###################################################################################
  
!*test_function:* just for testing the new module
  subroutine import_terminal_ventilation_c(FLOWFILE, filename_len) bind(C, name="import_terminal_ventilation_c")
  
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use imports, only: import_terminal_ventilation
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: FLOWFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, FLOWFILE, filename_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_import_terminal_ventilation(filename_f)
#else
    call import_terminal_ventilation(filename_f)
#endif

  end subroutine import_terminal_ventilation_c
  
end module imports_c

