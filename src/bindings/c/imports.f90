module imports_c
  implicit none
  private

contains
!
!###################################################################################
!
  subroutine import_ventilation_c(FLOWFILE, filename_len) bind(C, name="import_ventilation_c")
  
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use imports, only: import_ventilation
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: FLOWFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, FLOWFILE, filename_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_import_ventilation(filename_f)
#else
    call import_ventilation(filename_f)
#endif

  end subroutine import_ventilation_c
  
!
!###################################################################################
!
  subroutine import_perfusion_c(FLOWFILE, filename_len) bind(C, name="import_perfusion_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use imports, only: import_perfusion
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: FLOWFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, FLOWFILE, filename_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_import_perfusion(filename_f)
#else
    call import_perfusion(filename_f)
#endif

  end subroutine import_perfusion_c

end module imports_c

