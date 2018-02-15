module imports_c
  implicit none
  private

contains
!
!###################################################################################
  
!*test_function:* just for testing the new module
  subroutine import_terminal_ventilation_c() bind(C, name="import_terminal_ventilation_c")
  
    use imports, only: import_terminal_ventilation
    implicit none


#if defined _WIN32 && defined __INTEL_COMPILER
    call so_import_terminal_ventilation
#else
    call import_terminal_ventilation
#endif

  end subroutine import_terminal_ventilation_c
  
end module imports_c

