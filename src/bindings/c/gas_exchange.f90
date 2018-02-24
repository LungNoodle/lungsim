module gas_exchange_c
  implicit none
  private

contains
!
!###################################################################################
  
!*test_function:* just for testing the new module
  subroutine initial_gasexchange_c() bind(C, name="initial_gasexchange_c")
  
    use gas_exchange, only: initial_gasexchange
    implicit none


#if defined _WIN32 && defined __INTEL_COMPILER
    call so_initial_gasexchange
#else
    call initial_gasexchange
#endif

  end subroutine initial_gasexchange_c
  
end module gas_exchange_c

