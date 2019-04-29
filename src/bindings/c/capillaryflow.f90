module capillaryflow_c

  implicit none

  private

contains

  !!!######################################################################
  subroutine calc_cap_imped_c(ha, hv, omega) bind(C, name="calc_cap_imped_c")

    use capillaryflow, only: calc_cap_imped
    use arrays, only: dp

    implicit none

    real(dp),intent(in) :: ha,hv,omega

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_calc_cap_imped(ha,hv,omega)
#else
    call calc_cap_imped(ha,hv,omega)
#endif

  end subroutine calc_cap_imped_c
  !!!######################################################################

end module capillaryflow_c
