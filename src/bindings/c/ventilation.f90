module ventilation_c
  implicit none
  private

contains

!!!###################################################################################

  subroutine evaluate_vent_c(num_breaths, dt) bind(C, name="evaluate_vent_c")

    use arrays,only: dp
    use ventilation, only: evaluate_vent
    implicit none

    integer, intent(in) :: num_breaths
    real(dp), intent(in) :: dt

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_evaluate_vent(num_breaths, dt)
#else
    call evaluate_vent(num_breaths, dt)
#endif

  end subroutine evaluate_vent_c


  !###################################################################################

  subroutine evaluate_uniform_flow_c() bind(C, name="evaluate_uniform_flow_c")

    use ventilation, only: evaluate_uniform_flow
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_evaluate_uniform_flow
#else
    call evaluate_uniform_flow
#endif

  end subroutine evaluate_uniform_flow_c


!###################################################################################

  subroutine two_unit_test_c() bind(C, name="two_unit_test_c")
    use ventilation, only: two_unit_test
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_two_unit_test
#else
    call two_unit_test
#endif

  end subroutine two_unit_test_c

!###################################################################################
end module ventilation_c
