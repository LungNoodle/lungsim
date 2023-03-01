module ventilation_c
  implicit none
  private

contains

!!!###################################################################################

  subroutine evaluate_vent_c() bind(C, name="evaluate_vent_c")

    use arrays,only: dp
    use ventilation, only: evaluate_vent
    implicit none

    call evaluate_vent()

  end subroutine evaluate_vent_c


  !###################################################################################

  subroutine evaluate_uniform_flow_c() bind(C, name="evaluate_uniform_flow_c")

    use ventilation, only: evaluate_uniform_flow
    implicit none

    call evaluate_uniform_flow

  end subroutine evaluate_uniform_flow_c


!###################################################################################

  subroutine two_unit_test_c() bind(C, name="two_unit_test_c")
    use ventilation, only: two_unit_test
    implicit none

    call two_unit_test

  end subroutine two_unit_test_c

!###################################################################################
end module ventilation_c
