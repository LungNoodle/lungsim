module ventilation_c
  implicit none
  private

contains

!!!###################################################################################

  subroutine evaluate_flow_c() bind(C, name="evaluate_flow_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_EVALUATE_FLOW_C" :: EVALUATE_FLOW_C

    use ventilation, only: evaluate_flow
    implicit none

    call evaluate_flow

  end subroutine evaluate_flow_c


  !###################################################################################

  subroutine evaluate_uniform_flow_c() bind(C, name="evaluate_uniform_flow_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_EVALUATE_UNIFORM_FLOW_C" :: EVALUATE_UNIFORM_FLOW_C

    use ventilation, only: evaluate_uniform_flow
    implicit none

    call evaluate_uniform_flow

  end subroutine evaluate_uniform_flow_c


!###################################################################################

  subroutine two_unit_test_c() bind(C, name="two_unit_test_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_TWO_UNIT_TEST_C" :: TWO_UNIT_TEST_C
    use ventilation, only: two_unit_test
    implicit none

    call two_unit_test

  end subroutine two_unit_test_c

!###################################################################################
end module ventilation_c
