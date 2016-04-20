!> \file
!> \author Merryn Tawhai
!> \brief This module handles diagnostics. 
!>
!> \section LICENSE
!>
!>
!> Contributor(s):
!>
!>\Description
!> This module handles diagnostics
module diagnostics

  use other_consts

  implicit none

  private
  public enter_exit, set_diagnostics_on, set_test_cpu_time

contains

!!!######################################################################

  subroutine enter_exit(sub_name,type)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_ENTER_EXIT" :: ENTER_EXIT
    use other_consts, only: diagnostics_on
    implicit none

    !COMMON /DLL_OTHER_CONSTS/ diagnostics_on
    integer,intent(in) :: type
    character,intent(in) :: sub_name(*)

    integer :: iend

    if(diagnostics_on)then
       iend = len(sub_name)
       iend=60
       if(type.eq.1)then
          write(*,'('' Entering subroutine '',60A,'':'')') sub_name(1:iend)
       else
          write(*,'('' Exiting subroutine '',60A,'':'')') sub_name(1:iend)
       endif
    endif

  end subroutine enter_exit

  subroutine set_diagnostics_on(state)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_SET_DIAGNOSTICS_ON" :: SET_DIAGNOSTICS_ON
    use other_consts, only: diagnostics_on
    implicit none
    
    logical, intent(in) :: state
    
    diagnostics_on = state
    
  end subroutine set_diagnostics_on

  subroutine set_test_cpu_time(state)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_SET_TEST_CPU_TIME" :: SET_TEST_CPU_TIME
    use other_consts, only: test_cpu_time
    implicit none
    
    logical, intent(in) :: state
    
    test_cpu_time = state
    
  end subroutine set_test_cpu_time


end module diagnostics
