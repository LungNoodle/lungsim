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

  implicit none
  logical :: diagnostics_on

  private
  public enter_exit, set_diagnostics_on

contains

!!!######################################################################

  subroutine enter_exit(sub_name,type)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_ENTER_EXIT" :: ENTER_EXIT
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
    implicit none
    
    logical, intent(in) :: state
    
    diagnostics_on = state
    
  end subroutine set_diagnostics_on

end module diagnostics
