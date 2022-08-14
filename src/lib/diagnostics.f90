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

  logical,public :: diagnostics_on

  private
  public enter_exit, get_diagnostics_on, set_diagnostics_on

contains

!!!######################################################################

  subroutine enter_exit(sub_name, state)

    integer,intent(in) :: state
    character(len=MAX_SUBNAME_LEN), intent(in) :: sub_name

    if(diagnostics_on)then
      if(state.eq.1)then
        write(*,'('' Entering subroutine '',60A,'':'')') sub_name(1:MAX_SUBNAME_LEN)
      else
        write(*,'('' Exiting subroutine '',60A,'':'')') sub_name(1:MAX_SUBNAME_LEN)
      endif
    endif

  end subroutine enter_exit

!!!######################################################################

  subroutine set_diagnostics_on(state)

    logical, intent(in) :: state

    diagnostics_on = state

  end subroutine set_diagnostics_on

  subroutine get_diagnostics_on(state)
    implicit none

    logical :: state

    state = diagnostics_on

  end subroutine get_diagnostics_on

end module diagnostics
