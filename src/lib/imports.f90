module imports
!*Brief Description:* This module contains all the subroutines required to
!import fields, previous model results, etc.
!*LICENSE:*
!
!
!
!*Full Description:*
!More info on what the module does if necessary
!
  use other_consts
  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private 
  public import_terminal_ventilation

contains
!
!##############################################################################
!
 subroutine import_terminal_ventilation(FLOWFILE)
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_TEST_FUNCTION: TEST_FUNCTION
   use indices,only: ne_radius
   use arrays, only: dp
   use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
   use diagnostics, only: enter_exit

   character(len=MAX_FILENAME_LEN),intent(in) :: FLOWFILE
   !local variables
   integer :: other,stuff

   character(len=60) :: sub_name

   sub_name = 'import_terminal_ventilation'
   call enter_exit(sub_name,1)

   write(*,*) ne_radius

   call enter_exit(sub_name,2)
 end subroutine import_terminal_ventilation

!
!###########################################################################################
!
end module imports
