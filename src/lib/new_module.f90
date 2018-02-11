module new_module
!*Brief Description:* A one line descriptor of what the module does.
!
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
  public test_function

contains
!
!##############################################################################
!
 subroutine test_function()
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_TEST_FUNCTION: TEST_FUNCTION
   use indices,only: ne_radius
   use arrays, only: dp
   use diagnostics, only: enter_exit

   !local variables
   integer :: other,stuff

   character(len=60) :: sub_name

   sub_name = 'test_function'
   call enter_exit(sub_name,1)

   write(*,*) ne_radius

   call enter_exit(sub_name,2)
 end subroutine test_function

!
!###########################################################################################
!
end module new_module
