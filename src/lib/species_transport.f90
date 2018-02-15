module species_transport
!*Brief Description:* This module contains all the subroutines common
!to species transport models, this includes gas exchange, gas mixing,
!and particle transport models
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
  public initialise_transport

contains
!
!##############################################################################
!
 subroutine initialise_transport()
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_TEST_FUNCTION: TEST_FUNCTION
   use indices,only: ne_radius
   use arrays, only: dp
   use diagnostics, only: enter_exit

   !local variables
   integer :: other,stuff

   character(len=60) :: sub_name

   sub_name = 'initialise_transport'
   call enter_exit(sub_name,1)

   write(*,*) ne_radius

   call enter_exit(sub_name,2)
 end subroutine initialise_transport

!
!###########################################################################################
!
end module species_transport
