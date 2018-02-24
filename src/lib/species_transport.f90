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
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_INITIALISE_TRANSPORT: TEST_INITIALISE_TRANSPORT
   use indices
   use arrays, only: dp
   use diagnostics, only: enter_exit

   !local variables

   character(len=60) :: sub_name

   sub_name = 'initialise_transport'
   call enter_exit(sub_name,1)

   select case (model_type)
     case ('gas_exchange')
       print *, 'You are solving a gas exchange model'
       !Note that as V, Q are prerequisites something needs to be added here that checks
       !these have been read in and if not sets up linear gradient based on some default parameters
     case ('gas_mix')
       print *, 'You are solving a gas mixing model'
       !Note that as V is prerequisites something needs to be added here that checks
       !these have been read in and if not sets up linear gradient based on some default parameters
     case ('gas_transfer')
       print *, 'You are solving a gas transfer model'
       !Note that as V, Q are prerequisites something needs to be added here that checks
       !these have been read in and if not sets up linear gradient based on some default parameters
    end select
   call enter_exit(sub_name,2)
 end subroutine initialise_transport

!
!###########################################################################################
!
end module species_transport
