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
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_INITIALISE_TRANSPORT" :: INITIALISE_TRANSPORT
   use indices
   use arrays, only: dp
   use gas_exchange, only: initial_gasexchange
   use diagnostics, only: enter_exit

   !local variables

   character(len=60) :: sub_name

   sub_name = 'initialise_transport'
   call enter_exit(sub_name,1)

   call allocate_memory_speciestrans

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
       !note a linear q gradient should  be set up to scale for shunt fraction automatically
       call initial_gasexchange(149.0_dp)
       call solve_transport

    end select
   call enter_exit(sub_name,2)
 end subroutine initialise_transport

 !
!##############################################################################
!
 subroutine solve_transport()
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_INITIALISE_TRANSPORT" :: INITIALISE_TRANSPORT
   use indices
   use arrays, only: dp
   use gas_exchange, only: steadystate_gasexchange
   use diagnostics, only: enter_exit

   !local variables
   real(dp) c_art_o2, c_ven_o2,p_art_co2,p_art_o2, p_ven_co2,p_ven_o2

   character(len=60) :: sub_name

   sub_name = 'solve_transport'
   call enter_exit(sub_name,1)


   select case (model_type)
     case ('gas_exchange')
       print *, 'Nothing implemented'
       !Note that as V, Q are prerequisites something needs to be added here that checks
       !these have been read in and if not sets up linear gradient based on some default parameters
     case ('gas_mix')
       print *, 'Nothing implemented'
       !Note that as V is prerequisites something needs to be added here that checks
       !these have been read in and if not sets up linear gradient based on some default parameters
     case ('gas_transfer')
       print *, 'Calling gas transfer model '
       p_art_co2=40.0_dp
       p_ven_co2=45.0_dp
       p_art_o2=100.0_dp
       p_ven_o2=40.0_dp
       call steadystate_gasexchange(c_art_o2,c_ven_o2,&
       p_art_co2,p_art_o2,149.0_dp,p_ven_co2,p_ven_o2,0.03_dp,&
       0.8_dp*(260.0_dp*1.0e+3_dp/60.0_dp),260.0_dp*1.0e+3_dp/60.0_dp )

    end select
   call enter_exit(sub_name,2)
 end subroutine solve_transport

!
!###########################################################################################
!
 subroutine allocate_memory_speciestrans()
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_INITIALISE_TRANSPORT" :: INITIALISE_TRANSPORT
   use indices
   use arrays, only: dp,gasex_field,num_units
   use diagnostics, only: enter_exit

   character(len=60) :: sub_name
   sub_name = 'allocate_memory_speciestrans'
   call enter_exit(sub_name,1)

!!! allocate memory for the gasex_field array, if not already allocated
    if(.not.allocated(gasex_field)) allocate(gasex_field(num_gx,num_units))


   call enter_exit(sub_name,2)
 end subroutine allocate_memory_speciestrans

end module species_transport
