!> \file
!> \author Merryn Tawhai, Alys Clark
!> \brief This module handles all geometry read/write/generation. 
!>
!> \section LICENSE
!>
!>
!> Contributor(s):
!>
!>\Description
!> This module handles all all geometry read/write/generation. 

module indices
  implicit none
!parameters
  ! indices for elem_ordrs 
  integer :: num_ord=3,no_gen=1,no_hord=2,no_sord=3
  ! indices for node_field
  integer :: num_nj,nj_radius,nj_radius0,nj_press,nj_conc1
  ! indices for elem_field
  integer ::num_ne,ne_radius,ne_length,ne_vol,&
      ne_resist,ne_t_resist,ne_flow,ne_flow0,ne_a_A,&
       ne_dvdt,ne_radius_in, ne_radius_out
  ! indices for unit_field
  integer :: num_nu,nu_vol,nu_comp,nu_flow0,nu_flow1, &
       nu_flow2,nu_dpdt,nu_pe,nu_vt,nu_press,nu_conc1,nu_vent

public num_ord,no_gen,no_hord,no_sord,&
      num_nj,nj_radius,nj_radius0,nj_press,nj_conc1,&
       num_ne,ne_radius,ne_length,ne_vol,&
      ne_resist,ne_t_resist,ne_flow,ne_flow0,ne_a_A,&
       ne_dvdt,ne_radius_in,ne_radius_out,&
      num_nu,nu_vol,nu_comp,nu_flow0,nu_flow1, &
       nu_flow2,nu_dpdt,nu_pe,nu_vt,nu_press,nu_conc1,nu_vent

!Interfaces
private
public ventilation_indices, perfusion_indices, get_ne_radius, get_nj_conc1

contains

  !> Ventilation indices 
  subroutine ventilation_indices
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_VENTILATION_INDICES" :: VENTILATION_INDICES

    use diagnostics, only: enter_exit
    implicit none
    character(len=60) :: sub_name
    
    sub_name = 'ventilation_indices'
    call enter_exit(sub_name,1)
    ! indices for elem_ordrs. These dont usually change.
    ! indices for node_field
    num_nj=3
    nj_press=2 
    nj_conc1=3
    ! indices for elem_field
    num_ne=9
    ne_radius=1 
    ne_length=2 
    ne_vol=3
    ne_resist=4
    ne_t_resist=5
    ne_flow=6 
    ne_flow0=7
    ne_a_A=8 
    ne_dvdt=9
    ! indices for unit_field
    num_nu=11
    nu_vol=1 
    nu_comp=2
    nu_flow0=3 
    nu_flow1=4
    nu_flow2=5 
    nu_dpdt=6
    nu_pe=7
    nu_vt=8
    nu_press=9
    nu_conc1=10
    nu_vent=11
    call enter_exit(sub_name,2)
  end subroutine ventilation_indices
!
!######################################################################
!
!> Perfusion indices 
  subroutine perfusion_indices
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_PERFUSION_INDICES" :: PERFUSION_INDICES

    use diagnostics, only: enter_exit
    implicit none
    character(len=60) :: sub_name
    
    sub_name = 'perfusion_indices'
    call enter_exit(sub_name,1)

    ! indices for node_field
    num_nj=1
    nj_press=1 !pressure
    ! indices for elem_field
    num_ne=5
    ne_radius=1 !strained average radius over whole element
    ne_radius_in=2 !unstrained radius into an element
    ne_radius_out=3 !unstrained radius out of an element
    ne_length=4
    !ne_group=5!Groups vessels into arteries (field=0), capillaries (field=1) and veins(field=2)
    !ne_vol=3, ne_resist=4,ne_t_resist=5,ne_flow=6,ne_flow0=7,ne_radius0=8

     call enter_exit(sub_name,2)
  end subroutine perfusion_indices

  function get_ne_radius result(res)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_GET_NE_RADIUS" :: GET_NE_RADIUS

    use diagnostics, only: enter_exit
    implicit none
    character(len=60) :: sub_name
    integer :: res
    
    sub_name = 'get_ne_radius'
    call enter_exit(sub_name,1)
    
    res = ne_radius
    
    call enter_exit(sub_name,2)
  end function get_ne_radius

  function get_nj_conc1 result(res)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_GET_NJ_CONC1" :: GET_NJ_CONC1

    use diagnostics, only: enter_exit
    implicit none
    character(len=60) :: sub_name
    integer :: res
    
    sub_name = 'get_nj_conc1'
    call enter_exit(sub_name,1)
    
    res = nj_conc1
    
    call enter_exit(sub_name,2)
  end function get_nj_conc1

end module indices
