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
  ! indices for node_fields
  integer :: num_nj,nj_aw_press,nj_bv_press,nj_conc1
  ! indices for elem_field
  integer ::num_ne,ne_radius,ne_length,ne_vol,&
      ne_resist,ne_t_resist,ne_Vdot,ne_Vdot0,ne_a_A,&
       ne_dvdt,ne_radius_in,ne_radius_in0,&
       ne_radius_out,ne_radius_out0,ne_group,ne_Qdot
  ! indices for unit_field
  integer :: num_nu,nu_vol,nu_comp,nu_Vdot0,nu_Vdot1, &
       nu_Vdot2,nu_dpdt,nu_pe,nu_vt,nu_air_press,nu_conc1,nu_vent,&
       nu_perf,nu_blood_press
  character :: problem_type

public num_ord,no_gen,no_hord,no_sord

public num_nj,nj_aw_press,nj_bv_press, nj_conc1

public num_ne,ne_radius,ne_length,ne_vol,&
      ne_resist,ne_t_resist,ne_Vdot,ne_Vdot0,ne_a_A,&
      ne_dvdt,ne_radius_in,ne_radius_in0,ne_radius_out,&
      ne_radius_out0,ne_group,ne_Qdot

public num_nu,nu_vol,nu_comp,nu_Vdot0,nu_Vdot1, &
       nu_Vdot2,nu_dpdt,nu_pe,nu_vt,nu_air_press,&
       nu_conc1,nu_vent,&
       nu_perf,nu_blood_press

public problem_type

!Interfaces
private
public define_problem_type,ventilation_indices, perfusion_indices, get_ne_radius, get_nj_conc1

contains

  !> Define problem type
  subroutine define_problem_type(PROBLEM_TYPE)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_PROBLEM_TYPE" :: DEFINE_PROBLEM_TYPE
    use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
    use diagnostics, only: enter_exit

    character(len=MAX_FILENAME_LEN),intent(in) :: PROBLEM_TYPE
    character(len=60) :: sub_name
    write(*,*) PROBLEM_TYPE
    select case (PROBLEM_TYPE)
      case ('gas_exchange')
        print *, 'You are solving a gas exchange model, setting up indices'
      case ('gas_mix')
        print *, 'You are solving a gas mixing model, setting up indices'
        call gasmix_indices
      case ('gas_transfer')
        print *, 'You are solving a gas transfer model, setting up indices'
      case ('perfusion')
        print *, 'You are solving a static perfusion model, setting up indices'
        call perfusion_indices
      case ('ventilation')
        print *, 'You are solving a ventilation model, setting up indices'
        call ventilation_indices
    end select

    call enter_exit(sub_name,2)
  end subroutine define_problem_type

  !>Gas mixing indices
  subroutine gasmix_indices
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GASMIX_INDICES" :: GASMIX_INDICES

    use diagnostics, only: enter_exit
    implicit none
    character(len=60) :: sub_name

    sub_name = 'ventilation_indices'
    call enter_exit(sub_name,1)
    ! indices for elem_ordrs. These dont usually change.
    ! indices for node_field
    num_nj=3
    nj_aw_press=2
    nj_conc1=3
    ! indices for elem_field
    num_ne=9
    ne_radius=1
    ne_length=2
    ne_vol=3
    ne_resist=4
    ne_t_resist=5
    ne_Vdot=6
    ne_Vdot0=7
    ne_a_A=8
    ne_dvdt=9
    ! indices for unit_field
    num_nu=11
    nu_vol=1
    nu_comp=2
    nu_Vdot0=3
    nu_Vdot1=4
    nu_Vdot2=5
    nu_dpdt=6
    nu_pe=7
    nu_vt=8
    nu_air_press=9
    nu_conc1=10
    nu_vent=11
    call enter_exit(sub_name,2)
  end subroutine gasmix_indices

  !> Ventilation indices
  subroutine ventilation_indices
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_VENTILATION_INDICES" :: VENTILATION_INDICES

    use diagnostics, only: enter_exit
    implicit none
    character(len=60) :: sub_name

    sub_name = 'ventilation_indices'
    call enter_exit(sub_name,1)
    ! indices for elem_ordrs. These dont usually change.
    ! indices for node_field
    num_nj=2 !number of nodal fields
    nj_aw_press=2 !air pressure
    ! indices for elem_field
    num_ne=8 !number of element fields
    ne_radius=1 !radius of airway
    ne_length=2 !length of airway
    ne_vol=3 !volume
    ne_resist=4 !resistance of airway
    ne_t_resist=5
    ne_Vdot=6 !Air flow, current time step
    ne_Vdot0=7 !air flow, last timestep
    ne_dvdt=8
    ! indices for unit_field
    num_nu=10
    nu_vol=1
    nu_comp=2
    nu_Vdot0=3
    nu_Vdot1=4
    nu_Vdot2=5
    nu_dpdt=6
    nu_pe=7
    nu_vt=8
    nu_air_press=9
    nu_vent=10
    call enter_exit(sub_name,2)
  end subroutine ventilation_indices
!
!######################################################################
!
!> Perfusion indices
  subroutine perfusion_indices
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_PERFUSION_INDICES" :: PERFUSION_INDICES

    use diagnostics, only: enter_exit
    implicit none
    character(len=60) :: sub_name

    sub_name = 'perfusion_indices'
    call enter_exit(sub_name,1)

    ! indices for node_field
    num_nj=1
    nj_bv_press=1 !pressure in blood vessel
    ! indices for elem_field
    num_ne=9
    ne_radius=1 !strained average radius over whole element
    ne_radius_in=2 !strained radius into an element
    ne_radius_out=3 !strained radius out of an element
    ne_length=4!length of an elevent
    ne_radius_in0=5!unstrained radius into an element
    ne_radius_out0=6!unstrained radius out of an element
    ne_Qdot=7 !flow in an element
    ne_resist=8 !resistance of a blood vessel
    ne_group=9!Groups vessels into arteries (field=0), capillaries (field=1) and veins(field=2)
    !indices for units
    num_nu=2
    nu_perf=1
    nu_blood_press=2

     call enter_exit(sub_name,2)
  end subroutine perfusion_indices

  function get_ne_radius() result(res)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GET_NE_RADIUS" :: GET_NE_RADIUS

    use diagnostics, only: enter_exit

    implicit none
    character(len=60) :: sub_name
    integer :: res

    sub_name = 'get_ne_radius'
    call enter_exit(sub_name,1)

    res=ne_radius

    call enter_exit(sub_name,2)
  end function get_ne_radius

  function get_nj_conc1() result(res)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GET_NJ_CONC1" :: GET_NJ_CONC1

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
