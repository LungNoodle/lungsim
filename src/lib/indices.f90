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

  use arrays, only: dp, num_units, unit_field
  use diagnostics
  use other_consts
  
  implicit none
  
  !parameters
  ! indices for elem_ordrs
  integer :: num_ord=4,no_gen=1,no_hord=2,no_sord=3,no_type = 4
  ! indices for node_fields
  integer :: num_nj,nj_aw_press=0,nj_bv_press=0,nj_conc1=0,&
       nj_conc2=0,nj_pp=0
  ! indices for elem_field
  integer ::num_ne,ne_radius=0,ne_length=0,ne_mass=0,ne_vol=0,&
       ne_resist=0,ne_t_resist=0,ne_Vdot=0,ne_Vdot0=0,ne_a_A=0,&
       ne_dvdt=0,ne_radius_in=0,ne_radius_in0=0,&
       ne_radius_out=0,ne_radius_out0=0,ne_group=0,ne_Qdot=0, &
       ne_vd_bel=0, ne_vol_bel=0, ne_unit=0, ne_lymph_flux=0, ne_lymph_intsat=0
  ! indices for unit_field
  integer :: num_nu,nu_vol=0,nu_comp=0,nu_conc2=0,nu_Vdot0=0,nu_Vdot1=0, &
       nu_Vdot2=0,nu_dpdt=0,nu_pe=0,nu_vt=0,nu_air_press=0,nu_conc1=0,nu_vent=0,&
       nu_vd=0,nu_perf=0,nu_blood_press=0,nu_capillary_press=0,nu_sa=0,nu_tt=0,nu_Pe_max=0, &
       nu_Pe_min=0,nu_flux=0,nu_intsat=0,nu_av_flux=0,nu_lymphflow=0,nu_time=0
  !indices for gas exchange field
  ! indices for gasex_field
  integer,parameter :: num_gx = 12
  integer,parameter :: ng_p_alv_o2=1      ! index for alveolar partial pressure of O2
  integer,parameter :: ng_p_alv_co2=2     ! index for alveolar partial pressure of CO2
  integer,parameter :: ng_p_ven_o2=3      ! index for local venous partial pressure of O2
  integer,parameter :: ng_p_ven_co2=4     ! index for local venous partial pressure of CO2
  integer,parameter :: ng_p_cap_o2=5      ! index for local end capillary partial pressure of O2
  integer,parameter :: ng_p_cap_co2=6     ! index for local end capillary partial pressure of CO2
  integer,parameter :: ng_source_o2=7     ! index for source (flux) of O2
  integer,parameter :: ng_source_co2=8    ! index for source (flux) of CO2
  integer,parameter :: ng_Vc=9            ! index for unit's capillary blood volume
  integer,parameter :: ng_sa=10           ! index for unit's capillary surface area
  integer,parameter :: ng_tt=11           ! index for transit time in unit
  integer,parameter :: ng_time=12         ! index for time elapsed for RBC in capillaries
  
  !model type
  character(len=60) :: model_type
  
  public num_ord,no_gen,no_hord,no_sord,no_type
  
  public num_nj,nj_aw_press,nj_bv_press,nj_conc1,nj_conc2, nj_pp
  
  public num_ne,ne_radius,ne_length,ne_mass,ne_vol,&
       ne_resist,ne_t_resist,ne_Vdot,ne_Vdot0,ne_a_A,&
       ne_dvdt,ne_radius_in,ne_radius_in0,ne_radius_out,&
       ne_radius_out0,ne_group,ne_Qdot, &
       ne_vd_bel, ne_vol_bel, ne_unit, ne_lymph_flux, ne_lymph_intsat
  
  public num_nu,nu_vol,nu_comp, nu_conc2,nu_Vdot0,nu_Vdot1, &
       nu_Vdot2,nu_dpdt,nu_pe,nu_vt,nu_air_press,&
       nu_conc1,nu_vent,nu_vd,&
       nu_perf,nu_blood_press,nu_capillary_press,nu_sa,nu_tt,nu_Pe_max,nu_Pe_min,nu_flux,nu_intsat, &
       nu_av_flux,nu_lymphflow,nu_time
  
  public num_gx, ng_p_alv_o2,ng_p_alv_co2,ng_p_ven_o2,ng_p_ven_co2, &
       ng_p_cap_o2, ng_p_cap_co2,ng_source_o2,ng_source_co2, &
       ng_Vc, ng_sa, ng_tt, ng_time
  
  
  public model_type
  
  !Interfaces
  private
  public define_problem_type, lymphatic_indices, ventilation_indices, perfusion_indices, &
       get_ne_radius, get_ne_lymph_flux, get_ne_lymph_intsat, get_nj_conc1, growing_indices
  
contains
  
  !> Define problem type
  subroutine define_problem_type(PROBLEM_TYPE)
    
    character(len=MAX_FILENAME_LEN),intent(in) :: PROBLEM_TYPE
    
    character(len=60) :: sub_name
    
    sub_name = 'define_problem_type'
    call enter_exit(sub_name,1)
    select case (PROBLEM_TYPE)
    case ('gas_exchange')
       print *, 'You are solving a gas exchange model, setting up indices'
       call exchange_indices
    case ('gas_mix')
       print *, 'You are solving a gas mixing model, setting up indices'
       call gasmix_indices
    case ('gas_transfer')
       print *, 'You are solving a gas transfer model, setting up indices'
       call exchange_indices
    case ('lymphatic_transport')
       print *, 'You are solving a lymphatic transport model, setting up indices'
       call lymphatic_indices
       call update_units
    case ('perfusion')
       print *, 'You are solving a static perfusion model, setting up indices'
       call perfusion_indices
    case ('ventilation')
       print *, 'You are solving a ventilation model, setting up indices'
       call ventilation_indices
    case('grow_tree')
       print *, 'You are solving a growing problem, setting up indices'
       call growing_indices
    end select
    model_type=TRIM(PROBLEM_TYPE)
    call enter_exit(sub_name,2)
  end subroutine define_problem_type
  
  !>Gas mixing indices
  subroutine exchange_indices
    
    character(len=60) :: sub_name
    
    sub_name = 'exchange_indices'
    call enter_exit(sub_name,1)
    ne_unit = 0
    ne_lymph_flux = 0
    ne_lymph_intsat = 0
    ! indices for elem_ordrs. These dont usually change.
    ! indices for node_field
    num_nj=5
    nj_conc1=2
    nj_conc2=3
    nj_aw_press=4 !air pressure
    nj_pp = 5
    
    ! indices for elem_field
    num_ne = 15
    ne_radius = 1
    ne_length = 2
    ne_vol = 3
    ne_resist = 4
    ne_t_resist = 5
    ne_Vdot = 6 !Air flow, current time step
    ne_Vdot0 = 7 !air flow, last timestep
    ne_dvdt = 8
    ne_vd_bel = 9
    ne_vol_bel = 10
    ne_Qdot = 11
    ne_mass = 12
    ne_a_A = 13
    ne_radius_in = 14
    ne_radius_out = 15
    
    ! indices for unit_field
    num_nu=14
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
    nu_vd=11
    nu_perf=12
    nu_conc1=13
    nu_conc2=14
    
    call enter_exit(sub_name,2)
  end subroutine exchange_indices
  
  !>Gas mixing indices
  subroutine gasmix_indices
    
    character(len=60) :: sub_name
    
    sub_name = 'gasmix_indices'
    call enter_exit(sub_name,1)
    ne_unit = 0
    ne_lymph_flux = 0
    ne_lymph_intsat = 0
    ! indices for elem_ordrs. These dont usually change.
    ! indices for node_field
    num_nj=3
    nj_aw_press=2
    nj_conc1=3
    ! indices for elem_field
    num_ne = 11
    ne_radius = 1
    ne_length = 2
    ne_vol = 3
    ne_resist = 4
    ne_t_resist = 5
    ne_Vdot = 6
    ne_Vdot0 = 7
    ne_a_A = 8
    ne_dvdt = 9
    ne_vd_bel = 10
    ne_vol_bel = 11
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
    
    character(len=60) :: sub_name
    
    sub_name = 'ventilation_indices'
    call enter_exit(sub_name,1)
    ne_unit = 0
    ne_lymph_flux = 0
    ne_lymph_intsat = 0
    ! indices for elem_ordrs. These dont usually change.
    ! indices for node_field
    num_nj=2 !number of nodal fields
    nj_aw_press=2 !air pressure
    ! indices for elem_field
    num_ne = 10 !number of element fields
    ne_radius = 1 !radius of airway
    ne_length = 2 !length of airway
    ne_vol = 3 !volume
    ne_resist = 4 !resistance of airway
    ne_t_resist = 5
    ne_Vdot = 6 !Air flow, current time step
    ne_Vdot0 = 7 !air flow, last timestep
    ne_dvdt = 8
    ne_vd_bel = 9
    ne_vol_bel = 10
    ! indices for unit_field
    num_nu=12
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
    nu_Pe_max=11
    nu_Pe_min=12
    call enter_exit(sub_name,2)
  end subroutine ventilation_indices

!!!#############################################################################

  !> Indices used by the published pulmonary lymphatic transport model.
  subroutine lymphatic_indices

    character(len=60) :: sub_name

    sub_name = 'lymphatic_indices'
    call enter_exit(sub_name,1)

    num_nu = 12
    nu_perf = 1
    nu_blood_press = 2
    nu_sa = 3
    nu_tt = 4
    nu_capillary_press = 5
    nu_Pe_max = 6
    nu_Pe_min = 7
    nu_flux = 8
    nu_intsat = 9
    nu_av_flux = 10
    nu_lymphflow = 11
    nu_time = 12

    call enter_exit(sub_name,2)

  end subroutine lymphatic_indices

!!!#############################################################################

  subroutine growing_indices
    !* Growing indices:* set up indices for growing (1D tree) arrays
    
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'growing_indices'
    call enter_exit(sub_name,1)
    ne_unit = 0
    ne_lymph_flux = 0
    ne_lymph_intsat = 0
    
    ! indices for elem_ordrs. These dont usually change.
    ! indices for node_field
    num_nj = 0 !number of nodal fields
    ! indices for elem_field
    num_ne = 8 !number of element fields
    ne_radius = 1 !radius of branch
    ne_radius_in = 2
    ne_radius_out = 3
    ne_length = 4 !length of branch
    ne_vol = 5
    ne_a_A = 6 !ratio of duct to total cross-section (airway)
    ne_vd_bel = 7
    ne_vol_bel = 8
    ! indices for unit_field
    num_nu = 0
    
    call enter_exit(sub_name,2)
    
  end subroutine growing_indices
  !
  !######################################################################
  !
  !> Perfusion indices
  subroutine perfusion_indices
    
    character(len=60) :: sub_name
    
    sub_name = 'perfusion_indices'
    call enter_exit(sub_name,1)
    
    ! indices for node_field
    num_nj=1
    nj_bv_press=1 !pressure in blood vessel
    ! indices for elem_field
    num_ne=12
    ne_radius=1 !strained average radius over whole element
    ne_radius_in=2 !strained radius into an element
    ne_radius_out=3 !strained radius out of an element
    ne_length=4!length of an elevent
    ne_radius_in0=5!unstrained radius into an element
    ne_radius_out0=6!unstrained radius out of an element
    ne_Qdot=7 !flow in an element
    ne_resist=8 !resistance of a blood vessel
    ne_group=9!Groups vessels into arteries (field=0), capillaries (field=1) and veins(field=2)
    ne_unit=10 !maps terminal arterial elements to unit_field columns
    ne_lymph_flux=11 !averaged lymphatic flux mapped onto the arterial tree
    ne_lymph_intsat=12 !interstitial saturation mapped onto the arterial tree
    !indices for units
    num_nu=5
    nu_perf=1
    nu_blood_press=2
    nu_sa=3
    nu_tt=4
    nu_capillary_press=5
    
    call enter_exit(sub_name,2)
  end subroutine perfusion_indices

  subroutine update_units
    ! Reallocate unit_field when switching from perfusion to lymphatics while
    ! preserving the fields already calculated or imported.
    integer :: fields_to_copy, num_nu_old
    real(dp), allocatable :: unit_temp(:,:)

    if (allocated(unit_field)) then
       num_nu_old = size(unit_field,1)
       allocate(unit_temp(num_nu_old,num_units))
       unit_temp = unit_field
       deallocate(unit_field)
       allocate(unit_field(num_nu,num_units))
       unit_field = 0.0_dp
       fields_to_copy = min(num_nu_old,num_nu)
       unit_field(1:fields_to_copy,:) = unit_temp(1:fields_to_copy,:)
       deallocate(unit_temp)
    endif

  end subroutine update_units
  
  function get_ne_radius() result(res)
    
    implicit none
    character(len=60) :: sub_name
    integer :: res
    
    sub_name = 'get_ne_radius'
    call enter_exit(sub_name,1)
    
    res=ne_radius
    
    call enter_exit(sub_name,2)
  end function get_ne_radius

  function get_ne_lymph_flux() result(res)

    integer :: res

    res = ne_lymph_flux

  end function get_ne_lymph_flux

  function get_ne_lymph_intsat() result(res)

    integer :: res

    res = ne_lymph_intsat

  end function get_ne_lymph_intsat
  
  function get_nj_conc1() result(res)
    
    character(len=60) :: sub_name
    integer :: res
    
    sub_name = 'get_nj_conc1'
    call enter_exit(sub_name,1)
    
    res = nj_conc1
    
    call enter_exit(sub_name,2)
  end function get_nj_conc1
  
end module indices
