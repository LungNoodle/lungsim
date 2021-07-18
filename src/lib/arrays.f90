module arrays
!*Brief Description:* This module defines arrays.
!
!*LICENSE:*
!
!
!*Contributor(s):* Merryn Tawhai, Alys Clark
!
!*Full Description:*
!
!This module defines arrays

  use precision
  
  implicit none

  integer :: num_elems,num_elems_2d,num_nodes,num_data,num_nodes_2d,num_units,num_lines_2d,maxgen

  integer,allocatable :: nodes(:) !allocated in define_node_geometry
  integer,allocatable :: nodes_2d(:) !allocated in define_node_geometry_2d
  integer,allocatable :: node_versn_2d(:) !allocated in define_node_geometry_2d
  integer,allocatable :: elems(:) !allocated in define_1d_elements
  integer,allocatable :: lines_2d(:)
  integer,allocatable :: parentlist(:)
  integer,allocatable :: line_versn_2d(:,:,:)
  integer,allocatable :: lines_in_elem(:,:)
  integer,allocatable :: nodes_in_line(:,:,:)
  integer,allocatable :: elems_2d(:) !allocated in define_elem_geometry_2d
  integer,allocatable :: elem_cnct(:,:,:)  !NXI(-ni:ni,1,ne)
  integer,allocatable :: elem_cnct_2d(:,:,:)
  integer,allocatable :: elem_nodes(:,:)
  integer,allocatable :: elem_nodes_2d(:,:)
  integer,allocatable :: elem_versn_2d(:,:)
  integer,allocatable :: elem_lines_2d(:,:)
  integer,allocatable :: elem_ordrs(:,:)
  integer,allocatable :: elem_symmetry(:)
  integer,allocatable :: elem_units_below(:)
  integer,allocatable :: elems_at_node(:,:)
  integer,allocatable :: elems_at_node_2d(:,:)
  integer,allocatable :: units(:)

  real(dp),allocatable :: arclength(:,:)
  real(dp),allocatable :: elem_field(:,:) !properties of elements
  real(dp),allocatable :: elem_direction(:,:)
  real(dp),allocatable :: node_xyz(:,:)
  real(dp),allocatable :: data_xyz(:,:)
  real(dp),allocatable :: data_weight(:,:)
  real(dp),allocatable :: node_xyz_2d(:,:,:,:)
  real(dp),allocatable :: gasex_field(:,:) !gasexchange specific fields
  real(dp),allocatable :: unit_field(:,:) !properties of elastic units
  real(dp),allocatable :: node_field(:,:)
  real(dp),allocatable :: scale_factors_2d(:,:)

  logical,allocatable :: expansile(:)

  type capillary_bf_parameters
    integer :: num_symm_gen=9 !no units
    real(dp) :: total_cap_area=0.63000e02_dp !m
    real(dp) :: Palv=0.0_dp!Pa
    real(dp) :: H0=0.35000e-05_dp !m
    real(dp) :: K_cap=0.12000e02_dp
    real(dp) :: F_cap=0.18000e01_dp
    real(dp) :: F_sheet=0.10400e00_dp
    real(dp) :: sigma_cap=0.43637e03_dp !Pa
    real(dp) :: mu_c=0.19200e-02_dp !Pa.s
    real(dp) :: alpha_a=2.33e-08_dp !m/Pa
    real(dp) :: alpha_v=2.33e-08_dp !m/Pa
    real(dp) :: F_rec=0.64630e00_dp
    real(dp) :: sigma_rec=0.22300e04_dp
    real(dp) :: L_c=0.11880e-02_dp !m
    real(dp) :: Plb_c=0.0_dp !Pa
    real(dp) :: Pub_c=3138.24_dp !Pa
    real(dp) :: Pub_a_v=3138.24_dp !Pa
    real(dp) :: L_art_terminal=0.13000e-03_dp !m
    real(dp) :: L_vein_terminal=0.13000e-03_dp !m
    real(dp) :: R_art_terminal=0.10000e-04_dp !m
    real(dp) :: R_vein_terminal=0.90000e-05!m
  end type capillary_bf_parameters

  type admittance_param
    character (len=20) :: admittance_type
    character (len=20) :: bc_type
  end type admittance_param
  type, EXTENDS (admittance_param) :: two_parameter
     real(dp) :: admit_P1=1.0_dp
     real(dp) :: admit_P2=1.0_dp
  end type two_parameter
  type, EXTENDS (two_parameter) :: three_parameter
    real(dp) :: admit_P3=1.0_dp
  end type three_parameter
  type, EXTENDS (three_parameter) :: four_parameter
    real(dp) :: admit_P4=1.0_dp
  end type four_parameter
  type,EXTENDS (four_parameter) :: all_admit_param
  end type all_admit_param

  type elasticity_vessels
    character(len=20) ::vessel_type
  end type elasticity_vessels
  type, EXTENDS(elasticity_vessels) :: elasticity_param
    real(dp) :: elasticity_parameters(3)=0.0_dp
  end type elasticity_param

  type fluid_properties
     real(dp) :: blood_viscosity=0.33600e-02_dp !Pa.s
     real(dp) :: blood_density=0.10500e-02_dp !kg/cm3
     real(dp) :: air_viscosity = 1.8e-5_dp   ! Pa.s
     real(dp) :: air_density = 1.146e-6_dp ! g.mm^-3
  end type fluid_properties
  
  type default_lung_mechanics
     ! default values for Fung exponential, as per Tawhai et al (2009)
     real(dp) :: a = 0.433_dp 
     real(dp) :: b = -0.611_dp
     real(dp) :: c = 2500.0_dp
     real(dp) :: refvol_ratio = 0.5_dp
     real(dp) :: chest_wall_compliance = 2000.0_dp
  end type default_lung_mechanics
  
  type default_lung_volumes
     ! default values for the 'typical' upright lung
     real(dp) :: frc = 3.0e+6_dp  ! frc in mm3
     real(dp) :: Rmax = 0.79_dp   ! ratio of density in non-dependent tissue to mean density 
     real(dp) :: Rmin = 1.29_dp   ! ratio of density in dependent tissue to mean density
     real(dp) :: COV = 0.1_dp     ! coefficient of variation for density
  end type default_lung_volumes

  type default_ventilation
     ! default values for ventilation
     real(dp) :: tidal_volume = 4.0e+5_dp  ! mm^3
     real(dp) :: i_to_e_ratio = 1.0_dp     ! dim.
     real(dp) :: time_breath = 4.0_dp      ! sec
     real(dp) :: P_air_inlet = 0.0_dp      ! Pa
     real(dp) :: P_muscle_estimate = -98.0665_dp * 2.0_dp  ! 2 cmH2O converted to Pa
     real(dp) :: factor_P_muscle_insp = 1.0_dp  ! multiplier to scale inspiratory pressure
     real(dp) :: factor_P_muscle_expn = 1.0_dp  ! multiplier to scale expiratory pressure
     character(len=7) :: expiration_type = 'active' ! or passive
  end type default_ventilation

  type default_ventilation_solver
     ! default values for the iterative solution in ventilation code
     integer :: num_iterations = 200
     real(dp) :: error_tolerance = 1.0e-08_dp
  end type default_ventilation_solver

!!! arrays that start with default values, updated during simulations
  type(default_lung_mechanics) :: lung_mechanics
  type(default_lung_volumes) :: lung_volumes
  type(default_ventilation) :: ventilation_values
  type(default_ventilation_solver) :: ventilation_solver

! temporary, for debugging:
  real(dp) :: unit_before

  private

  public set_node_field_value, elem_field, num_elems, num_elems_2d, elem_nodes, node_xyz, &
         nodes,nodes_2d, elems, num_nodes, num_nodes_2d, num_data, data_xyz, data_weight, &
         node_xyz_2d, node_versn_2d, units, num_units, unit_field, node_field, dp, &
         elem_cnct, elem_ordrs, elem_direction, elems_at_node, elem_symmetry, expansile, &
         elem_units_below, maxgen,capillary_bf_parameters, zero_tol,loose_tol,gasex_field, &
         num_lines_2d, lines_2d, line_versn_2d, lines_in_elem, nodes_in_line, elems_2d, &
         elem_cnct_2d, elem_nodes_2d, elem_versn_2d, elem_lines_2d, elems_at_node_2d, arclength, &
         scale_factors_2d, parentlist, fluid_properties, elasticity_vessels, admittance_param, &
         elasticity_param, all_admit_param, lung_mechanics, lung_volumes, ventilation_values, &
         ventilation_solver, update_parameter

contains
  subroutine set_node_field_value(row, col, value)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SET_NODE_FIELD_VALUE" :: SET_NODE_FIELD_VALUE
    implicit none

    integer, intent(in) :: row, col
    real(dp), intent(in) :: value

    node_field(row, col) = value

  end subroutine set_node_field_value

  subroutine update_parameter(parameter_name, parameter_value)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_UPDATE_PARAMETER" :: UPDATE_PARAMETER
    implicit none
    real(dp), intent(in) :: parameter_value
    character(len=*), intent(in) :: parameter_name

    select case(parameter_name)
       
!!! lung_volumes
    case('COV')
       lung_volumes%COV = parameter_value
    case('FRC')
       lung_volumes%FRC = parameter_value
    case('Rmax')
       lung_volumes%Rmax = parameter_value
    case('Rmin')
       lung_volumes%Rmin = parameter_value

!!! lung_mechanics
    case('chest_wall_compliance')
       lung_mechanics%chest_wall_compliance = parameter_value
    case('mech_a')
       lung_mechanics%a = parameter_value
    case('mech_b')
       lung_mechanics%b = parameter_value
    case('mech_c')
       lung_mechanics%c = parameter_value
    case('refvol_ratio')
       lung_mechanics%refvol_ratio = parameter_value

!!! ventilation_values
    case('i_to_e_ratio')
       ventilation_values%i_to_e_ratio = parameter_value
    case('tidal_volume')
       ventilation_values%tidal_volume = parameter_value
    case('time_breath')
       ventilation_values%time_breath = parameter_value
    case('P_muscle_estimate')
       ventilation_values%P_muscle_estimate = parameter_value
    case('P_air_inlet')
       ventilation_values%P_air_inlet = parameter_value

!!! ventilation_solver
    case('vent_error_tol')
       ventilation_solver%error_tolerance = parameter_value
    case('vent_num_iterations')
       ventilation_solver%num_iterations = int(parameter_value)
       
    end select
    
  end subroutine update_parameter


end module arrays
