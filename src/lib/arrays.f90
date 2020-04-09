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
  real(dp),allocatable :: field_xyz(:)
  real(dp),allocatable :: data_field(:)

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
    real(dp) :: air_viscosity
    real(dp) :: air_density
  end type fluid_properties

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
         elasticity_param, all_admit_param, field_xyz, data_field

contains
  subroutine set_node_field_value(row, col, value)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SET_NODE_FIELD_VALUE" :: SET_NODE_FIELD_VALUE
    implicit none

    integer, intent(in) :: row, col
    real(dp), intent(in) :: value

    node_field(row, col) = value

  end subroutine set_node_field_value


end module arrays
