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

  integer :: num_elems,num_elems_2d,num_groups,num_nodes,num_data, &
       num_nodes_2d,num_triangles,num_units,num_vertices,num_lines_2d,maxgen

  integer,allocatable :: nodes(:) !allocated in define_node_geometry
  integer,allocatable :: nodes_2d(:) !allocated in define_node_geometry_2d
  integer,allocatable :: node_versn_2d(:) !allocated in define_node_geometry_2d
  integer :: ndata_groups(20,2)
  integer,allocatable :: nelem_groups(:,:)
  integer,allocatable :: elems(:) !allocated in define_1d_elements
  integer,allocatable :: lines_2d(:)
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
  integer,allocatable :: triangle(:,:)
  integer,allocatable :: units(:)

  ! from p-r-f
  integer,allocatable :: mesh_from_depvar(:,:,:)
  integer, allocatable :: depvar_at_node(:,:,:)
  integer, allocatable :: depvar_at_elem(:,:,:)
  integer, allocatable :: SparseCol(:)
  integer, allocatable :: SparseRow(:)
  integer, allocatable :: update_resistance_entries(:)
  real(dp), allocatable :: SparseVal(:)
  real(dp), allocatable :: RHS(:)
  real(dp), allocatable :: prq_solution(:,:),solver_solution(:)
  logical, allocatable :: FIX(:)
  
  real(dp),allocatable :: arclength(:)
  real(dp),allocatable :: elem_field(:,:) !properties of elements
  real(dp),allocatable :: elem_direction(:,:)
  real(dp),allocatable :: node_xyz(:,:)
  real(dp),allocatable :: data_field(:,:)
  real(dp),allocatable :: data_xyz(:,:)
  real(dp),allocatable :: data_weight(:,:)
  real(dp),allocatable :: node_xyz_2d(:,:,:,:)
  real(dp),allocatable :: gasex_field(:,:) !gasexchange specific fields
  real(dp),allocatable :: unit_field(:,:) !properties of elastic units
  real(dp),allocatable :: vertex_xyz(:,:)
  real(dp),allocatable :: node_field(:,:)
  real(dp),allocatable :: scale_factors_2d(:,:)

  character(len=20),dimension(20) :: data_group_names,elem_group_names
  
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
     real(dp) :: blood_viscosity = 0.33600e-02_dp ! Pa.s
     real(dp) :: blood_density = 0.10500e-02_dp   ! kg/cm3
     real(dp) :: air_viscosity = 1.8e-5_dp        ! Pa.s
     real(dp) :: air_density = 1.146e-6_dp        ! g.mm^-3
  end type fluid_properties
  
! temporary, for debugging:
  real(dp) :: unit_before

  private

  public set_node_field_value, elem_field, num_elems, num_elems_2d, num_groups, elem_nodes, node_xyz, &
       nodes,nodes_2d, elems, num_nodes, num_nodes_2d, num_data, num_triangles, num_vertices, &
       data_field, data_xyz, data_weight, &
       node_xyz_2d, node_versn_2d, units, num_units, unit_field, node_field, dp, &
       data_group_names, elem_group_names, ndata_groups, nelem_groups, &
       elem_cnct, elem_ordrs, elem_direction, elems_at_node, elem_symmetry, expansile, &
       elem_units_below, maxgen,capillary_bf_parameters, zero_tol,loose_tol,gasex_field, &
       num_lines_2d, lines_2d, line_versn_2d, lines_in_elem, nodes_in_line, elems_2d, &
       elem_cnct_2d, elem_nodes_2d, elem_versn_2d, elem_lines_2d, elems_at_node_2d, arclength, &
       scale_factors_2d, fluid_properties, elasticity_vessels, admittance_param, &
       elasticity_param, two_parameter, three_parameter, four_parameter, all_admit_param, &
       mesh_from_depvar, depvar_at_node, depvar_at_elem, SparseCol, SparseRow, triangle, &
       update_resistance_entries, vertex_xyz, &
       SparseVal, RHS, prq_solution, solver_solution, FIX

contains
  subroutine set_node_field_value(row, col, value)
    implicit none

    integer, intent(in) :: row, col
    real(dp), intent(in) :: value

    node_field(row, col) = value

  end subroutine set_node_field_value

end module arrays
