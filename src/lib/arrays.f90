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

  implicit none

  integer :: num_elems,num_elems_2d,num_nodes,num_data,num_nodes_2d,num_units,num_lines_2d,maxgen

  integer, parameter :: dp=kind(0.d0) !  for double precision
  real(dp),parameter :: zero_tol = 1.0e-12_dp
  real(dp),parameter :: loose_tol = 1.0e-6_dp

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
    real(dp) :: total_cap_area=0.63000e+02_dp !m
    real(dp) :: Palv=0.0_dp!Pa
    real(dp) :: H0=0.35000e-05_dp !m
    real(dp) :: K_cap=0.12000e+02_dp
    real(dp) :: F_cap=0.18000e+01_dp
    real(dp) :: F_sheet=0.10400e+00_dp
    real(dp) :: sigma_cap=0.43637e+03_dp !Pa
    real(dp) :: mu_c=0.19200e-02_dp !Pa.s
    real(dp) :: alpha_a=2.33e-08_dp !/Pa
    real(dp) :: alpha_v=2.33e-08_dp !/Pa
    real(dp) :: F_rec=0.64630e+00_dp
    real(dp) :: sigma_rec=0.22300e+04_dp
    real(dp) :: L_c=0.11880e-02_dp !m
    real(dp) :: Plb_c=0.0_dp !Pa
    real(dp) :: Pub_c=3138.24_dp !Pa
    real(dp) :: Pub_a_v=3138.24_dp !Pa
    real(dp) :: L_art_terminal=0.13000e-03_dp !m
    real(dp) :: L_vein_terminal=0.13000e-03_dp !m
    real(dp) :: R_art_terminal=0.10000e-04_dp !m
    real(dp) :: R_vein_terminal=0.90000e-05!m
  end type capillary_bf_parameters

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
         scale_factors_2d, parentlist

contains
  subroutine set_node_field_value(row, col, value)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SET_NODE_FIELD_VALUE" :: SET_NODE_FIELD_VALUE
    implicit none

    integer, intent(in) :: row, col
    real(dp), intent(in) :: value

    node_field(row, col) = value

  end subroutine set_node_field_value


end module arrays
