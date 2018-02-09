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

  integer :: num_elems,num_nodes,num_units,maxgen

  integer, parameter :: dp=kind(0.d0) !  for double precision

  integer,allocatable :: nodes(:) !allocated in define_node_geometry, the number index of each node, the node number may not be continuous
  integer,allocatable :: elems(:) !allocated in define_1d_elements, the number of each element, the element number may not be continuous
  !! Assume each elem has two nodes np1, np2, for a random elem ne, elem_cnct(-1,0,ne) represents the number of elements which connect element "ne" through node "np1",
  !! elem_cnct(-1,1:2,ne) represents the element number which connect this element "ne" through node "np1", similarly, elem_cnct(1,0:2,ne) represents the connection
  !! between this element "ne" and other element through node "np2".
  !! elem_cnct(-1,0:2,ne) represents the parent element connections, elem_cnct(1,0:2,ne) represents the children element connections.
  integer,allocatable :: elem_cnct(:,:,:)  !NXI(-ni:ni,1,ne)
  integer,allocatable :: elem_nodes(:,:) !the node numbers for each element, the size is elem_nodes(2,num_elems)
  integer,allocatable :: elem_ordrs(:,:)
  integer,allocatable :: elem_symmetry(:)
  integer,allocatable :: elem_units_below(:) ! the effective number of elements below each element
  integer,allocatable :: elems_at_node(:,:) !! For each line of elems_at_node(number_of_nodes,0:3), index 0 represents the total number of element this node connects, index 1-3 represnts the element number this node connects
  integer,allocatable :: units(:) ! the element number index of each unit element

  real(dp),allocatable :: elem_field(:,:) !properties of elements
  real(dp),allocatable :: elem_direction(:,:) !the sin value of each dircion, the size is elem_direction(3,num_elems)
  real(dp),allocatable :: node_xyz(:,:) ! the node coords, size is node_xyz(3,num_nodes)
  real(dp),allocatable :: unit_field(:,:) !properties of elastic units
  real(dp),allocatable :: node_field(:,:)
  real(dp),allocatable :: node_label(:) !the label index of each node
  real(dp),allocatable :: elem_label(:) !the label index of each element
  real(dp),allocatable :: unit_elem_label(:) !the label index of each terminal element unit

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
  public set_node_field_value, elem_field, num_elems, elem_nodes, node_xyz, nodes, elems, &
    num_nodes, units, num_units, unit_field, node_field, dp, elem_cnct, elem_ordrs, elem_direction, &
    elems_at_node, elem_symmetry, expansile, elem_units_below, maxgen,capillary_bf_parameters,node_label,unit_elem_label, &
    elem_label

contains
  subroutine set_node_field_value(row, col, value)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SET_NODE_FIELD_VALUE" :: SET_NODE_FIELD_VALUE
    implicit none

    integer, intent(in) :: row, col
    real(dp), intent(in) :: value

    node_field(row, col) = value

  end subroutine set_node_field_value


end module arrays
