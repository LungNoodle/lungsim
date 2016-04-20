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

  integer,allocatable :: nodes(:) !allocated in define_node_geometry
  integer,allocatable :: elems(:) !allocated in define_1d_elements
  integer,allocatable :: elem_cnct(:,:,:)  !NXI(-ni:ni,1,ne)
  integer,allocatable :: elem_nodes(:,:)
  integer,allocatable :: elem_ordrs(:,:)
  integer,allocatable :: elem_symmetry(:)
  integer,allocatable :: elem_units_below(:)
  integer,allocatable :: elems_at_node(:,:)
  integer,allocatable :: units(:)

  real(dp),allocatable :: elem_field(:,:) !properties of elements
  real(dp),allocatable :: elem_direction(:,:)
  real(dp),allocatable :: node_xyz(:,:)
  real(dp),allocatable :: unit_field(:,:) !properties of elastic units
  real(dp),allocatable :: node_field(:,:)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_NODE_FIELD" :: NODE_FIELD

  logical,allocatable :: expansile(:)

! temporary, for debugging:
  real(dp) :: unit_before

end module arrays
