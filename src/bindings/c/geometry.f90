module geometry_c
  use arrays
  use diagnostics
  use indices
  !use mesh_functions
  !use precision ! sets dp for precision
  !use math_constants !pi  

implicit none
  private

contains
!
!###################################################################################
!
!*add_mesh:* Reads in an ipmesh file and adds this mesh to the terminal branches of an existing tree geometry
  subroutine add_mesh_c(AIRWAY_MESHFILE, filename_len) bind(C, name="add_mesh_c")
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: add_mesh
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: AIRWAY_MESHFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, AIRWAY_MESHFILE, filename_len)

    call add_mesh(filename_f)

  end subroutine add_mesh_c
!
!###################################################################################
!
!*add_matching_mesh:* Replicates an existing mesh, continuing node and element numbers
  subroutine add_matching_mesh_c() bind(C, name="add_matching_mesh_c")
    use geometry, only: add_matching_mesh
    implicit none

    call add_matching_mesh

  end subroutine add_matching_mesh_c

!
!###################################################################################
!
!*append_units:* Appends terminal units at the end of a tree structure
  subroutine append_units_c() bind(C, name="append_units_c")
    use geometry, only: append_units
    implicit none

    call append_units

  end subroutine append_units_c

!
!###################################################################################
!
  subroutine define_1d_elements_c(ELEMFILE, filename_len) bind(C, name="define_1d_elements_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: define_1d_elements
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: ELEMFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, ELEMFILE, filename_len)

    call define_1d_elements(filename_f)

  end subroutine define_1d_elements_c

!
!###################################################################################
!
  subroutine define_elem_geometry_2d_c(ELEMFILE, filename_len, SF_OPTION, sf_option_len) bind(C, name="define_elem_geometry_2d_c")
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: define_elem_geometry_2d
    implicit none

    integer,intent(in) :: filename_len, sf_option_len
    type(c_ptr), value, intent(in) :: ELEMFILE, SF_OPTION
    character(len=MAX_FILENAME_LEN) :: filename_f, sf_option_f

    call strncpy(filename_f, ELEMFILE, filename_len)
    call strncpy(sf_option_f, SF_OPTION, sf_option_len)

    call define_elem_geometry_2d(filename_f,sf_option_f)

  end subroutine define_elem_geometry_2d_c
!
!###################################################################################
!
!*define_mesh_geometry_test:*
  subroutine define_mesh_geometry_test_c() bind(C, name="define_mesh_geometry_test_c")
    use geometry, only: define_mesh_geometry_test
    implicit none

    call define_mesh_geometry_test

  end subroutine define_mesh_geometry_test_c
!
!###################################################################################
!
  subroutine define_node_geometry_c(NODEFILE, filename_len) bind(C, name="define_node_geometry_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: define_node_geometry
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: NODEFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, NODEFILE, filename_len)

    call define_node_geometry(filename_f)

  end subroutine define_node_geometry_c

!
!###################################################################################
!
  subroutine import_node_geometry_2d_c(NODEFILE, filename_len) bind(C, name="import_node_geometry_2d_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: import_node_geometry_2d
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: NODEFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, NODEFILE, filename_len)

    call import_node_geometry_2d(filename_f)

  end subroutine import_node_geometry_2d_c

!
!###################################################################################
!
  subroutine import_ply_triangles_c(ply_file, filename_len) bind(C, name="import_ply_triangles_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: import_ply_triangles
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: ply_file
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, ply_file, filename_len)

    call import_ply_triangles(filename_f)

  end subroutine import_ply_triangles_c

!
!###################################################################################
!
  subroutine make_data_grid_c(surface_elems_len, surface_elems, offset, spacing, &
       filename, filename_len, groupname, groupname_len)&
 bind(C, name="make_data_grid_c")
    
    use arrays,only: dp
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
    use geometry, only: make_data_grid
    implicit none

    integer,intent(in) :: surface_elems_len
    integer,intent(in) :: surface_elems(surface_elems_len)
    real(dp),intent(in) :: offset, spacing
    integer,intent(in) :: filename_len, groupname_len
    type(c_ptr), value, intent(in) :: filename, groupname
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: groupname_f

    call strncpy(filename_f, filename, filename_len)
    call strncpy(groupname_f, groupname, groupname_len)

    call make_data_grid(surface_elems, offset, spacing, filename_f, groupname_f)

  end subroutine make_data_grid_c

!!!###################################################################################

  subroutine make_2d_vessel_from_1d_c(elemlist, elemlist_len) bind(C, name="make_2d_vessel_from_1d_c")
    use geometry, only: make_2d_vessel_from_1d
    implicit none

    integer,intent(in) :: elemlist_len
    integer,intent(in) :: elemlist(elemlist_len)

    call make_2d_vessel_from_1d(elemlist)

  end subroutine make_2d_vessel_from_1d_c
  
!
!###################################################################################
!
  subroutine define_data_geometry_c(DATAFILE, filename_len) bind(C, name="define_data_geometry_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: define_data_geometry
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: DATAFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, DATAFILE, filename_len)

    call define_data_geometry(filename_f)

  end subroutine define_data_geometry_c

!
!###################################################################################
!
  subroutine define_node_geometry_2d_c(NODEFILE, filename_len) bind(C, name="define_node_geometry_2d_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: define_node_geometry_2d
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: NODEFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, NODEFILE, filename_len)

    call define_node_geometry_2d(filename_f)

  end subroutine define_node_geometry_2d_c

!
!###################################################################################
!

  subroutine define_rad_from_file_c(FIELDFILE, filename_len, radius_type, radius_type_len) bind(C, name="define_rad_from_file_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
    use geometry, only: define_rad_from_file
    implicit none

    integer,intent(in) :: filename_len, radius_type_len
    type(c_ptr), value, intent(in) :: FIELDFILE, radius_type
    character(len=MAX_FILENAME_LEN) :: filename_f, radius_type_f

    call strncpy(filename_f, FIELDFILE, filename_len)
    call strncpy(radius_type_f, radius_type, radius_type_len)

    call define_rad_from_file(filename_f, radius_type_f)

    end subroutine define_rad_from_file_c
!
!##################################################################################
!
!*define_rad_from_geom:* Defines vessel or airway radius based on their geometric structure
  subroutine define_rad_from_geom_c(order_system, order_system_len, control_param, &
        start_from, start_from_len, start_rad, group_type, group_type_len, group_options, group_options_len) &
        bind(C, name="define_rad_from_geom_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_STRING_LEN
    use arrays, only: dp
    use geometry, only: define_rad_from_geom
    implicit none

    real(dp),intent(in) :: control_param, start_rad
    integer,intent(in) :: order_system_len, start_from_len, group_type_len, group_options_len
    type(c_ptr), value, intent(in) :: order_system, start_from, group_type, group_options
    character(len=MAX_STRING_LEN) :: order_system_f, start_from_f, group_type_f, group_options_f

    call strncpy(order_system_f, order_system, order_system_len)
    call strncpy(start_from_f, start_from, start_from_len)
    call strncpy(group_options_f, group_options, group_options_len)
    call strncpy(group_type_f, group_type, group_type_len)

    call define_rad_from_geom(order_system_f, control_param, start_from_f, start_rad, group_type_f, group_options_f)

  end subroutine define_rad_from_geom_c
!
!###########################################################################
!
!*element_connectivity_1d:*  Calculates element connectivity in 1D and stores in elelem_cnct
  subroutine element_connectivity_1d_c() bind(C, name="element_connectivity_1d_c")
    use geometry, only: element_connectivity_1d
    implicit none

    call element_connectivity_1d

  end subroutine element_connectivity_1d_c

!
!###################################################################################
!
!*evaluate_ordering:* calculates generations, Horsfield orders, Strahler orders for a given tree
  subroutine evaluate_ordering_c() bind(C, name="evaluate_ordering_c")
    use geometry, only: evaluate_ordering
    implicit none

    call evaluate_ordering

  end subroutine evaluate_ordering_c

!
!###################################################################################
!
!*volume_of_mesh:* calculates the volume of an airway mesh including conducting and respiratory airways
  subroutine volume_of_mesh_c(volume_model,volume_tree) bind(C, name="volume_of_mesh_c")
    use arrays, only: dp
    use geometry, only: volume_of_mesh
    implicit none

    real(dp) :: volume_model,volume_tree

    call volume_of_mesh(volume_model, volume_tree)

  end subroutine volume_of_mesh_c


  function get_local_node_f_c(ndimension,np_global) result(get_local_node) bind(C, name="get_local_node_f_c")
    use arrays, only: dp
    use geometry, only: get_local_node_f
    implicit none
    
    integer :: ndimension,np_global
    integer :: get_local_node
    
    get_local_node=get_local_node_f(ndimension,np_global)

  end function get_local_node_f_c

!
!###################################################################################
!
  subroutine write_elem_geometry_2d_c(elemfile, filename_len) bind(C, name="write_elem_geometry_2d_c")
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: write_elem_geometry_2d
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: elemfile
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, elemfile, filename_len)

    call write_elem_geometry_2d(filename_f)

  end subroutine write_elem_geometry_2d_c
!
!###################################################################################
!
  subroutine write_geo_file_c(ntype, geofile, filename_len) bind(C, name="write_geo_file_c")
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: write_geo_file
    implicit none

    integer,intent(in) :: ntype, filename_len
    type(c_ptr), value, intent(in) :: geofile
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, geofile, filename_len)

    call write_geo_file(ntype, filename_f)

  end subroutine write_geo_file_c
!
!###################################################################################
!
  subroutine write_node_geometry_2d_c(nodefile, filename_len) bind(C, name="write_node_geometry_2d_c")
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: write_node_geometry_2d
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: nodefile
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, nodefile, filename_len)

    call write_node_geometry_2d(filename_f)

  end subroutine write_node_geometry_2d_c
!
!###################################################################################
!

end module geometry_c


