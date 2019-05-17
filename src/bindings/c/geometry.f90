module geometry_c
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
#if defined _WIN32 && defined __INTEL_COMPILER
    call so_add_mesh(filename_f)
#else
    call add_mesh(filename_f)
#endif

  end subroutine add_mesh_c
!
!###################################################################################
!
!*add_matching_mesh:* Replicates an existing mesh, continuing node and element numbers
  subroutine add_matching_mesh_c() bind(C, name="add_matching_mesh_c")
    use geometry, only: add_matching_mesh
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_add_matching_mesh
#else
    call add_matching_mesh
#endif

  end subroutine add_matching_mesh_c

!
!###################################################################################
!
!*append_units:* Appends terminal units at the end of a tree structure
  subroutine append_units_c() bind(C, name="append_units_c")
    use geometry, only: append_units
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_append_units
#else
    call append_units
#endif

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

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_define_1d_elements(filename_f)
#else
    call define_1d_elements(filename_f)
#endif

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

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_define_elem_geometry_2d(filename_f,sf_option_f)
#else
    call define_elem_geometry_2d(filename_f,sf_option_f)
#endif

  end subroutine define_elem_geometry_2d_c
!
!###################################################################################
!
!*define_mesh_geometry_test:*
  subroutine define_mesh_geometry_test_c() bind(C, name="define_mesh_geometry_test_c")
    use geometry, only: define_mesh_geometry_test
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_define_mesh_geometry_test
#else
    call define_mesh_geometry_test
#endif

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

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_define_node_geometry(filename_f)
#else
    call define_node_geometry(filename_f)
#endif

  end subroutine define_node_geometry_c

!
!###################################################################################
!
  subroutine make_data_grid_c(surface_elems, spacing, to_export, filename, filename_len, groupname, groupname_len)&
 bind(C, name="make_data_grid_c")
    
    use arrays,only: dp
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
    use geometry, only: make_data_grid
    implicit none

    integer,intent(in) :: surface_elems(:)
    real(dp),intent(in) :: spacing
    logical,intent(in) :: to_export
    integer,intent(in) :: filename_len, groupname_len
    type(c_ptr), value, intent(in) :: filename, groupname
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: groupname_f

    call strncpy(filename_f, filename, filename_len)
    call strncpy(groupname_f, groupname, groupname_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_make_data_grid(surface_elems, spacing, to_export, filename_f, groupname_f)
#else
    call make_data_grid(surface_elems, spacing, to_export, filename_f, groupname_f)
#endif

  end subroutine make_data_grid_c

!
!###################################################################################
!
  subroutine group_elem_parent_term_c(ne_parent) bind(C, name="group_elem_parent_term_c")

    use iso_c_binding, only: c_ptr
    use geometry, only: group_elem_parent_term
    implicit none

    integer,intent(in) :: ne_parent

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_group_elem_parent_term(ne_parent)
#else
    call group_elem_parent_term(ne_parent)
#endif

  end subroutine group_elem_parent_term_c

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

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_define_data_geometry(filename_f)
#else
    call define_data_geometry(filename_f)
#endif

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

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_define_node_geometry_2d(filename_f)
#else
    call define_node_geometry_2d(filename_f)
#endif

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

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_define_rad_from_file(filename_f, radius_type_f)
#else
    call define_rad_from_file(filename_f, radius_type_f)
#endif

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

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_define_rad_from_geom(order_system_f, control_param, start_from_f, start_rad, group_type_f, group_options_f)
#else
    call define_rad_from_geom(order_system_f, control_param, start_from_f, start_rad, group_type_f, group_options_f)
#endif

  end subroutine define_rad_from_geom_c
!
!###########################################################################
!
!*element_connectivity_1d:*  Calculates element connectivity in 1D and stores in elelem_cnct
  subroutine element_connectivity_1d_c() bind(C, name="element_connectivity_1d_c")
    use geometry, only: element_connectivity_1d
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_element_connectivity_1d
#else
    call element_connectivity_1d
#endif

  end subroutine element_connectivity_1d_c

!
!###################################################################################
!
!*evaluate_ordering:* calculates generations, Horsfield orders, Strahler orders for a given tree
  subroutine evaluate_ordering_c() bind(C, name="evaluate_ordering_c")
    use geometry, only: evaluate_ordering
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_evaluate_ordering
#else
    call evaluate_ordering
#endif

  end subroutine evaluate_ordering_c
!
!###################################################################################
!
!>*set_initial_volume:* assigns a volume to terminal units appended on a tree structure
!>based on an assumption of a linear gradient in the gravitational direction with max
!> min and COV values defined.
  subroutine set_initial_volume_c(Gdirn, COV, total_volume, Rmax, Rmin) bind(C, name="set_initial_volume_c")

    use geometry, only: set_initial_volume
    use arrays, only: dp
    implicit none

    !     Parameter List
    integer,intent(in) :: Gdirn
    real(dp),intent(in) :: COV, total_volume, Rmax, Rmin

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_set_initial_volume(Gdirn, COV, total_volume, Rmax, Rmin)
#else
    call set_initial_volume(Gdirn, COV, total_volume, Rmax, Rmin)
#endif

  end subroutine set_initial_volume_c

!
!###################################################################################
!
!*volume_of_mesh:* calculates the volume of an airway mesh including conducting and respiratory airways
  subroutine volume_of_mesh_c(volume_model,volume_tree) bind(C, name="volume_of_mesh_c")
    use arrays, only: dp
    use geometry, only: volume_of_mesh
    implicit none

    real(dp) :: volume_model,volume_tree

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_volume_of_mesh(volume_model, volume_tree)
#else
    call volume_of_mesh(volume_model, volume_tree)
#endif

  end subroutine volume_of_mesh_c

!
!###################################################################################
!
  function get_local_node_f_c(ndimension,np_global) result(get_local_node) bind(C, name="get_local_node_f_c")
    use arrays, only: dp
    use geometry, only: get_local_node_f
    implicit none
    
    integer :: ndimension,np_global
    integer :: get_local_node
    
    get_local_node=get_local_node_f(ndimension,np_global)

  end function get_local_node_f_c

end module geometry_c

