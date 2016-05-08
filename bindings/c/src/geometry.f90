module geometry_c
  implicit none
  private

contains
!
!###################################################################################
!
!*add_mesh:* Reads in an ipmesh file and adds this mesh to the terminal branches of an existing tree geometry
  subroutine add_mesh_c(AIRWAY_MESHFILE, filename_len) bind(C, name="add_mesh_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_ADD_MESH_C" :: ADD_MESH_C
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: add_mesh
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: AIRWAY_MESHFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, AIRWAY_MESHFILE, filename_len)
    call so_add_mesh(filename_f)

  end subroutine add_mesh_c
!
!###################################################################################
!
!*append_units:* Appends terminal units at the end of a tree structure
  subroutine append_units_c() bind(C, name="append_units_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_APPEND_UNITS_c" :: APPEND_UNITS_C
    use geometry, only: append_units
    implicit none

    call so_append_units

  end subroutine append_units_c

!
!###################################################################################
!
  subroutine define_1d_elements_c(ELEMFILE, filename_len) bind(C, name="define_1d_elements_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_DEFINE_1D_ELEMENTS" :: DEFINE_1D_ELEMENTS

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use geometry, only: define_1d_elements
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: ELEMFILE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, ELEMFILE, filename_len)
    call so_define_1d_elements(filename_f)

  end subroutine define_1d_elements_c
!
!###################################################################################
!
!*define_mesh_geometry_test:*
  subroutine define_mesh_geometry_test_c() bind(C, name="define_mesh_geometry_test_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_DEFINE_MESH_GEOMETRY_TEST_C" :: DEFINE_MESH_GEOMETRY_TEST_C
    use geometry, only: define_mesh_geometry_test
    implicit none

    call define_mesh_geometry_test

  end subroutine define_mesh_geometry_test_c
!
!###################################################################################
!
  subroutine define_node_geometry_c(NODEFILE, filename_len) bind(C, name="define_node_geometry_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_DEFINE_NODE_GEOMETRY_C" :: DEFINE_NODE_GEOMETRY_C

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
  subroutine define_rad_from_file_c(FIELDFILE, filename_len, radius_type, radius_type_len) bind(C, name="define_rad_from_file_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_DEFINE_RAD_FROM_FILE_C" :: DEFINE_RAD_FROM_FILE_C

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
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_DEFINE_RAD_FROM_GEOM_C" :: DEFINE_RAD_FROM_GEOM_C

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_STRING_LEN, dp
    use geometry, only: define_rad_from_geom
    implicit none

    real(dp),intent(in) :: control_param, start_rad
    integer,intent(in) :: order_system_len, start_from_len, group_type_len, group_options_len
    type(c_ptr), value, intent(in) :: order_system, start_from, group_type, group_options
    character(len=MAX_STRING_LEN) :: order_system_f, start_from_f, group_type_f, group_options_f

    call strncpy(order_system_f, order_system, order_system_len)
    call strncpy(start_from_f, start_from, start_from_len)
    call strncpy(group_options_f, group_options, group_options_len)
    call strncpy(group_options_f, group_type, group_type_len)
    call define_rad_from_geom(order_system_f, control_param, start_from_f, start_rad, group_type_f, group_options_f)

  end subroutine define_rad_from_geom_c
!
!###########################################################################
!
!*element_connectivity_1d:*  Calculates element connectivity in 1D and stores in elelem_cnct
  subroutine element_connectivity_1d_c() bind(C, name="element_connectivity_1d_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_ELEMENT_CONNECTIVITY_1D_C" :: ELEMENT_CONNECTIVITY_1D_C
    use geometry, only: element_connectivity_1d
    implicit none

    call element_connectivity_1d

  end subroutine element_connectivity_1d_c

!
!###################################################################################
!
!*evaluate_ordering:* calculates generations, Horsfield orders, Strahler orders for a given tree
  subroutine evaluate_ordering_c() bind(C, name="evaluate_ordering_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_EVALUATE_ORDERING_C" :: EVALUATE_ORDERING_C
    use geometry, only: evaluate_ordering
    implicit none

    call evaluate_ordering

  end subroutine evaluate_ordering_c
!
!###################################################################################
!
!>*set_initial_volume:* assigns a volume to terminal units appended on a tree structure
!>based on an assumption of a linear gradient in the gravitational direction with max
!> min and COV values defined.
  subroutine set_initial_volume_c(Gdirn, COV, total_volume, Rmax, Rmin) bind(C, name="set_initial_volume_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_SET_INITIAL_VOLUME_C" :: SET_INITIAL_VOLUME_C

    use geometry, only: set_initial_volume
    use other_consts, only: dp
    implicit none

    !     Parameter List
    integer,intent(in) :: Gdirn
    real(dp),intent(in) :: COV, total_volume, Rmax, Rmin

    call set_initial_volume(Gdirn, COV, total_volume, Rmax, Rmin)

  end subroutine set_initial_volume_c

!
!###################################################################################
!
!*volume_of_mesh:* calculates the volume of an airway mesh including conducting and respiratory airways
  subroutine volume_of_mesh_c(volume_model,volume_tree) bind(C, name="volume_of_mesh_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_VOLUME_OF_MESH_C" :: VOLUME_OF_MESH_C
    use other_consts,only: dp
    use geometry, only: volume_of_mesh
    implicit none

    real(dp) :: volume_model,volume_tree

    call volume_of_mesh(volume_model, volume_tree)

  end subroutine volume_of_mesh_c

end module geometry_c

