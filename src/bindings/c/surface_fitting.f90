module surface_fitting_c
  implicit none

  private

contains
!!!################################################################

  subroutine fit_surface_geometry_c(niterations, fitting_file, filename_len) &
    bind(C, name="fit_surface_geometry_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    use surface_fitting, only: fit_surface_geometry
    implicit none

    integer,intent(in) :: niterations,filename_len
    type(c_ptr), value, intent(in) :: fitting_file
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, fitting_file, filename_len)

    call fit_surface_geometry(niterations, filename_f)

  end subroutine fit_surface_geometry_c

!!!################################################################

  subroutine reset_fitting_c() &
    bind(C, name="reset_fitting_c")

    use surface_fitting, only: reset_fitting
    implicit none

    call reset_fitting()

  end subroutine reset_fitting_c

!!!############################################################################

  subroutine define_data_fit_group_c(datafile, datafile_len, groupname, groupname_len) &
       bind(C, name="define_data_fit_group_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    use surface_fitting, only: define_data_fit_group
    implicit none
    
    integer,intent(in) :: datafile_len, groupname_len
    type(c_ptr), value, intent(in) :: datafile, groupname
    character(len=MAX_FILENAME_LEN) :: datafile_f, groupname_f

    call strncpy(datafile_f, datafile, datafile_len)
    call strncpy(groupname_f, groupname, groupname_len)

    call define_data_fit_group(datafile_f, groupname_f)

  end subroutine define_data_fit_group_c
    
!!!############################################################################

  subroutine initialise_fit_mesh_c() bind(C, name="initialise_fit_mesh_c")
    use surface_fitting, only: initialise_fit_mesh
    implicit none
    
    call initialise_fit_mesh

  end subroutine initialise_fit_mesh_c
    
!!!############################################################################
end module surface_fitting_c
