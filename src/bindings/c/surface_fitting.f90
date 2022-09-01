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

!!!############################################################################

  subroutine initialise_fit_mesh_c() bind(C, name="initialise_fit_mesh_c")
    use surface_fitting, only: initialise_fit_mesh
    implicit none
    
    call initialise_fit_mesh

  end subroutine initialise_fit_mesh_c
    
!!!############################################################################
end module surface_fitting_c
