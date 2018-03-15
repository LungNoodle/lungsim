module exports_c
  implicit none

  private

contains
!!!################################################################

  subroutine export_1d_elem_field_c(ne_field, EXELEMFILE, filename_len, group_name, group_name_len, field_name, field_name_len) &
    bind(C, name="export_1d_elem_field_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use exports, only: export_1d_elem_field
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: ne_field, filename_len, group_name_len, field_name_len
    type(c_ptr), value, intent(in) :: EXELEMFILE, group_name, field_name
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: group_name_f, field_name_f

    call strncpy(filename_f, EXELEMFILE, filename_len)
    call strncpy(group_name_f, group_name, group_name_len)
    call strncpy(field_name_f, field_name, field_name_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_export_1d_elem_field(ne_field, filename_f, group_name_f, field_name_f)
#else
    call export_1d_elem_field(ne_field, filename_f, group_name_f, field_name_f)
#endif

  end subroutine export_1d_elem_field_c

!!!############################################################################

  subroutine export_1d_elem_geometry_c(EXELEMFILE, filename_len, name, name_len) bind(C, name="export_1d_elem_geometry_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use exports, only: export_1d_elem_geometry
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: filename_len, name_len
    type(c_ptr), value, intent(in) :: EXELEMFILE, name
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: name_f

    call strncpy(filename_f, EXELEMFILE, filename_len)
    call strncpy(name_f, name, name_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_export_1d_elem_geometry(filename_f, name_f)
#else
    call export_1d_elem_geometry(filename_f, name_f)
#endif

  end subroutine export_1d_elem_geometry_c


!!!##########################################################################

  subroutine export_node_geometry_c(EXNODEFILE, filename_len, name, name_len) bind(C, name="export_node_geometry_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    use exports, only: export_node_geometry
    implicit none

    integer,intent(in) :: filename_len, name_len
    type(c_ptr), value, intent(in) :: EXNODEFILE, name
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: name_f

    call strncpy(filename_f, EXNODEFILE, filename_len)
    call strncpy(name_f, name, name_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_export_node_geometry(filename_f, name_f)
#else
    call export_node_geometry(filename_f, name_f)
#endif

  end subroutine export_node_geometry_c


!!!########################################################################

  subroutine export_terminal_solution_c(EXNODEFILE, filename_len, name, name_len) bind(C, name="export_terminal_solution_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use exports, only: export_terminal_solution
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: filename_len, name_len
    type(c_ptr), value, intent(in) :: EXNODEFILE, name
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: name_f

    call strncpy(filename_f, EXNODEFILE, filename_len)
    call strncpy(name_f, name, name_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_export_terminal_solution(filename_f, name_f)
#else
    call export_terminal_solution(filename_f, name_f)
#endif

  end subroutine export_terminal_solution_c

  !!!########################################################################

  subroutine export_terminal_perfusion_c(EXNODEFILE, filename_len, name, name_len) bind(C, name="export_terminal_perfusion_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use exports, only: export_terminal_perfusion
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: filename_len, name_len
    type(c_ptr), value, intent(in) :: EXNODEFILE, name
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: name_f

    call strncpy(filename_f, EXNODEFILE, filename_len)
    call strncpy(name_f, name, name_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_export_terminal_perfusion(filename_f, name_f)
#else
    call export_terminal_perfusion(filename_f, name_f)
#endif

  end subroutine export_terminal_perfusion_c

   !!!########################################################################

  subroutine export_terminal_ssgexch_c(EXNODEFILE, filename_len, name, name_len) bind(C, name="export_terminal_ssgexch_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use exports, only: export_terminal_ssgexch
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: filename_len, name_len
    type(c_ptr), value, intent(in) :: EXNODEFILE, name
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: name_f

    call strncpy(filename_f, EXNODEFILE, filename_len)
    call strncpy(name_f, name, name_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_export_terminal_ssgexch(filename_f, name_f)
#else
    call export_terminal_ssgexch(filename_f, name_f)
#endif

  end subroutine export_terminal_ssgexch_c


!!! #################################################################

  subroutine export_node_field_c(nj_field, EXNODEFIELD, filename_len, name, name_len, field_name, field_name_len) &
    bind(C, name="export_node_field_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use exports, only: export_node_field
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: nj_field, filename_len, name_len, field_name_len
    type(c_ptr), value, intent(in) :: EXNODEFIELD, name, field_name
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: name_f, field_name_f

    call strncpy(filename_f, EXNODEFIELD, filename_len)
    call strncpy(name_f, name, name_len)
    call strncpy(field_name_f, field_name, field_name_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_export_node_field(nj_field, filename_f, name_f, field_name_f)
#else
    call export_node_field(nj_field, filename_f, name_f, field_name_f)
#endif

  end subroutine export_node_field_c


!!! ###########################################################

  subroutine export_elem_field_c(EXELEMFIELD, filename_len, name, name_len, field_name, field_name_len) &
    bind(C, name="export_elem_field_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use exports, only: export_elem_field
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none

    integer,intent(in) :: filename_len, name_len, field_name_len
    type(c_ptr), value, intent(in) :: EXELEMFIELD, name, field_name
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: name_f, field_name_f

    call strncpy(filename_f, EXELEMFIELD, filename_len)
    call strncpy(name_f, name, name_len)
    call strncpy(field_name_f, field_name, field_name_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_export_elem_field(filename_f, name_f, field_name_f)
#else
    call export_elem_field(filename_f, name_f, field_name_f)
#endif

  end subroutine export_elem_field_c

end module exports_c
