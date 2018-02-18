
module indices_c
  implicit none
!Interfaces
private

contains
  subroutine define_problem_type_c(PROBLEMTYPE, filename_len) bind(C, name="define_problem_type_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use indices, only: define_problem_type
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    implicit none
    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: PROBLEMTYPE
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, PROBLEMTYPE, filename_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_define_problem_type(filename_f)
#else
    call define_problem_type(filename_f)
#endif

  end subroutine define_problem_type_c


  !> Ventilation indices
  subroutine ventilation_indices_c() bind(C, name="ventilation_indices_c")

    use indices, only: ventilation_indices
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_ventilation_indices()
#else
    call ventilation_indices()
#endif

  end subroutine ventilation_indices_c
!
!######################################################################
!
!> Perfusion indices
  subroutine perfusion_indices_c() bind(C, name="perfusion_indices_c")

    use indices, only: perfusion_indices
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_perfusion_indices()
#else
    call perfusion_indices()
#endif

  end subroutine perfusion_indices_c

  function get_ne_radius_c() result(res) bind(C, name="get_ne_radius_c")

    use indices, only: get_ne_radius
    implicit none
    integer :: res

    res = get_ne_radius()

  end function get_ne_radius_c

  function get_nj_conc1_c() result(res) bind(C, name="get_nj_conc1_c")

    use indices, only: get_nj_conc1
    implicit none
    integer :: res

    res = get_nj_conc1()

  end function get_nj_conc1_c

end module indices_c
