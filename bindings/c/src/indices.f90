
module indices_c
  implicit none
!Interfaces
private

contains

  !> Ventilation indices
  subroutine ventilation_indices_c() bind(C, name="ventilation_indices_c")

    use indices, only: ventilation_indices
    implicit none

#if defined _WIN32 .and. defined __INTEL_COMPILER
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

#if defined _WIN32 .and. defined __INTEL_COMPILER
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
