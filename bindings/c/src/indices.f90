
module indices_c
  implicit none
!Interfaces
private

contains

  !> Ventilation indices
  subroutine ventilation_indices_c() bind(C, name="ventilation_indices_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_VENTILATION_INDICES_C" :: VENTILATION_INDICES_C

    use indices, only: ventilation_indices
    implicit none

    call ventilation_indices()

  end subroutine ventilation_indices_c
!
!######################################################################
!
!> Perfusion indices
  subroutine perfusion_indices_c() bind(C, name="perfusion_indices_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_PERFUSION_INDICES" :: PERFUSION_INDICES

    use indices, only: perfusion_indices
    implicit none

    call perfusion_indices()

  end subroutine perfusion_indices_c

  function get_ne_radius_c() result(res) bind(C, name="get_ne_radius_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_GET_NE_RADIUS_C" :: GET_NE_RADIUS_C

    use indices, only: get_ne_radius
    implicit none
    integer :: res

    res = get_ne_radius()

  end function get_ne_radius_c

  function get_nj_conc1_c() result(res) bind(C, name="get_nj_conc1_c")
  !!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_GET_NJ_CONC1_C" :: GET_NJ_CONC1_C

    use indices, only: get_nj_conc1
    implicit none
    integer :: res

    res = get_nj_conc1()

  end function get_nj_conc1_c

end module indices_c
