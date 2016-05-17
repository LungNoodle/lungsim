
module gasmix_c

contains

!!!########################################################################

  subroutine initial_gasmix_c(initial_concentration,inlet_concentration) bind(C, name="initial_gasmix_c")

    use gasmix, only: initial_gasmix
    use arrays, only: dp
    implicit none

    real(dp),intent(in) :: initial_concentration,inlet_concentration

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_initial_gasmix(initial_concentration, inlet_concentration)
#else
    call initial_gasmix(initial_concentration, inlet_concentration)
#endif

  end subroutine initial_gasmix_c

!!!#######################################################################

  subroutine solve_gasmix_c(fileid,inr_itr_max,out_itr_max,diffusion_coeff,&
       dt,initial_volume,inlet_concentration,inlet_flow,solve_tolerance,time_end,&
       time_start,inspiration) bind(C, name="solve_gasmix_c")

    use gasmix, only: solve_gasmix
    use arrays, only: dp
    implicit none

    integer,intent(in) :: fileid,inr_itr_max,out_itr_max
    real(dp),intent(in) :: diffusion_coeff,dt,initial_volume,&
         inlet_concentration,inlet_flow,solve_tolerance,time_end,time_start
    logical,intent(in) :: inspiration

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_solve_gasmix(fileid,inr_itr_max,out_itr_max,diffusion_coeff,&
      dt,initial_volume,inlet_concentration,inlet_flow,solve_tolerance,time_end,&
      time_start,inspiration)
#else
    call solve_gasmix(fileid,inr_itr_max,out_itr_max,diffusion_coeff,&
      dt,initial_volume,inlet_concentration,inlet_flow,solve_tolerance,time_end,&
      time_start,inspiration)
#endif

  end subroutine solve_gasmix_c

!!! ##################################################################

  subroutine transfer_flow_vol_from_units_c() bind(C, name="transfer_flow_vol_from_units_c")

    use gasmix, only: transfer_flow_vol_from_units
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_transfer_flow_vol_from_units()
#else
    call transfer_flow_vol_from_units()
#endif

  end subroutine transfer_flow_vol_from_units_c

end module gasmix_c
