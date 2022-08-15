
module gasmix_c

contains

!!!########################################################################

  subroutine initial_gasmix_c(initial_concentration,inlet_concentration) bind(C, name="initial_gasmix_c")

    use gasmix, only: initial_gasmix
    use arrays, only: dp
    implicit none

    real(dp),intent(in) :: initial_concentration,inlet_concentration

    call initial_gasmix(initial_concentration, inlet_concentration)

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

    call solve_gasmix(fileid,inr_itr_max,out_itr_max,diffusion_coeff,&
      dt,initial_volume,inlet_concentration,inlet_flow,solve_tolerance,time_end,&
      time_start,inspiration)

  end subroutine solve_gasmix_c

!!! ##################################################################

  subroutine transfer_flow_vol_from_units_c() bind(C, name="transfer_flow_vol_from_units_c")

    use gasmix, only: transfer_flow_vol_from_units
    implicit none

    call transfer_flow_vol_from_units()

  end subroutine transfer_flow_vol_from_units_c

end module gasmix_c
