module gas_exchange_c
  implicit none
  private

contains
  !!!######################################################################
  subroutine initial_gasexchange_c(initial_concentration,surface_area,V_cap) &
    bind(C, name="initial_gasexchange_c")
    use gas_exchange,only: initial_gasexchange
    use arrays,only: dp
    implicit none
    !!! Parameter List
    real(dp),intent(in) :: initial_concentration,surface_area,V_cap

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_initial_gasexchange(initial_concentration,surface_area,V_cap)
#else
    call initial_gasexchange(initial_concentration,surface_area,V_cap)
#endif
    
  end subroutine initial_gasexchange_c

!!!######################################################################
  subroutine steadystate_gasexchange_c(deadspace,p_i_o2,shunt_fraction, &
       target_p_art_co2,target_p_ven_o2,VCO2,VO2) bind(C, name="steadystate_gasexchange_c")
    use gas_exchange, only: steadystate_gasexchange
    use arrays,only: dp
    implicit none

    !!! Parameter List
    real(dp),intent(in) :: deadspace,p_i_o2,shunt_fraction,target_p_art_co2, &
         target_p_ven_o2,VCO2,VO2

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_steadystate_gasexchange(deadspace,p_i_o2,shunt_fraction, &
       target_p_art_co2,target_p_ven_o2,VCO2,VO2)
#else
    call steadystate_gasexchange(deadspace,p_i_o2,shunt_fraction, &
       target_p_art_co2,target_p_ven_o2,VCO2,VO2)
#endif

  end subroutine steadystate_gasexchange_c



end module gas_exchange_c

