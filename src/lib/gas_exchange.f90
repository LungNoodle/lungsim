module gas_exchange
!*Brief Description:* This module is for simulating lung gas exchange.
!
!*LICENSE:*
!
!
!
!*Full Description:*
!More info on what the module does if necessary
!
  use other_consts
  use arrays, only: dp,zero_tol
  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private 
  public initial_gasexchange,steadystate_gasexchange

  real(dp),parameter :: standard_molar_vol = 22.4136e+3_dp ! at STP; mm^3/mmol
  real(dp),parameter :: p_water = 47.0_dp !mmHg
  real(dp),parameter :: Kr = 3.6e+6_dp !L/mol
  real(dp),parameter :: Kt = 10.0e+3_dp !L/mol
  real(dp),parameter :: L = 171.2e+6_dp
  real(dp),parameter :: mcv = 90.0e-15_dp !L
  real(dp),parameter :: mch = 30.0e-12_dp !grams
  real(dp),parameter :: mw = 64458_dp !molecular weight of Hb, g/mol
  real(dp),parameter :: pH=7.4_dp ! pH of plasma
  real(dp),parameter :: temperature=37.0_dp !blood temperature,degrees
  real(dp),parameter :: press_atm=760.0_dp !atmospheric pressure, mmHg
  real(dp),parameter :: o2molvol = 25.44e+3_dp ! mm^3/mmol (converted from 22.41e3 at STP using V2=T2*V1/T1)
  real(dp),parameter :: Hb = 2.33e-3_dp !haemoglobin
  real(dp),parameter :: alphaO2 = 1.46e-6_dp ! O2 solubilitiy in water at T=37, mol/mmHg
  real(dp),parameter :: Wbl=0.81_dp !fractional water content of blood
  !!! The O2 and CO2 concentrations are stored in node_field(nj_conc1/nj_conc2,np)

!!! The O2 and CO2 partial pressures are stored in gasex_field(ng_p_x_y,nunit),
!!! using the nomenclature below for the indices.

!!! uses the following nomenclature:
  !  p_art_o2 ! arterial partial pressure of oxygen (PaO2)
  !  p_alv_o2 ! alveolar partial pressure of oxygen (PAO2)
  !  p_ven_o2 ! venous partial pressure of oxygen (PvO2)

  !  p_art_co2 ! arterial partial pressure of CO2 (PaCO2)
  !  p_alv_co2 ! alveolar partial pressure of CO2 (PACO2)
  !  p_ven_co2 ! venous partial pressure of CO2 (PvCO2)

  !  c_art_o2 ! arterial content of oxygen (CaO2)
  !  c_ven_o2 ! mixed venous content of oxygen (CvO2)

  !  c_art_co2 ! arterial content of CO2 (CaCO2)
  !  c_ven_co2 ! mixed venous content of CO2 (CvCO2)

contains
!
!##############################################################################
!
 subroutine initial_gasexchange(initial_concentration,surface_area,V_cap)
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_INITIAL_GASEXCHANGE" :: INITIAL_GASEXCHANGE
   use indices
   use arrays, only: dp,elem_units_below,gasex_field,node_field,num_nodes,&
         num_units,unit_field
   use diagnostics, only: enter_exit

   !local variables
   real(dp),intent(in) :: initial_concentration
   real(dp), optional ::  surface_area,V_cap

    integer :: nunit
    real(dp) :: Vcap_unit
    real(dp),parameter :: p_water = 47.0_dp
    real(dp),parameter :: press_atm=760.0_dp !atmospheric pressure, mmHg


   character(len=60) :: sub_name

   sub_name = 'initial_gasexchange'
   call enter_exit(sub_name,1)

!!! allocate memory for the gasex_field array, if not already allocated
    if(.not.allocated(gasex_field)) allocate(gasex_field(num_gx,num_units))

!!! initialiase nj_conc2 (for CO2 concentration); currently hardcoded to 40 mmHg
    node_field(nj_conc2,1:num_nodes) = 40.0_dp/(o2molvol*(press_atm-p_water))
    write(*,'('' Initialising Palv_CO2 to 40 mmHg'')')

!!! initialise the gas exchange field for o2 partial pressures
    gasex_field(ng_p_alv_o2,1:num_units) = initial_concentration* &
         o2molvol*(press_atm-p_water)
    gasex_field(ng_p_cap_o2,1:num_units) = initial_concentration*&
         o2molvol*(press_atm-p_water)

    gasex_field(ng_p_alv_co2,1:num_units) = 40.0_dp ! mmHg; should make this user defined
    gasex_field(ng_p_ven_o2,1:num_units) = 40.0_dp ! mmHg; should make this user defined

    unit_field(nu_conc1,1:num_units) = gasex_field(ng_p_alv_o2,1:num_units)/&
         (o2molvol*(press_atm-p_water)) ! from mmHg to mmol/mm^3
    unit_field(nu_conc2,1:num_units) = gasex_field(ng_p_alv_co2,1:num_units)/&
         (o2molvol*(press_atm-p_water)) ! from mmHg to mmol/mm^3

!!! initialise the gas exchange field for co2 partial pressures
    gasex_field(ng_p_alv_co2,1:num_units) = 40.0_dp ! mmHg; should make this user defined
    gasex_field(ng_p_cap_co2,1:num_units) = 40.0_dp ! mmHg; should make this user defined
    gasex_field(ng_p_ven_co2,1:num_units) = 45.0_dp ! mmHg; should make this user defined
    if(present(surface_area))then
!!! initialise the time blood has been in capillaries
      gasex_field(ng_time,1:num_units) = 0.0_dp

!!! capillary volume per gas exchange unit = transit time * flow
    ! elem_units_below is the EFFECTIVE number of units, so this is correct
    !Note that these are calculated on a per unit basis in the perfusion model so can be read in for future iterations
      Vcap_unit = V_cap/elem_units_below(1) ! the capillary volume per gas exchange unit
      gasex_field(ng_Vc,1:num_units) = Vcap_unit
      gasex_field(ng_sa,1:num_units) = surface_area/elem_units_below(1)

!!! transit time through the gas exchange unit = capillary volume/flow
      forall (nunit=1:num_units) gasex_field(ng_tt,nunit) = &
           Vcap_unit/unit_field(nu_perf,nunit)
    endif

   call enter_exit(sub_name,2)
 end subroutine initial_gasexchange

!
!###########################################################################################
!
 subroutine steadystate_gasexchange(c_art_o2,c_ven_o2,&
       p_art_co2,p_art_o2,p_i_o2,p_ven_co2,p_ven_o2,shunt_fraction,&
       VCO2,VO2)
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_STEADYSTATE_GASEXCHANGE" :: STEADYSTATE_GASEXCHANGE
    use arrays,only: dp,elem_field,num_units,gasex_field,node_field,units,unit_field,&
         elem_units_below,elem_nodes
    use indices
    use diagnostics, only:enter_exit

!!! Parameter List
    real(dp),intent(in) :: p_i_o2,shunt_fraction,VCO2,VO2
    real(dp), intent(inout) :: c_art_o2,c_ven_o2,p_art_co2,p_art_o2,p_ven_o2,p_ven_co2
!!! Local Variables
    integer :: counter,k,ne,np,nunit
    real(dp) :: c_art_co2,c_cap_co2,c_cap_o2,c_ven_co2,fun_o2, &
         fun_co2,fdash,p_cap_co2,p_cap_o2,p_art_co2_last, &
         p_art_o2_last,p_ven_co2_last,p_ven_o2_last,Q_total,V_total, &
         target_c_ven_co2,target_c_ven_o2,v_q,p_alv_o2,p_alv_co2

    real(dp),parameter :: m = 0.02386_dp, tol = 1.0e-6_dp
    logical :: continue
    character(len=60) :: sub_name

    sub_name = 'steadystate_gasexchange'
    call enter_exit(sub_name,1)


!!! Calculate steady state gas exchange for CO2
    p_ven_co2_last = p_ven_co2 ! updates at each iteration, until converged
    counter = 1                ! count the number of iterations
    continue = .true.
    do while(continue)
       Q_total = 0.0_dp        ! sum the blood flows; should be same as cardiac output!
       V_total = 0.0_dp        !sum the ventilations
       c_art_co2 = 0.0_dp      ! initialise the content in arterial blood
       do nunit = 1,num_units  ! for each elastic/gas exchange unit
          ne = units(nunit)    ! terminal element number
          p_cap_co2 = gasex_field(ng_p_cap_co2,nunit)      ! initialise capillary CO2
          v_q = unit_field(nu_Vdot0,nunit) &
               /unit_field(nu_perf,nunit)             ! the unit v/q
          if(dabs(v_q) .le. 1.0e-3_dp)then ! no ventilation; cap CO2 == venous CO2
             p_cap_co2 = p_ven_co2
          else                             ! calculate the steady-state PCO2
             fun_co2 = function_co2(v_q,p_cap_co2,p_ven_co2)
             fdash = fdash_co2(v_q,p_cap_co2)
             K=0
             do while(dabs(fun_co2).ge.1.0e-4_dp.and.(k.LT.200))
                K=K+1
                p_cap_co2 = p_cap_co2 - fun_CO2/fdash
                fun_co2 = function_co2(v_q,p_cap_co2,p_ven_co2)
                fdash = fdash_co2(v_q,p_cap_co2)
             enddo
          endif

          Q_total = Q_total + elem_units_below(ne) * dabs(unit_field(nu_perf,nunit)) !mm3/s
          V_total = V_total + elem_units_below(ne) * dabs(unit_field(nu_Vdot0,nunit))


!!! including a limitation that p_cap_co2 cannot be less than zero
          p_cap_co2 = max(p_cap_co2,0.0_dp)


          gasex_field(ng_p_cap_co2,nunit) = p_cap_co2 ! store the capillary/alveolar CO2
          gasex_field(ng_p_alv_co2,nunit) = p_cap_co2 ! store the capillary/alveolar CO2

!!! calculate the content of CO2 in capillary blood, given the partial pressure
          c_cap_co2 = m*p_cap_co2/(1 + m*p_cap_co2)
!!! sum the content in arterial blood (flow weighted sum)
          c_art_co2 = c_art_co2 + elem_units_below(ne)* &
               (c_cap_co2*dabs(unit_field(nu_perf,nunit))) !flow-weighted
!! sum the alveolar co2
          p_alv_co2=p_alv_co2 + elem_units_below(ne)* &
               (p_cap_co2*dabs(unit_field(nu_Vdot0,nunit))) !flow-weighted

       enddo !nunit
!!! update the arterial content of CO2
       c_art_co2 = c_art_co2/Q_total !ml CO2 / ml blood
!!! include the shunt fraction in total arterial CO2
       c_ven_co2 = m * p_ven_co2/(1+m*p_ven_co2)
       c_art_co2 = c_art_co2*(1-SHUNT_FRACTION)+c_ven_co2*SHUNT_FRACTION
 !!  summed alveolar pco2
       p_alv_co2=p_alv_co2/elem_field(ne_Vdot,1)

!!! calculate the partial pressure of pulmonary arterial CO2:
       p_art_co2 = 1/(m*(1-c_art_co2)) ! initialise p_art_co2
       K=0 !counter
       fun_co2 = m*p_art_co2/(1+m*p_art_co2)-c_art_co2
       do while (dabs(fun_co2).ge.1.0e-4_dp.and.(k.lt.200))
          K=K+1
          fdash=m/(1+m*p_art_co2)**2
          p_art_co2 = p_art_co2 - fun_co2/fdash
          fun_co2 = m*p_art_co2/(1+m*p_art_co2)-c_art_co2
       enddo !while
!!! find the p_ven_co2 for the new (target) content of venous CO2
       target_c_ven_co2 = c_art_co2 + VCO2/(elem_field(ne_Qdot,1)*(1+SHUNT_FRACTION))
       !mL(CO2)/mL(blood)   mL/mL  [mm^3/s]/[mm^3/s]
       p_ven_co2 = 1/(m*(1-target_c_ven_co2))
       K=0
       fun_co2=m*p_ven_co2/(1+m*p_ven_co2)-target_c_ven_CO2
       do while (dabs(fun_co2).ge.1.0e-4_dp.and.(k.lt.200))
          K=K+1
          fdash=m/(1+m*p_ven_co2)**2
          p_ven_co2 = p_ven_co2-fun_co2/fdash
          fun_co2 = m*p_ven_co2/(1+m*p_ven_co2)-target_c_ven_co2
       enddo !while
!!! now have updated values for p_art_co2 and p_ven_co2
       write(*,'('' Interim PPs:'',4(f8.3))') p_art_o2,p_ven_o2,p_art_co2,p_ven_co2
!!! check whether p_ven_co2 and p_art_co2 have converged
       if(counter.gt.1)then
          if(abs(p_ven_co2-p_ven_co2_last)/p_ven_co2_last.lt.tol.and. &
               abs(p_art_co2-p_art_co2_last)/p_art_co2_last.lt.tol) then
             continue = .false.
          else
             if(counter.gt.200) continue = .false.
             counter=counter+1
             p_ven_co2_last = p_ven_co2
             p_art_co2_last = p_art_co2
          endif !convergence check
       else
          counter = counter+1
          p_ven_co2_last = p_ven_co2
          p_art_co2_last = p_art_co2
       endif

    enddo !while continue
!    read(*,*)

    write(*,'('' Total blood flow ='',F10.1,'' mm3/s,&
         & alveolar ventilation='',F10.1,'' mm3/s'')') Q_total,V_total
    write(*,'('' Steady-state P_art_CO2 ='',F6.1,'' mmHg,&
         & P_ven_CO2='',F6.1,'' mmHg'')') p_art_co2,p_ven_co2
    write(*,'(''               P_alv_CO2 ='',F6.1,'' mmHg,&
         &  P(A-a)CO2='',F6.1,'' mmHg'')') p_alv_co2,p_alv_co2-p_art_co2

!!! Calculate steady state gas exchange for O2
    p_ven_o2_last = p_ven_o2
    counter = 1
    continue = .true.
    do while (continue)
       c_art_o2 = 0.0_dp
       p_alv_o2=0.0_dp
       do nunit=1,num_units
          ne = units(nunit)
          p_cap_o2 = gasex_field(ng_p_cap_o2,nunit) !initialise
          v_q = unit_field(nu_Vdot0,nunit) &
               /unit_field(nu_perf,nunit)             ! the unit v/q
          if(abs(v_q) .le. 1.0e-3_dp)then
             p_cap_o2 = p_ven_o2
          else ! calculate the steady-state p_cap_o2
             p_cap_co2 = gasex_field(ng_p_cap_co2,nunit)
             fun_o2 = function_o2(p_cap_co2,p_cap_o2,p_i_o2,&
                  p_ven_co2,p_ven_o2,v_q)
             K=0
             do while (abs(fun_o2).ge.1.0e-4_dp.and.(k.lt.200))
                K=K+1
                fdash = fdash_o2(p_cap_co2,p_cap_o2,v_q)
                p_cap_o2 = p_cap_o2 - fun_o2/fdash
                fun_o2 = function_o2(p_cap_co2,p_cap_o2,p_i_o2,&
                     p_ven_co2,p_ven_o2,v_q)
             enddo
          endif
!!! including a limitation that p_cap_o2 cannot be less than p_ven_o2
          p_cap_o2 = max(p_cap_o2,p_ven_o2)

          gasex_field(ng_p_cap_o2,nunit) = p_cap_o2
          gasex_field(ng_p_alv_o2,nunit) = p_cap_o2

!!! calculate the content of O2 in capillary blood, given the partial pressure
          c_cap_o2 = content_from_po2(p_cap_co2,p_cap_o2)
!!! sum the content in arterial blood (flow weighted sum)
          c_art_o2 = c_art_o2 + elem_units_below(ne)* &
               (c_cap_o2*dabs(unit_field(nu_perf,nunit))) !flow-weighted
!! sum the alveolar o2
          p_alv_o2=p_alv_o2 + elem_units_below(ne)* &
               (p_cap_o2*dabs(unit_field(nu_Vdot0,nunit))) !flow-weighted

         ! write(*,*) 'V/Q=',v_q,' pO2=',p_cap_o2,c_cap_o2,c_art_o2
       enddo !nunit

!!! update the arterial content of O2
       c_art_o2 = c_art_o2/Q_total !ml O2 / ml blood
!!! include the shunt fraction in total arterial O2
       c_ven_o2 = content_from_po2(p_ven_co2,p_ven_o2)
       c_art_o2 = c_art_o2*(1.0_dp-SHUNT_FRACTION)+c_ven_o2*SHUNT_FRACTION
!!! calculate the partial pressure of pulmonary arterial O2:
       p_art_o2 = po2_from_content(c_art_o2,p_art_co2)
!!  summed alveolar po2
       p_alv_o2=p_alv_o2/elem_field(ne_Vdot,1)
!!! find the p_ven_o2 for the new (target) content of venous O2
       target_c_ven_o2 = c_art_o2 - VO2/(elem_field(ne_Qdot,1)*(1+SHUNT_FRACTION))
       !mL(O2)/mL(blood)   mL/mL   [mm^3/s]/[mm^3/s]
       p_ven_o2 = po2_from_content(target_c_ven_o2,p_ven_co2)


!!! now have updated values for p_art_o2 and p_ven_o2
!!! check whether p_ven_o2 and p_art_o2 have converged
       if(counter.gt.1)then
          if(abs(p_ven_o2-p_ven_o2_last)/p_ven_o2_last.lt.tol.and. &
               abs(p_art_o2-p_art_o2_last)/p_art_o2_last.lt.tol) then
             continue = .false.
          else
             if(counter.gt.200) continue = .false. !ARC made this one
             counter=counter+1
             p_ven_o2_last = p_ven_o2
             p_art_o2_last = p_art_o2
          endif !convergence check
       else
          counter=counter+1
          p_ven_o2_last = p_ven_o2
          p_art_o2_last = p_art_o2
       endif

    enddo !while continue

    write(*,'('' Steady-state  P_art_O2 ='',F6.1,'' mmHg,&
         &  P_ven_O2='',F6.1,'' mmHg'')') p_art_o2,p_ven_o2
    write(*,'(''               P_alv_O2 ='',F6.1,'' mmHg,&
         &  P(A-a)O2='',F6.1,'' mmHg'')') p_alv_o2,p_alv_o2-p_art_o2

    do nunit=1,num_units
       ne=units(nunit)
       np=elem_nodes(2,ne)
       node_field(nj_conc1,np) = p_cap_o2/(o2molvol*(press_atm-p_water))
       node_field(nj_conc2,np) = p_cap_co2/(o2molvol*(press_atm-p_water))
!!! note: and update the pco2
    enddo

!!! calculate concentrations in the gas exchange units from partial pressures.
!!! this information is used in 'track_back' during expiration
    unit_field(nu_conc1,1:num_units) = gasex_field(ng_p_alv_o2,1:num_units)/&
         (o2molvol*(press_atm-p_water)) ! from mmHg to mmol/mm^3
    unit_field(nu_conc2,1:num_units) = gasex_field(ng_p_alv_co2,1:num_units)/&
         (o2molvol*(press_atm-p_water)) ! from mmHg to mmol/mm^3

    call enter_exit(sub_name,2)

  end subroutine steadystate_gasexchange

  !!! ####################################################

  function function_co2 ( v_q, p_cap_co2, p_ven_co2)

    real(dp),intent(in) :: v_q, p_cap_co2, p_ven_co2
    real(dp),parameter :: m = 0.02386_dp
    real(dp) :: function_co2

    function_co2 = v_q * p_cap_co2 + (press_atm - p_water) * &
         (m*p_cap_co2/(1 + m*p_cap_co2) - m*p_ven_co2/(1 + m*p_ven_co2))

  end function function_co2

!!! ####################################################

  function function_o2(p_cap_co2,p_cap_o2,p_i_o2,&
       p_ven_co2,p_ven_o2,v_q)

!!! Parameters
    real(dp),intent (in) :: p_cap_co2,p_cap_o2,p_i_o2,p_ven_co2,&
         p_ven_o2,v_q
!!! Local variables
    real(dp) :: c_cap_o2,c_ven_o2,function_o2

    c_cap_o2 = content_from_po2(p_cap_co2,p_cap_o2)
    c_ven_o2 = content_from_po2(p_ven_co2,p_ven_o2)

    function_o2 = v_q * (p_i_o2 - p_cap_o2) - 713.0_dp*(c_cap_o2 - c_ven_o2)

  end function function_o2

!!! ####################################################

  function fdash_co2 (v_q,p_cap_co2)

!!! Parameters
    real(dp),intent(in) :: v_q,p_cap_co2
!!! Local variables
    real(dp) :: fdash_co2
    real(dp),parameter :: m = 0.02386_dp

    fdash_co2 = v_q + 713.0_dp * m/(1 + m * p_cap_co2)**2

  end function fdash_co2

!!! ####################################################

  function fdash_o2 (p_x_co2,p_x_o2,v_q)

!!! Parameters
    real(dp),intent(in) :: p_x_co2, v_q
    real(dp) :: p_x_o2
!!! Local variables
    real(dp),parameter :: A1=-8.538889e+3_dp, A2=2.121401e+3_dp, A3=-6.707399e+1_dp,&
         A4=9.359609e+5_dp, A5=-3.134626e+4_dp, A6=2.396167e+3_dp, A7=-6.710441e+1_dp
    real(dp) :: aa,bb,aa_dash,bb_dash,C,gamma,X,fdash_o2

    gamma = 10.0_dp**(0.024_dp*(37.0_dp-temperature)+0.4_dp*(pH-7.4_dp)+ &
         0.06_dp*(DLOG10(DBLE(40.0_dp))-DLOG10(DBLE(p_x_co2))))
    X = p_x_o2*gamma

    aa = (X*(X*(X*(X+A3)+A2)+A1))
    bb = (X*(X*(X*(X+A7)+A6)+A5)+A4)
    aa_dash = gamma*(4.0_dp*X**3 + 3.0_dp*A3*X**2 + 2.0_dp*A2*X+A1)
    bb_dash = gamma*(4.0_dp*X**3 + 3.0_dp*A7*X**2 + 2.0_dp*A6*X+A5)
    C = (Wbl*alphaO2 + 4.0_dp*HB*(aa_dash*bb-aa*bb_dash)/bb**2)*(O2molVol*1.0e-3_dp)
    FDASH_O2 = -v_q - 713.0_dp * C

    RETURN
  END function fdash_o2


!!! ####################################################

  function content_from_po2 (PCO2,po2)

!!! Kelman method for calculating the content of O2 from partial pressure

!!! Parameters
    real(dp) :: PCo2,po2
!!! Local variables
    real(dp) :: content_from_po2,ShbO2

    if(dabs(po2).lt.zero_tol)then
       SHbO2 = 0.0_dp
       content_from_po2 = 0.0_dp
    else

       SHbO2 = saturation_of_o2(pco2,po2)

!!! Calculate O2 content (convert from molar to ml O2 per ml blood)
!!! o2molvol is in units of mm^3/mmol; alphaO2 is mol/mmHg; content should be ml/ml
       content_from_po2 = (Wbl * alphaO2 * PO2 + 4.0_dp * Hb * SHbO2) * (o2molvol*1.0e-3_dp)

    endif

    if(content_from_po2.LT.0.0_dp) content_from_po2=0.0_dp !curve fit behaves poorly at low PO2

  end function content_from_po2


!!! ####################################################

  function saturation_of_o2 (PCO2,po2)

!!! Kelman method for calculating the saturation of O2 from partial pressure

!!! Parameters
    real(dp),intent(in) :: PCo2,po2
!!! Local variables
    real(dp),parameter :: A1=-8.538889e+3_dp, A2=2.121401e+3_dp, A3=-6.707399e+1_dp,&
         A4=9.359609e+5_dp, A5=-3.134626e+4_dp, A6=2.396167e+3_dp, A7=-6.710441e+1_dp
    real(dp) :: saturation_of_o2,X,ShbO2

    if(dabs(po2).lt.zero_tol)then
       SHbO2 = 0.0_dp
    else

!!! Calculate Hb-O2 saturation
       X=PO2*10.0_dp**(0.024_dp*(37.0_dp-temperature)+0.4_dp*(pH-7.4_dp)+ &
            0.06_dp*(DLOG10(DBLE(40.0_dp))-DLOG10(DBLE(PCO2))))
       SHbO2=(X*(X*(X*(X+A3)+A2)+A1))/(X*(X*(X*(X+A7)+A6)+A5)+A4)
    endif
    if(SHbO2.LT.0.0_dp) SHbO2 = 0.0_dp

    saturation_of_o2 = SHbO2

  end function saturation_of_o2

   function po2_from_content(c_o2,p_co2)

!!! Parameter List
    real(dp),intent(in) :: c_o2,p_co2
!!! Local Variables
    integer :: i
    integer,parameter :: max_iterations = 100
    real(dp) :: c_o2_new,c_o2_old,diff_new,diff_old,diff_step,&
         inc,p_o2_new,p_o2_old,po2_from_content
    real(dp),parameter :: tolerance=1.0e-5_dp
    logical :: converged

    if(dabs(c_o2).lt.tolerance)then
       po2_from_content = 0.0_dp
    else
       converged = .false.
       i = 1
       inc = 10.0_dp
       ! initial guess for p_x_o2
       p_o2_new = 50.0_dp  ! mmHg
       c_o2_old = 0.0_dp   ! updated after each iteration from c_o2_new
       c_o2_new = content_from_po2(p_co2,p_o2_new)
       ! Check convergence
       if(abs((c_o2_new - c_o2)/c_o2).lt.tolerance*c_o2) converged =.true.
       ! Loop to find PO2 value
       do while (.not.converged.and.(i.lt.max_iterations))
          ! Modify increment size
          if(c_o2_new.gt.c_o2)then
             inc = -dabs(inc)
          elseif(c_o2_new.lt.c_o2)then
             inc = dabs(inc)
          endif
          if(i.gt.1)then
             diff_new = c_o2_new - c_o2
             diff_old = c_o2_old - c_o2
             diff_step = dabs(c_o2_new-c_o2_old)
             if((diff_old.gt.0.0_dp.and.diff_new.lt.0.0_dp).or. &
                  (diff_old.lt.0.0_dp.and.diff_new.gt.0.0_dp))then ! the last 2 steps straddle point
                inc=inc/2.0_dp
             elseif(dabs(diff_new).gt.diff_step)THEN
                inc=inc*2.0_dp
             endif
          endif

          ! Increment to find new PO2
          p_o2_old = p_o2_new
          c_o2_old = c_o2_new
          p_o2_new = p_o2_new + inc
          c_o2_new = content_from_po2(p_co2,p_o2_new)
          ! Check convergence
          if(dabs((c_o2_new-c_o2)/c_o2).LT.tolerance*c_o2) converged = .true.

          i=i+1

       enddo !while

       if(.not.converged) write(*,'(''>>Error: PO2 value not found'')')

       po2_from_content = p_o2_new

    endif

  end function po2_from_content

end module gas_exchange
