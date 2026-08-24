module lymphatics
  !*Brief Description:* This module contains all lymphatic-specific subroutines
  !
  !*LICENSE:*
  !
  !Test Test Test
  !
  !*Full Description:*
  !
  !This module contains code for pulmonary fluid flux within the alveolo-capillary network,
  !and lymph transport through lymphatic collecting vessels.

  use arrays
  use diagnostics
  use indices
  use other_consts
  use parameter_types, only: lymphatic_params
  use precision ! sets dp for precision

  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
!  private
  public alveolar_capillary_flux
  public lymphatic_transport
contains

!!!#############################################################################

  subroutine alveolar_capillary_flux(ne,write_out)
    !*alveolar_capillary_flux:* calculate fluid flux from blood to interstitium

    use other_consts,only: pi

    integer,intent(in) :: ne
    logical,intent(in) :: write_out

    ! Osmotic pressure parameters
    ! This value is currently g/L and likely needs revision before the osmotic
    ! pathway is activated.
    real(dp),parameter :: capillary_molar_conc = 0.0010250_dp
    real(dp),parameter :: IGC = 62.3630_dp ! ideal gas constant in mmHg.L.mol-1.K-1
    real(dp),parameter :: T = 310.0_dp ! temperature in Kelvin - based on 37C blood temp

    ! Local variables
    integer :: i,liflowcount,nunit,n_timesteps,printcount
    real(dp) :: alveolar_volume,breathing_function,capillary_flow,capillary_vps,capillary_osmotic,capillary_osm_n, &
         capillary_volume,capillary_volume_raw,cap_osm_conc,diffusion,dt,excess,flux_a,flux_b,flux_c,gas_diffusion_restriction, &
         initial_lymphatic_flow,initial_lymphatic_pressure,initial_lymphatic_surface_area,initial_lymphatic_volume, &
         initial_lymph_conc,initial_osm_n,interstitial_capacity,interstitial_capacity_a,interstitial_capacity_b, &
         interstitial_osmotic,interstitial_pressure_a,interstitial_pressure_b,interstitial_saturation,interstitial_volume, &
         int_osm_conc,int_osm_n,interstitial_volume_a,interstitial_volume_b,lung_mass,lymphatic_conductivity,max_Pe, &
         min_Pe,net_flux,fluctuation,mx_pe,mn_pe,intPmax,intPmin,lymphPmax,lymphPmin, &
         open_capillaries,osm_flux,osm_n_flux,overflow,sumflux,sumuptake,test_time,time,time_period,time_sum, &
         time_variable,total_flux,total_hydro_flux,total_osm_flux,capillary_pressure,transit_time,capillary_SA
    real(dp) :: sat1,sat2,sat3,sat4,sat5

    logical :: continue


    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'alveolar_capillary_flux'
    call enter_exit(sub_name,1)

    ! get information for the unit fron unit_field
    ! ne is the 'linker' element in the artery-capillary-vein model, so nunit is for the parent element
    nunit = int(elem_field(ne_unit,elem_cnct(-1,1,ne)))
    capillary_pressure = unit_field(nu_capillary_press,nunit)/133.32239_dp !converted to mmHg
    transit_time = unit_field(nu_tt,nunit)
    capillary_SA = unit_field(nu_sa,nunit)
    max_Pe = unit_field(nu_Pe_max,nunit)
    min_Pe = unit_field(nu_Pe_min,nunit)


    if(write_out)then
       write(*,'('' Unit'',i8,'': Pblood='',f7.2,'' mmHg; TT='',f7.2,'' s; SA='',f8.2,'' mm^2; Pe range='',f6.2,'' mmHg'')') &
            nunit,capillary_pressure,transit_time,capillary_SA,(max_Pe-min_Pe)/133.32239_dp
    endif

    ! Calculated values
    lung_mass = lymphatic_params%lung_mass_g
    open_capillaries = 1.0_dp/6.0_dp !based on open capillaries at rest. Should be solved for by perfusion model?
    capillary_volume_raw = 213.00_dp  !Gehr 1978 based on having a body mass of 74 kg - should be unique to each person
    capillary_volume = (capillary_volume_raw*open_capillaries)/real(num_units)
    !Only used for the osmotic model at the moment (which isn't operational) can volume be obtained elsewhere?

    ! interstitial values
    ! interstitial_capacity == maximal volume before spillover into alveolar in mm^3 - based on 30ml.100g of fluid (Drake 2002)
    interstitial_capacity = ((lymphatic_params%interstitial_capacity_ml_per_100g * &
         (lung_mass/100.0_dp))/real(num_units))*1000.0_dp
    interstitial_capacity_a = lymphatic_params%interstitial_compartment_a_fraction * &
         interstitial_capacity
    interstitial_capacity_b = (1.0_dp-lymphatic_params%interstitial_compartment_a_fraction) * &
         interstitial_capacity
    interstitial_volume_a = 0.0_dp
    interstitial_volume_b = lymphatic_params%initial_interstitial_saturation * interstitial_capacity
    alveolar_volume = 0.0_dp !alveolar volume likely greater at rest, but is lost to respiration - further information needed to put in model
    liflowcount = 0
    initial_lymphatic_surface_area = lymphatic_params%lymphatic_density * capillary_SA
    !both initial lymphatic SA and initial lymphatic hydraulic conductivity make up the filtration coefficient and currently it is the
    !hydraulic conductivity that is adjusted to compensate

    ! initial lymphatic values
    initial_lymphatic_volume = 0.0_dp

    ! Osmotic pressures
    i = 1
    capillary_osmotic = real(i)*capillary_molar_conc*IGC*T
    capillary_osm_n = 0.0_dp
    interstitial_osmotic = 0.0_dp
    int_osm_n = 0.0_dp
    initial_osm_n = 0.0_dp
    total_osm_flux = 0.0_dp

    time = 0.0_dp
    printcount = 0
    total_hydro_flux = 0.0_dp

    intPmax = lymphatic_params%interstitial_pressure_max_mmhg
    intPmin = lymphatic_params%interstitial_pressure_min_mmhg
    lymphPmax = lymphatic_params%lymphatic_pressure_max_mmhg
    lymphPmin = lymphatic_params%lymphatic_pressure_min_mmhg
    mx_pe = max_Pe/133.32239_dp
    mn_pe = min_Pe/133.32239_dp
    fluctuation = ((mx_pe-mn_pe)/2.0_dp)

    if(write_out) write(*,'(''fluctuation '',f8.4)')fluctuation
    ! dt or n_timesteps should be controlled by the user
    n_timesteps = lymphatic_params%integration_steps_per_transit
    dt = transit_time/real(n_timesteps)

    breathing_function = (2.0_dp*pi)/(60.0_dp/lymphatic_params%breathing_rate_bpm)

    continue = .true.

    sat1 = 1.0_dp
    sat2 = 2.0_dp
    sat3 = 3.0_dp
    sat4 = 4.0_dp
    sat5 = 5.0_dp

    !do while(time < lymphatic_params%test_time)
    do while (continue)
       time_sum = dt
       if(transit_time.le.0.001)then
          continue =.false.
       endif
       do while(time_sum < transit_time)
          interstitial_volume = interstitial_volume_a + interstitial_volume_b
          interstitial_saturation = interstitial_volume / interstitial_capacity  ! saturation as a proportion of 0-100%
          time_variable = time + time_sum

          ! calculating flux from capillary into interstitium
          interstitial_pressure_a = fluctuation * sin(time_variable * breathing_function) + &
               (((intPmin-intPmax+(fluctuation*2.0_dp)) * (interstitial_volume_a / interstitial_capacity_a)**2.0_dp) + &
               ((intPmin-intPmax+(fluctuation*2.0_dp))*(-2.0_dp)) * (interstitial_volume_a / interstitial_capacity_a) + &
               (intPmin + fluctuation))
          !arbitrarily defined mathematical relationship between interstitial volume and pressure: (same for a and b)
          !interstitial pressure changes a lot at low volumes with a small volume change, but at high volumes
          !a large volume change is needed to cause a small change in pressure
          interstitial_pressure_b = fluctuation * sin(time_variable * breathing_function) + &
               (((intPmin-intPmax+(fluctuation*2.0_dp)) * (interstitial_volume_b / interstitial_capacity_b)**2.0_dp) + &
               ((intPmin-intPmax+(fluctuation*2.0_dp))*(-2.0_dp)) * (interstitial_volume_b / interstitial_capacity_b) + &
               (intPmin + fluctuation))
          !write(*,'(''Pint: '',f8.4)')interstitial_pressure_b
          ! pressure determined from saturation equation based of literature (currently linear, but likely not)
          if(capillary_pressure > interstitial_pressure_a)then
             flux_a = 0.5_dp * (lymphatic_params%capillary_hydraulic_conductivity * &
                  capillary_SA * (capillary_pressure - &
                  interstitial_pressure_a)) * dt
          else
             flux_a = 0.0_dp
          endif
          if(capillary_pressure > interstitial_pressure_b)then
             flux_b = 0.5_dp * (lymphatic_params%capillary_hydraulic_conductivity * &
                  capillary_SA * (capillary_pressure - &
                  interstitial_pressure_b)) * dt
          else
             flux_b = 0.0_dp
          endif
          flux_c = flux_a + flux_b
          total_hydro_flux = total_hydro_flux + flux_c

          if(interstitial_volume_a + flux_a > interstitial_capacity_a)then
             excess = flux_a - (interstitial_capacity_a - interstitial_volume_a)
             interstitial_volume_a = interstitial_capacity_a
             alveolar_volume = alveolar_volume + 0.5_dp*excess

             if((interstitial_volume_b + 0.5_dp*excess) > interstitial_capacity_b)then
                overflow = 0.5_dp * excess - (interstitial_capacity_b - interstitial_volume_b)
                alveolar_volume = alveolar_volume + overflow
                interstitial_volume_b = interstitial_capacity_b
             else
                interstitial_volume_b = interstitial_volume_b + 0.5_dp*excess
             endif
          else
             interstitial_volume_a = interstitial_volume_a + flux_a
          endif

          interstitial_volume_b = interstitial_volume_b + flux_b
          diffusion = (((interstitial_volume_a/interstitial_capacity_a)-(interstitial_volume_b/interstitial_capacity_b))/ &
               (200_dp)) * dt
          interstitial_volume_b = interstitial_volume_b + diffusion
          interstitial_volume_a = interstitial_volume_a - diffusion

          ! Osmotic
          int_osm_conc = int_osm_n/interstitial_volume
          interstitial_osmotic = real(i)*int_osm_conc*IGC*T
          osm_flux = (lymphatic_params%reflection_coefficient * capillary_SA * &
               (capillary_osmotic - interstitial_osmotic))* dt
          cap_osm_conc = capillary_osm_n / capillary_volume
          int_osm_conc = int_osm_n / interstitial_volume

          if(capillary_osmotic > interstitial_osmotic)then
             capillary_volume = capillary_volume - osm_flux
             interstitial_volume_b = interstitial_volume_b + osm_flux
          else
             capillary_volume = capillary_volume + osm_flux
             interstitial_volume_b = interstitial_volume_b - osm_flux
          endif

          net_flux = osm_flux + flux_a + flux_b
          if(net_flux > 0.0_dp)then
             osm_n_flux = net_flux * cap_osm_conc
             int_osm_n = int_osm_n + osm_n_flux
             !assuming capillary is constant and does not need updating
          else
             osm_n_flux = net_flux * int_osm_conc
             int_osm_n = int_osm_n - osm_n_flux
          endif

          total_osm_flux = total_osm_flux + osm_flux

          !calculating flux from interstitium to initial lymphatics
          if(interstitial_volume_b/interstitial_capacity_b < lymphatic_params%lymphatic_saturation_threshold)then
             lymphatic_conductivity = lymphatic_params%lymphatic_baseline_conductivity_ratio * &
                  lymphatic_params%capillary_hydraulic_conductivity
          !no information on the size of pores or similar for lympatic conductivity so assumed to be similar to capillary.
          else

             lymphatic_conductivity = ((lymphatic_params%lymphatic_conductivity_coefficient_1 * &
                  (interstitial_volume_b / interstitial_capacity_b)**5.0_dp) + &
                  (lymphatic_params%lymphatic_conductivity_coefficient_2 * &
                  (interstitial_volume_b / interstitial_capacity_b)**4.0_dp) + &
                  (lymphatic_params%lymphatic_conductivity_coefficient_3 * &
                  (interstitial_volume_b / interstitial_capacity_b)**3.0_dp) + &
                  (lymphatic_params%lymphatic_conductivity_coefficient_4 * &
                  (interstitial_volume_b / interstitial_capacity_b)**2.0_dp) + &
                  (lymphatic_params%lymphatic_conductivity_coefficient_5 * &
                  (interstitial_volume_b / interstitial_capacity_b)) + &
                  lymphatic_params%lymphatic_conductivity_coefficient_6) * &
                  lymphatic_params%capillary_hydraulic_conductivity
          endif

          initial_lymphatic_pressure = fluctuation * sin((time_variable * breathing_function) + &
               lymphatic_params%pressure_phase_offset_radians) + &
               ((((lymphPmax-lymphPmin-(fluctuation*2.0_dp))* ((interstitial_volume_b / interstitial_capacity_b)**2.0_dp)) + &
               (lymphPmin + fluctuation)))
          !arbitrarily defined mathematical relationship to show that lymphatic pressure does not change much at low volumes with a
          !large volume change, but at high volumes only a small volume change is needed to cause a large change in pressure
          !write(*,'(''Plym: '',f8.4)')initial_lymphatic_pressure
          if(interstitial_volume.le.0.0_dp)then
             initial_lymphatic_flow = 0.0_dp
             interstitial_volume = 0.0_dp
          elseif (interstitial_pressure_b > initial_lymphatic_pressure)then
             initial_lymphatic_flow = (lymphatic_conductivity * initial_lymphatic_surface_area * (interstitial_pressure_b - &
                     initial_lymphatic_pressure)) * dt
          else
             initial_lymphatic_flow = 0.0_dp
          endif

          interstitial_volume_b = interstitial_volume_b - initial_lymphatic_flow
          initial_lymphatic_volume = initial_lymphatic_volume + initial_lymphatic_flow

          int_osm_conc = int_osm_n/interstitial_volume_b
          liflowcount = liflowcount + initial_lymphatic_flow
          initial_osm_n = initial_osm_n + (initial_lymphatic_flow*int_osm_conc)
          int_osm_n = int_osm_n - (initial_lymphatic_flow*int_osm_conc)
          if (initial_lymphatic_volume > 0.0_dp)then
             initial_lymph_conc = initial_osm_n/initial_lymphatic_volume
          else
             initial_lymph_conc = 0.0_dp
          endif

          time_sum = time_sum + dt

          total_flux = total_hydro_flux ! +total_osm_flux
          sumuptake = sumuptake + initial_lymphatic_flow

       enddo

       time = time + transit_time

       printcount = printcount + 1
       if (printcount.eq.100)then
          sat5 = sat4
          sat4 = sat3
          sat3 = sat2
          sat2 = sat1
          sat1 = interstitial_saturation
          if(write_out)then
             write(*,'(f12.2, e12.4, f9.4, e12.4, f11.4, f9.4, e12.4, 2(f11.4), f12.8,f12.4)') time_variable, &
                  flux_c/transit_time*1000.0_dp,interstitial_volume,interstitial_volume_a,interstitial_volume_b, &
                  100.0_dp*interstitial_saturation,initial_lymphatic_flow*1000.0_dp, &
                  initial_lymphatic_volume,total_flux,alveolar_volume,interstitial_pressure_b
          endif
          printcount = 0
          if(time.gt.200.0_dp*transit_time)then
             if((abs(((sat1 + sat2 + sat3 + sat4 +sat5)/5.0_dp)-sat1)).le. &
                  lymphatic_params%convergence_tolerance)then
                continue = .false.
             endif
          endif
       endif

       if(time.gt.10000.0_dp*transit_time) continue = .false.

    enddo !while

    unit_field(nu_intsat,nunit) = interstitial_saturation
    unit_field(nu_time,nunit) = time
    unit_field(nu_av_flux,nunit) = total_flux/time
    unit_field(nu_lymphflow,nunit) = initial_lymphatic_volume/time
    !write(*,'('' T='',e12.3,'': intsat='',e12.6,'' %; flux='',e12.3,'' ul/s; avFlux='',e12.3,'' ul/s; lyFlo='',e12.6,'' ul/s'')') &
    !     unit_field(nu_time,nunit),unit_field(nu_intsat,nunit),&
    !     unit_field(nu_flux,nunit),unit_field(nu_av_flux,nunit),unit_field(nu_lymphflow,nunit)

    call enter_exit(sub_name,2)

  end subroutine alveolar_capillary_flux

!!!#############################################################################

  subroutine lymphatic_transport(filename)
    !*lymphatic_transport:* whole system transport

    character(len=MAX_FILENAME_LEN), intent(in) :: filename
    ! Local parameters
    integer :: i,ne,ne_child,np,nunit
    real(dp) :: time_to_run,time_0
    character(len=300) :: writefile
    character(len=60) :: sub_name
    !real(dp) :: interstitial_saturation,interstitial_pressure_b,nu_av_flux,nu_lymphflow,nu_time
    ! --------------------------------------------------------------------------

    sub_name = 'lymphatic_transport'
    call enter_exit(sub_name,1)

    if(index(filename, ".oplymph")> 0) then !full filename is given
       writefile = filename
    else ! need to append the correct filename extension
       writefile = trim(filename)//'.oplymph'
    endif

    open(10, file=writefile, status='replace')

    call cpu_time(time_0)
    do ne = 1,num_elems
       if(elem_field(ne_group,ne).eq.1.0_dp)then!(elem_field(ne_group,ne)-1.0_dp).lt.TOLERANCE)then
          nunit = int(elem_field(ne_unit,elem_cnct(-1,1,ne)))
          call alveolar_capillary_flux(ne,.false.)
          np = elem_nodes(2,ne)
          write(10,'(i8,11(e14.5))') ne,node_xyz(1:3,np),unit_field(nu_av_flux,nunit),unit_field(nu_intsat,nunit), &
               unit_field(nu_lymphflow,nunit),unit_field(nu_capillary_press,nunit),unit_field(nu_tt,nunit), &
               unit_field(nu_sa,nunit),unit_field(nu_Pe_max,nunit),unit_field(nu_Pe_min,nunit)
       endif
    enddo
    time_to_run = time_0
    call cpu_time(time_to_run)
    time_to_run = time_to_run - time_0
    write(*,*) 'run time=',time_to_run
    close(10)

    do ne = num_elems,1,-1
       if(elem_field(ne_group,ne).eq.1.0_dp)then
          nunit = int(elem_field(ne_unit,elem_cnct(-1,1,ne)))
          elem_field(ne_lymph_flux,ne) = unit_field(nu_av_flux,nunit)
          elem_field(ne_lymph_intsat,ne) = unit_field(nu_intsat,nunit)
       else if(elem_field(ne_group,ne).eq.0.0_dp)then ! artery
          elem_field(ne_lymph_flux,ne) = 0.0_dp
          elem_field(ne_lymph_intsat,ne) = 0.0_dp
          do i = 1,elem_cnct(1,0,ne) ! each child branch
             ne_child = elem_cnct(1,i,ne)
             elem_field(ne_lymph_flux,ne) = elem_field(ne_lymph_flux,ne) + elem_field(ne_lymph_flux,ne_child)
             elem_field(ne_lymph_intsat,ne) = elem_field(ne_lymph_intsat,ne) + elem_field(ne_lymph_intsat,ne_child)
          enddo
          elem_field(ne_lymph_flux,ne) = elem_field(ne_lymph_flux,ne)/real(elem_cnct(1,0,ne))
          elem_field(ne_lymph_intsat,ne) = elem_field(ne_lymph_intsat,ne)/real(elem_cnct(1,0,ne))
       endif
    enddo

    call enter_exit(sub_name,2)

  end subroutine lymphatic_transport

!!!#############################################################################

end module lymphatics


!FUTURE DIRECTIONS
!input a constant to account for difference between current values and expected values
     !Model appeared to be working within the range of the literature but is likely off by a factor of 1000 due to nl to ul conversion error.
     !Need to check that outputted units are correct - most things are in ml and mmHg
!Lymphatic network tree
     !currently all lymph is returned to the circulation immediately, in reality it moves up a tree of lymphatics against a pressure gradient
     !would require excessive modelling perhaps
!impairment of gas diffusion caused by high interstitial saturation
     !unclear at what level this would occur
!alveolar flooding changes
     !alveolar flooding should be able to move between adjacent compartments
     !alveolar fluid should be removed via respiration naturally and therefore should always occur naturally at some low level
!individuality needs to be added in line with the other modules
     !currently operates only on preset male/female values
