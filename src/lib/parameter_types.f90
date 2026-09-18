module parameter_types

  use precision

  implicit none

  ! make the defined types available across modules
  public :: lung_params
  public :: gx_params
  public :: V_params
  public :: Q_params
  public :: solve_gx_params
  public :: species_params
  public :: lymphatic_params

  ! make the 'update' subroutines accessible via python bindings
  public :: update_lung
  public :: update_mechs
  public :: update_gasexchange
  public :: update_ventilation
  public :: update_cardiac
  public :: update_solve_gx
  public :: update_solve_V
  public :: update_species
  public :: update_lymphatics

  type :: fundamental_constants
     ! fixed constants; no update option
     real(dp) :: o2molvol_37deg = 27.128e+3_dp          ! mm^3/mmol, O2 molecular volume @BTPD = R.T/P_dry
     real(dp) :: o2molvol_btps  = 25.452e+3_dp          ! mm^3/mmol, O2 molecular volume @BTPS = R.T/P_atm
     real(dp) :: o2molvol_stpd  = 22.414e+3_dp          ! mm^3/mmol, O2 molecular volume @STPD
     real(dp) :: max_o2_concentration = 3.93236e-5_dp   ! mmol/mm^3, maximum concentration (at 100% O2)
     real(dp) :: mw = 64458.0_dp                        ! g/mol, molecular weight of Hb
     real(dp) :: alphaO2 = 1.19e-9_dp                   ! mmol/mm^3/mmHg, O2 solubility in water at T=37 (converted from 0.0031 mL/dL/mmHg)
     real(dp) :: alphaCO2 = 3.07e-8_dp                  ! mmol/mm^3/mmHg, CO2 solubility in plasma at T=37 (T dependent)
     real(dp) :: R = 6.2364e4_dp                        ! mm^3.mmHg/mmol/K
     real(dp) :: Wbl = 0.81_dp                          ! fractional water content of blood
     real(dp) :: kappa_o2 = 3.85_dp                     ! mol(O2)/mol(blood); O2 carrying capacity of haemoglobin [pg.26]
     real(dp) :: kc_O2 = 4.4e8_dp                       ! mm^3/mmol/s; forward reaction velocity for O2 with Hb (Weibel 1997)
     real(dp) :: sigma_o2 = 1.4e-9_dp                   ! mmol/mm^3/mmHg; solubility of O2 in blood (Hill et al., 1973a)
     real(dp) :: K = 5.5e-8_dp                          ! mm^2/s/mmHg. Krogh's permeation coefficient for O2. (==3.3e-8 cm2/min/mmHg); Weibel 1993
     real(dp) :: K_CO = 4.47e-8_dp                      ! mm^2/s/mmHg. Krogh's permeation coefficient for CO. (==2.68e-8 cm2/min/mmHg)
  end type fundamental_constants
  
  type :: lung_parameters
     ! parameters for lung orientation and sizing
     integer  :: gravity_dirn = 3                       ! gravity direction, 1== on side, 2==supine, 3==upright          
     real(dp) :: surface_area = 3.0e3_dp * 32.0e3_dp    ! mm^2, gas exchange surface area == 30 mm^2/acinus * 32K acini
     real(dp) :: capillary_volume = 80.0e3_dp           ! mm^3, capillary blood volume
     real(dp) :: FRC = 3.0e6_dp                         ! mm^3, functional residual capacity
     real(dp) :: TLC = 6.0e6_dp                         ! mm^3, total lung capacity
     real(dp) :: anatomical_deadspace = 150.0e3_dp      ! mm^3, volume of airways
     real(dp) :: chest_wall_compliance = 2039.4324_dp   ! mm^3/Pa, == 0.2 L/cmH2O * 1e6 mm^3/L / (98.0665_dp Pa/cmH2O)
     real(dp) :: cov = 0.1_dp                           ! dim, COV for 'randomness' in distal unit sizing
     real(dp) :: Rmax = 1.29_dp                         ! dim, ratio of maximum to average distal unit volume
     real(dp) :: Rmin = 0.79_dp                         ! dim, ratio of minimum to average distal unit volume
  end type lung_parameters

  type :: mechanics_parameters
     ! parameters for 3D soft tissue mechanics and elastic tissue units
     real(dp) :: ref_vol_ratio = 0.5_dp                 ! dim, ratio of the reference volume to the initialised volume (e.g. 0.5 of FRC)
     real(dp) :: a = 0.433_dp                           ! dim, coefficient 'a' in Fung SEDF
     real(dp) :: b = -0.611_dp                          ! dim, coefficient 'b' in Fung SEDF
     real(dp) :: cc = 2500.0_dp                         ! dim, coefficient 'c' in Fung SEDF
  end type mechanics_parameters

  type :: gasexchange_parameters
     ! parameters related to gas properties and gas exchange
     real(dp) :: press_atm = 760.0_dp                   ! mmHg, atmospheric pressure
     real(dp) :: press_H2O = 47.0_dp                    ! mmHg, water vapour pressure
     real(dp) :: diffusion_coeff = 22.5_dp              ! mm^2/s, binary diffusion coefficient
     real(dp) :: FiO2 = 0.21_dp                         ! fractional inspired O2
     real(dp) :: init_p_alv_o2 = 100.0_dp               ! mmHg, initial PO2 in alveoli
     real(dp) :: target_p_art_co2 = 40.0_dp             ! mmHg, target PaCO2
     real(dp) :: target_p_ven_o2 = 38.0_dp              ! mmHg, target PvO2
     real(dp) :: VO2 = 250.0e3_dp /60.0_dp              ! mm^3/s, metabolic consumption of O2
     real(dp) :: VCO2 = 0.8_dp * 250.0e3_dp /60.0_dp    ! mm^3/s, metabolic production of CO2
     real(dp) :: Hb = 150.0_dp * 0.1_dp                 ! g/dL, haemoglobin concentration [M 13.5 - 17.5; F 12.0 - 15.5]
     real(dp) :: pHa = 7.4_dp                           ! pH of arterial blood
     real(dp) :: body_temp = 37.0_dp                    ! degC, body temperature
     character(len=9) :: sat_model = 'dash'             ! the model of O2 saturation. options: dash, kelman, valsecchi
  end type gasexchange_parameters

  type :: ventilation_parameters
     integer  :: breaths_per_minute = 15                ! breaths per minute
     real(dp) :: tidal_volume = 400.0e3_dp              ! mm^3, tidal volume
     real(dp) :: i_to_e_ratio = 1.0_dp                  ! dim, ratio of inspiration to expiration time
     real(dp) :: T_interval = 4.0_dp                    ! s, the total length of the breath
     real(dp) :: press_in = 0.0_dp                      ! Pa, constant pressure applied at the model inlet
     real(dp) :: insp_press_muscle = -196.133_dp        ! Pa, total driving pressure over an insp, == 2 cmH2O * 98.0665 Pa/cmH2O
     character(len=8) :: expiration_type = 'active'     ! type, either active or passive or pressure
  end type ventilation_parameters
  
  type :: cardiac_parameters
     ! parameters for cardiac output
     real(dp) :: cardiac_output = 6.0e6_dp /60.0_dp     ! mm^3/s, cardiac output
     real(dp) :: shunt_fraction = 0.02_dp    !proportion of cardiac output that is shunt
  end type cardiac_parameters

  type :: lymphatic_parameters
     ! Published human defaults for the pulmonary lymphatic transport model.
     real(dp) :: lung_mass_g = 639.0_dp
     real(dp) :: breathing_rate_bpm = 15.0_dp
     real(dp) :: capillary_hydraulic_conductivity = 4.41335e-8_dp
     real(dp) :: interstitial_capacity_ml_per_100g = 30.0_dp
     real(dp) :: initial_interstitial_saturation = 0.48_dp
     real(dp) :: interstitial_compartment_a_fraction = 0.005_dp
     real(dp) :: interstitial_pressure_min_mmhg = -8.0_dp
     real(dp) :: interstitial_pressure_max_mmhg = -1.0_dp
     real(dp) :: lymphatic_pressure_min_mmhg = -8.0_dp
     real(dp) :: lymphatic_pressure_max_mmhg = 1.0_dp
     real(dp) :: lymphatic_density = 1.0_dp
     real(dp) :: lymphatic_saturation_threshold = 0.3_dp
     real(dp) :: lymphatic_baseline_conductivity_ratio = 1.48_dp
     real(dp) :: lymphatic_conductivity_coefficient_1 = 845.87_dp
     real(dp) :: lymphatic_conductivity_coefficient_2 = -2416.7_dp
     real(dp) :: lymphatic_conductivity_coefficient_3 = 2388.5_dp
     real(dp) :: lymphatic_conductivity_coefficient_4 = -922.24_dp
     real(dp) :: lymphatic_conductivity_coefficient_5 = 125.85_dp
     real(dp) :: lymphatic_conductivity_coefficient_6 = -0.0067_dp
     real(dp) :: pressure_phase_offset_radians = 1.570796326794895_dp
     integer :: integration_steps_per_transit = 96
     real(dp) :: convergence_tolerance = 0.000005_dp
     ! Retained legacy fields. Integrity and test_time are not used by the
     ! published equations; reflection_coefficient only affects the inactive
     ! osmotic pathway.
     real(dp) :: lymphatic_integrity = 1.0_dp
     real(dp) :: reflection_coefficient = 0.0_dp
     real(dp) :: test_time = 86400.0_dp
  end type lymphatic_parameters

  type :: solve_gx_parameters
     ! parameters to control gas exchange and gas mixing solutions and solver
     integer  :: num_breaths = 20                       ! max # breaths to solve for
     integer  :: out_itr_max = 200                      ! max # (outer) iterations using GMRES solver.
     integer  :: inr_itr_max = 100                      ! max # (inner) iterations using GMRES solver.
     real(dp) :: theta = 2.0_dp/3.0_dp                  ! weighting for matrices in reduced system: A = M+K*dt*theta; B = -K*c^(n)*dt
     real(dp) :: solve_tolerance = 1.0e-8_dp            ! tolerance for comparing residuals
     real(dp) :: dt = 0.025_dp                          ! time step for PDE solution
     real(dp) :: dt_gx = 0.0025_dp                      ! time step for gas exchange model solution
  end type solve_gx_parameters
    
  type :: solve_vent_parameters
     ! parameters to control ventilation solutions and solver
     integer  :: num_breaths = 10                       ! max # breaths to solve for
     integer  :: max_iterations = 200                   ! max # iterations using GMRES solver.
     real(dp) :: err_tolerance = 1.0e-8_dp              ! tolerance for comparing residuals
     real(dp) :: dt = 0.05_dp                           ! time step for ventilation model solution
  end type solve_vent_parameters
    
  type :: species_parameters
     ! define species-specific parameters (non-geometric or functional)

     ! gas exchange - default human values. use update_species to set new values
     real(dp) :: mcv = 90e-3_dp                         ! pico-L, mean RBC volume for human (ref range 80-96)
     real(dp) :: mch = 30_dp                            ! pico-grams, mean mass Hb/RBC for human (ref 27-31 picograms/cell)
     real(dp) :: Hct = 0.45_dp                          ! hematocrit (Dietel & Kampmann 1971)
     real(dp) :: tau_h = 1.11e-3_dp                     ! mm (1.11 um). Thickness of tissue barrier plus plasma. Weibel (1993)
     real(dp) :: S2 = 2.34e4_dp                         ! coefficient in Severinghaus 
  end type species_parameters
  
  ! retain between modules
  type(fundamental_constants)  :: constants
  type(lung_parameters)        :: lung_params
  type(mechanics_parameters)   :: mech_params
  type(gasexchange_parameters) :: gx_params
  type(ventilation_parameters) :: V_params
  type(cardiac_parameters)     :: Q_params
  type(lymphatic_parameters)   :: lymphatic_params
  type(solve_gx_parameters)    :: solve_gx_params
  type(solve_vent_parameters)  :: solve_V_params
  type(species_parameters)     :: species_params

  
  contains

    subroutine update_lymphatics(param_name, param_value)
      character(len=*), intent(in) :: param_name
      real(dp), intent(in) :: param_value

      select case (trim(param_name))
      case ('lung_mass_g')
         lymphatic_params%lung_mass_g = param_value
      case ('breathing_rate_bpm')
         lymphatic_params%breathing_rate_bpm = param_value
      case ('capillary_hydraulic_conductivity')
         lymphatic_params%capillary_hydraulic_conductivity = param_value
      case ('interstitial_capacity_ml_per_100g')
         lymphatic_params%interstitial_capacity_ml_per_100g = param_value
      case ('initial_interstitial_saturation')
         lymphatic_params%initial_interstitial_saturation = param_value
      case ('interstitial_compartment_a_fraction')
         lymphatic_params%interstitial_compartment_a_fraction = param_value
      case ('interstitial_pressure_min_mmhg')
         lymphatic_params%interstitial_pressure_min_mmhg = param_value
      case ('interstitial_pressure_max_mmhg')
         lymphatic_params%interstitial_pressure_max_mmhg = param_value
      case ('lymphatic_pressure_min_mmhg')
         lymphatic_params%lymphatic_pressure_min_mmhg = param_value
      case ('lymphatic_pressure_max_mmhg')
         lymphatic_params%lymphatic_pressure_max_mmhg = param_value
      case ('lymphatic_surface_area_ratio')
         ! Descriptive alias for the historical lymphatic_density name.
         lymphatic_params%lymphatic_density = param_value
      case ('lymphatic_density')
         lymphatic_params%lymphatic_density = param_value
      case ('lymphatic_saturation_threshold')
         lymphatic_params%lymphatic_saturation_threshold = param_value
      case ('lymphatic_baseline_conductivity_ratio')
         lymphatic_params%lymphatic_baseline_conductivity_ratio = param_value
      case ('lymphatic_conductivity_coefficient_1')
         lymphatic_params%lymphatic_conductivity_coefficient_1 = param_value
      case ('lymphatic_conductivity_coefficient_2')
         lymphatic_params%lymphatic_conductivity_coefficient_2 = param_value
      case ('lymphatic_conductivity_coefficient_3')
         lymphatic_params%lymphatic_conductivity_coefficient_3 = param_value
      case ('lymphatic_conductivity_coefficient_4')
         lymphatic_params%lymphatic_conductivity_coefficient_4 = param_value
      case ('lymphatic_conductivity_coefficient_5')
         lymphatic_params%lymphatic_conductivity_coefficient_5 = param_value
      case ('lymphatic_conductivity_coefficient_6')
         lymphatic_params%lymphatic_conductivity_coefficient_6 = param_value
      case ('pressure_phase_offset_radians')
         lymphatic_params%pressure_phase_offset_radians = param_value
      case ('integration_steps_per_transit')
         lymphatic_params%integration_steps_per_transit = nint(param_value)
      case ('convergence_tolerance')
         lymphatic_params%convergence_tolerance = param_value
      case ('lymphatic_integrity')
         lymphatic_params%lymphatic_integrity = param_value
      case ('reflection_coefficient')
         lymphatic_params%reflection_coefficient = param_value
      case ('test_time')
         lymphatic_params%test_time = param_value
      case ('help')
         call print_lymphatic_parameters()
      case default
         write(*,*) 'WARNING: unknown lymphatic parameter name: ', trim(param_name)
         write(*,*) '         parameters are case sensitive: use all lowercase'
      end select

    end subroutine update_lymphatics

    subroutine print_lymphatic_parameters()
      write(*,'('' Current values for update_lymphatics:'')')
      write(*,'(''    - lung_mass_g = '',es12.5)') lymphatic_params%lung_mass_g
      write(*,'(''    - breathing_rate_bpm = '',es12.5)') lymphatic_params%breathing_rate_bpm
      write(*,'(''    - capillary_hydraulic_conductivity = '',es12.5)') &
           lymphatic_params%capillary_hydraulic_conductivity
      write(*,'(''    - interstitial_capacity_ml_per_100g = '',es12.5)') &
           lymphatic_params%interstitial_capacity_ml_per_100g
      write(*,'(''    - initial_interstitial_saturation = '',es12.5)') &
           lymphatic_params%initial_interstitial_saturation
      write(*,'(''    - interstitial_compartment_a_fraction = '',es12.5)') &
           lymphatic_params%interstitial_compartment_a_fraction
      write(*,'(''    - interstitial_pressure_min_mmhg = '',es12.5)') &
           lymphatic_params%interstitial_pressure_min_mmhg
      write(*,'(''    - interstitial_pressure_max_mmhg = '',es12.5)') &
           lymphatic_params%interstitial_pressure_max_mmhg
      write(*,'(''    - lymphatic_pressure_min_mmhg = '',es12.5)') &
           lymphatic_params%lymphatic_pressure_min_mmhg
      write(*,'(''    - lymphatic_pressure_max_mmhg = '',es12.5)') &
           lymphatic_params%lymphatic_pressure_max_mmhg
      write(*,'(''    - lymphatic_density = '',es12.5)') lymphatic_params%lymphatic_density
      write(*,'(''    - lymphatic_saturation_threshold = '',es12.5)') &
           lymphatic_params%lymphatic_saturation_threshold
      write(*,'(''    - lymphatic_baseline_conductivity_ratio = '',es12.5)') &
           lymphatic_params%lymphatic_baseline_conductivity_ratio
      write(*,'(''    - lymphatic_conductivity_coefficient_1 = '',es12.5)') &
           lymphatic_params%lymphatic_conductivity_coefficient_1
      write(*,'(''    - lymphatic_conductivity_coefficient_2 = '',es12.5)') &
           lymphatic_params%lymphatic_conductivity_coefficient_2
      write(*,'(''    - lymphatic_conductivity_coefficient_3 = '',es12.5)') &
           lymphatic_params%lymphatic_conductivity_coefficient_3
      write(*,'(''    - lymphatic_conductivity_coefficient_4 = '',es12.5)') &
           lymphatic_params%lymphatic_conductivity_coefficient_4
      write(*,'(''    - lymphatic_conductivity_coefficient_5 = '',es12.5)') &
           lymphatic_params%lymphatic_conductivity_coefficient_5
      write(*,'(''    - lymphatic_conductivity_coefficient_6 = '',es12.5)') &
           lymphatic_params%lymphatic_conductivity_coefficient_6
      write(*,'(''    - pressure_phase_offset_radians = '',es12.5)') &
           lymphatic_params%pressure_phase_offset_radians
      write(*,'(''    - integration_steps_per_transit = '',i0)') &
           lymphatic_params%integration_steps_per_transit
      write(*,'(''    - convergence_tolerance = '',es12.5)') &
           lymphatic_params%convergence_tolerance
    end subroutine print_lymphatic_parameters

    
    subroutine update_lung(param_name, param_value)

      character(len=*), intent(in) :: param_name
      real(dp), intent(in) :: param_value

      select case (trim(param_name))
      case ('gravity_dirn')
         lung_params%gravity_dirn = int(param_value)
      case ('surface_area')
         lung_params%surface_area = param_value
      case ('capillary_volume')
         lung_params%capillary_volume = param_value
      case('FRC')
         lung_params%FRC = param_value
      case('TLC')
         lung_params%TLC = param_value
      case('anatomical_deadspace')
         lung_params%anatomical_deadspace = param_value
      case('chest_wall_compliance')
         lung_params%chest_wall_compliance = param_value
      case('cov')
         lung_params%cov = param_value
      case('rmax')
         lung_params%Rmax = param_value
      case('rmin')
         lung_params%Rmin = param_value
      case ('help')
         write(*,'('' Current values for update_lung:'')') 
         write(*,'(''    - gravity_dirn  = '', i6)') lung_params%gravity_dirn
         write(*,'(''    - surface_area  = '', d8.3, '' mm2'')') lung_params%surface_area
         write(*,'(''    - capillary_volume  = '', d8.3, '' mm3'')') lung_params%capillary_volume
         write(*,'(''    - FRC  = '', d8.2, '' mm3'')') lung_params%FRC
         write(*,'(''    - TLC  = '', d8.2, '' mm3'')') lung_params%TLC
         write(*,'(''    - anatomical_deadspace  = '', d8.2, '' mm3'')') lung_params%anatomical_deadspace
         write(*,'(''    - chest_wall_compliance  = '', d8.2, '' mm3/Pa'')') lung_params%chest_wall_compliance
         write(*,'(''    - cov  = '', f8.2)') lung_params%cov
         write(*,'(''    - Rmax  = '', f8.2)') lung_params%Rmax
         write(*,'(''    - Rmin  = '', f8.2)') lung_params%Rmin
      case default
         write(*,*) 'WARNING: unknown parameter name: ', trim(param_name)
         write(*,*) '         parameters are case sensitive: use all lowercase'
      end select

    end subroutine update_lung

    subroutine update_mechs(param_name, param_value)
      character(len=*), intent(in) :: param_name
      real(dp), intent(in) :: param_value

      select case (trim(param_name))
      case('ref_vol_ratio')
         mech_params%ref_vol_ratio = param_value
      case('a')
         mech_params%a = param_value
      case('b')
         mech_params%b = param_value
      case('c')
         mech_params%cc = param_value
      case ('help')
         write(*,'('' Current values for update_mechs:'')') 
         write(*,'(''    -  ref_vol_ratio = '', f8.2)') mech_params%ref_vol_ratio
         write(*,'(''    -  a = '', f8.2)') mech_params%a
         write(*,'(''    -  b = '', f8.2)') mech_params%b
         write(*,'(''    -  c = '', f8.2)') mech_params%cc
      case default
         write(*,*) 'WARNING: unknown parameter name: ', trim(param_name)
         write(*,*) '         parameters are case sensitive: use all lowercase'
      end select
    end subroutine update_mechs

    subroutine update_gasexchange(param_name, param_value)
      
      character(len=*), intent(in) :: param_name
      real(dp), intent(in) :: param_value

      select case (trim(param_name))
      case ('press_atm')
         gx_params%press_atm = param_value
      case ('press_h2o')
         gx_params%press_H2O = param_value
      case ('diffusion_coefficient')
         gx_params%diffusion_coeff = param_value
      case ('fio2')
         gx_params%FiO2 = param_value
      case('init_p_alv_o2')
         gx_params%init_p_alv_o2 = param_value
      case ('target_p_art_co2')
         gx_params%target_p_art_co2 = param_value
      case ('target_p_ven_o2')
         gx_params%target_p_ven_o2 = param_value
      case ('vo2')
         gx_params%VO2 = param_value
      case ('vco2')
         gx_params%VCO2 = param_value
      case ('hb')
         gx_params%Hb = param_value
      case ('pha')
         gx_params%pHa = param_value
      case ('body_temp')
         gx_params%body_temp = param_value
      case ('help')
         write(*,'('' Current values for update_gasexchange:'')') 
         write(*,'(''    -  press_atm = '', f8.1, '' mmHg'')') gx_params%press_atm
         write(*,'(''    -  press_h2o = '', f8.1, '' mmHg'')') gx_params%press_H2O
         write(*,'(''    -  diffusion_coefficient = '', f8.2, '' mm2/s'')') gx_params%diffusion_coeff
         write(*,'(''    -  fio2 = '', f8.2)') gx_params%FiO2
         write(*,'(''    -  init_p_alv_o2 = '', f8.1, '' mmHg'')') gx_params%init_p_alv_o2
         write(*,'(''    -  target_p_art_co2 = '', f8.1, '' mmHg'')') gx_params%target_p_art_co2
         write(*,'(''    -  target_p_ven_o2 = '', f8.1, '' mmHg'')') gx_params%target_p_ven_o2
         write(*,'(''    -  vo2 = '', f8.2, '' mm3/s'')') gx_params%VO2
         write(*,'(''    -  vco2 = '', f8.2, '' mm3/s'')') gx_params%VCO2
         write(*,'(''    -  hb = '', f8.2, '' g/dL'')') gx_params%Hb
         write(*,'(''    -  pha = '', f8.2)') gx_params%pHa
         write(*,'(''    -  body_temp = '', f8.2, '' degC'')') gx_params%body_temp
      case default
         write(*,*) 'WARNING: unknown parameter name: ', trim(param_name)
         write(*,*) '         parameters are case sensitive: use all lowercase'
      end select
    end subroutine update_gasexchange

    
    subroutine update_ventilation(param_name, param_value)
      
      character(len=*), intent(in) :: param_name
      real(dp), intent(in) :: param_value
      
      select case (trim(param_name))
      case ('breaths_per_min')
         V_params%breaths_per_minute = param_value
      case ('tidal_volume')
         V_params%tidal_volume = param_value
      case ('i_to_e_ratio')
         V_params%i_to_e_ratio = param_value
      case ('t_interval')
         V_params%t_interval = param_value
      case ('press_in')
         V_params%press_in = param_value
      case ('insp_press_muscle')
         V_params%insp_press_muscle = param_value
      case ('help')
         write(*,'('' Current values for update_ventilation:'')') 
         write(*,'(''    -  breaths_per_min = '', i6)')  V_params%breaths_per_minute
         write(*,'(''    -  tidal_volume = '', f8.1, '' mm3/s'')')  V_params%tidal_volume
         write(*,'(''    -  i_to_e_ratio = '', f6.2)')  V_params%i_to_e_ratio
         write(*,'(''    -  t_interval = '', f6.2, '' s'')')  V_params%t_interval
         write(*,'(''    -  press_in = '', f6.2, '' Pa'')')  V_params%press_in
         write(*,'(''    -  insp_press_muscle = '', f8.3, '' Pa'')')  V_params%insp_press_muscle
      case default
         write(*,*) 'WARNING: unknown parameter name: ', trim(param_name)
         write(*,*) '         parameters are case sensitive: use all lowercase'
      end select
    end subroutine update_ventilation


    subroutine update_cardiac(param_name, param_value)
      
      character(len=*), intent(in) :: param_name
      real(dp), intent(in) :: param_value

      select case (trim(param_name))
      case ('cardiac_output')
         Q_params%cardiac_output = param_value
      case ('shunt_fraction')
         Q_params%shunt_fraction = param_value
      case ('help')
         write(*,'('' Current values for update_cardiac:'')') 
         write(*,'(''    -  cardiac_output = '', f8.1, '' mm^3/s'')')  Q_params%cardiac_output
         write(*,'(''    -  shunt_fraction = '', f6.3)')  Q_params%shunt_fraction
      case default
         write(*,*) 'WARNING: unknown parameter name: ', trim(param_name)
         write(*,*) '         parameters are case sensitive: use all lowercase'
      end select
    end subroutine update_cardiac
    

    subroutine update_solve_gx(param_name, param_value)
      
      character(len=*), intent(in) :: param_name
      real(dp), intent(in) :: param_value

      select case (trim(param_name))
      case ('number_of_breaths')
         solve_gx_params%num_breaths = param_value
      case ('max_outer_iterations')
         solve_gx_params%out_itr_max = param_value
      case ('max_inner_iterations')
         solve_gx_params%inr_itr_max = param_value
      case ('solver_tolerance')
         solve_gx_params%solve_tolerance = param_value
      case ('dt_solve')
         solve_gx_params%dt = param_value
      case ('dt_gx')
         solve_gx_params%dt_gx = param_value
      case ('help')
         write(*,'('' Current values for update_solve_gx:'')') 
         write(*,'(''    - number_of_breaths  = '', i6, '' '')')  solve_gx_params%num_breaths
         write(*,'(''    - max_outer_iterations  = '', i6, '' '')')  solve_gx_params%out_itr_max
         write(*,'(''    - max_inner_iterations  = '', i6, '' '')')  solve_gx_params%inr_itr_max
         write(*,'(''    - solver_tolerance  = '', d8.3, '' '')')  solve_gx_params%solve_tolerance
         write(*,'(''    - dt_solve  = '', d8.3, '' '')')  solve_gx_params%dt
         write(*,'(''    - dt_gx  = '', d8.3, '' '')')  solve_gx_params%dt_gx
      case default
         write(*,*) 'WARNING: unknown parameter name: ', trim(param_name)
         write(*,*) '         parameters are case sensitive: use all lowercase'
      end select
    end subroutine update_solve_gx

    subroutine update_solve_V(param_name, param_value)
      
      character(len=*), intent(in) :: param_name
      real(dp), intent(in) :: param_value

      select case (trim(param_name))
      case ('number_of_breaths')
         solve_V_params%num_breaths = param_value
      case ('max_iterations')
         solve_V_params%max_iterations = param_value
      case ('err_tolerance')
         solve_V_params%err_tolerance = param_value
      case ('dt')
         solve_V_params%dt = param_value
      case ('help')
         write(*,'('' Current values for update_solve_V:'')') 
         write(*,'(''    - number_of_breaths  = '', i6, '' '')')  solve_V_params%num_breaths
         write(*,'(''    - max_iterations  = '', i6, '' '')')  solve_V_params%max_iterations
         write(*,'(''    - err_tolerance  = '', d8.3, '' '')')  solve_V_params%err_tolerance
         write(*,'(''    - dt  = '', d8.3, '' '')')  solve_V_params%dt
       case default
         write(*,*) 'WARNING: unknown parameter name: ', trim(param_name)
         write(*,*) '         parameters are case sensitive: use all lowercase'
      end select
    end subroutine update_solve_V
    
    subroutine update_species(param_name)
      
      character(len=*), intent(in) :: param_name    ! human, rabbit, rat, mouse

      select case (trim(param_name))
      case ('Human')
         species_params%mcv = 90e-3_dp 
         species_params%mch = 30_dp    
         species_params%Hct = 0.45_dp  
         species_params%tau_h = 1.11e-3_dp
         species_params%S2 = 2.34e4_dp    
      case ('Rabbit')
         species_params%mcv = 66.7e-3_dp  
         species_params%mch = 20.95_dp    
         species_params%Hct = 0.436_dp    
         species_params%tau_h = 0.8e-3_dp 
         species_params%S2 = 3.5e4_dp     
      case ('Rat')
         species_params%mcv = 59.35e-3_dp 
         species_params%mch = 17.9_dp     
         species_params%Hct = 0.5182_dp   
         species_params%tau_h = 0.754e-3_dp
         species_params%S2 = 5.0e4_dp      
      case ('Mouse')
         species_params%mcv = 55.1e-3_dp     
         species_params%mch = 15.95_dp       
         species_params%Hct = 0.523_dp       
         species_params%tau_h = 0.7e-3_dp   
         species_params%S2 = 8.0e4_dp        
      case default
         write(*,*) 'WARNING: unknown parameter name: ', trim(param_name)
         write(*,*) '         parameters are case sensitive: use all lowercase'
      end select
    end subroutine update_species


  end module parameter_types
