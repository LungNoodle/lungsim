module ventilation
  !*Brief Description:* This module handles all code specific to
  ! simulating ventilation
  !
  !*LICENSE:*
  !TBC
  !
  !
  !*Full Description:*
  !
  ! This module handles all code specific to simulating ventilation 
  
  use arrays
  use diagnostics
  use exports
  use geometry
  use indices
  use other_consts
  use precision
  
  implicit none
  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private
  public evaluate_vent
  public evaluate_uniform_flow
  public two_unit_test
  public sum_elem_field_from_periphery

  real(dp),parameter,private :: gravity = 9.81e3_dp         ! mm/s2
!!! for air
  real(dp),parameter,private :: gas_density =   1.146e-6_dp ! g.mm^-3
  real(dp),parameter,private :: gas_viscosity = 1.8e-5_dp   ! Pa.s

contains

!!!#############################################################################

  subroutine evaluate_vent
    !*evaluate_vent:* Sets up and solves dynamic ventilation model
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_VENT" :: EVALUATE_VENT

    ! Local variables
    integer :: gdirn                  ! 1(x), 2(y), 3(z); upright lung (for our
    !                                   models) is z, supine is y.
    integer :: iter_step,n,ne,num_brths,num_itns,nunit
    real(dp) :: chestwall_restvol     ! resting volume of chest wall
    real(dp) :: chest_wall_compliance ! constant compliance of chest wall
    real(dp) :: constrict             ! for applying uniform constriction
    real(dp) :: COV                   ! COV of tissue compliance
    real(dp) :: i_to_e_ratio          ! ratio inspiration to expiration time
    real(dp) :: p_mus                 ! muscle (driving) pressure
    real(dp) :: pmus_factor_ex        ! pmus_factor (_in and _ex) used to scale 
    real(dp) :: pmus_factor_in        ! modifies driving pressures to converge 
    !                                   tidal volume and expired volume to the 
    !                                   target volume.
    real(dp) :: pmus_step             ! change in Ppl for driving flow (Pa)
    real(dp) :: press_in              ! constant pressure at entry to model (Pa)
    real(dp) :: press_in_total        ! dynamic pressure at entry to model (Pa)
    real(dp) :: refvol                ! proportion of model for 'zero stress'
    real(dp) :: RMaxMean              ! ratio max to mean volume
    real(dp) :: RMinMean              ! ratio min to mean volume
    real(dp) :: sum_expid             ! sum of expired volume  (mm^3)
    real(dp) :: sum_tidal             ! sum of inspired volume  (mm^3)
    real(dp) :: Texpn                 ! time for expiration (s)
    real(dp) :: T_interval            ! the total length of the breath (s)
    real(dp) :: Tinsp                 ! time for inspiration (s)
    real(dp) :: undef                 ! the zero stress volume. undef < RV 
    real(dp) :: volume_target         ! the target tidal volume (mm^3)

    real(dp) :: dpmus,dt,endtime,err_est,err_tol,FRC,init_vol,last_vol, &
         current_vol,Pcw,ppl_current,pptrans,prev_flow,ptrans_frc, &
         sum_dpmus,sum_dpmus_ei,time,totalc,Tpass,ttime,volume_tree,WOBe,WOBr, &
         WOBe_insp,WOBr_insp,WOB_insp
    character :: expiration_type*(10) ! active (sine wave), passive, pressure
    logical :: CONTINUE,converged

    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'evaluate_vent'
    call enter_exit(sub_name,1)
    
!!! Initialise variables:
    pmus_factor_in = 1.0_dp
    pmus_factor_ex = 1.0_dp
    time = 0.0_dp !initialise the simulation time.
    n = 0 !initialise the 'breath number'. incremented at start of each breath.
    sum_tidal = 0.0_dp ! initialise the inspired and expired volumes
    sum_expid = 0.0_dp
    last_vol = 0.0_dp

!!! set default values for the parameters that control the breathing simulation
!!! these should be controlled by user input (showing hard-coded for now)

    call read_params_evaluate_flow(gdirn, chest_wall_compliance, &
       constrict, COV, FRC, i_to_e_ratio, pmus_step, press_in,&
       refvol, RMaxMean, RMinMean, T_interval, volume_target, expiration_type)
    call read_params_main(num_brths, num_itns, dt, err_tol)

!!! set dynamic pressure at entry. only changes for the 'pressure' option
    press_in_total = press_in
    
!!! calculate key variables from the boundary conditions/problem parameters
    Texpn = T_interval / (1.0_dp+i_to_e_ratio)
    Tinsp = T_interval - Texpn

!!! store initial branch lengths, radii, resistance etc. in array 'elem_field'
    call update_elem_field(1.0_dp)
    call update_resistance
    call volume_of_mesh(init_vol,volume_tree)
    
!!! distribute the initial tissue unit volumes along the gravitational axis.
    call set_initial_volume(gdirn,COV,FRC*1.0e+6_dp,RMaxMean,RMinMean)
    undef = refvol * (FRC*1.0e+6_dp-volume_tree)/dble(elem_units_below(1))

!!! calculate the total model volume
    call volume_of_mesh(init_vol,volume_tree)

    write(*,'('' Anatomical deadspace = '',F8.3,'' ml'')') &
         volume_tree/1.0e+3_dp ! in mL
    write(*,'('' Respiratory volume   = '',F8.3,'' L'')') &
         (init_vol-volume_tree)/1.0e+6_dp !in L
    write(*,'('' Total lung volume    = '',F8.3,'' L'')') &
         init_vol/1.0e+6_dp !in L

    unit_field(nu_dpdt,1:num_units) = 0.0_dp

!!! calculate the compliance of each tissue unit
    call tissue_compliance(chest_wall_compliance,undef)
    totalc = SUM(unit_field(nu_comp,1:num_units)) !the total model compliance
    call update_pleural_pressure(ppl_current) !calculate new pleural pressure
    pptrans=SUM(unit_field(nu_pe,1:num_units))/num_units

    chestwall_restvol = init_vol + chest_wall_compliance * (-ppl_current)
    Pcw = (chestwall_restvol - init_vol)/chest_wall_compliance
    write(*,'('' Chest wall RV = '',F8.3,'' L'')') chestwall_restvol/1.0e+6_dp
        
    call write_flow_step_results(chest_wall_compliance,init_vol, &
         current_vol,ppl_current,pptrans,Pcw,p_mus,0.0_dp,0.0_dp)
    
    continue = .true.
    do while (continue)
       n = n + 1 ! increment the breath number
       ttime = 0.0_dp ! each breath starts with ttime=0
       endtime = T_interval * n - 0.5_dp * dt ! the end time of this breath
       p_mus = 0.0_dp 
       ptrans_frc = SUM(unit_field(nu_pe,1:num_units))/num_units !ptrans at frc

       if(n.gt.1)then !write out 'end of breath' information
          call write_end_of_breath(init_vol,current_vol,pmus_factor_in, &
               pmus_step,sum_expid,sum_tidal,volume_target,WOBe_insp, &
               WOBr_insp,WOB_insp)
          
          if(abs(volume_target).gt.1.0e-5_dp)THEN
             ! modify driving muscle pressure by volume_target/sum_tidal
             ! this increases p_mus for volume_target>sum_tidal, and
             ! decreases p_mus for volume_target<sum_tidal
             pmus_factor_in = pmus_factor_in * abs(volume_target/sum_tidal)
             pmus_factor_ex = pmus_factor_ex * abs(volume_target/sum_expid)
          endif
          sum_tidal = 0.0_dp !reset the tidal volume
          sum_expid = 0.0_dp !reset the expired volume
          unit_field(nu_vt,1:num_units) = 0.0_dp !reset acinar tidal volume
          sum_dpmus = 0.0_dp
          sum_dpmus_ei = 0.0_dp
       endif

!!! solve for a single breath (for time up to endtime)
       do while (time.lt.endtime) 
          ttime = ttime + dt ! increment the breath time
          time = time + dt ! increment the whole simulation time
!!!.......calculate the flow and pressure distribution for one time-step
          call evaluate_vent_step(num_itns,chest_wall_compliance, &
               chestwall_restvol,dt,err_tol,init_vol,last_vol,current_vol, &
               Pcw,pmus_factor_ex,pmus_factor_in,pmus_step,p_mus,ppl_current, &
               pptrans,press_in_total,prev_flow,ptrans_frc,sum_dpmus,sum_dpmus_ei, &
               sum_expid,sum_tidal,texpn,time,tinsp,ttime,undef,WOBe,WOBr, &
               WOBe_insp,WOBr_insp,WOB_insp,expiration_type, &
               dpmus,converged,iter_step)
!!!.......update the estimate of pleural pressure
          call update_pleural_pressure(ppl_current) ! new pleural pressure
           
          call write_flow_step_results(chest_wall_compliance,init_vol, &
               current_vol,ppl_current,pptrans,Pcw,p_mus,time,ttime)

       enddo !while time<endtime
       
!!!....check whether simulation continues
       continue = ventilation_continue(n,num_brths,sum_tidal,volume_target)

    enddo !...WHILE(CONTINUE)

    call write_end_of_breath(init_vol,current_vol,pmus_factor_in,pmus_step, &
         sum_expid,sum_tidal,volume_target,WOBe_insp,WOBr_insp,WOB_insp)

!!! Transfer the tidal volume for each elastic unit to the terminal branches,
!!! and sum up the tree. Divide by inlet flow. This gives the time-averaged and
!!! normalised flow field for the tree.
    do nunit = 1,num_units 
       ne = units(nunit) !local element number
       elem_field(ne_Vdot,ne) = unit_field(nu_vt,nunit)
    enddo
    unit_field(nu_vent,:) = unit_field(nu_vt,:)/(Tinsp+Texpn)
    call sum_elem_field_from_periphery(ne_Vdot)
    elem_field(ne_Vdot,1:num_elems) = &
         elem_field(ne_Vdot,1:num_elems)/elem_field(ne_Vdot,1)

!    call export_terminal_solution(TERMINAL_EXNODEFILE,'terminals')

    call enter_exit(sub_name,2)

  end subroutine evaluate_vent

!!!#############################################################################

  subroutine evaluate_vent_step(num_itns,chest_wall_compliance, &
       chestwall_restvol,dt,err_tol,init_vol,last_vol,current_vol,Pcw, &
       pmus_factor_ex,pmus_factor_in,pmus_step,p_mus,ppl_current,pptrans, &
       press_in_total,prev_flow,ptrans_frc,sum_dpmus,sum_dpmus_ei,sum_expid, &
       sum_tidal,texpn,time,tinsp,ttime,undef,WOBe,WOBr,WOBe_insp,WOBr_insp, &
       WOB_insp,expiration_type,dpmus,converged,iter_step)

    integer,intent(in) :: num_itns
    real(dp),intent(in) :: chest_wall_compliance,chestwall_restvol,dt, &
         err_tol,init_vol,pmus_factor_ex,pmus_factor_in,pmus_step,pptrans, &
         press_in_total,ptrans_frc,texpn,time,tinsp,ttime,undef
    real(dp) :: last_vol,current_vol,Pcw,ppl_current,prev_flow,p_mus, &
         sum_dpmus,sum_dpmus_ei,sum_expid,sum_tidal,WOBe,WOB_insp,WOBe_insp, &
         WOBr,WOBr_insp
    character,intent(in) :: expiration_type*(*)
    ! Local variables
    integer :: iter_step
    real(dp) :: dpmus,err_est,totalC,Tpass,volume_tree
    logical :: converged
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'evaluate_vent_step'
    call enter_exit(sub_name,1)

!!! Solve for a new flow and pressure field
!!! We will estimate the flow into each terminal lumped
!!! parameter unit (assumed to be an acinus), so we can calculate flow
!!! throughout the rest of the tree simply by summation. After summing
!!! the flows we can use the resistance equation (P0-P1=R1*Q1) to update
!!! the pressures throughout the tree.

    ! set the increment in driving (muscle) pressure
    call set_driving_pressures(dpmus,dt,pmus_factor_ex,pmus_factor_in, &
         pmus_step,p_mus,Texpn,Tinsp,ttime,expiration_type)
    prev_flow = elem_field(ne_Vdot,1)
    
!!! Solve for a new flow and pressure field
!!! We will estimate the flow into each terminal lumped
!!! parameter unit (assumed to be an acinus), so we can calculate flow
!!! throughout the rest of the tree simply by summation. After summing
!!! the flows we can use the resistance equation (P0-P1=R1*Q1) to update
!!! the pressures throughout the tree.
    
    !initialise Qinit to the previous flow
    elem_field(ne_Vdot0,1:num_elems) = elem_field(ne_Vdot,1:num_elems)
    converged = .FALSE.
    iter_step=0
    do while (.not.converged)
       iter_step = iter_step+1 !count the iterative steps
       call estimate_flow(dpmus,dt,err_est) !analytic solution for Q
       if(iter_step.gt.1.and.err_est.lt.err_tol)then
          converged = .TRUE.
       else if(iter_step.gt.num_itns)then
          converged = .TRUE.
          write(*,'('' Warning: lower convergence '// &
               'tolerance and time step - check values, Error='',D10.3)') &
               err_est
       endif
       call sum_elem_field_from_periphery(ne_Vdot) !sum flows UP tree
       call update_elem_field(1.0_dp)
       call update_resistance ! updates resistances
       call update_node_pressures(press_in_total) ! updates the pressures at nodes
       call update_unit_dpdt(dt) ! update dP/dt at the terminal units
    enddo !converged
    
    call update_unit_volume(dt) ! Update tissue unit volumes, unit tidal vols
    call volume_of_mesh(current_vol,volume_tree) ! calculate mesh volume
    call update_elem_field(1.0_dp)
    call update_resistance  !update element lengths, volumes, resistances
    call tissue_compliance(chest_wall_compliance,undef) ! unit compliances
    totalc = SUM(unit_field(nu_comp,1:num_units)) !the total model compliance
    call update_proximal_pressure ! pressure at proximal nodes of end branches
    call calculate_work(current_vol-init_vol,current_vol-last_vol,WOBe,WOBr, &
         pptrans)!calculate work of breathing
    last_vol=current_vol
    Pcw = (chestwall_restvol - current_vol)/chest_wall_compliance
    
    ! increment the tidal volume, or the volume expired
    if(elem_field(ne_Vdot,1).gt.0.0_dp)then
       sum_tidal = sum_tidal+elem_field(ne_Vdot,1)*dt
    else
       sum_expid = sum_expid-elem_field(ne_Vdot,1)*dt
       if(prev_flow.gt.0.0_dp)then
          WOBe_insp = (WOBe+sum_tidal*ptrans_frc*1.0e-9_dp)*(30.0_dp/Tinsp)
          WOBr_insp = WOBr*(30.0_dp/Tinsp)
          WOB_insp = WOBe_insp+WOBr_insp
          WOBe = 0.0_dp
          WOBr = 0.0_dp
       endif
    endif

  end subroutine evaluate_vent_step

!!!#############################################################################

  subroutine evaluate_uniform_flow
    !*evaluate_uniform_flow:* Sets up and solves uniform ventilation model
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_UNIFORM_FLOW" :: EVALUATE_UNIFORM_FLOW
  
    ! Local variables
    integer :: ne,nunit
    real(dp) :: init_vol,volume_tree
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'evaluate_uniform_flow'
    call enter_exit(sub_name,1)

!!! calculate the total model volume
    call volume_of_mesh(init_vol,volume_tree)

!!! initialise the flow field to zero
    elem_field(ne_Vdot,1:num_elems) = 0.0_dp

!!! For each elastic unit, calculate uniform ventilation
    do nunit = 1,num_units
       ne = units(nunit) !local element number
       unit_field(nu_Vdot0,nunit) = unit_field(nu_vol,nunit)/ &
            (init_vol-volume_tree)
       elem_field(ne_Vdot,ne) = unit_field(nu_Vdot0,nunit)
    enddo

    call sum_elem_field_from_periphery(ne_Vdot)

    call enter_exit(sub_name,2)

  end subroutine evaluate_uniform_flow


!!!#############################################################################

  subroutine set_driving_pressures(dpmus,dt,pmus_factor_ex,pmus_factor_in, &
       pmus_step,p_mus,Texpn,Tinsp,ttime,expiration_type)

    real(dp),intent(in) :: dt,pmus_factor_ex,pmus_factor_in,pmus_step,Texpn, &
         Tinsp,ttime
    real(dp) :: dpmus,p_mus
    character(len=*),intent(in) :: expiration_type
    ! Local variables
    real(dp) :: sum_dpmus,sum_dpmus_ei,Tpass
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'set_driving_pressures'
    call enter_exit(sub_name,1)

    select case(expiration_type)
       
    case("active")
       if(ttime.lt.Tinsp)then
          dpmus = pmus_step*pmus_factor_in*PI* &
               sin(pi/Tinsp*ttime)/(2.0_dp*Tinsp)*dt
       elseif(ttime.le.Tinsp+Texpn)then
          dpmus = pmus_step*pmus_factor_ex*PI* &
               sin(2.0_dp*pi*(0.5_dp+(ttime-Tinsp)/(2.0_dp*Texpn)))/ &
               (2.0_dp*Texpn)*dt
       endif
       
    case("passive")
       if(ttime.le.Tinsp+0.5_dp*dt)then
          dpmus = pmus_step*pmus_factor_in*PI*dt* &
               sin(pi*ttime/Tinsp)/(2.0_dp*Tinsp)
          sum_dpmus = sum_dpmus+dpmus
          sum_dpmus_ei = sum_dpmus
       else
          Tpass = 0.1_dp
          dpmus = MIN(-sum_dpmus_ei/(Tpass*Texpn)*dt,-sum_dpmus)
          sum_dpmus = sum_dpmus+dpmus
       endif
       
    end select
    
    p_mus = p_mus + dpmus !current value for muscle pressure

    call enter_exit(sub_name,2)

  end subroutine set_driving_pressures

!!!#############################################################################

  subroutine update_unit_dpdt(dt)
    !*update_unit_dpdt:* updates the rate of change of pressure at the proximal
    ! end of element that supplies tissue unit. i.e. not the rate of change of
    ! pressure within the unit.
    real(dp), intent(in) :: dt
    ! Local variables
    integer :: ne,np1,nunit
    real(dp) :: est
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'update_unit_dpdt'
    call enter_exit(sub_name,1)

    do nunit = 1,num_units
       ne = units(nunit)
       np1 = elem_nodes(1,ne)
       ! linear estimate
       est = (node_field(nj_aw_press,np1) &
            - unit_field(nu_air_press,nunit))/dt
!!!    For stability, weight new estimate with the previous dP/dt
       unit_field(nu_dpdt,nunit) = 0.5_dp*(est+unit_field(nu_dpdt,nunit))
    enddo !nunit

    call enter_exit(sub_name,2)

  end subroutine update_unit_dpdt


!!!#############################################################################

  subroutine update_proximal_pressure
    !*update_proximal_pressure:* Update the pressure at the proximal node of
    ! the element that feeds an elastic unit

    ! Local variables
    integer :: ne,np1,nunit
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'update_proximal_pressure'
    call enter_exit(sub_name,1)

    do nunit = 1,num_units
       ne = units(nunit)
       np1 = elem_nodes(1,ne)
!!!    store the entry node pressure as an elastic unit air pressure
       unit_field(nu_air_press,nunit) = node_field(nj_aw_press,np1) 
    enddo !noelem

    call enter_exit(sub_name,2)

  end subroutine update_proximal_pressure


!!!#############################################################################

  subroutine update_pleural_pressure(ppl_current)
    !*update_pleural_pressure:* Update the mean pleural pressure based on
    ! current Pel (=Ptp) and Palv, i.e. Ppl(unit) = -Pel(unit)+Palv(unit)

    real(dp),intent(out) :: ppl_current
    ! Local variables
    integer :: ne,np2,nunit
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'update_pleural_pressure'
    call enter_exit(sub_name,1)

    ppl_current = 0.0_dp
    do nunit = 1,num_units
       ne = units(nunit)
       np2 = elem_nodes(2,ne)
       ppl_current = ppl_current - unit_field(nu_pe,nunit) + &
            node_field(nj_aw_press,np2)
    enddo !noelem
    ppl_current = ppl_current/num_units

    call enter_exit(sub_name,2)

  end subroutine update_pleural_pressure


!!!#############################################################################

  subroutine update_node_pressures(press_in)
    !*update_node_pressures:* Use the known resistances and flows to calculate
    ! nodal pressures through whole tree

    real(dp),intent(in) :: press_in
    !Local parameters
    integer :: ne,np1,np2
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'update_node_pressures'
    call enter_exit(sub_name,1)

    ! set the initial node pressure to be the input pressure (usually zero)
    ne = 1 !element number at top of tree, usually = 1
    np1 = elem_nodes(1,ne) !first node in element
    node_field(nj_aw_press,np1) = press_in !set pressure at top of tree

    do ne = 1,num_elems !for each element
       np1 = elem_nodes(1,ne) !start node number
       np2 = elem_nodes(2,ne) !end node number
       !P(np2) = P(np1) - Resistance(ne)*Flow(ne)
       node_field(nj_aw_press,np2) = node_field(nj_aw_press,np1) &
            - (elem_field(ne_resist,ne)*elem_field(ne_Vdot,ne))* &
            dble(elem_ordrs(no_type,ne))
    enddo !noelem

    call enter_exit(sub_name,2)

  end subroutine update_node_pressures


!!!#############################################################################

  subroutine tissue_compliance(chest_wall_compliance,undef)

    real(dp), intent(in) :: chest_wall_compliance,undef
    ! Local variables
    integer :: ne,nunit
    real(dp),parameter :: a = 0.433_dp, b = -0.611_dp, cc = 2500.0_dp
    real(dp) :: exp_term,lambda,ratio
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'update_tissue_compliance'
    call enter_exit(sub_name,1)

    !.....dV/dP=1/[(1/2h^2).c/2.(3a+b)exp().(4h(h^2-1)^2)+(h^2+1)/h^2)]

    do nunit = 1,num_units
       ne = units(nunit)
       !calculate a compliance for the tissue unit
       ratio = unit_field(nu_vol,nunit)/undef
       lambda = ratio**(1.0_dp/3.0_dp) !uniform extension ratio
       exp_term = exp(0.75_dp*(3.0_dp*a+b)*(lambda**2-1.0_dp)**2)

       unit_field(nu_comp,nunit) = cc*exp_term/6.0_dp*(3.0_dp*(3.0_dp*a+b)**2 &
            *(lambda**2-1.0_dp)**2/lambda**2+(3.0_dp*a+b) &
            *(lambda**2+1.0_dp)/lambda**4)
       unit_field(nu_comp,nunit) = undef/unit_field(nu_comp,nunit) ! V/P
       ! add the chest wall (proportionately) in parallel
       unit_field(nu_comp,nunit) = 1.0_dp/(1.0_dp/unit_field(nu_comp,nunit)&
            +1.0_dp/(chest_wall_compliance/dble(num_units)))
       !estimate an elastic recoil pressure for the unit
       unit_field(nu_pe,nunit) = cc/2.0_dp*(3.0_dp*a+b)*(lambda**2.0_dp &
            -1.0_dp)*exp_term/lambda
    enddo !nunit

    call enter_exit(sub_name,2)

  end subroutine tissue_compliance


!!!#############################################################################

  subroutine sum_elem_field_from_periphery(ne_field)

    integer,intent(in) :: ne_field
    !Local parameters
    real(dp) :: field_value
    integer :: i,ne,ne2
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'sum_elem_field_from_periphery'
    call enter_exit(sub_name,1)

    do ne = num_elems,1,-1
       if(elem_cnct(1,0,ne).gt.0)then !not terminal
          field_value = 0.0_dp
          do i = 1,elem_cnct(1,0,ne) !for each possible daughter branch (max 2)
             ne2 = elem_cnct(1,i,ne) !the daughter element number
             field_value = field_value+dble(elem_symmetry(ne2))* &
                  elem_field(ne_field,ne2) !sum daughter fields
          enddo !noelem2
          elem_field(ne_field,ne) = field_value
       endif
    enddo !noelem

    call enter_exit(sub_name,2)

  end subroutine sum_elem_field_from_periphery

!!!#############################################################################

  subroutine update_unit_volume(dt)

    real(dp),intent(in) :: dt
    ! Local variables
    integer :: ne,np,nunit
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'update_unit_volume'
    call enter_exit(sub_name,1)

    do nunit = 1,num_units
       ne = units(nunit)
       np = elem_nodes(2,ne)
       ! update the volume of the lumped tissue unit
       unit_field(nu_vol,nunit) = unit_field(nu_vol,nunit)+dt* &
            elem_field(ne_Vdot,ne) !in mm^3
       if(elem_field(ne_Vdot,1).gt.0.0_dp)then  !only store inspired volume
          unit_field(nu_vt,nunit) = unit_field(nu_vt,nunit)+dt* &
               elem_field(ne_Vdot,ne)
       endif
    enddo !nunit

    call enter_exit(sub_name,2)

  end subroutine update_unit_volume

!!!#############################################################################

  subroutine update_elem_field(alpha)

    real(dp),intent(in) :: alpha   ! the factor by which the radius changes
    ! Local variables
    integer :: ne,np1,np2
    real(dp) :: gamma,resistance,reynolds,zeta
    real(dp) :: rad,le
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'update_elem_field'
    call enter_exit(sub_name,1)

    do ne = 1,num_elems
       np1 = elem_nodes(1,ne)
       np2 = elem_nodes(2,ne)

       ! element length
       elem_field(ne_length,ne) = sqrt((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)

       ! element radius
       elem_field(ne_radius,ne) = sqrt(alpha) * elem_field(ne_radius,ne)

       ! element volume
       elem_field(ne_vol,ne) = PI * elem_field(ne_radius,ne)**2 * &
            elem_field(ne_length,ne)
    enddo ! ne
    
    call enter_exit(sub_name,2)
    
  end subroutine update_elem_field

!!!#############################################################################

  subroutine update_resistance

    ! Local variables
    integer :: i,ne,ne2,np1,np2,nunit
    real(dp) :: ett_resistance,gamma,le,rad,resistance,reynolds,sum,zeta
    real(dp) :: tissue_resistance
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'update_resistance'
    call enter_exit(sub_name,1)

    elem_field(ne_t_resist,1:num_elems) = 0.0_dp

    tissue_resistance = 0.0_dp  ! 0.35_dp * 98.0665_dp/1.0e6_dp 

    do nunit = 1,num_units
       ne = units(nunit)
       elem_field(ne_t_resist,ne) = tissue_resistance * dble(elem_units_below(1))
    enddo
    
    do ne = 1,num_elems
       np1 = elem_nodes(1,ne)
       np2 = elem_nodes(2,ne)
       
       le = elem_field(ne_length,ne)
       rad = elem_field(ne_radius,ne)

       ! element Poiseuille (laminar) resistance in units of Pa.s.mm-3   
       resistance = 8.0_dp*GAS_VISCOSITY*elem_field(ne_length,ne)/ &
            (PI*elem_field(ne_radius,ne)**4) !laminar resistance
       
       ! element turbulent resistance (flow in bifurcating tubes)
       gamma = 0.357_dp !inspiration
       if(elem_field(ne_Vdot,ne).lt.0.0_dp) gamma = 0.46_dp !expiration
       
       reynolds = abs(elem_field(ne_Vdot,ne)*2.0_dp*GAS_DENSITY/ &
            (pi*elem_field(ne_radius,ne)*GAS_VISCOSITY))
       zeta = MAX(1.0_dp,dsqrt(2.0_dp*elem_field(ne_radius,ne)* &
            reynolds/elem_field(ne_length,ne))*gamma)
       elem_field(ne_resist,ne) = resistance * zeta
       elem_field(ne_t_resist,ne) = elem_field(ne_resist,ne) + &
            elem_field(ne_t_resist,ne)
    enddo !noelem
    
    do ne = num_elems,1,-1
       sum = 0.0_dp
       if(elem_cnct(1,0,ne).gt.0)then !not terminal
          do i = 1,elem_cnct(1,0,ne) !for each possible daughter branch (max 2)
             ne2 = elem_cnct(1,i,ne) !the daughter element number
             ! line below is sum = sum + 1/R, where 1/R is multiplied by
             !  2 if this is a symmetric child branch
             sum = sum + dble(elem_symmetry(ne2))* &
                  dble(elem_ordrs(no_type,ne2))/elem_field(ne_t_resist,ne2)
          enddo
          if(sum.gt.zero_tol) elem_field(ne_t_resist,ne) = &
               elem_field(ne_t_resist,ne) + 1.0_dp/sum
       endif
    enddo

    call enter_exit(sub_name,2)

  end subroutine update_resistance

!!!#############################################################################

  subroutine estimate_flow(dp_external,dt,err_est)

    real(dp),intent(in) :: dp_external,dt
    real(dp),intent(out) :: err_est
    ! Local variables
    integer :: ne,nunit
    real(dp) :: alpha,beta,flow_diff,flow_sum,Q,Qinit
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'estimate_flow'
    call enter_exit(sub_name,1)

    err_est = 0.0_dp
    flow_sum = 0.0_dp

!!! For each elastic unit, calculate Qbar (equation 4.13 from Swan thesis)
    do nunit = 1,num_units !for each terminal only (with tissue units attached)
       ne = units(nunit) !local element number
       ! Calculate the mean flow into the unit in the time step
       ! alpha is rate of change of pressure at start node of terminal element
       alpha = unit_field(nu_dpdt,nunit) !dPaw/dt, updated each iter
       Qinit = elem_field(ne_Vdot0,ne) !terminal element flow, updated each dt
       ! beta is rate of change of 'external' pressure, incl muscle and entrance
       beta = dp_external/dt ! == dPmus/dt (-ve for insp), updated each dt

!!!    Q = C*(alpha-beta)+(Qinit-C*(alpha-beta))*exp(-dt/(C*R))
       Q = unit_field(nu_comp,nunit)*(alpha-beta)+ &
            (Qinit-unit_field(nu_comp,nunit)*(alpha-beta))* &
            exp(-dt/(unit_field(nu_comp,nunit)*elem_field(ne_t_resist,ne)))

       unit_field(nu_Vdot2,nunit) = unit_field(nu_Vdot1,nunit) !flow at iter-2
       unit_field(nu_Vdot1,nunit) = unit_field(nu_Vdot0,nunit) !flow at iter-1

!!!    for stability the flow estimate for current iteration
!!!    includes flow estimates from previous two iterations
       unit_field(nu_Vdot0,nunit) = 0.75_dp*unit_field(nu_Vdot2,nunit)+ &
            0.25_dp*(Q+unit_field(nu_Vdot1,nunit))*0.5_dp

       flow_diff = unit_field(nu_Vdot0,nunit) - elem_field(ne_Vdot,ne)
       if(abs(flow_diff).gt.zero_tol) &
            err_est = err_est+flow_diff**2 !sum up the error for all elements
       if(abs(unit_field(nu_Vdot0,nunit)).gt.zero_tol) &
            flow_sum = flow_sum+unit_field(nu_Vdot0,nunit)**2
       

!!! ARC: DO NOT CHANGE BELOW. THIS IS NEEDED FOR THE ITERATIVE STEP
!!! - SIMPLER OPTIONS JUST FORCE IT TO CONVERGE WHEN ITS NOT
       elem_field(ne_Vdot,ne) = (unit_field(nu_Vdot0,nunit)&
            +unit_field(nu_Vdot1,nunit))/2.0_dp
       unit_field(nu_Vdot0,nunit) = elem_field(ne_Vdot,ne)
    enddo !nunit

    ! the estimate of error for the iterative solution
    if(abs(flow_sum*dble(num_units)).gt.zero_tol) then
       err_est = err_est/(flow_sum*dble(num_units))
    else
       err_est = err_est/dble(num_units)
    endif

    call enter_exit(sub_name,2)

  end subroutine estimate_flow

!!!#############################################################################

  subroutine calculate_work(breath_vol,dt_vol,WOBe,WOBr,pptrans)

    real(dp) :: breath_vol,dt_vol,WOBe,WOBr,pptrans
    ! Local variables
    integer :: ne,np1,nunit
    real(dp) :: p_resis,p_trans
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'calculate_work'
    call enter_exit(sub_name,1)

    p_resis = 0.0_dp
    !estimate elastic and resistive WOB for each dt (sum dP.V)
    p_trans = SUM(unit_field(nu_pe,1:num_units))/num_units
    do nunit = 1,num_units
       ne = units(nunit)
       np1 = elem_nodes(2,ne)
       p_resis = p_resis+node_field(nj_aw_press,1)-node_field(nj_aw_press,np1)
    enddo
    p_resis=p_resis/num_units
    ! vol in mm3 *1e-9=m3, pressure in Pa, hence *1d-9 = P.m3 (Joules)
    WOBe = WOBe+(p_trans-pptrans)*breath_vol*1.0e-9_dp
    WOBr = WOBr+p_resis*dt_vol*1.0e-9_dp

    pptrans = p_trans

    call enter_exit(sub_name,2)

  end subroutine calculate_work

!!!#############################################################################

  subroutine read_params_main(num_brths, num_itns, dt, err_tol)

    integer,intent(out) :: num_brths, num_itns
    real(dp) :: dt,err_tol

    ! Local variables
    character(len=100) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios
    integer :: line
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'read_params_main'
    call enter_exit(sub_name,1)

    ios = 0
    line = 0
    open(fh, file='Parameters/params_main.txt')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, '    ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          select case (label)
          case ('num_brths')
             read(buffer, *, iostat=ios) num_brths
             print *, 'Read num_brths: ', num_brths
          case ('num_itns')
             read(buffer, *, iostat=ios) num_itns
             print *, 'Read num_itns: ', num_itns
          case ('dt')
             read(buffer, *, iostat=ios) dt
             print *, 'Read dt: ', dt
          case ('err_tol')
             read(buffer, *, iostat=ios) err_tol
             print *, 'Read err_tol: ', err_tol
          case default
             print *, 'Skipping invalid label at line', line
          end select
       end if
    end do

    close(fh)
    call enter_exit(sub_name,2)

  end subroutine read_params_main

!!!#############################################################################

  subroutine read_params_evaluate_flow (gdirn, chest_wall_compliance, &
       constrict, COV, FRC, i_to_e_ratio, pmus_step, press_in,&
       refvol, RMaxMean, RMinMean, T_interval, volume_target, expiration_type)

    integer,intent(out) :: gdirn
    real(dp),intent(out) :: chest_wall_compliance, constrict, COV,&
       FRC, i_to_e_ratio, pmus_step, press_in,&
       refvol, RMaxMean, RMinMean, T_interval, volume_target
    character,intent(out) :: expiration_type*(*)

    ! Local variables
    character(len=100) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios
    integer :: line
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    ios = 0
    line = 0
    sub_name = 'read_params_evaluate_flow'
    call enter_exit(sub_name,1)

    ! following values are examples from control.txt
    !    T_interval = 4.0_dp !s
    !    gdirn = 3
    !    press_in = 0.0_dp !Pa
    !    COV = 0.2_dp
    !    RMaxMean = 1.29_dp
    !    RMinMean = 0.78_dp
    !    i_to_e_ratio = 0.5_dp !dimensionless
    !    refvol = 0.6_dp !dimensionless
    !    volume_target = 8.0e5_dp !mm^3  800 ml
    !    pmus_step = -5.4_dp * 98.0665_dp !-5.4 cmH2O converted to Pa
    !    expiration_type = 'passive' ! or 'active'
    !    chest_wall_compliance = 0.2e6_dp/98.0665_dp !(0.2 L/cmH2O --> mm^3/Pa)

    open(fh, file='Parameters/params_evaluate_flow.txt')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, '    ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          select case (label)
          case ('FRC')
             read(buffer, *, iostat=ios) FRC
             print *, 'Read FRC: ', FRC
          case ('constrict')
             read(buffer, *, iostat=ios) constrict
             print *, 'Read constrict: ', constrict
          case ('T_interval')
             read(buffer, *, iostat=ios) T_interval
             print *, 'Read T_interval: ', T_interval
          case ('Gdirn')
             read(buffer, *, iostat=ios) gdirn
             print *, 'Read Gdirn: ', gdirn
          case ('press_in')
             read(buffer, *, iostat=ios) press_in
             print *, 'Read press_in: ', press_in
          case ('COV')
             read(buffer, *, iostat=ios) COV
             print *, 'Read COV: ', COV
          case ('RMaxMean')
             read(buffer, *, iostat=ios) RMaxMean
             print *, 'Read RMaxMean: ', RMaxMean
          case ('RMinMean')
             read(buffer, *, iostat=ios) RMinMean
             print *, 'Read RMinMean: ', RMinMean
          case ('i_to_e_ratio')
             read(buffer, *, iostat=ios) i_to_e_ratio
             print *, 'Read i_to_e_ratio: ', i_to_e_ratio
          case ('refvol')
             read(buffer, *, iostat=ios) refvol
             print *, 'Read refvol: ', refvol
          case ('volume_target')
             read(buffer, *, iostat=ios) volume_target
             print *, 'Read volume_target: ', volume_target
          case ('pmus_step')
             read(buffer, *, iostat=ios) pmus_step
             print *, 'Read pmus_step_coeff: ', pmus_step
          case ('expiration_type')
             read(buffer, *, iostat=ios) expiration_type
             print *, 'Read expiration_type: ', expiration_type
          case ('chest_wall_compliance')
             read(buffer, *, iostat=ios) chest_wall_compliance
             print *, 'Read chest_wall_compliance: ', chest_wall_compliance
          case default
             print *, 'Skipping invalid label at line', line
          end select
       end if
    end do

    close(fh)
    call enter_exit(sub_name,2)

  end subroutine read_params_evaluate_flow

!!!#############################################################################

  subroutine two_unit_test
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_TWO_UNIT_TEST" :: TWO_UNIT_TEST

    ! Local variables
    integer ne,noelem,nonode,np
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'two_unit_test'
    call enter_exit(sub_name,1)

    ! set up a test geometry. this has only three branches and two elastic units

    num_nodes=4 !four nodes (branch junctions)
    num_elems=3 !three elements (branches)
    num_units = 2

    allocate (nodes(num_nodes))
    allocate (node_xyz(3,num_nodes))
    allocate (node_field(num_nj,num_nodes))
    allocate(elems(num_elems))
    allocate(elem_cnct(-1:1,0:2,0:num_elems))
    allocate(elem_nodes(2,num_elems))
    allocate(elem_ordrs(num_ord,num_elems))
    allocate(elem_symmetry(num_elems))
    allocate(elem_units_below(num_elems))
    allocate(elems_at_node(num_nodes,0:3))
    allocate(elem_field(num_ne,num_elems))
    allocate(elem_direction(3,num_elems))
    allocate(units(num_units))
    allocate(unit_field(num_nu,num_units))

    nodes=0
    node_xyz = 0.0_dp !initialise all values to 0
    node_field = 0.0_dp
    elem_field = 0.0_dp
    unit_field=0.0_dp
    elems=0
    units=0
    elem_cnct = 0

    do nonode=1,num_nodes !loop over all of the nodes
       np=nonode
       nodes(nonode)=np !set local node number same as order in node list
    enddo !nonode

    node_xyz(3,2) = -100.0_dp !setting the z coordinate of node 2
    node_xyz(2,3) = -50.0_dp !setting the y coordinate of node 3
    node_xyz(2,4) = 50.0_dp !setting the y coordinate of node 4
    node_xyz(3,3) = -150.0_dp !setting the z coordinate of node 3
    node_xyz(3,4) = -150.0_dp !setting the z coordinate of node 4

    elem_field(ne_radius,1) = 8.0_dp
    elem_field(ne_radius,2) = 6.0_dp
    elem_field(ne_radius,3) = 5.0_dp

    ! set up elems:
    do noelem=1,num_elems !loop over all of the elements
       ne=noelem
       elems(noelem)=ne !set local elem number same as order in elem list
       elem_nodes(2,noelem)=ne+1
    enddo !noelem
    elem_nodes(1,1)=1
    elem_nodes(1,2)=2
    elem_nodes(1,3)=2

    elem_cnct(-1,0,1:num_elems) = 1 !initialise all branches to have 1 parent
    elem_cnct(-1,0,1) = 0 !element 1 has 0 adjacent branches in -Xi1 direction
    elem_cnct(1,0,1)=2 ! element 1 has 2 adjacent branches in +Xi1 direction
    elem_cnct(1,1,1)=2 !element number of 1st adjacent branch
    elem_cnct(1,2,1)=3 !element number of 2nd adjacent branch
    elem_cnct(-1,1,2)=1 !element number of parent branch
    elem_cnct(-1,1,3)=1 !element number of parent branch

    elem_ordrs(no_gen,1)=1 !branch generation
    elem_ordrs(no_gen,2)=2 !branch generation
    elem_ordrs(no_gen,3)=3 !branch generation
    elem_ordrs(no_Hord,1)=2 !branch Horsfield order
    elem_ordrs(no_Hord,2)=1 !branch Horsfield order
    elem_ordrs(no_Hord,3)=1 !branch Horsfield order
    elem_ordrs(no_Sord,1)=2 !branch Strahler order
    elem_ordrs(no_Sord,2)=1 !branch Strahler order
    elem_ordrs(no_Sord,3)=1 !branch Strahler order

    call append_units

    unit_field(nu_vol,1) = 1.5d6 !arbitrary volume for element 2
    unit_field(nu_vol,2) = 1.5d6 !arbitrary volume for element 3

    elem_units_below(1)=2
    elem_units_below(2)=1
    elem_units_below(3)=1

    elem_symmetry(1:num_elems) = 1

    call enter_exit(sub_name,2)

  end subroutine two_unit_test

!!!#############################################################################

  subroutine write_end_of_breath(init_vol,current_vol,pmus_factor_in, &
       pmus_step,sum_expid,sum_tidal,volume_target,WOBe_insp,WOBr_insp,WOB_insp)

    real(dp),intent(in) :: init_vol,current_vol,pmus_factor_in,pmus_step, &
         sum_expid,sum_tidal,volume_target,WOBe_insp,WOBr_insp,WOB_insp
    ! Local variables
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'write_end_of_breath'
    call enter_exit(sub_name,1)

    write(*,'('' End of breath, inspired = '',F10.2,'' L'')') &
         sum_tidal/1.0e+6_dp
    write(*,'('' End of breath, expired  = '',F10.2,'' L'')') &
         sum_expid/1.0e+6_dp
    write(*,'('' Peak muscle pressure    = '',F10.2,'' cmH2O'')') &
         pmus_step*pmus_factor_in/98.0665_dp
    write(*,'('' Drift in FRC from start = '',F10.2,'' %'')') &
         100*(current_vol-init_vol)/init_vol
    write(*,'('' Difference from target Vt = '',F8.2,'' %'')') &
         100*(volume_target-sum_tidal)/volume_target
    write(*,'('' Total Work of Breathing ='',F7.3,''J/min'')')WOB_insp
    write(*,'('' elastic WOB ='',F7.3,''J/min'')')WOBe_insp
    write(*,'('' resistive WOB='',F7.3,''J/min'')')WOBr_insp
          
    call enter_exit(sub_name,2)

  end subroutine write_end_of_breath

!!!#############################################################################

  subroutine write_flow_step_results(chest_wall_compliance,init_vol, &
       current_vol,ppl_current,pptrans,Pcw,p_mus,time,ttime)

    real(dp),intent(in) :: chest_wall_compliance,init_vol,current_vol, &
         ppl_current,pptrans,Pcw,p_mus,time,ttime
    ! Local variables
    real(dp) :: totalC,Precoil
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'write_flow_step_results'
    call enter_exit(sub_name,1)

    !the total model compliance
    totalC = 1.0_dp/(1.0_dp/sum(unit_field(nu_comp,1:num_units))+ &
         1.0_dp/chest_wall_compliance)
    Precoil = sum(unit_field(nu_pe,1:num_units))/num_units
    
    if(abs(time).lt.zero_tol)then
!!! write out the header information for run-time output
       write(*,'(2X,''Time'',3X,''Inflow'',4X,''V_t'',5X,''Raw'',5X,&
            &''Comp'',4X,''Ppl'',5X,''Ptp'',5X,''VolL'',4X,''Pmus'',&
            &4X,''Pcw'',2X,''Pmus-Pcw'')')
       write(*,'(3X,''(s)'',4X,''(mL/s)'',3X,''(mL)'',1X,''(cmH/L.s)'',&
            &1X,''(L/cmH)'',1X,''(...cmH2O...)'',&
            &4X,''(L)'',5X,''(......cmH2O.......)'')')
       
       write(*,'(F7.3,2(F8.1),8(F8.2))') &
            0.0_dp,0.0_dp,0.0_dp, &  !time, flow, tidal
            elem_field(ne_t_resist,1)*1.0e+6_dp/98.0665_dp, & !res (cmH2O/L.s)
            totalC*98.0665_dp/1.0e+6_dp, & !total model compliance
            ppl_current/98.0665_dp, & !Ppl (cmH2O)
            -ppl_current/98.0665_dp, & !mean Ptp (cmH2O)
            init_vol/1.0e+6_dp, & !total model volume (L)
            0.0_dp, & !Pmuscle (cmH2O)
            Pcw/98.0665_dp, & !Pchest_wall (cmH2O)
            (-Pcw)/98.0665_dp !Pmuscle - Pchest_wall (cmH2O)
    else
       write(*,'(F7.3,2(F8.1),8(F8.2))') &
            time, & !time through breath (s)
            elem_field(ne_Vdot,1)/1.0e+3_dp, & !flow at the inlet (mL/s)
            (current_vol - init_vol)/1.0e+3_dp, & !current tidal volume (mL)
            elem_field(ne_t_resist,1)*1.0e+6_dp/98.0665_dp, & !res (cmH2O/L.s)
            totalC*98.0665_dp/1.0e+6_dp, & !total model compliance
            ppl_current/98.0665_dp, & !Ppl (cmH2O)
            pptrans/98.0665_dp, & !mean Ptp (cmH2O)
            current_vol/1.0e+6_dp, & !total model volume (L)
            p_mus/98.0665_dp, & !Pmuscle (cmH2O)
            -Pcw/98.0665_dp, & !Pchest_wall (cmH2O)
            (p_mus+Pcw)/98.0665_dp !Pmuscle - Pchest_wall (cmH2O)
       
    endif

    call enter_exit(sub_name,2)

  end subroutine write_flow_step_results

!!!#############################################################################

  function ventilation_continue(n,num_brths,sum_tidal,volume_target)

    integer,intent(in) :: n,num_brths
    real(dp),intent(in) :: sum_tidal,volume_target
    ! Local variables
    logical :: ventilation_continue

    ! --------------------------------------------------------------------------

    ventilation_continue = .true.
    if(n.ge.num_brths)then
       ventilation_continue = .false.
    elseif(abs(volume_target).gt.1.0e-3_dp)then
       if(abs(100.0_dp*(volume_target-sum_tidal) &
            /volume_target).gt.0.1_dp.or.(n.lt.2))then
          ventilation_continue = .true.
       else
          ventilation_continue = .false.
       endif
    endif

  end function ventilation_continue

!!!#############################################################################

end module ventilation
