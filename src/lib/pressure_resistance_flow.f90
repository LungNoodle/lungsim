module pressure_resistance_flow
!*Brief Description:* This module contains tools that are used to solve systems of equations representing steady pressure, resistance and flow problems in any branching geometry. The subroutines in this module are core subroutines that are used in many problem types and are applicable beyond lung modelling
!
!*LICENSE:*
!TBC
!
!
!*Full Description:*
!
!This module contains tools that are used to solve systems of equations representing steady pressure, resistance and flow problems in any branching geometry. The subroutines in this module are core subroutines that are used in many problem types and are applicable beyond lung modelling

  use arrays
  use capillaryflow
  use diagnostics
  use indices
  use other_consts
  use solve, only: BICGSTAB_LinSolv,pmgmres_ilu_cr

  implicit none
  integer,allocatable :: occ_list(:),compl_list(:)
  integer,allocatable :: perfusion_flow_elements(:)
  real(dp),allocatable :: perfusion_flow_values(:)
  ! Prescribed perfusion flows are imposed by adjusting vessel resistances.
  ! ctrl_group(ne): -1 = not controlled, 0 = independently controlled element,
  ! >0 = controlled child of junction node ctrl_group(ne) (a 'junction group').
  character(len=20) :: perfusion_flow_mode = 'geometric'
  integer,allocatable :: ctrl_group(:)
  logical,allocatable :: ctrl_prescribed(:)
  real(dp),allocatable :: ctrl_target(:),ctrl_logscale(:),ctrl_mult(:),resist_base(:)
  integer :: num_ctrl = 0
  real(dp),parameter :: ctrl_flow_tol = 1.0e-5_dp
  integer,parameter :: ctrl_max_iterations = 500

  !Module parameters

  !Module types

  !Module depvar

  !Interfaces
  private
  public evaluate_prq,calculate_ppl,occlusion_list,compliance_list
  public set_perfusion_flows,clear_perfusion_flows,set_perfusion_flow_mode
contains
!###################################################################################
!
!*evaluate_PRQ:* Solves for pressure and flow in a rigid or compliant tree structure
  subroutine evaluate_prq(mesh_type,vessel_type,grav_dirn,grav_factor,bc_type,inlet_bc,outlet_bc,remodeling_grade)

    !local variables
    integer :: mesh_dof,depvar_types
    integer, allocatable :: mesh_from_depvar(:,:,:)
    integer, allocatable :: depvar_at_node(:,:,:)
    integer, allocatable :: depvar_at_elem(:,:,:)
    integer, dimension(0:2,2) :: depvar_totals
    integer, allocatable :: SparseCol(:)
    integer, allocatable :: SparseRow(:)
    integer, allocatable :: update_resistance_entries(:)
    integer, allocatable :: update_fixed_flow_rhs_rows(:)
    real(dp), allocatable :: fixed_flow_rhs_base(:)
    real(dp), allocatable :: SparseVal(:)
    real(dp), allocatable :: RHS(:)
    integer :: num_vars,NonZeros,MaxNonZeros,MatrixSize
    integer :: AllocateStatus

    real(dp), intent(in) :: remodeling_grade ! Remodeling if applicable, 0 stands for healthy and anything from 10 to 100 is hypertension
    real(dp), allocatable :: prq_solution(:,:),solver_solution(:)
    real(dp) :: viscosity,density,inlet_bc,outlet_bc,inletbc,outletbc,grav_vect(3),gamma,total_resistance,ERR
    logical, allocatable :: FIX(:)
    logical :: ADD=.FALSE.,CONVERGED=.FALSE.,flow_control,flow_ok
    real(dp) :: flow_err
    character(len=60) :: sub_name,mesh_type,vessel_type,mechanics_type,bc_type
    integer :: grav_dirn,no,depvar,KOUNT,nz,ne,SOLVER_FLAG,ne0,ne1,nj
    real(dp) :: MIN_ERR,N_MIN_ERR,elasticity_parameters(3),mechanics_parameters(2),grav_factor,P1
    real(dp) :: P2,Q01,Rin,Rout,x_cap,y_cap,z_cap,Ppl,LPM_R,Lin,Lout

    sub_name = 'evaluate_prq'
    call enter_exit(sub_name,1)
    ! ADD and CONVERGED are initialised in their declarations (implicit SAVE),
    ! so reset them to allow repeated solves in the same process.
    ADD=.FALSE.
    CONVERGED=.FALSE.
!!---------DESCRIPTION OF MODEL Types -----------
!mesh_type: can be simple_tree, full_plus_ladder, full_sheet, full_tube The first can be airways, arteries, veins but no special features at the terminal level, the last one has arteries and veins connected by capillary units of some type (lung ladder acinus, lung sheet capillary bed, capillaries are just tubes represented by an element)

!vessel_type:
  !rigid, no elasticity, no parameters required
  !elastic_g0_beta, R=R0*((Ptm/G0)+1.d0)^(1.d0/elasticity_parameters(2)),with an optional maximum pressure beyond which the vessel radius is constant three parameters, g0, elasticity_parameters(2), elasticity_parameters(3)
  !elastic alpha,  R=R0*(alpha*Ptm+1.d0), up to a limit elasticity_parameters(3) two parameters alpha, elasticity_parameters(3)
  !elastic_hooke, two parameters E and h,R=R0+3.0_dp*R0**2*Ptm/(4.0_dp*E*h*R0)

!mechanics type:
  !linear two parmeters, transpulmonary pressure (average) and pleural density (gradient)
  !mechanics, two parameters, pressure and stretch fields

!bc_type:
    !pressure (at inlet and outlets)
    !flow (flow at inlet pressure at outlet).


mechanics_type='linear'

if (vessel_type.eq.'rigid') then
    elasticity_parameters=0.0_dp
elseif (vessel_type.eq.'elastic_g0_beta') then
    elasticity_parameters(1)=6.67e3_dp!G0 (Pa)
    elasticity_parameters(2)=1.0_dp!elasticity_parameters(2)
    elasticity_parameters(3)=32.0_dp*98.07_dp !elasticity_parameters(3) (Pa)
elseif (vessel_type.eq.'elastic_alpha') then
    elasticity_parameters(1)=1.503e-4_dp!alpha (1/Pa)
    elasticity_parameters(2)=32.0_dp*98.07_dp !elasticity_parameters(3) (Pa)
    elasticity_parameters(3)=0.0_dp !Not used
elseif (vessel_type.eq.'elastic_hooke') then
     elasticity_parameters(1)=1.5e6_dp !Pa
     elasticity_parameters(2)=0.1_dp!this is a fraction of the radius so is unitless
     elasticity_parameters(3)=0.0_dp !Not used
else
    print *, 'WARNING: Your chosen vessel type does not seem to be implemented assuming rigid'
    vessel_type='rigid'
    elasticity_parameters=0.0_dp
endif

if (mechanics_type.eq.'linear') then
    mechanics_parameters(1)=5.0_dp*98.07_dp !average pleural pressure (Pa)
    mechanics_parameters(2)=0.25_dp*0.1e-2_dp !pleural density, defines gradient in pleural pressure
else
    print *, 'ERROR: Only linear mechanics models have been implemented to date,assuming default parameters'
     call exit(0)
endif

grav_vect=0.d0
if (grav_dirn.eq.1) then
    grav_vect(1)=1.0_dp
elseif (grav_dirn.eq.2) then
    grav_vect(2)=1.0_dp
elseif (grav_dirn.eq.3) then
    grav_vect(3)=1.0_dp
else
     print *, "ERROR: Posture not recognised (currently only x=1,y=2,z=3))"
     call exit(0)
endif
grav_vect=grav_vect*grav_factor

    if(bc_type.eq.'pressure')then
    inletbc=inlet_bc
    outletbc=outlet_bc
elseif(bc_type.eq.'flow')then
    inletbc=inlet_bc
    outletbc=outlet_bc
    elseif((bc_type.NE.'pressure').AND.(bc_type.NE.'flow'))then
      print *,"unsupported bc_type",bc_type
      call exit(1)
    endif

!!---------PHYSICAL PARAMETERS-----------
!viscosity: fluid viscosity
!density:fluid density
!gamma:Pedley correction factor
density=0.10500e-02_dp !kg/cm3
viscosity=0.33600e-02_dp !Pa.s
gamma = 0.327_dp !=1.85/(4*sqrt(2))

!! Allocate memory to depvar arrays
    mesh_dof=num_elems+num_nodes
    depvar_types=2 !pressure/flow
    allocate (mesh_from_depvar(0:2,mesh_dof,0:2), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for mesh_from_depvar array ***"
    allocate (depvar_at_elem(0:2,depvar_types,num_elems), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for depvar_at_elem array ***"
    allocate (depvar_at_node(num_nodes,0:2,depvar_types), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for depvar_at_node array ***"
    allocate (prq_solution(mesh_dof,2), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for prq_solution array ***"
    prq_solution=0.0_dp !initialise
    allocate (FIX(mesh_dof), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for FIX array ***"



!! Setting up mappings between nodes, elements and solution depvar
    call calc_depvar_maps(mesh_from_depvar,depvar_at_elem,&
                depvar_totals,depvar_at_node,mesh_dof,num_vars)

!! Define boundary conditions
    !first call to define inlet boundary conditions
    call boundary_conditions(ADD,FIX,bc_type,grav_vect,density,inletbc,outletbc,&
      depvar_at_node,depvar_at_elem,prq_solution,mesh_dof,mesh_type)
    !second call if simple tree need to define pressure bcs at all terminal branches
    if(mesh_type.eq.'simple_tree')then
        ADD=.TRUE.
        call boundary_conditions(ADD,FIX,bc_type,grav_vect,density,inletbc,outletbc,&
            depvar_at_node,depvar_at_elem,prq_solution,mesh_dof,mesh_type)
    elseif(mesh_type.eq.'full_plus_ladder')then
        ADD=.TRUE.
        call boundary_conditions(ADD,FIX,bc_type,grav_vect,density,inletbc,outletbc,&
            depvar_at_node,depvar_at_elem,prq_solution,mesh_dof,mesh_type)
    endif

    ! Prescribed flows (if any) are matched by scaling vessel resistances.
    call setup_flow_control(bc_type,inletbc,flow_control)
    flow_ok=.TRUE.
    flow_err=0.0_dp

    KOUNT=0
!! Calculate resistance of each element
   call calculate_resistance(viscosity,KOUNT)
   if(flow_control)call store_and_scale_resistance()

!! Calculate sparsity structure for solution matrices
    !Determine size of and allocate solution vectors/matrices
    call calc_sparse_size(mesh_dof,depvar_at_elem,depvar_at_node,FIX,MaxNonZeros,MatrixSize)
    allocate (SparseCol(MaxNonZeros), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for SparseCol array ***"
    allocate (SparseRow(MatrixSize+1), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for SparseRow array ***"
    allocate (SparseVal(MaxNonZeros), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for SparseVal array ***"
    allocate (RHS(MatrixSize), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for RHS array ***"
    allocate (solver_solution(MatrixSize), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (update_resistance_entries(num_elems), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    update_resistance_entries = 0
    allocate (update_fixed_flow_rhs_rows(num_elems), fixed_flow_rhs_base(num_elems), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for fixed-flow RHS maps ***"
    update_fixed_flow_rhs_rows = 0
    fixed_flow_rhs_base = 0.0_dp
    !calculate the sparsity structure
    call calc_sparse_1dtree(bc_type,density,FIX,grav_vect,mesh_dof,depvar_at_elem, &
        depvar_at_node,MaxNonZeros,MatrixSize,SparseCol,SparseRow,SparseVal,RHS, &
        prq_solution,update_resistance_entries,update_fixed_flow_rhs_rows, &
        fixed_flow_rhs_base,NonZeros)
!!! --ITERATIVE LOOP--
    MIN_ERR=1.d10
    N_MIN_ERR=0
    do while(.NOT.CONVERGED)
      KOUNT=KOUNT+1
      print*, 'Outer loop iterations:',KOUNT
!!! Initialise solution vector based on bcs and rigid vessel resistance
      if(KOUNT.eq.1)then!set up boundary conditions
        if(bc_type.eq.'pressure')then
          if(mesh_type.eq.'full_plus_ladder')then
            total_resistance=1000.0_dp
          else
            call tree_resistance(total_resistance)
          endif
          call initialise_solution(inletbc,outletbc,(inletbc-outletbc)/total_resistance, &
              mesh_dof,prq_solution,depvar_at_node,depvar_at_elem,FIX)
          !move initialisation to solver solution (skipping BCs).
          no=0
          do depvar=1,mesh_dof !loop over mesh dofs
            if(.NOT.FIX(depvar))then
              no=no+1
              solver_solution(no)=prq_solution(depvar,1)
            endif
          enddo !mesh_dof
        else!flow BCs to be implemented
        endif
      else!Need to update just the resistance values in the solution matrix
        do ne=1,num_elems !update for all ne
          if(update_resistance_entries(ne).gt.0)then
            nz=update_resistance_entries(ne)
            SparseVal(nz)=-elem_field(ne_resist,ne) !Just updating resistance
          endif
        enddo
        if(bc_type.eq.'flow')then
          do ne=1,num_elems
            if(update_fixed_flow_rhs_rows(ne).gt.0)then
              depvar = depvar_at_elem(1,1,ne)
              RHS(update_fixed_flow_rhs_rows(ne)) = fixed_flow_rhs_base(ne) + &
                prq_solution(depvar,1)*elem_field(ne_resist,ne)
            endif
          enddo
          !SparseVal(nz)=-elem_field(ne_resist,ne) !Just updating resistance
        endif
      endif!first or subsequent iteration
!! ----CALL SOLVER----
      call pmgmres_ilu_cr(MatrixSize, NonZeros, SparseRow, SparseCol, SparseVal, &
         solver_solution, RHS, 500, 500,1.d-5,1.d-4,SOLVER_FLAG)
       if(SOLVER_FLAG == 0)then
          print *, 'Warning: pmgmres has reached max iterations. Solution may not be valid if this warning persists'
       elseif(SOLVER_FLAG ==2)then
          print *, 'ERROR: pmgmres has failed to converge'
          deallocate (mesh_from_depvar, STAT = AllocateStatus)
          deallocate (depvar_at_elem, STAT = AllocateStatus)
          deallocate (depvar_at_node, STAT = AllocateStatus)
          deallocate (prq_solution, STAT = AllocateStatus)
          deallocate (FIX, STAT = AllocateStatus)
          deallocate (solver_solution, STAT = AllocateStatus)
          deallocate (SparseCol, STAT = AllocateStatus)
          deallocate (SparseVal, STAT = AllocateStatus)
          deallocate (SparseRow, STAT = AllocateStatus)
          deallocate (RHS, STAT = AllocateStatus)
          deallocate (update_resistance_entries, STAT=AllocateStatus)
          exit
       endif
!!--TRANSFER SOLVER SOLUTIONS TO FULL SOLUTIONS
      ERR=0.0_dp
      no=0
      do depvar=1,mesh_dof
         if(.NOT.FIX(depvar)) THEN
            no=no+1
            prq_solution(depvar,2)=prq_solution(depvar,1) !temp storage of previous solution
            prq_solution(depvar,1)=solver_solution(no) !new pressure & flow solutions
            if(DABS(prq_solution(depvar,1)).GT.0.d-6)THEN
               ERR=ERR+(prq_solution(depvar,2)-prq_solution(depvar,1))**2.d0/prq_solution(depvar,1)**2
            endif
         endif
      enddo !no2
!rigid vessels no need to update - tag as converged and exit
      if(vessel_type.eq.'rigid')then
        ERR=0.0_dp
        CONVERGED=.TRUE.
        if(flow_control)then
          ! Resistances only change through the flow-control multipliers.
          call update_flow_control(prq_solution,depvar_at_elem,mesh_dof,flow_err)
          flow_ok=(flow_err.le.ctrl_flow_tol)
          print *,"Prescribed flow mismatch (max relative) =",flow_err
          if(.not.flow_ok)then
            call rescale_resistance()
            CONVERGED=.FALSE.
          endif
          if((.not.flow_ok).and.(KOUNT.ge.ctrl_max_iterations))then
            print *,'WARNING: prescribed flows not matched after',KOUNT,'iterations'
            CONVERGED=.TRUE.
          endif
        endif
      else
!Update vessel radii based on predicted pressures and then update resistance through tree
        call calc_press_area(grav_vect,KOUNT,depvar_at_node,prq_solution,&
           mesh_dof,vessel_type,elasticity_parameters,mechanics_parameters,remodeling_grade)
        call calculate_resistance(viscosity,KOUNT)

!Put the ladder stuff here --> See solve11.f
         if(mesh_type.eq.'full_plus_ladder')then
           do ne=1,num_elems
              if(elem_field(ne_group,ne).eq.1.0_dp)then!(elem_field(ne_group,ne)-1.0_dp).lt.TOLERANCE)then
                ne0=elem_cnct(-1,1,ne)!upstream element number
                ne1=elem_cnct(1,1,ne)
                P1=prq_solution(depvar_at_node(elem_nodes(2,ne0),0,1),1) !pressure at start node of capillary element
                P2=prq_solution(depvar_at_node(elem_nodes(1,ne1),0,1),1)!pressure at end node of capillary element
                Q01=prq_solution(depvar_at_elem(1,1,ne0),1) !flow in element upstream of capillary element !mm^3/s
                Rin=elem_field(ne_radius_out0,ne0)!radius of upstream element
                Rout=elem_field(ne_radius_out0,ne1) !radius of downstream element
                x_cap=node_xyz(1,elem_nodes(1,ne))
                y_cap=node_xyz(2,elem_nodes(1,ne))
                z_cap=node_xyz(3,elem_nodes(1,ne))
                call calculate_ppl(elem_nodes(1,ne),grav_vect,mechanics_parameters,Ppl)
                Lin=elem_field(ne_length,ne0)
                Lout=elem_field(ne_length,ne1)
                 call cap_flow_ladder(ne,LPM_R,Lin,Lout,P1,P2,&
                        Ppl,Q01,Rin,Rout,x_cap,y_cap,z_cap,&
                        .FALSE.)
                 elem_field(ne_resist,ne)=LPM_R
              endif
           enddo
         endif
         if(flow_control)then
           call update_flow_control(prq_solution,depvar_at_elem,mesh_dof,flow_err)
           flow_ok=(flow_err.le.ctrl_flow_tol)
           call store_and_scale_resistance()
           print *,"Prescribed flow mismatch (max relative) =",flow_err
           if((.not.flow_ok).and.(KOUNT.ge.ctrl_max_iterations))then
             print *,'WARNING: prescribed flows not matched after',KOUNT,'iterations'
             flow_ok=.TRUE.
           endif
         endif

         ERR=ERR/MatrixSize !sum of error divided by no of unknown depvar
         if(ERR.LE.1.d-6.AND.(KOUNT.NE.1).AND.flow_ok)then
           CONVERGED=.TRUE.
            print *,"Convergence achieved after",KOUNT,"iterations",ERR
         else !if error not converged
            if(ERR.GE.MIN_ERR) then
              N_MIN_ERR=N_MIN_ERR+1
            else
              MIN_ERR=ERR
            endif
            print *,"Not converged, error =",ERR
         endif !ERR not converged
      endif!vessel type
    enddo !notconverged

!need to write solution to element/nodal fields for export
    call map_solution_to_mesh(prq_solution,depvar_at_elem,depvar_at_node,mesh_dof)
    if(flow_control)then
      call report_flow_control()
      call free_flow_control()
    endif
    !NEED TO UPDATE TERMINAL SOLUTION HERE. LOOP THO' UNITS AND TAKE FLOW AND PRESSURE AT TERMINALS
    call map_flow_to_terminals
    !EXPORT LADDER SOLUTION
    if(mesh_type.eq.'full_plus_ladder')then
      open(10, file='micro_flow_ladder.out', status='replace')
      open(20, file='micro_flow_unit.out', status='replace')
      do ne=1,num_elems
        if(elem_field(ne_group,ne).eq.1.0_dp)then!(elem_field(ne_group,ne)-1.0_dp).lt.TOLERANCE)then
          ne0=elem_cnct(-1,1,ne)!upstream element number
          ne1=elem_cnct(1,1,ne)
          P1=prq_solution(depvar_at_node(elem_nodes(2,ne0),0,1),1) !pressure at start node of capillary element
          P2=prq_solution(depvar_at_node(elem_nodes(1,ne1),0,1),1)!pressure at end node of capillary element
          Q01=prq_solution(depvar_at_elem(1,1,ne0),1) !flow in element upstream of capillary element !mm^3/s
          Rin=elem_field(ne_radius_out0,ne0)!radius of upstream element
          Rout=elem_field(ne_radius_out0,ne1) !radius of downstream element
          x_cap=node_xyz(1,elem_nodes(1,ne))
          y_cap=node_xyz(2,elem_nodes(1,ne))
          z_cap=node_xyz(3,elem_nodes(1,ne))
          call calculate_ppl(elem_nodes(1,ne),grav_vect,mechanics_parameters,Ppl)
          Lin=elem_field(ne_length,ne0)
          Lout=elem_field(ne_length,ne1)
          call cap_flow_ladder(ne,LPM_R,Lin,Lout,P1,P2,&
            Ppl,Q01,Rin,Rout,x_cap,y_cap,z_cap,&
            .TRUE.)
        endif
      enddo
      close(10)
      close(20)
    endif

    deallocate (mesh_from_depvar, STAT = AllocateStatus)
    deallocate (depvar_at_elem, STAT = AllocateStatus)
    deallocate (depvar_at_node, STAT = AllocateStatus)
    deallocate (prq_solution, STAT = AllocateStatus)
    deallocate (FIX, STAT = AllocateStatus)
    deallocate (solver_solution, STAT = AllocateStatus)
    deallocate (SparseCol, STAT = AllocateStatus)
    deallocate (SparseVal, STAT = AllocateStatus)
    deallocate (SparseRow, STAT = AllocateStatus)
    deallocate (RHS, STAT = AllocateStatus)
    deallocate (update_resistance_entries, STAT=AllocateStatus)
    deallocate (update_fixed_flow_rhs_rows, fixed_flow_rhs_base, STAT=AllocateStatus)
    call enter_exit(sub_name,2)
  end subroutine evaluate_prq
!
!###################################################################################
!
!*boundary_conditions:* Defines boundary conditions for prq problems
 subroutine boundary_conditions(ADD,FIX,bc_type,grav_vect,density,inletbc,outletbc,&
       depvar_at_node,depvar_at_elem,prq_solution,mesh_dof,mesh_type)

    integer :: mesh_dof
    integer :: depvar_at_elem(0:2,2,num_elems)
    integer :: depvar_at_node(num_nodes,0:2,2)
    real(dp) :: prq_solution(mesh_dof,2),inletbc,outletbc,density,grav_vect(3)
    logical:: ADD
    logical :: FIX(mesh_dof)
    character(len=60) ::bc_type,mesh_type

  ! local variables
    integer :: nonode,np,ne,ny1,nj,np_in
    real(dp) :: grav
    character(len=60) :: sub_name

  sub_name = 'boundary_conditions'
  call enter_exit(sub_name,1)
  if(.NOT.ADD)THEN
     ! Initial values
     FIX(1:mesh_dof)=.FALSE.
     prq_solution = 0
     ! Fixed boundary conditions
     ! These are inlet BCs, apply to all inlet BCs (there should only be one)
     do ne=1,num_elems
        !ne=elems(noelem)
        if (elem_cnct(-1,0,ne) == 0) THEN !Entry element
           if(BC_TYPE == 'pressure')THEN
              np=elem_nodes(1,ne)
              ny1=depvar_at_node(np,1,1) !for fixed pressure BC
              FIX(ny1)=.TRUE. !set fixed
              prq_solution(ny1,1)=inletbc !Putting BC value into solution array
              np_in=np !inlet node set here, gravity reference to this point
           else if(BC_TYPE == 'flow')THEN
              ny1=depvar_at_elem(0,1,ne) !fixed
              FIX(ny1)=.TRUE. !set fixed
              prq_solution(ny1,1)=inletbc !Putting BC value into solution array
              np_in=elem_nodes(1,ne)
           endif
        endif
     enddo
  else !Add terminal pressure BC for all terminal branches
    if(mesh_type.eq.'simple_tree')then !NEED TO SET GRAVITY IN THIS CASE
      do nonode=1,num_units
         np=elem_nodes(2,units(nonode)) !Second node in element = terminal node
         ny1=depvar_at_node(np,1,1) !for fixed pressure BC
         FIX(ny1)=.TRUE. !set fixed
!! NB// Add gravitational factor in here
         if(np_in.eq.0) then
            print *,"Warning --> np_in is not set yet, setting to first node as default"
            np_in=1 !Setting to first node as default
         endif
         grav=0.d0
         do nj=1,3
            grav=grav+density*grav_vect(nj)*9810.d0*(node_xyz(nj,np_in)-node_xyz(nj,np))
         enddo
         prq_solution(ny1,1)=outletbc-grav !Putting BC value into solution array
         !print *,"BC----",ny1,outletbc,grav,prq_solution(ny1,1)
       enddo
     else
       do ne=1,num_elems
        !ne=elems(noelem)
        if (elem_cnct(1,0,ne) == 0) THEN !EXIT ELEMENT
              np=elem_nodes(2,ne)
              ny1=depvar_at_node(np,1,1) !for fixed pressure BC
              FIX(ny1)=.TRUE. !set fixed
              prq_solution(ny1,1)=outletbc !Putting BC value into solution array
              np_in=np !inlet node set here, gravity reference to this point
          endif
        enddo

     endif
  endif
    call enter_exit(sub_name,2)
  end subroutine boundary_conditions

!
!###################################################################################
!
  !* occlusion_list gets an element list from user which includes the element that were
  ! occluded for CTEPH and now the occlusions are removed and we want to make sure there
  ! is no remodelling downstream of those elements
  subroutine occlusion_list(elem_list)

    integer,intent(in) :: elem_list(:) ! list of surface elements defining the host region
  !local variables
    integer :: i
    character(len=60) :: sub_name

    sub_name = 'occlusion_list'
    call enter_exit(sub_name,1)

    if(count(elem_list.ne.0).gt.0)then ! a surface element list is given for converting to
       !                                 create a list of occluded elements
       allocate(occ_list(count(elem_list.ne.0)))
       !!! get the list of occlusion list from the given list
       do i = 1,count(elem_list.ne.0)
          occ_list(i) = elem_list(i)
       enddo
    endif
    !
    ! write(*,*) 'occlusion_list:', occ_list
    ! pause


    call enter_exit(sub_name,2)
    end subroutine occlusion_list
!
!###################################################################################
!
  !* compliance_list gets an element list from user which includes the element that were
  ! occluded for CTEPH the compliance of these vessels is reduced to zero - basically make these
  ! vessels rigid
  subroutine compliance_list(elem_list2)

    integer,intent(in) :: elem_list2(:) ! list of surface elements defining the host region
  !local variables
    integer :: i
    character(len=60) :: sub_name

    sub_name = 'compliance_list'
    call enter_exit(sub_name,1)

    if(count(elem_list2.ne.0).gt.0)then ! a surface element list is given for converting to
       !                                 create a list of occluded elements
       allocate(compl_list(count(elem_list2.ne.0)))
       !!! get the list of occlusion list from the given list
       do i = 1,count(elem_list2.ne.0)
          compl_list(i) = elem_list2(i)
       enddo
    endif

    call enter_exit(sub_name,2)
  end subroutine compliance_list

!###################################################################################
!*set_perfusion_flows:* Opt in to prescribed flows for selected elements.  A
! zero-length list disables the option and restores the original PRQ behaviour.
  subroutine set_perfusion_flows(element_ids,flow_values)
    integer,intent(in) :: element_ids(:)
    real(dp),intent(in) :: flow_values(:)
    integer :: allocate_status

    call clear_perfusion_flows()
    if(size(element_ids).ne.size(flow_values))then
      print *, 'ERROR: prescribed perfusion element and flow lists must have the same length'
      call exit(1)
    endif
    if(size(element_ids).eq.0)return

    allocate(perfusion_flow_elements(size(element_ids)), &
      perfusion_flow_values(size(flow_values)),stat=allocate_status)
    if(allocate_status.ne.0)stop '*** Not enough memory for prescribed perfusion flows ***'
    perfusion_flow_elements=element_ids
    perfusion_flow_values=flow_values
  end subroutine set_perfusion_flows

!###################################################################################
  subroutine clear_perfusion_flows()
    if(allocated(perfusion_flow_elements))deallocate(perfusion_flow_elements)
    if(allocated(perfusion_flow_values))deallocate(perfusion_flow_values)
  end subroutine clear_perfusion_flows

!###################################################################################
!*set_perfusion_flow_mode:* Chooses how resistances are shared between the
! children of a junction whose flows are all known (see setup_flow_control).
!   'geometric' (default): the children's resistance multipliers have a
!                geometric mean of one (smallest symmetric change).
!   'constrict': resistances only increase (smallest multiplier is one).
!   'dilate':    resistances only decrease (largest multiplier is one).
  subroutine set_perfusion_flow_mode(mode)
    character(len=*),intent(in) :: mode

    select case(trim(adjustl(mode)))
    case('geometric','constrict','dilate')
      perfusion_flow_mode=trim(adjustl(mode))
    case default
      print *, 'ERROR: unknown prescribed perfusion flow mode: ',trim(mode)
      print *, '       expected geometric, constrict or dilate'
      call exit(1)
    end select
  end subroutine set_perfusion_flow_mode

!###################################################################################
!*elem_root_node:* The node at which an element joins its parent (trunk) vessel:
! the proximal node for arteries/capillaries, the distal node for veins.
  integer function elem_root_node(ne)
    integer,intent(in) :: ne

    if(nint(elem_field(ne_group,ne)).eq.2)then
      elem_root_node=elem_nodes(2,ne)
    else
      elem_root_node=elem_nodes(1,ne)
    endif
  end function elem_root_node

!###################################################################################
!*elem_trunk:* The parent element of ne (towards the inlet for arteries, towards
! the outlet for veins) within the same vessel group, or 0 if there is none.
  integer function elem_trunk(ne)
    integer,intent(in) :: ne
    integer :: ne0

    elem_trunk=0
    if(nint(elem_field(ne_group,ne)).eq.2)then
      if(elem_cnct(1,0,ne).ge.1)then
        ne0=elem_cnct(1,1,ne)
      else
        return
      endif
    else
      if(elem_cnct(-1,0,ne).ge.1)then
        ne0=elem_cnct(-1,1,ne)
      else
        return
      endif
    endif
    if(nint(elem_field(ne_group,ne0)).eq.nint(elem_field(ne_group,ne)))elem_trunk=ne0
  end function elem_trunk

!###################################################################################
!*setup_flow_control:* Prepares resistance control for prescribed flows.
!
! A prescribed flow is imposed by multiplying the resistance of every vessel in
! the subtree fed by the prescribed element (down to the next controlled element)
! by a factor that is found iteratively; Poiseuille's law is retained in every
! element.  Multipliers are nested, so a controlled element inside a controlled
! subtree redistributes flow within it.
!
! Flows that follow from conservation are found first, propagating from the inlet
! (flow BCs) towards the periphery.  At a junction where the parent flow is known
! and every child flow is known (prescribed or implied), all children become
! controlled: this makes the result independent of which children were named
! (e.g. prescribing RPA=70% or LPA=30% of a known MPA flow is the same problem).
! Such a junction has one more multiplier than independent constraints; the
! extra degree of freedom is fixed by perfusion_flow_mode.  Any other prescribed
! element is controlled independently by its own multiplier.
  subroutine setup_flow_control(bc_type,inletbc,flow_control)
    character(len=60),intent(in) :: bc_type
    real(dp),intent(in) :: inletbc
    logical,intent(out) :: flow_control

    integer :: no,ne,ne2,np,k,n_unknown,ne_unknown,num_inlets,num_outlets,allocate_status
    integer :: n_trunk,n_child
    logical :: changed,all_known,any_prescribed
    logical,allocatable :: known(:)
    real(dp),allocatable :: flow_value(:)
    real(dp) :: trunk_sum,child_sum,scale

    flow_control=.FALSE.
    num_ctrl=0
    if(.not.allocated(perfusion_flow_elements))return
    flow_control=.TRUE.

    call free_flow_control()
    allocate(ctrl_group(num_elems),ctrl_prescribed(num_elems),ctrl_target(num_elems), &
      ctrl_logscale(num_elems),ctrl_mult(num_elems),resist_base(num_elems), &
      known(num_elems),flow_value(num_elems),stat=allocate_status)
    if(allocate_status.ne.0)stop '*** Not enough memory for prescribed perfusion flows ***'
    ctrl_group=-1
    ctrl_prescribed=.FALSE.
    ctrl_target=0.0_dp
    ctrl_logscale=0.0_dp
    ctrl_mult=1.0_dp
    resist_base=0.0_dp
    known=.FALSE.
    flow_value=0.0_dp

    !--- validate the prescribed elements
    do no=1,size(perfusion_flow_elements)
      ne=perfusion_flow_elements(no)
      if(ne.lt.1.or.ne.gt.num_elems)then
        print *, 'ERROR: prescribed perfusion element is outside the mesh:',ne
        call exit(1)
      endif
      if(ctrl_prescribed(ne))then
        print *, 'ERROR: prescribed perfusion element occurs more than once:',ne
        call exit(1)
      endif
      if(elem_cnct(-1,0,ne).eq.0)then
        print *, 'ERROR: prescribed perfusion flow is on an inlet element:',ne
        call exit(1)
      endif
      if(nint(elem_field(ne_group,ne)).eq.1)then
        print *, 'ERROR: prescribed perfusion flow is on a capillary element:',ne
        call exit(1)
      endif
      if(perfusion_flow_values(no).le.0.0_dp)then
        print *, 'ERROR: prescribed perfusion flows must be positive (element',ne,')'
        call exit(1)
      endif
      ctrl_prescribed(ne)=.TRUE.
      known(ne)=.TRUE.
      flow_value(ne)=perfusion_flow_values(no)
    enddo

    !--- flows known from the boundary conditions
    if(bc_type.eq.'flow')then
      num_inlets=0
      num_outlets=0
      do ne=1,num_elems
        if(elem_cnct(-1,0,ne).eq.0)then
          num_inlets=num_inlets+1
          known(ne)=.TRUE.
          flow_value(ne)=inletbc
        endif
      enddo
      do ne=1,num_elems
        if(elem_cnct(1,0,ne).eq.0.and.nint(elem_field(ne_group,ne)).eq.2)num_outlets=num_outlets+1
      enddo
      if(num_outlets.eq.1)then
        do ne=1,num_elems
          if(elem_cnct(1,0,ne).eq.0.and.nint(elem_field(ne_group,ne)).eq.2)then
            known(ne)=.TRUE.
            flow_value(ne)=inletbc*num_inlets
          endif
        enddo
      endif
    endif

    !--- propagate known flows from trunk to children by conservation
    changed=.TRUE.
    do while(changed)
      changed=.FALSE.
      do np=1,num_nodes
        if(elems_at_node(np,0).lt.2)cycle
        all_known=.TRUE.
        n_trunk=0
        n_unknown=0
        ne_unknown=0
        trunk_sum=0.0_dp
        child_sum=0.0_dp
        do k=1,elems_at_node(np,0)
          ne=elems_at_node(np,k)
          if(elem_root_node(ne).eq.np)then
            if(known(ne))then
              child_sum=child_sum+flow_value(ne)
            else
              n_unknown=n_unknown+1
              ne_unknown=ne
            endif
          else
            n_trunk=n_trunk+1
            if(known(ne))then
              trunk_sum=trunk_sum+flow_value(ne)
            else
              all_known=.FALSE.
            endif
          endif
        enddo
        if(.not.all_known.or.n_trunk.eq.0)cycle
        if(n_unknown.eq.1)then
          known(ne_unknown)=.TRUE.
          flow_value(ne_unknown)=trunk_sum-child_sum
          changed=.TRUE.
        endif
      enddo
    enddo

    !--- identify junction groups and independently controlled elements
    do np=1,num_nodes
      if(elems_at_node(np,0).lt.2)cycle
      all_known=.TRUE.
      any_prescribed=.FALSE.
      n_trunk=0
      n_child=0
      trunk_sum=0.0_dp
      child_sum=0.0_dp
      do k=1,elems_at_node(np,0)
        ne=elems_at_node(np,k)
        if(.not.known(ne))all_known=.FALSE.
        if(elem_root_node(ne).eq.np)then
          n_child=n_child+1
          child_sum=child_sum+flow_value(ne)
          if(ctrl_prescribed(ne))any_prescribed=.TRUE.
        else
          n_trunk=n_trunk+1
          trunk_sum=trunk_sum+flow_value(ne)
        endif
      enddo
      if(.not.all_known.or.n_trunk.eq.0.or.n_child.eq.0)cycle
      scale=max(abs(trunk_sum),abs(child_sum),1.0e-30_dp)
      if(abs(trunk_sum-child_sum).gt.1.0e-6_dp*scale)then
        print *, 'ERROR: prescribed perfusion flows violate conservation at node',np
        print *, '       flow in =',trunk_sum,' flow out =',child_sum
        call exit(1)
      endif
      if(.not.any_prescribed)cycle
      do k=1,elems_at_node(np,0)
        ne=elems_at_node(np,k)
        if(elem_root_node(ne).ne.np)cycle
        if(nint(elem_field(ne_group,ne)).eq.1)then
          print *, 'ERROR: prescribed perfusion flows imply a capillary flow at node',np
          call exit(1)
        endif
        if(flow_value(ne).le.0.0_dp)then
          print *, 'ERROR: prescribed perfusion flows imply a non-positive flow in element',ne
          print *, '       implied flow =',flow_value(ne)
          call exit(1)
        endif
        ctrl_group(ne)=np
        ctrl_target(ne)=flow_value(ne)
      enddo
    enddo
    do no=1,size(perfusion_flow_elements)
      ne=perfusion_flow_elements(no)
      if(ctrl_group(ne).lt.0)then
        ctrl_group(ne)=0
        ctrl_target(ne)=flow_value(ne)
      endif
    enddo
    num_ctrl=count(ctrl_group.ge.0)

    deallocate(known,flow_value)
  end subroutine setup_flow_control

!###################################################################################
!*compute_ctrl_multipliers:* Resistance multiplier of every element: product of
! the scale factors of all controlled elements on its path to the trunk.
  subroutine compute_ctrl_multipliers()
    integer :: ne,ne0
    real(dp) :: logsum

    do ne=1,num_elems
      logsum=0.0_dp
      if(nint(elem_field(ne_group,ne)).ne.1)then
        ne0=ne
        do while(ne0.gt.0)
          if(ctrl_group(ne0).ge.0)logsum=logsum+ctrl_logscale(ne0)
          ne0=elem_trunk(ne0)
        enddo
      endif
      ctrl_mult(ne)=exp(logsum)
    enddo
  end subroutine compute_ctrl_multipliers

!###################################################################################
!*store_and_scale_resistance:* Call after resistances have been freshly computed
! from vessel radii: keeps them as the base and applies the control multipliers.
  subroutine store_and_scale_resistance()
    resist_base(1:num_elems)=elem_field(ne_resist,1:num_elems)
    call rescale_resistance()
  end subroutine store_and_scale_resistance

!###################################################################################
  subroutine rescale_resistance()
    call compute_ctrl_multipliers()
    elem_field(ne_resist,1:num_elems)=resist_base(1:num_elems)*ctrl_mult(1:num_elems)
  end subroutine rescale_resistance

!###################################################################################
!*update_flow_control:* Compares solved flows with the targets and updates the
! scale factors (multiplicative update, exact for parallel branches), then
! applies the junction normalisation chosen by perfusion_flow_mode.
  subroutine update_flow_control(prq_solution,depvar_at_elem,mesh_dof,flow_err)
    integer,intent(in) :: mesh_dof
    real(dp),intent(in) :: prq_solution(mesh_dof,2)
    integer,intent(in) :: depvar_at_elem(0:2,2,num_elems)
    real(dp),intent(out) :: flow_err

    integer :: ne,ne2,np,k,n
    real(dp) :: q,step,shift
    real(dp),parameter :: max_step=log(4.0_dp),max_log=log(1.0e8_dp)

    flow_err=0.0_dp
    do ne=1,num_elems
      if(ctrl_group(ne).lt.0)cycle
      q=prq_solution(depvar_at_elem(0,1,ne),1)
      flow_err=max(flow_err,abs(q-ctrl_target(ne))/ctrl_target(ne))
      if(q.gt.0.0_dp)then
        step=log(q/ctrl_target(ne))
      else
        step=-max_step !reversed flow: lower the resistance
      endif
      ctrl_logscale(ne)=ctrl_logscale(ne)+max(-max_step,min(max_step,step))
    enddo

    ! Junction groups: remove the redundant common factor.
    do ne=1,num_elems
      np=ctrl_group(ne)
      if(np.le.0)cycle
      if(elem_root_node(ne).ne.np)cycle
      ! process each group once, from its first child at node np
      if(first_group_child(np).ne.ne)cycle
      n=0
      select case(trim(perfusion_flow_mode))
      case('constrict')
        shift=huge(1.0_dp)
      case('dilate')
        shift=-huge(1.0_dp)
      case default
        shift=0.0_dp
      end select
      do k=1,elems_at_node(np,0)
        ne2=elems_at_node(np,k)
        if(ctrl_group(ne2).ne.np)cycle
        n=n+1
        select case(trim(perfusion_flow_mode))
        case('constrict')
          shift=min(shift,ctrl_logscale(ne2))
        case('dilate')
          shift=max(shift,ctrl_logscale(ne2))
        case default
          shift=shift+ctrl_logscale(ne2)
        end select
      enddo
      if(trim(perfusion_flow_mode).eq.'geometric')shift=shift/n
      do k=1,elems_at_node(np,0)
        ne2=elems_at_node(np,k)
        if(ctrl_group(ne2).eq.np)ctrl_logscale(ne2)=ctrl_logscale(ne2)-shift
      enddo
    enddo

    do ne=1,num_elems
      if(ctrl_group(ne).lt.0)cycle
      if(abs(ctrl_logscale(ne)).gt.max_log)then
        print *,'WARNING: resistance multiplier for prescribed flow in element',ne, &
          'reached its limit; the requested flow may not be achievable'
        ctrl_logscale(ne)=sign(max_log,ctrl_logscale(ne))
      endif
    enddo
  end subroutine update_flow_control

!###################################################################################
  integer function first_group_child(np)
    integer,intent(in) :: np
    integer :: k

    first_group_child=0
    do k=1,elems_at_node(np,0)
      if(ctrl_group(elems_at_node(np,k)).eq.np)then
        first_group_child=elems_at_node(np,k)
        return
      endif
    enddo
  end function first_group_child

!###################################################################################
  subroutine report_flow_control()
    integer :: ne
    character(len=10) :: kind

    print *,'Prescribed perfusion flows (mode = ',trim(perfusion_flow_mode),'):'
    print '(a)','     element       kind   junction        target      achieved  R multiplier'
    do ne=1,num_elems
      if(ctrl_group(ne).lt.0)cycle
      if(ctrl_prescribed(ne))then
        kind='prescribed'
      else
        kind='   implied'
      endif
      print '(i12,1x,a10,i11,2es14.5,es14.5)',ne,kind,ctrl_group(ne),ctrl_target(ne), &
        elem_field(ne_Qdot,ne),exp(ctrl_logscale(ne))
    enddo
  end subroutine report_flow_control

!###################################################################################
  subroutine free_flow_control()
    if(allocated(ctrl_group))deallocate(ctrl_group)
    if(allocated(ctrl_prescribed))deallocate(ctrl_prescribed)
    if(allocated(ctrl_target))deallocate(ctrl_target)
    if(allocated(ctrl_logscale))deallocate(ctrl_logscale)
    if(allocated(ctrl_mult))deallocate(ctrl_mult)
    if(allocated(resist_base))deallocate(resist_base)
  end subroutine free_flow_control
!
!###################################################################################
!
subroutine calculate_resistance(viscosity,KOUNT)

  real(dp):: viscosity
    integer :: KOUNT
!local variables
    integer :: ne,np1,np2
    real(dp) :: resistance,zeta
    character(len=60) :: sub_name

  sub_name = 'calculate_resistance'
  call enter_exit(sub_name,1)

!Loop over all elements in model and define resistance for that branch.
    do ne=1,num_elems
       !ne=elems(noelem)
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       ! element length
       elem_field(ne_length,ne) = DSQRT((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)
       elem_field(ne_radius,ne)=(elem_field(ne_radius_in,ne)+elem_field(ne_radius_out,ne))/2.0_dp
       ! element Poiseuille (laminar) resistance in units of Pa.s.mm-3
       resistance = 8.d0*viscosity*elem_field(ne_length,ne)/ &
            (PI*elem_field(ne_radius,ne)**4) !laminar resistance
       ! element turbulent resistance (flow in bifurcating tubes)
       !reynolds=DABS(elem_field(ne_Qdot,ne)*2.d0*density/ &
         !   (PI*elem_field(ne_radius,ne)*viscosity))
       zeta = 1.0_dp!MAX(1.d0,dsqrt(2.d0*elem_field(ne_radius,ne)* &
            !reynolds/elem_field(ne_length,ne))*gamma)
       if(elem_field(ne_group,ne).eq.1.0_dp)then
         elem_field(ne_resist,ne) = 1000.0_dp !initialises resistance for first iteration
       else
        elem_field(ne_resist,ne) = resistance * zeta
       endif
       !print *,"TESTING RESISTANCE",ne,elem_field(ne_resist,ne),elem_field(ne_radius,ne),elem_ordrs(2,ne)
      !endif
    enddo

    call enter_exit(sub_name,2)
  end subroutine calculate_resistance
!
!##################################################################
!
!*calc_depvar_maps:* calculates the mapping between the nodes and elements and
!the problem depvar that are needed for matrix setup and solution
  subroutine calc_depvar_maps(mesh_from_depvar,depvar_at_elem,&
                depvar_totals,depvar_at_node,mesh_dof,num_vars)

    character(len=60) :: sub_name

     integer :: mesh_dof
     integer :: mesh_from_depvar(0:2,mesh_dof,0:2)
     integer :: depvar_at_elem(0:2,2,num_elems)
     integer :: depvar_at_node(num_nodes,0:2,2)
     integer :: depvar_totals(0:2,2)
     integer :: num_vars
!     local variables
    integer :: ny_start=0
    integer :: nc,ne,nn,np,nrc,ny

    sub_name = 'calc_depvar_maps'
    call enter_exit(sub_name,1)

     depvar_totals = 0
     mesh_from_depvar = 0
     depvar_at_elem = 0
     depvar_at_node = 0

!nrc = loops from 0,1,2

!  Set up mapping arrays for current region:
!  includes nested loop creating depvar_at_node & depvar_at_elem consecutively
     do nrc=0,2 !row or column no
        ny=0 !depvar tag
        do nc=1,2 !no of dependent depvar types
          ! if(nrc.NE.0) ny=ny_start !--> resets to ny_start only for nrc=1,2. Maybe should also do for nrc=1??
           ny=ny_start
           do ne=1,num_elems
              !ne=elems(noelem)
              do nn=1,2 !Loop over number of nodes in each element: always=2 in 1D element
                 np=elem_nodes(nn,ne)
                 if(depvar_at_node(np,nrc,nc).EQ.0)THEN
                    !Checks if this node already been done
                    ny=ny+1
                     if(ny.GT.mesh_dof) THEN
                       print *,"1. Need to increase mesh_dof! ny=",ny
                       call exit(1)
                    endif
                    if(nrc.NE.0) THEN
                       if(ny.GT.depvar_totals(nrc,nc)) depvar_totals(nrc,nc)=ny
                    endif
                    depvar_at_node(np,nrc,nc)=ny
                    if(nrc.NE.1.OR.nc.EQ.1) THEN
                       ! don't set mesh_from_depvar for nrc=1(rows), nc<>1(LHS) cases
                       mesh_from_depvar(0,ny,nrc)=1 !mesh dof is nodal based
                       mesh_from_depvar(1,ny,nrc)=np
                       mesh_from_depvar(2,ny,nrc)=nc
                    endif !nrc.NE.1
                 endif !depvar_at_node.EQ.0
              enddo !nn (np)
              ny=ny+1
              if(ny.GT.mesh_dof) THEN
                 print *,"2. Need to increase mesh_dof! ny=",ny
                 call exit(1)
              endif
              if(nrc.NE.0) THEN
                 if(ny.GT.depvar_totals(nrc,nc)) depvar_totals(nrc,nc)=ny
              endif
              depvar_at_elem(nrc,nc,ne)=ny
              if(nrc.NE.1.OR.nc.EQ.1) THEN
                 !                 don't set mesh_from_depvar for nrc=1(rows), nc<>1(LHS) cases
                 mesh_from_depvar(0,ny,nrc)=2 !mesh dof is element based
                 mesh_from_depvar(1,ny,nrc)=nc
                 mesh_from_depvar(2,ny,nrc)=ne
              endif
           enddo !noelem (ne)
        enddo !nc
     enddo !nrc
     num_vars=ny
     !print *,"Max ny number=",ny


    call enter_exit(sub_name,2)
  end subroutine calc_depvar_maps
!
!##################################################################
!
!*tree_resistance:* Calculates the total resistance of a tree (arterial tree only)

  subroutine tree_resistance(resistance)

    character(len=60) :: sub_name
!local variables
    real(dp), intent(out) :: resistance
    real(dp) :: invres,elem_res(num_elems)
    integer :: num2,ne,ne2

    sub_name = 'tree_resistance'
    call enter_exit(sub_name,1)

    elem_res(1:num_elems)=elem_field(ne_resist,1:num_elems)
    do ne=num_elems,1,-1
       !ne=elems(num)
      invres=0.0_dp
      do num2=1,elem_cnct(1,0,ne)
         ne2=elem_cnct(1,num2,ne)
         invres=invres+1.0_dp/elem_res(ne2)
      enddo
      if(elem_cnct(1,0,ne).gt.0)then
        elem_res(ne)=elem_res(ne)+1.0_dp/invres
       endif
    enddo
    resistance=elem_res(1)

    call enter_exit(sub_name,2)
  end subroutine tree_resistance
!
!##################################################################
!
!*initialise_solution:* Calculates an estimate for initial solution to a prq problem based on cardiact output and pressure BCs

subroutine initialise_solution(pressure_in,pressure_out,cardiac_output,mesh_dof,prq_solution,&
    depvar_at_node,depvar_at_elem,FIX)

    integer, intent(in) :: mesh_dof
    integer,intent(in) :: depvar_at_elem(0:2,2,num_elems)
    integer,intent(in) :: depvar_at_node(num_nodes,0:2,2)
    real(dp), intent(in) :: pressure_in, pressure_out,cardiac_output
    real(dp) :: prq_solution(mesh_dof,2)
    logical, intent(in) :: FIX(mesh_dof)
!local variables
    integer :: nn,ne,np,n_depvar
    character(len=60) :: sub_name
    sub_name = 'intialise_solution'
    call enter_exit(sub_name,1)
    do ne=1,num_elems
       !ne=elems(noelem)
       do nn=1,2 !Loop over number of nodes in each element: always=2 in 1D element
         np=elem_nodes(nn,ne) !Node number
         n_depvar=depvar_at_node(np,0,1) !--> This will be a pressure because it's the depvar at a node
         if(.NOT.FIX(n_depvar))then
           prq_solution(n_depvar,1)=(pressure_in+pressure_out)/2.0
         endif
         n_depvar=depvar_at_elem(0,1,ne) !--> This will be a flow because it's the depvar at the element
         if(.NOT.FIX(n_depvar))then
           prq_solution(n_depvar,1)=cardiac_output/(2**(elem_ordrs(1,ne)-1)) !Here you can use the generation number to split flow
         endif
       enddo
    enddo
    call enter_exit(sub_name,2)
  end subroutine initialise_solution
!
!##################################################################
!
!*calc_sparse_1d_tree:* Calculates sparsity structure for 1d tree problems

subroutine calc_sparse_1dtree(bc_type,density,FIX,grav_vect,mesh_dof,depvar_at_elem,&
        depvar_at_node,MaxNonZeros,MatrixSize,SparseCol,SparseRow,SparseVal,RHS,&
        prq_solution,update_resistance_entries,update_fixed_flow_rhs_rows,&
        fixed_flow_rhs_base,NonZeros)

    character(len=60), intent(in) :: bc_type
    real(dp), intent(in) :: density
    real(dp), intent(in) :: grav_vect(3)
    integer, intent(in) :: mesh_dof
    logical, intent(in) :: FIX(mesh_dof)
    integer,intent(in) :: depvar_at_elem(0:2,2,num_elems)
    integer,intent(in) :: depvar_at_node(num_nodes,0:2,2)

    integer, intent(in) :: MaxNonZeros,MatrixSize
    integer, intent(inout) :: SparseCol(MaxNonZeros)
    integer, intent(inout) :: SparseRow(MatrixSize+1)
    real(dp), intent(inout) :: SparseVal(MaxNonZeros)
    real(dp), intent(inout) :: RHS(MatrixSize)
    real(dp), intent(inout) :: prq_solution(mesh_dof,2)
    integer, intent(inout) :: update_resistance_entries(num_elems)
    integer, intent(inout) :: update_fixed_flow_rhs_rows(num_elems)
    real(dp), intent(inout) :: fixed_flow_rhs_base(num_elems)
    integer, intent(out) :: NonZeros
!local variables
    integer :: ne,nn,np,np1,np2,depvar,depvar1,depvar2,depvar3,flow_var,offset,nzz,&
      nzz_row,ne2,noelem2,ne3,row_nz
    logical :: FlowBalancedNodes(num_nodes)
    logical :: NodePressureDone(num_nodes)
    logical :: ElementPressureEquationDone(num_elems)
    logical :: elem_found,one_node_balanced,has_diagonal
    real(dp) :: flow_term
    real(dp) :: grav
    integer :: nj

    character(len=60) :: sub_name
    sub_name = 'calc_sparse_1dtree'
    call enter_exit(sub_name,1)
    !Initialise matrices and indices
    SparseCol=0
    SparseRow=1
    SparseVal=0.0_dp
    RHS=0.0_dp

    nzz=1 !position in SparseCol and SparseVal
    nzz_row=1 !position in SparseRow
    FlowBalancedNodes = .FALSE. !.TRUE. for nodes which have had a conservation of flow equation done
    NodePressureDone = .FALSE.  !.TRUE. for nodes which have been processed
    ElementPressureEquationDone = .FALSE.
    offset=0!variable position offset

    ! Keep the coupling-style variable-driven equation ordering so the reduced
    ! matrix stays structurally aligned with the ILU diagonal.
    do ne=1,num_elems
      !look at pressure variables at each node
      do nn=1,2 !2 nodes in 1D element
        np=elem_nodes(nn,ne)
        depvar = depvar_at_node(np,1,1)
        if((.NOT.NodePressureDone(np)).AND.(.NOT.FIX(depvar)))then
            ne2=0
            if(nn.EQ.1)then
                ne2=ne
            elseif(nn.EQ.2)then
                if((bc_type.EQ.'pressure').OR.(.NOT.ElementPressureEquationDone(ne)))then
                    ne2=ne
                endif
            endif
            if((ne2.EQ.0).AND.(elems_at_node(np,0).GT.1))then
                elem_found=.FALSE.
                noelem2=1
                do while ((.NOT.elem_found).AND.(noelem2.LE.elems_at_node(np,0)))
                    ne3=elems_at_node(np,noelem2)
                    if((ne3.NE.ne).AND.(.NOT.ElementPressureEquationDone(ne3)))then
                        ne2=ne3
                        elem_found=.TRUE.
                    endif
                    noelem2=noelem2+1
                enddo
            endif
            if(ne2.GT.0)then
                np1=elem_nodes(1,ne2)
                depvar1=depvar_at_node(np1,1,1)
                np2=elem_nodes(2,ne2)
                depvar2=depvar_at_node(np2,1,1)
                depvar3=depvar_at_elem(0,1,ne2)
                grav=0.d0
                if(elem_field(ne_group,ne2).eq.1.0_dp)then
                elseif(elem_ordrs(no_gen,ne2).eq.1)then
                else
                  do nj=1,3
                    grav=grav+density*grav_vect(nj)*9810.0_dp* &
                      (node_xyz(nj,elem_nodes(1,ne2))-node_xyz(nj,elem_nodes(2,ne2)))
                  enddo
                endif
                RHS(nzz_row)=grav
                if(FIX(depvar1))then
                    RHS(nzz_row)=RHS(nzz_row)-prq_solution(depvar1,1)
                else
                    call get_variable_offset(depvar1,mesh_dof,FIX,offset)
                    SparseCol(nzz)=depvar1-offset
                    SparseVal(nzz)=1.0_dp
                    nzz=nzz+1
                endif
                if(FIX(depvar2))then
                    RHS(nzz_row)=RHS(nzz_row)+prq_solution(depvar2,1)
                else
                    call get_variable_offset(depvar2,mesh_dof,FIX,offset)
                    SparseCol(nzz)=depvar2-offset
                    SparseVal(nzz)=-1.0_dp
                    nzz=nzz+1
                endif
                if(FIX(depvar3))then
                    fixed_flow_rhs_base(ne2)=RHS(nzz_row)
                    RHS(nzz_row)=RHS(nzz_row)+ &
                      prq_solution(depvar3,1)*elem_field(ne_resist,ne2)
                    update_fixed_flow_rhs_rows(ne2)=nzz_row
                else
                    call get_variable_offset(depvar3,mesh_dof,FIX,offset)
                    SparseCol(nzz)=depvar3-offset
                    SparseVal(nzz)=-elem_field(ne_resist,ne2)
                    update_resistance_entries(ne2)=nzz
                    nzz=nzz+1
                endif
                nzz_row=nzz_row+1
                SparseRow(nzz_row)=nzz
                NodePressureDone(np)=.TRUE.
                ElementPressureEquationDone(ne2)=.TRUE.
            endif
        endif
      enddo

      ! Select the equation associated with each remaining flow variable.
      ! A fixed flow is already represented in pressure and conservation RHS
      ! terms, so its variable row is omitted from the reduced system.
      flow_var=depvar_at_elem(0,1,ne)
      if(.NOT.FIX(flow_var))then
        one_node_balanced=.FALSE.
        do nn=1,2
          np=elem_nodes(nn,ne)
          if((elems_at_node(np,0).GT.1).AND.(.NOT.FlowBalancedNodes(np)))then
            if((bc_type.EQ.'pressure').OR.((bc_type.EQ.'flow').AND.(.NOT.one_node_balanced)))then
                do noelem2=1,elems_at_node(np,0)
                    ne2=elems_at_node(np,noelem2)
                    depvar=depvar_at_elem(1,1,ne2)
                    flow_term=0.0_dp
                    if(np.EQ.elem_nodes(2,ne2))then
                        flow_term=1.0_dp
                    elseif(np.EQ.elem_nodes(1,ne2))then
                        flow_term=-1.0_dp
                    endif
                    if(FIX(depvar))then
                        RHS(nzz_row)=RHS(nzz_row)-prq_solution(depvar,1)*flow_term
                    else
                        call get_variable_offset(depvar,mesh_dof,FIX,offset)
                        SparseCol(nzz)=depvar-offset
                        SparseVal(nzz)=flow_term
                        nzz=nzz+1
                    endif
                enddo
                FlowBalancedNodes(np)=.TRUE.
                nzz_row=nzz_row+1
                SparseRow(nzz_row)=nzz
                one_node_balanced=.TRUE.
            endif
          endif
        enddo

        if((.NOT.one_node_balanced).AND.(.NOT.ElementPressureEquationDone(ne)))then
            np1=elem_nodes(1,ne)
            depvar1=depvar_at_node(np1,1,1)
            np2=elem_nodes(2,ne)
            depvar2=depvar_at_node(np2,1,1)

            call get_variable_offset(depvar1,mesh_dof,FIX,offset)
            SparseCol(nzz)=depvar1-offset
            SparseVal(nzz)=1.0_dp
            nzz=nzz+1

            if(FIX(depvar2))then
                RHS(nzz_row)=prq_solution(depvar2,1)
            else
                call get_variable_offset(depvar2,mesh_dof,FIX,offset)
                SparseCol(nzz)=depvar2-offset
                SparseVal(nzz)=-1.0_dp
                nzz=nzz+1
            endif

            call get_variable_offset(flow_var,mesh_dof,FIX,offset)
            SparseCol(nzz)=flow_var-offset
            SparseVal(nzz)=-elem_field(ne_resist,ne)
            update_resistance_entries(ne)=nzz
            nzz=nzz+1

            nzz_row=nzz_row+1
            SparseRow(nzz_row)=nzz
            ElementPressureEquationDone(ne)=.TRUE.
        endif
      endif
    enddo

    NonZeros=nzz-1
    if(nzz_row-1.ne.MatrixSize)then
      print *, 'ERROR: sparse PRQ assembly row mismatch:',nzz_row-1,MatrixSize
      call exit(1)
    elseif(NonZeros.gt.MaxNonZeros)then
      print *, 'ERROR: sparse PRQ assembly exceeds allocated storage:',NonZeros,MaxNonZeros
      call exit(1)
    endif

    do np=1,MatrixSize
      has_diagonal=.FALSE.
      do row_nz=SparseRow(np),SparseRow(np+1)-1
        if(SparseCol(row_nz).eq.np)then
          has_diagonal=.TRUE.
          exit
        endif
      enddo
      if(.NOT.has_diagonal)then
        print *, 'ERROR: sparse PRQ assembly has no ILU diagonal entry in row:',np
        call exit(1)
      endif
    enddo

    call enter_exit(sub_name,2)
  end subroutine calc_sparse_1dtree


!
!##################################################################
!
!*calc_sparse_size:* Calculates sparsity sizes

subroutine calc_sparse_size(mesh_dof,depvar_at_elem,depvar_at_node,FIX,MaxNonZeros,MatrixSize)

  integer, intent(in) :: mesh_dof
    integer,intent(in) :: depvar_at_elem(0:2,2,num_elems)
    integer,intent(in) :: depvar_at_node(num_nodes,0:2,2)
    logical, intent(in) :: FIX(mesh_dof)
    integer :: MaxNonZeros,MatrixSize
!local variables
    integer :: i,np,fixed_variables
    character(len=60) :: sub_name
    sub_name = 'calc_sparse_size'
    call enter_exit(sub_name,1)


    fixed_variables = 0
    !count fixed variables
    do i=1,mesh_dof
        if(FIX(i))then
            fixed_variables = fixed_variables + 1
        endif
    enddo
    MatrixSize = mesh_dof - fixed_variables

    ! Allocate an upper bound and let calc_sparse_1dtree return the exact
    ! number used.  Fixed internal flows can remove whole conservation rows,
    ! so subtracting a fixed count here is not generally exact.
    MaxNonZeros=num_elems*3
    do np=1,num_nodes
      if(elems_at_node(np,0).gt.1)then
        MaxNonZeros=MaxNonZeros+elems_at_node(np,0)
      endif
    enddo
    call enter_exit(sub_name,2)
  end subroutine calc_sparse_size

!##################################################################
!
!*calc_press_area:* Calculates new radii based on pressure area relnships

subroutine calc_press_area(grav_vect,KOUNT,depvar_at_node,prq_solution,&
    mesh_dof,vessel_type,elasticity_parameters,mechanics_parameters,remodeling_grade)

    character(len=60), intent(in) :: vessel_type
    real(dp), intent(in) :: grav_vect(3)
    integer,intent(in) :: KOUNT,mesh_dof
    integer,intent(in) :: depvar_at_node(num_nodes,0:2,2)
    real(dp),intent(in) ::  prq_solution(mesh_dof,2),remodeling_grade
    real(dp),intent(in) :: elasticity_parameters(3),mechanics_parameters(2)

!local variables
    integer :: nj,np,ne,ny,nn,k,no_compl_ne
    real(dp) :: h,Ptm,R0,Pblood,Ppl,counter,cc1,cc2,cc3
    real(dp) :: alt_hyp,alt_fib,prox_fib,narrow_rad_one,narrow_rad_two,narrow_factor,prune_rad,prune_fraction,counter1,counter2
    integer,allocatable :: templss(:)
    logical:: check, FOUND

    character(len=60) :: sub_name
    sub_name = 'calc_press_area'
    call enter_exit(sub_name,1)
    if(KOUNT.eq.1)then !store initial, unstressed radius values
      do  ne=1,num_elems
        elem_field(ne_radius_in0,ne)=elem_field(ne_radius_in,ne)
        elem_field(ne_radius_out0,ne)=elem_field(ne_radius_out,ne)
      enddo !elems
    endif

  if(remodeling_grade.eq.0.0_dp) then  ! Solving for Healthy
    do ne=1,num_elems
      do nn=1,2
        if(nn.eq.1) np=elem_nodes(1,ne)
        if(nn.eq.2) np=elem_nodes(2,ne)
        ny=depvar_at_node(np,0,1)
        call calculate_ppl(np,grav_vect,mechanics_parameters,Ppl)
        Pblood=prq_solution(ny,1) !Pa
        Ptm=Pblood+Ppl     ! Pa
        if(nn.eq.1)R0=elem_field(ne_radius_in0,ne)
        if(nn.eq.2)R0=elem_field(ne_radius_out0,ne)
        if(vessel_type.eq.'elastic_g0_beta')then
          if(Ptm.lt.elasticity_parameters(3).and.elasticity_parameters(1).gt.0.0_dp)then
            if(nn.eq.1) elem_field(ne_radius_in,ne)=R0*((Ptm/elasticity_parameters(1))+1.d0)**(1.d0/elasticity_parameters(2))
            if(nn.eq.2) elem_field(ne_radius_out,ne)=R0*((Ptm/elasticity_parameters(1))+1.d0)**(1.d0/elasticity_parameters(2))
          elseif(Ptm.lt.0.0_dp.or.elasticity_parameters(1).LT.TOLERANCE)then
            if(Ptm.lt.0)write(*,*) 'Transmural pressure < zero',ne,Ptm,Pblood,Ppl
            if(nn.eq.1) elem_field(ne_radius_in,ne)=R0
            if(nn.eq.2) elem_field(ne_radius_out,ne)=R0
          else !ptm>ptmmax
            if(nn.eq.1)then
               elem_field(ne_radius_in,ne)=R0*((elasticity_parameters(3)/elasticity_parameters(1))+1.d0) &
                 **(1.d0/elasticity_parameters(2))
            endif
            if(nn.eq.2)then
              elem_field(ne_radius_out,ne)=R0*((elasticity_parameters(3)/elasticity_parameters(1))+1.d0) &
                **(1.d0/elasticity_parameters(2))
            endif
          endif
        elseif(vessel_type.eq.'elastic_alpha')then
           if(Ptm.LT.elasticity_parameters(2))then
             if(nn.eq.1) elem_field(ne_radius_in,ne)=R0*((Ptm*elasticity_parameters(1))+1.d0)
             if(nn.eq.2) elem_field(ne_radius_out,ne)=R0*((Ptm*elasticity_parameters(1))+1.d0)
           elseif(Ptm.lt.0.0_dp)then
            if(Ptm.lt.0)write(*,*) 'Transmural pressure < zero',ne,Ptm,Pblood,Ppl
            if(nn.eq.1) elem_field(ne_radius_in,ne)=R0
            if(nn.eq.2) elem_field(ne_radius_out,ne)=R0
           else !ptm>ptmmax
              if(nn.eq.1)then
               elem_field(ne_radius_in,ne)=R0*((elasticity_parameters(2)*elasticity_parameters(1))+1.d0)
              endif
              if(nn.eq.2)then
                elem_field(ne_radius_out,ne)=R0*((elasticity_parameters(2)*elasticity_parameters(1))+1.d0)
              endif
          endif
        elseif(vessel_type.eq.'elastic_hooke')then
          h=elasticity_parameters(2)*R0
          if(nn.eq.1) elem_field(ne_radius_in,ne)=R0+3.0_dp*R0**2*Ptm/(4.0_dp*elasticity_parameters(1)*h)
          if(nn.eq.2) elem_field(ne_radius_out,ne)=R0+3.0_dp*R0**2*Ptm/(4.0_dp*elasticity_parameters(1)*h)
        else
        print *, 'no vessel type defined, assuming rigid'
        if(nn.eq.1) elem_field(ne_radius_in,ne)=R0
        if(nn.eq.2) elem_field(ne_radius_out,ne)=R0
      endif
      enddo!nn
    enddo!ne
  else ! Solving for remodeling case - only implemented for elastic_alpha

    if(remodeling_grade.ge.60)then  ! hypertophy effect
      alt_hyp = 1.0_dp/6
    else
      alt_hyp = (-1.0_dp/60)*remodeling_grade + 7.0_dp/6
    endif

    ! Narrowing effect
    narrow_rad_one=0.015_dp
    narrow_rad_two=0.15_dp
    if(remodeling_grade.le.20.0_dp)then
      narrow_factor = 1.0_dp
    elseif(remodeling_grade.ge.60.0_dp)then
      narrow_factor = 0.55_dp
    else
      narrow_factor = (-9.0_dp/800) * remodeling_grade + 1.225_dp
    endif

    ! pruning fraction
    if(remodeling_grade.le.20.0_dp)then
      prune_fraction = 0
      prune_rad = 0.16_dp
    elseif(remodeling_grade.ge.50.0_dp)then
      prune_rad = 0.25_dp
      prune_fraction = (1.0_dp/160) * remodeling_grade - 1.0_dp/8
    else
      prune_rad = 0.16_dp
      prune_fraction = (1.0_dp/160) * remodeling_grade - 1.0_dp/8
    endif

    ! fibrosis effect
    if(remodeling_grade.le.50.0_dp)then
      alt_fib = 1.0_dp
    else
      alt_fib = (-1.0_dp/60) * remodeling_grade + 11.0_dp/6
    endif
    counter=0.0_dp
    cc1=0.0_dp
    cc2=0.0_dp
    cc3=0.0_dp
    counter1 = 1.0_dp
    counter2 = 1.0_dp
    do ne=1,num_elems
      do nn=1,2
        if(nn.eq.1) np=elem_nodes(1,ne)
        if(nn.eq.2) np=elem_nodes(2,ne)
        ny=depvar_at_node(np,0,1)
        call calculate_ppl(np,grav_vect,mechanics_parameters,Ppl)
        Pblood=prq_solution(ny,1) !Pa
        Ptm=Pblood+Ppl     ! Pa
        if(nn.eq.1) R0=elem_field(ne_radius_in0,ne)
        if(nn.eq.2) R0=elem_field(ne_radius_out0,ne)
        if(elem_field(ne_group,ne).eq.0.0_dp) then !only applying on arteries
          FOUND = .False.
          if (allocated(occ_list))then
            do k=1,count(occ_list.ne.0)
              check = is_downstream(ne,occ_list(k))
              if(check)then
                FOUND = .True.
              endif
            end do
          endif
          if(.NOT.FOUND)then ! FOUND
            if(nn.eq.1) then ! nn.eq.1
              if(R0.lt.prune_rad.and.elem_ordrs(no_sord,ne).eq.1) then
                if(counter1/100.le.prune_fraction) then ! pruning the right percentage based on the fraction defined
                  cc1 = cc1+1.0_dp
                  R0=0.005_dp ! Setting the radius to a small value
                else ! the remaining of the canditates that are not pruned because of the fraction
                  R0=elem_field(ne_radius_in0,ne)
                endif
                counter1 = counter1 + 1.0_dp ! since a canditate was found, one is added to the counter1
                if(counter1.ge.101.0_dp) counter1=1.0_dp ! now that the fraction out of hundred was blocked set the counter back to start
              else ! R0 greater than prune_rad
                cc3=cc3+1.0_dp
                R0=elem_field(ne_radius_in0,ne) ! treating the artery as normal unstrained radius (no constraints)
              endif
            endif ! nn.eq.1
            if(nn.eq.2) then  ! same thing as nn=1
              if(R0.lt.prune_rad.and.elem_ordrs(no_sord,ne).eq.1) then
                if(counter2/100.le.prune_fraction) then
                  cc2=cc2+1.0_dp
                  R0=0.005_dp
                else
                  R0=elem_field(ne_radius_out0,ne)
                endif
                counter2 = counter2 + 1.0_dp
                if(counter2.ge.101.0_dp) counter2=1.0_dp
              else
                R0=elem_field(ne_radius_out0,ne)
              endif
            endif ! same thing as nn=1
          else ! treat as healthy
            if(nn.eq.1) R0=elem_field(ne_radius_in0,ne)
            if(nn.eq.2) R0=elem_field(ne_radius_out0,ne)
          endif ! FOUND
          if(vessel_type.eq.'elastic_g0_beta') then ! elastic_g0_beta
            if(elem_field(ne_group,ne).eq.0.0_dp) then !only applying on arteries
              if(Ptm.LT.elasticity_parameters(3).and.elasticity_parameters(1).gt.0.0_dp)then
                if(nn.eq.1) then
                  if((R0.gt.0.015).and.(R0.lt.0.15)) then
                    elem_field(ne_radius_in,ne)=0.55_dp*R0*((Ptm/(0.16_dp*elasticity_parameters(1)))+1.d0) &
                    **(1.d0/elasticity_parameters(2))
                  else
                    elem_field(ne_radius_in,ne)=R0*((Ptm/elasticity_parameters(1))+1.d0)**(1.d0/elasticity_parameters(2))
                  endif
                endif
                if(nn.eq.2) then
                  if((R0.gt.0.015).and.(R0.lt.0.15)) then
                    elem_field(ne_radius_out,ne)=0.55_dp*R0*((Ptm/(0.16_dp*elasticity_parameters(1)))+1.d0) &
                    **(1.d0/elasticity_parameters(2))
                  else
                    elem_field(ne_radius_out,ne)=R0*((Ptm/elasticity_parameters(1))+1.d0)**(1.d0/elasticity_parameters(2))
                  endif
                endif
              elseif(Ptm.lt.0.0_dp.or.elasticity_parameters(1).LT.TOLERANCE)then
                if(Ptm.lt.0)write(*,*) 'Transmural pressure < zero',ne,Ptm,Pblood,Ppl
                if(nn.eq.1) elem_field(ne_radius_in,ne)=R0
                if(nn.eq.2) elem_field(ne_radius_out,ne)=R0
              else !ptm>ptmmax
                if(nn.eq.1) then
                  if((R0.gt.0.015).and.(R0.lt.0.15)) then
                    elem_field(ne_radius_in,ne)=0.55_dp*R0*((elasticity_parameters(3)/(0.16_dp*elasticity_parameters(1)))+1.d0) &
                   **(1.d0/elasticity_parameters(2))
                  else
                    elem_field(ne_radius_in,ne)=R0*((elasticity_parameters(3)/elasticity_parameters(1))+1.d0) &
                   **(1.d0/elasticity_parameters(2))
                  endif
                endif
                if(nn.eq.2) then
                  if((R0.gt.0.015).and.(R0.lt.0.15)) then
                    elem_field(ne_radius_out,ne)=0.55_dp*R0*((elasticity_parameters(3)/(0.16_dp*elasticity_parameters(1)))+1.d0) &
                    **(1.d0/elasticity_parameters(2))
                  else
                    elem_field(ne_radius_out,ne)=R0*((elasticity_parameters(3)/elasticity_parameters(1))+1.d0) &
                    **(1.d0/elasticity_parameters(2))
                  endif
                endif
              endif
            else !other than arteries
              if(Ptm.LT.elasticity_parameters(3).and.elasticity_parameters(1).gt.0.0_dp)then
                if(nn.eq.1) elem_field(ne_radius_in,ne)=R0*((Ptm/elasticity_parameters(1))+1.d0)**(1.d0/elasticity_parameters(2))
                if(nn.eq.2) elem_field(ne_radius_out,ne)=R0*((Ptm/elasticity_parameters(1))+1.d0)**(1.d0/elasticity_parameters(2))
              elseif(Ptm.lt.0.0_dp.or.elasticity_parameters(1).LT.TOLERANCE)then
                if(Ptm.lt.0)write(*,*) 'Transmural pressure < zero',ne,Ptm,Pblood,Ppl
                if(nn.eq.1) elem_field(ne_radius_in,ne)=R0
                if(nn.eq.2) elem_field(ne_radius_out,ne)=R0
              else !ptm>ptmmax
                if(nn.eq.1) elem_field(ne_radius_in,ne)=R0*((elasticity_parameters(3)/elasticity_parameters(1))+1.d0) &
                **(1.d0/elasticity_parameters(2))
                if(nn.eq.2)then
                  elem_field(ne_radius_out,ne)=R0*((elasticity_parameters(3)/elasticity_parameters(1))+1.d0) &
                  **(1.d0/elasticity_parameters(2))
                endif
              endif
            endif !artery or vein
          elseif(vessel_type.eq.'elastic_alpha') then ! elastic alpha
            !if(elem_field(ne_group,ne).eq.0.0_dp) then !only applying on arteries  BEN WAS HERE
              FOUND = .False.
              if (allocated(occ_list))then ! occ allocated
                do k=1,count(occ_list.ne.0)
                  check = is_downstream(ne,occ_list(k))
                  if(check)then
                    FOUND = .True.
                  endif
                end do
              endif ! occ allocated
              if(.NOT.FOUND)then !FOUND
                if(Ptm.lt.elasticity_parameters(2))then
                  if(nn.eq.1) then ! nn.eq.1
                    if((R0.gt.narrow_rad_one).and.(R0.lt.0.5_dp)) then ! Hypertorphy+narrow factor effect
                      if(R0.lt.0.05_dp) then ! only Narrow_factor
                        elem_field(ne_radius_in,ne) = narrow_factor*R0*((Ptm*elasticity_parameters(1))+1.d0)
                      elseif(R0.gt.narrow_rad_two) then ! only Hypertophy
                        elem_field(ne_radius_in,ne) = R0*((Ptm*alt_hyp*elasticity_parameters(1))+1.d0)
                      else ! Both hypertophy and narrowing
                        elem_field(ne_radius_in,ne) = narrow_factor*R0*((Ptm*alt_hyp*alt_fib*elasticity_parameters(1))+1.d0)
                      endif
                    elseif(R0.gt.0.5_dp) then !Not within the range of our target radii, hence, no remodeling for this element
                      elem_field(ne_radius_in,ne) = R0*((Ptm*elasticity_parameters(1))+1.d0)
                    else !Pruning
                      elem_field(ne_radius_in,ne) = R0
                    endif
                  endif ! nn.eq.1
                  if(nn.eq.2) then ! nn.eq.2
                    if((R0.gt.narrow_rad_one).and.(R0.lt.0.5_dp)) then ! Hypertorphy+narrow factor effect
                      if(R0.lt.0.05_dp) then ! only Narrow_factor
                        elem_field(ne_radius_out,ne) = narrow_factor*R0*((Ptm*elasticity_parameters(1))+1.d0)
                      elseif(R0.gt.narrow_rad_two) then ! hypertophy only
                        elem_field(ne_radius_out,ne) = R0*((Ptm*alt_hyp*elasticity_parameters(1))+1.d0)
                      else ! Both hypertophy and narrowing
                        elem_field(ne_radius_out,ne) = narrow_factor*R0*((Ptm*alt_hyp*alt_fib*elasticity_parameters(1))+1.d0)
                      endif
                    elseif(R0.gt.0.5_dp) then !Not within the range of our target radii, hence, no remodeling
                        elem_field(ne_radius_out,ne) = R0*((Ptm*elasticity_parameters(1))+1.d0)
                    else !Pruning
                      elem_field(ne_radius_out,ne) = R0
                    endif
                  endif ! nn.eq.2
                  if (allocated(compl_list)) then ! compl allocated
                    do no_compl_ne=1,count(compl_list.ne.0)
                      !pause
                      if(ne.eq.compl_list(no_compl_ne))then
                        if(nn.eq.1) elem_field(ne_radius_in,compl_list(no_compl_ne))=R0!*((Ptm*elasticity_parameters(1)*0.0)+1.d0)
                        if(nn.eq.2) elem_field(ne_radius_out,compl_list(no_compl_ne))=R0!*((Ptm*elasticity_parameters(1)*0.0)+1.d0)
                      endif
                    end do
                  endif ! compl allocated
                elseif(Ptm.lt.0.0_dp)then !Ptm
                    if(Ptm.lt.0) write(*,*) 'Transmural pressure < zero',ne,Ptm,Pblood,Ppl
                    if(nn.eq.1) then ! nn.eq.1
                      if((R0.gt.narrow_rad_one).and.(R0.lt.0.5)) then ! Hypertophy+narrowing effect
                        if(R0.lt.narrow_rad_two) then !only narrowing
                          elem_field(ne_radius_in,ne)=narrow_factor*R0
                        else
                          elem_field(ne_radius_in,ne)=R0
                        endif
                      else  ! Not within the target range, hence, no remodeling
                        elem_field(ne_radius_in,ne)=R0
                      endif
                    endif ! nn.eq.1
                    if(nn.eq.2) then ! nn.eq.2
                      if((R0.gt.narrow_rad_one).and.(R0.lt.0.5_dp)) then
                        if(R0.lt.narrow_rad_two) then
                          elem_field(ne_radius_out,ne)=narrow_factor*R0
                        else
                          elem_field(ne_radius_out,ne)=R0
                        endif
                      else ! Not within the target range, hence, no remodeling
                        elem_field(ne_radius_out,ne)=R0
                      endif
                    endif ! nn.eq.2
                else !ptm>ptmmax
                    if(nn.eq.1) then ! nn.eq.1
                      if((R0.gt.narrow_rad_one).and.(R0.lt.0.5_dp)) then ! Hypertophy+narrowing effect
                        if(R0.lt.0.05_dp) then ! only Narrow_factor
                          elem_field(ne_radius_in,ne) = narrow_factor*R0*((Ptm*alt_hyp*alt_fib*elasticity_parameters(1))+1.d0)
                        elseif(R0.gt.narrow_rad_two) then ! hypertophy only
                          elem_field(ne_radius_in,ne)=R0*((elasticity_parameters(2)*alt_hyp*elasticity_parameters(1))+1.d0)
                        else ! Both hypertophy and narrowing
                          elem_field(ne_radius_in,ne)=narrow_factor*R0*(elasticity_parameters(2)* &
                          alt_hyp*alt_fib*elasticity_parameters(1)+1.d0)
                        endif
                      elseif(R0.gt.0.5_dp) then ! Not within the target range, hence, no remodeling
                          elem_field(ne_radius_in,ne)=R0*((elasticity_parameters(2)*(elasticity_parameters(1)))+1.d0)
                      else
                        elem_field(ne_radius_in,ne)=R0
                      endif
                    endif ! nn.eq.1
                    if(nn.eq.2) then ! nn.eq.2
                      if((R0.gt.narrow_rad_one).and.(R0.lt.0.5_dp)) then ! Hypertrophy+narrowing effect
                        if(R0.lt.0.05_dp) then ! only Narrow_factor
                          elem_field(ne_radius_out,ne)=narrow_factor*R0*((elasticity_parameters(2)*elasticity_parameters(1))+1.d0)
                        elseif(R0.gt.narrow_rad_two) then ! hypertophy only
                          elem_field(ne_radius_out,ne)=R0*((elasticity_parameters(2)*alt_hyp*elasticity_parameters(1))+1.d0)
                        else ! Both hypertophy and narrowing
                          elem_field(ne_radius_out,ne)=narrow_factor*R0*((elasticity_parameters(2)* &
                          alt_hyp*alt_fib*elasticity_parameters(1))+1.d0)
                        endif
                      elseif(R0.gt.0.5_dp) then ! Not within the target range, hence, no remodeling
                          elem_field(ne_radius_out,ne)=R0*((elasticity_parameters(2)*(elasticity_parameters(1)))+1.d0)
                      else
                        elem_field(ne_radius_out,ne)=R0
                      endif ! radius criteria
                    endif ! nn.eq.2
                    if (allocated(compl_list)) then ! compl allocated
                      do no_compl_ne=1,count(compl_list.ne.0)
                        !pause
                        if(ne.eq.compl_list(no_compl_ne))then
                          if(nn.eq.1) elem_field(ne_radius_in,compl_list(no_compl_ne))=R0!*((Ptm*elasticity_parameters(1)*0.0)+1.d0)
                          if(nn.eq.2) elem_field(ne_radius_out,compl_list(no_compl_ne))=R0!*((Ptm*elasticity_parameters(1)*0.0)+1.d0)
                        endif
                      end do
                    endif ! compl allocated
                endif !ptm
              else ! FOUND
                if(Ptm.lt.elasticity_parameters(2))then
                  if(nn.eq.1) elem_field(ne_radius_in,ne)=R0*((Ptm*elasticity_parameters(1))+1.d0)
                  if(nn.eq.2) elem_field(ne_radius_out,ne)=R0*((Ptm*elasticity_parameters(1))+1.d0)
                elseif(Ptm.lt.0.0_dp)then !Ptm
                  if(Ptm.lt.0)write(*,*) 'Transmural pressure < zero',ne,Ptm,Pblood,Ppl
                  if(nn.eq.1) elem_field(ne_radius_in,ne)=R0
                  if(nn.eq.2) elem_field(ne_radius_out,ne)=R0
                else !ptm>ptmmax
                  if(nn.eq.1) elem_field(ne_radius_in,ne)=R0*((elasticity_parameters(2)*elasticity_parameters(1))+1.d0)
                  if(nn.eq.2) elem_field(ne_radius_out,ne)=R0*((elasticity_parameters(2)*elasticity_parameters(1))+1.d0)
                endif
              endif ! FOUND
              if (allocated(compl_list)) then ! compl allocated
                do no_compl_ne=1,count(compl_list.ne.0)
                  ! pause
                  if(ne.eq.compl_list(no_compl_ne))then
                    if(nn.eq.1) elem_field(ne_radius_in,compl_list(no_compl_ne))=R0*((Ptm*elasticity_parameters(1)*0.0)+1.d0)
                    if(nn.eq.2) elem_field(ne_radius_out,compl_list(no_compl_ne))=R0*((Ptm*elasticity_parameters(1)*0.0)+1.d0)
                  endif
                end do
              endif ! compl allocated
          elseif(vessel_type.eq.'elastic_hooke')then
            h=elasticity_parameters(2)*R0
            if(nn.eq.1) elem_field(ne_radius_in,ne)=R0+3.0_dp*R0**2*Ptm/(4.0_dp*elasticity_parameters(1)*h)
            if(nn.eq.2) elem_field(ne_radius_out,ne)=R0+3.0_dp*R0**2*Ptm/(4.0_dp*elasticity_parameters(1)*h) !vessel type
          else
            print *, 'no vessel type defined, assuming rigid'
            if(nn.eq.1) elem_field(ne_radius_in,ne)=R0
            if(nn.eq.2) elem_field(ne_radius_out,ne)=R0 !vessel type
          endif !vessel type
        else !other than arteries
          if(Ptm.lt.elasticity_parameters(2))then
            if(nn.eq.1) elem_field(ne_radius_in,ne)=R0*((Ptm*elasticity_parameters(1))+1.d0)
            if(nn.eq.2) elem_field(ne_radius_out,ne)=R0*((Ptm*elasticity_parameters(1))+1.d0)
          elseif(Ptm.lt.0.0_dp)then !Ptm
            if(Ptm.lt.0)write(*,*) 'Transmural pressure < zero',ne,Ptm,Pblood,Ppl
            if(nn.eq.1) elem_field(ne_radius_in,ne)=R0
            if(nn.eq.2) elem_field(ne_radius_out,ne)=R0
          else !ptm>ptmmax
            if(nn.eq.1) elem_field(ne_radius_in,ne)=R0*((elasticity_parameters(2)*elasticity_parameters(1))+1.d0)
            if(nn.eq.2) elem_field(ne_radius_out,ne)=R0*((elasticity_parameters(2)*elasticity_parameters(1))+1.d0)
          endif
        endif ! if artery
      enddo!nn
    enddo!ne
  endif

call enter_exit(sub_name,2)
end subroutine calc_press_area
!##############################################################################
!
!*map_solution_to_mesh* maps the solution array to appropriate nodal and element fields
subroutine map_solution_to_mesh(prq_solution,depvar_at_elem,depvar_at_node,mesh_dof)

    integer, intent(in) :: mesh_dof
    real(dp),intent(in) ::  prq_solution(mesh_dof,2)
    integer,intent(in) :: depvar_at_elem(0:2,2,num_elems)
    integer,intent(in) :: depvar_at_node(num_nodes,0:2,2)
    !local variables
    integer :: np,ne,ny


    character(len=60) :: sub_name
    sub_name = 'map_solution_to_mesh'
    call enter_exit(sub_name,1)
      do  ne=1,num_elems
        ny=depvar_at_elem(1,1,ne)
        elem_field(ne_Qdot,ne)=prq_solution(ny,1)
      enddo !elems
      do np=1,num_nodes
        ny=depvar_at_node(np,0,1)
        node_field(nj_bv_press,np)=prq_solution(ny,1)
      enddo

    call enter_exit(sub_name,2)
end subroutine map_solution_to_mesh

!##############################################################################!##############################################################################
!
!*is_downstream* checks to see if a certain element is downstream of another element
!!!!!!!
! This function is useful for doing regional remodeling where you do not want to
! remodel vessels downstream of an occlusion (needed for post-PEA modelling)
!!!!!!!
recursive function is_downstream(ne,occ_ne) result(FOUND)

    integer, intent(in) :: ne,occ_ne
    !local variables
    logical :: FOUND
    integer :: ne_temp,np
    character(len=60) :: sub_name

    FOUND = .False.
    np=elem_nodes(1,ne)

    if (ne == occ_ne) then
        FOUND = .True.
    else if (ne == 1) then
        FOUND = .false.
    else
        ne_temp = elems_at_node(np,1)
        FOUND = is_downstream(ne_temp, occ_ne)
    end if

    call enter_exit(sub_name,2)
end function is_downstream

!##############################################################################
!
!*map_flow_to_terminals* maps the solution array to appropriate nodal and element fields
subroutine map_flow_to_terminals

  integer :: nu,ne,np
    character(len=60) :: sub_name
    sub_name = 'map_flow_to_terminals'
    call enter_exit(sub_name,1)

    do nu=1,num_units
      ne=units(nu)
      np=elem_nodes(2,ne)
      unit_field(nu_perf,nu)=elem_field(ne_Qdot,ne)
      unit_field(nu_blood_press,nu)=node_field(nj_bv_press,np)
    enddo

    call enter_exit(sub_name,2)
end subroutine map_flow_to_terminals
!
!
!
!*calculate_ppl* calculates pleural pressure at a node
!
subroutine calculate_ppl(np,grav_vect,mechanics_parameters,Ppl)

  integer, intent(in) :: np
    real(dp), intent(in) :: mechanics_parameters(2)
    real(dp), intent(in) :: grav_vect(3)
    real(dp), intent(out) :: Ppl
    !Local variables
    integer :: nj
    real(dp) :: G_PLEURAL,HEIGHT(3)
    character(len=60) :: sub_name

    sub_name = 'calculate_ppl'

    call enter_exit(sub_name,1)
        G_PLEURAL=0.0_dp    !gravitational force
        do nj=1,3
          HEIGHT(nj)=node_xyz(nj,np)-node_xyz(nj,1) !ARC - where to put grav reference height?
          G_PLEURAL=G_PLEURAL+mechanics_parameters(2)*grav_vect(nj)*9810.d0*HEIGHT(nj) !kg
        enddo
        Ppl=mechanics_parameters(1)-G_PLEURAL !Pa

    call enter_exit(sub_name,2)
end subroutine calculate_ppl

!
!##################################################################
!
subroutine get_variable_offset(depvar,mesh_dof,FIX,offset)
!*Description*: This subroutine returns the number of fixed variables in the depvar array that came before the input depvar
    integer, intent(in) :: depvar, mesh_dof
    logical, intent(in) :: FIX(mesh_dof)
    integer, intent(inout) :: offset

    integer :: depvar1

    offset = 0
    do depvar1 = 1,depvar
        if(FIX(depvar1))then
            offset = offset + 1
        endif
    enddo

end subroutine get_variable_offset

end module pressure_resistance_flow
