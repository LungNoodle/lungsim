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

  !Module parameters

  !Module types

  !Module depvar

  !Interfaces
  private
  public evaluate_prq,calculate_ppl
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
    real(dp), allocatable :: SparseVal(:)
    real(dp), allocatable :: RHS(:)
    integer :: num_vars,NonZeros,MatrixSize
    integer :: AllocateStatus

    real(dp), intent(in) :: remodeling_grade ! Remodeling if applicable, 0 stands for healthy and anything from 10 to 100 is hypertension
    real(dp), allocatable :: prq_solution(:,:),solver_solution(:)
    real(dp) :: viscosity,density,inlet_bc,outlet_bc,inletbc,outletbc,grav_vect(3),gamma,total_resistance,ERR
    logical, allocatable :: FIX(:)
    logical :: ADD=.FALSE.,CONVERGED=.FALSE.
    character(len=60) :: sub_name,mesh_type,vessel_type,mechanics_type,bc_type
    integer :: grav_dirn,no,depvar,KOUNT,nz,ne,SOLVER_FLAG,ne0,ne1,nj
    real(dp) :: MIN_ERR,N_MIN_ERR,elasticity_parameters(3),mechanics_parameters(2),grav_factor,P1
    real(dp) :: P2,Q01,Rin,Rout,x_cap,y_cap,z_cap,Ppl,LPM_R,Lin,Lout
    integer :: update_flow_nzz_row

    sub_name = 'evaluate_prq'
    call enter_exit(sub_name,1)
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

    KOUNT=0
!! Calculate resistance of each element
   call calculate_resistance(viscosity,KOUNT)

!! Calculate sparsity structure for solution matrices
    !Determine size of and allocate solution vectors/matrices
    call calc_sparse_size(mesh_dof,depvar_at_elem,depvar_at_node,FIX,NonZeros,MatrixSize)
    allocate (SparseCol(NonZeros), STAT = AllocateStatus)!Note we should be able to calculate the nonzeros and matrix size analtyically then we wont need this.
    if (AllocateStatus /= 0) STOP "*** Not enough memory for SparseCol array ***"
    allocate (SparseRow(MatrixSize+1), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for SparseRow array ***"
    allocate (SparseVal(NonZeros), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for SparseVal array ***"
    allocate (RHS(MatrixSize), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for RHS array ***"
    allocate (solver_solution(MatrixSize), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (update_resistance_entries(num_elems), STAT = AllocateStatus)
    if (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    update_resistance_entries = 0
    !calculate the sparsity structure
    call calc_sparse_1dtree(bc_type,density,FIX,grav_vect,mesh_dof,depvar_at_elem, &
        depvar_at_node,NonZeros,MatrixSize,SparseCol,SparseRow,SparseVal,RHS, &
        prq_solution,update_resistance_entries,update_flow_nzz_row)
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
        if(bc_type.eq.'flow')then !update RHS to account for element resistance
          do ne=1,num_elems
            depvar = depvar_at_elem(1,1,ne)
            if(FIX(depvar))then
              RHS(update_flow_nzz_row) = prq_solution(depvar,1)*elem_field(ne_resist,ne)
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

         ERR=ERR/MatrixSize !sum of error divided by no of unknown depvar
         if(ERR.LE.1.d-6.AND.(KOUNT.NE.1))then
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
        depvar_at_node,NonZeros,MatrixSize,SparseCol,SparseRow,SparseVal,RHS,&
        prq_solution,update_resistance_entries,update_flow_nzz_row)

    character(len=60) :: bc_type
    real(dp), intent(in) :: density
    real(dp), intent(in) :: grav_vect(3)
    integer, intent(in) :: mesh_dof
    logical, intent(in) :: FIX(mesh_dof)
    integer,intent(in) :: depvar_at_elem(0:2,2,num_elems)
    integer,intent(in) :: depvar_at_node(num_nodes,0:2,2)

    integer, intent(in) :: NonZeros,MatrixSize
    integer, intent(inout) :: SparseCol(NonZeros)
    integer, intent(inout) :: SparseRow(MatrixSize+1)
    real(dp), intent(inout) :: SparseVal(NonZeros)
    real(dp), intent(inout) :: RHS(MatrixSize)
    real(dp), intent(inout) :: prq_solution(mesh_dof,2)
    integer, intent(inout) :: update_resistance_entries(num_elems)
    integer, intent(inout) :: update_flow_nzz_row
!local variables
    integer :: ne,nn,np,np1,np2,depvar,depvar1,depvar2,depvar3,flow_var,fixed_var_index,offset,nzz,&
      nzz_row,ne2,noelem2,ne3
    logical :: FlowBalancedNodes(num_nodes)
    logical :: NodePressureDone(num_nodes)
    logical :: ElementPressureEquationDone(num_elems)
    logical :: elem_found,one_node_balanced
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

    do ne=1,num_elems
      !look at pressure variables at each node
      do nn=1,2 !2 nodes in 1D element
        np=elem_nodes(nn,ne)
        depvar = depvar_at_node(np,1,1)
        if((.NOT.NodePressureDone(np)).AND.(.NOT.FIX(depvar)))then !check if this node is not fixed and hasn't already been processed (as nodes are shared between elements)
            ne2=0
            if(nn.EQ.1)then !first node of the element
                ne2=ne! use the current element
            elseif(nn.EQ.2)then !second node of the element
                if((bc_type.EQ.'pressure').OR.(.NOT.ElementPressureEquationDone(ne)))then !if bc_type is pressure or element pressure equation for the current element hasn't been used
                    ne2=ne! use the current element
                else
                    !look for another element connected to this node with pressure equation that hasn't been used
                    if (elems_at_node(np,0).GT.1)then
                        elem_found=.FALSE.
                        noelem2 = 1
                        do while ((.NOT.elem_found).AND.(noelem2.LE.elems_at_node(np,0)))
                            ne3=elems_at_node(np,noelem2)
                            if((ne3.NE.ne).AND.(.NOT.ElementPressureEquationDone(ne3)))then
                                ne2 = ne3
                                elem_found=.TRUE.
                            endif
                            noelem2 = noelem2 + 1
                        end do
                    endif
                endif
            endif
            if(ne2.GT.0)then
                !do the pressure equation for element ne2
                !pressure for node 1 - pressure for node 2 - resistance * flow at element ne2 = 0
                np1=elem_nodes(1,ne2)
                depvar1=depvar_at_node(np1,1,1) !pressure variable for first node
                np2=elem_nodes(2,ne2) !second node
                depvar2=depvar_at_node(np2,1,1) !pressure variable for second node
                depvar3=depvar_at_elem(0,1,ne2) !flow variable for element
                grav=0.d0
                if(elem_field(ne_group,ne2).eq.1.0_dp)then
                elseif(elem_ordrs(no_gen,ne2).eq.1)then !gravitational head not applied in inlets
                else
                  do nj=1,3
                    grav=grav+density*grav_vect(nj)*9810.0_dp*(node_xyz(nj,elem_nodes(1,ne2))-node_xyz(nj,elem_nodes(2,ne2)))!rho g L cos theta (Pa)
                  enddo
                endif
                if(FIX(depvar1))then !checking if pressure at 1st node is fixed
                    !store known variable - inlet pressure
                    RHS(nzz_row) = -prq_solution(depvar1,1) + grav
                else
                    !unknown variable -pressure for node 1
                    call get_variable_offset(depvar1,mesh_dof,FIX,offset)
                    SparseCol(nzz) = depvar1 - offset !variable number
                    SparseVal(nzz)=1.0_dp !variable coefficient
                    nzz=nzz+1 !next column
                    RHS(nzz_row) = grav
                endif
                if(FIX(depvar2))then !checking if pressure at 2nd node is fixed
                    !store known variable - outlet pressure
                    RHS(nzz_row) = prq_solution(depvar2,1) + grav
                else
                    !unknown variable - pressure for node 2
                    call get_variable_offset(depvar2,mesh_dof,FIX,offset)
                    SparseCol(nzz) = depvar2 - offset !variable number
                    SparseVal(nzz)=-1.0_dp !variable coefficient
                    nzz=nzz+1 !next column
                endif
                if(FIX(depvar3))then !checking if flow at element ne2 is fixed
                    !store known variable - inlet flow * resistance for element ne
                    RHS(nzz_row) = prq_solution(depvar3,1)*elem_field(ne_resist,ne2)
                    update_flow_nzz_row = nzz_row
                else
                    !unknown flow
                    call get_variable_offset(depvar3,mesh_dof,FIX,offset)
                    SparseCol(nzz) = depvar3-offset !variable position in the unknown variable vector
                    SparseVal(nzz)=-elem_field(ne_resist,ne2) !variable coefficient = resistance for element ne2
                    update_resistance_entries(ne2) = nzz
                    nzz=nzz+1 !next column
                endif
                nzz_row=nzz_row+1 !store next row position
                    SparseRow(nzz_row)=nzz
                NodePressureDone(np) = .TRUE.
                ElementPressureEquationDone(ne2) = .TRUE.
            endif
        endif
      enddo !nn

      !look at flow variable for the element
      flow_var = depvar_at_elem(0,1,ne)
      if(.NOT.FIX(flow_var))then !don't do anything if flow is fixed
        one_node_balanced = .FALSE.
        !check if node 1 or node 2 are unbalanced
        do nn=1,2 !do flow balance for each element node
          np = elem_nodes(nn,ne)
          if((elems_at_node(np,0).GT.1).AND.(.NOT.FlowBalancedNodes(np)))then !if there is more than one element at a node and the node is not already flow balanced
            if((bc_type.EQ.'pressure').OR.((bc_type.EQ.'flow').AND.(.NOT.one_node_balanced)))then !do just one flow balance equation for bc_type flow
                !go through each element connected to node np and add the conservation of flow equation for the elements
                do noelem2=1,elems_at_node(np,0)
                    ne2=elems_at_node(np,noelem2)
                    depvar=depvar_at_elem(1,1,ne2)
                    flow_term = 0
                    if(np.EQ.elem_nodes(2,ne2))then !end node
                        flow_term = 1.0_dp
                    elseif(np.EQ.elem_nodes(1,ne2))then !start node
                        flow_term = -1.0_dp
                    endif
                    if(FIX(depvar))then
                        RHS(nzz_row)=-prq_solution(depvar,1)*flow_term
                    else
                        !populate SparseCol and SparseVal
                        call get_variable_offset(depvar,mesh_dof,FIX,offset)
                        SparseCol(nzz) = depvar - offset
                        SparseVal(nzz) = flow_term
                        nzz = nzz + 1
                    endif
                enddo
                FlowBalancedNodes(np) = .TRUE.
                nzz_row=nzz_row+1 !store next row position
                    SparseRow(nzz_row)=nzz
                    one_node_balanced = .TRUE.
                endif !checking bc_type
          endif !flow at node np is unbalanced
        enddo !nn

        !if flow balancing hasn't been done for any node for element ne and pressure equation hasn't already been done, do the pressure equation for the element
        if((.NOT.one_node_balanced).AND.(.NOT.ElementPressureEquationDone(ne)))then

            !do the pressure equation for element ne
            !pressure for node 1 - pressure for node 2 - resistance * flow at element ne = 0
            np1=elem_nodes(1,ne)
            depvar1=depvar_at_node(np1,1,1) !pressure variable for first node
            np2=elem_nodes(2,ne) !second node
            depvar2=depvar_at_node(np2,1,1) !pressure variable for second node

            !unknown variable -pressure for node 1
            call get_variable_offset(depvar1,mesh_dof,FIX,offset)
            SparseCol(nzz) = depvar1 - offset !variable number
            SparseVal(nzz)=1.0_dp !variable coefficient
            nzz=nzz+1 !next column

            if(FIX(depvar2))then !checking if pressure at 2nd node is fixed
                !store known variable - outlet pressure
                RHS(nzz_row) = prq_solution(depvar2,1)
            else
                !unknown variable - pressure for node 2
                call get_variable_offset(depvar2,mesh_dof,FIX,offset)
                SparseCol(nzz) = depvar2 - offset !variable number
                SparseVal(nzz)=-1.0_dp !variable coefficient
                nzz=nzz+1 !next column
            endif

            !unknown flow
            call get_variable_offset(flow_var,mesh_dof,FIX,offset)
            SparseCol(nzz) = flow_var-offset !variable position in the unknown variable vector
            SparseVal(nzz)=-elem_field(ne_resist,ne) !variable coefficient = resistance for element ne
            update_resistance_entries(ne) = nzz
            nzz=nzz+1 !next column

            nzz_row=nzz_row+1 !store next row position
            SparseRow(nzz_row)=nzz

            ElementPressureEquationDone(ne) = .TRUE.
        endif
      endif
    enddo !ne

    call enter_exit(sub_name,2)
  end subroutine calc_sparse_1dtree


!
!##################################################################
!
!*calc_sparse_size:* Calculates sparsity sizes

subroutine calc_sparse_size(mesh_dof,depvar_at_elem,depvar_at_node,FIX,NonZeros,MatrixSize)

  integer, intent(in) :: mesh_dof
    integer,intent(in) :: depvar_at_elem(0:2,2,num_elems)
    integer,intent(in) :: depvar_at_node(num_nodes,0:2,2)
    logical, intent(in) :: FIX(mesh_dof)
    integer :: NonZeros,MatrixSize
!local variables
    integer :: i,ne,np,fixed_variables, fixed_flows, fixed_pressures
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

    !get count of fixed flows
    fixed_flows = 0
    do ne=1,num_elems
        if(FIX(depvar_at_elem(1,1,ne)))then
            fixed_flows = fixed_flows + 1
        endif
    enddo

    fixed_pressures = fixed_variables - fixed_flows

    !count of pressure equations = (number of elements * 3 variables in each equation) - fixed pressures - fixed flows
    NonZeros = num_elems*3 - fixed_pressures - fixed_flows
    !count of conservation of flow equations = sum of elements connected to nodes which have at least 2 connected elements - fixed flows
    do np=1, num_nodes
        if(elems_at_node(np,0).gt.1)then
            NonZeros = NonZeros + elems_at_node(np,0)
        endif
    enddo
    NonZeros = NonZeros - fixed_flows
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
    integer :: nj,np,ne,ny,nn
    real(dp) :: h,Ptm,R0,Pblood,Ppl,counter,cc1,cc2,cc3
    real(dp) :: alt_hyp,alt_fib,prox_fib,narrow_rad_one,narrow_rad_two,narrow_factor,prune_rad,prune_fraction,counter1,counter2
    integer,allocatable :: templss(:)

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
          if(nn.eq.1) then
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
          endif
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
          endif
        endif
      if(vessel_type.eq.'elastic_g0_beta') then
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
        endif
      elseif(vessel_type.eq.'elastic_alpha') then
        if(elem_field(ne_group,ne).eq.0.0_dp) then !only applying on arteries
          if(Ptm.lt.elasticity_parameters(2))then
            if(nn.eq.1) then
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
            endif
            if(nn.eq.2) then
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
            endif
          elseif(Ptm.lt.0.0_dp)then !Ptm
            if(Ptm.lt.0) write(*,*) 'Transmural pressure < zero',ne,Ptm,Pblood,Ppl
            if(nn.eq.1) then
              if((R0.gt.narrow_rad_one).and.(R0.lt.0.5)) then ! Hypertophy+narrowing effect
                if(R0.lt.narrow_rad_two) then !only narrowing
                  elem_field(ne_radius_in,ne)=narrow_factor*R0
                else
                  elem_field(ne_radius_in,ne)=R0
                endif
              else  ! Not within the target range, hence, no remodeling
                elem_field(ne_radius_in,ne)=R0
              endif
            endif
            if(nn.eq.2) then
              if((R0.gt.narrow_rad_one).and.(R0.lt.0.5_dp)) then
                if(R0.lt.narrow_rad_two) then
                  elem_field(ne_radius_out,ne)=narrow_factor*R0
                else
                  elem_field(ne_radius_out,ne)=R0
                endif
              else ! Not within the target range, hence, no remodeling
                elem_field(ne_radius_out,ne)=R0
              endif
            endif
          else !ptm>ptmmax
            if(nn.eq.1) then
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
            endif
            if(nn.eq.2) then
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
              endif
            endif
          endif
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
