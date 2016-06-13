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
  use solve, only: BICGSTAB_LinSolv,pmgmres_ilu_cr
  implicit none
  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private
  public evaluate_prq
contains
!###################################################################################
!
!*evaluate_PRQ:* Solves for pressure and flow in a rigid or compliant tree structure
  subroutine evaluate_prq
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_PRQ" :: EVALUATE_PRQ
    use indices
    use arrays,only: dp,num_elems,num_nodes,elem_field
    use diagnostics, only: enter_exit
    !local variables
    integer :: mesh_dof,variable_types
    integer, allocatable :: mesh_from_variables(:,:,:)
    integer, allocatable :: variables_at_node(:,:,:)
    integer, allocatable :: variables_at_elem(:,:,:)
    integer, dimension(0:2,2) :: variable_totals
    integer, allocatable :: SparseCol(:)
    integer, allocatable :: SparseRow(:)
    integer, allocatable :: update_resistance_entries(:)
    real(dp), allocatable :: SparseVal(:)
    real(dp), allocatable :: RHS(:)
    integer :: num_vars,NonZeros,MatrixSize
    integer :: AllocateStatus

    real(dp), allocatable :: prq_solution(:,:),solver_solution(:)
    real(dp) :: viscosity,density,inletbc,outletbc,grav_vect(3),gamma,total_resistance,ERR
    logical, allocatable :: FIX(:)
    logical :: ADD=.FALSE.,CONVERGED=.FALSE.
    character(len=60) :: sub_name,mesh_type,vessel_type,mechanics_type,bc_type,g_type
    integer :: g_dirn,no,variable,KOUNT,nz,ne,SOLVER_FLAG
    real(dp) :: MIN_ERR,N_MIN_ERR,ptrans,beta,Go_artery,Ptm_max

    sub_name = 'evaluate_prq'
    call enter_exit(sub_name,1)
!!---------DESCRIPTION OF MODEL TypeS -----------
!mesh_type: can be simple_tree, full_ladder, full_sheet, full_tube The first can be airways, arteries, veins but no special features at the terminal level, the last one has arteries and veins connected by capillary units of some type (lung ladder acinus, lung sheet capillary bed, capillaries are just tubes represented by an element)
!
!
!vessel_type: rigid (no pressure area relationship), elastic_g_beta (non linear pressure area relationship controlled by g0,beta), elastic_alpha (linear pressure area relationship controlled by alpha.
!
!mechanics type: const (constant pressure external to vessels), linear (assumed gradient in Ppl along gravitational direction), mechanics (Ppl determined by solution to a mechanics model)
!bc_type: can be pressure (at inlet and outlets) or flow (flow at inlet pressure at outlet).
mesh_type='simple_tree'
vessel_type='linear_compliance'!rigid'
mechanics_type='linear'
bc_type='pressure'
g_type='off'
g_dirn=0

if (g_type.eq.'off') then
           grav_vect=0.d0
elseif (g_dirn.eq.1) then
           grav_vect(1)=1.d0
elseif (g_dirn.eq.-1) then
           grav_vect(1)=-1.d0
elseif (g_dirn.eq.2) then
           grav_vect(2)=1.d0
elseif (g_dirn.eq.-2) then
           grav_vect(2)=-1.d0
elseif (g_dirn.eq.3) then
           grav_vect(3)=1.d0
elseif (g_dirn.eq.-3) then
           grav_vect(3)=-1.d0
else
           write(*,*),"Posture not recognised (nogravity, upright, inverted, prone, supine?)"
endif

inletbc=2266.0_dp
outletbc=666.7_dp
ptrans=5.33_dp
beta=1.0_dp
Go_artery=6.67_dp
Ptm_max=32.0_dp

!!---------DESCRIPTION OF IMPORTANT PARAMETERS-----------
!viscosity: fluid viscosity
!density:fluid density
!

!!!set the default values for the parameters that control the prq simulation, these should be controlled by user input
    call read_params_evaluate_prq(viscosity,density,gamma)


!! Allocate memory to variable arrays
    mesh_dof=num_elems+num_nodes
    variable_types=2 !pressure/flow
    allocate (mesh_from_variables(0:2,mesh_dof,0:2), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for mesh_from_variables array ***"
    allocate (variables_at_elem(0:2,variable_types,num_elems), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for variables_at_elem array ***"
    allocate (variables_at_node(num_nodes,0:2,variable_types), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for variables_at_node array ***"
    allocate (prq_solution(mesh_dof,2), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for prq_solution array ***"
    prq_solution=0.0_dp !initialise
    allocate (FIX(mesh_dof), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for FIX array ***"



!! Setting up mappings between nodes, elements and solution variables
    call calc_variable_maps(mesh_from_variables,variables_at_elem,&
                variable_totals,variables_at_node,mesh_dof,num_vars)

!! Define boundary conditions
    !first call to define inlet boundary conditions
    call boundary_conditions(ADD,FIX,bc_type,g_type,grav_vect,density,inletbc,outletbc,&
      variables_at_node,variables_at_elem,prq_solution,mesh_dof)
    !second call if simple tree need to define pressure bcs at all terminal branches
    if(mesh_type.eq.'simple_tree')then
        ADD=.TRUE.
        call boundary_conditions(ADD,FIX,bc_type,g_type,grav_vect,density,inletbc,outletbc,&
            variables_at_node,variables_at_elem,prq_solution,mesh_dof)   
    endif

 
!! Calculate resistance of each element
   call calculate_resistance(density,gamma,viscosity)
        
!! Calculate sparsity structure for solution matrices
    !Determine size of and allocate solution vectors/matrices
    call calc_sparse_size(mesh_dof,variables_at_elem,variables_at_node,FIX,NonZeros,MatrixSize)
    allocate (SparseCol(NonZeros), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for SparseCol array ***"
    allocate (SparseRow(MatrixSize+1), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for SparseRow array ***"
    allocate (SparseVal(NonZeros), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for SparseVal array ***"
    allocate (RHS(MatrixSize), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for RHS array ***"
    allocate (solver_solution(MatrixSize), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"
    allocate (update_resistance_entries(num_elems), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for solver_solution array ***"

    call calc_sparse_1dtree(mesh_dof,variables_at_elem,variables_at_node,FIX,NonZeros,MatrixSize,SparseCol,&
	SparseRow,SparseVal,RHS,prq_solution,update_resistance_entries)

!!! --ITERATIVE LOOP--
    MIN_ERR=1.d10
    N_MIN_ERR=0
    KOUNT=0
    do while(.NOT.CONVERGED)
      KOUNT=KOUNT+1
      print*, 'Outer loop iterations:',KOUNT
!!! Initialise solution vector based on bcs and rigid vessel resistance
      if(KOUNT.eq.1)then!set up boundary conditions
        if(bc_type.eq.'pressure')then
          call tree_resistance(total_resistance)
          call initialise_solution(inletbc,outletbc,(inletbc-outletbc)/total_resistance, &
              mesh_dof,prq_solution,variables_at_node,variables_at_elem,FIX)
          !move initialisation to solver solution (skipping BCs).
          no=0
     	  do variable=1,mesh_dof !loop over mesh dofs
            if(.NOT.FIX(variable))then
              no=no+1
              solver_solution(no)=prq_solution(variable,1)
            endif
          enddo !mesh_dof
        else!flow BCs to be implemented
        endif
      else!Need to update just the resistance values in the solution matrix
        do ne=1,num_elems !update for all ne
          nz=update_resistance_entries(ne)
          SparseVal(nz)=-elem_field(ne_resist,ne) !Just updating resistance
        enddo
      endif!first or subsequent iteration
!! ----CALL SOLVER----

      call pmgmres_ilu_cr(MatrixSize, NonZeros, SparseRow, SparseCol, SparseVal, &
         solver_solution, RHS, 500, 500,1.d-5,1.d-4,SOLVER_FLAG)
       if(SOLVER_FLAG == 0)then 
          print *, 'Warning: pmgmres has reached max iterations. Solution may not be valid if this warning persists'
       elseif(SOLVER_FLAG ==2)then
          print *, 'ERROR: pmgmres has failed to converge'
          deallocate (mesh_from_variables, STAT = AllocateStatus)
          deallocate (variables_at_elem, STAT = AllocateStatus)
          deallocate (variables_at_node, STAT = AllocateStatus)
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
      do variable=1,mesh_dof
         if(.NOT.FIX(variable)) THEN
            no=no+1
            prq_solution(variable,2)=prq_solution(variable,1) !temp storage of previous solution
            prq_solution(variable,1)=solver_solution(no) !new pressure & flow solutions
            if(DABS(prq_solution(variable,1)).GT.0.d-6)THEN
               ERR=ERR+(prq_solution(variable,2)-prq_solution(variable,1))**2.d0/prq_solution(variable,1)**2
            endif
         endif
      enddo !no2
!rigid vessels no need to update - tag as converged and exit
      if(vessel_type.eq.'rigid')then
        ERR=0.0_dp
        CONVERGED=.TRUE.
      else
!Update vessel radii based on predicted pressures and then update resistance through tree
        call calc_press_area(KOUNT,variables_at_node,prq_solution,mesh_dof,ptrans,beta,Go_artery,Ptm_max)
        call calculate_resistance(viscosity,density,gamma)

!Put the ladder stuff here --> See solve11.f

         ERR=ERR/MatrixSize !sum of error divided by no of unknown variables
         if(ERR.LE.1.d-6.AND.(KOUNT.NE.1))then
           CONVERGED=.TRUE.
            write(*,*),"Convergence achieved after",KOUNT,"iterations",ERR
         else !if error not converged
            if(ERR.GE.MIN_ERR) then
              N_MIN_ERR=N_MIN_ERR+1
            else
              MIN_ERR=ERR
            endif
            write(*,*),"Not converged, error =",ERR
         endif !ERR not converged
      endif!vessel type
    enddo !notconverged

!need to write solution to element/nodal fields for export

    deallocate (mesh_from_variables, STAT = AllocateStatus)
    deallocate (variables_at_elem, STAT = AllocateStatus)
    deallocate (variables_at_node, STAT = AllocateStatus)
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
!*read_params_evaluate_prq:* Reads in important parameters for PRQ problems
  subroutine read_params_evaluate_prq(viscosity,density,gamma)
    use arrays,only: dp
    use diagnostics, only: enter_exit
    !local parameters
    real(dp) :: viscosity,density,gamma
    character(len=60) :: sub_name

    sub_name = 'read_params_evaluate_prq'
    call enter_exit(sub_name,1)
      density=0.10500e-02_dp
      viscosity=0.33600e-02_dp
      gamma = 0.327_dp !=1.85/(4*sqrt(2))


    call enter_exit(sub_name,2)
  end subroutine read_params_evaluate_prq
!
!###################################################################################
!
!*boundary_conditions:* Defines boundary conditions for prq problems
 subroutine boundary_conditions(ADD,FIX,bc_type,g_type,grav_vect,density,inletbc,outletbc,&
       variables_at_node,variables_at_elem,prq_solution,mesh_dof)
 use arrays,only: dp,num_elems,num_nodes,elems,elem_nodes,elem_cnct,node_xyz,units,&
        num_units
    use diagnostics, only: enter_exit

    integer :: mesh_dof
    integer :: variables_at_elem(0:2,2,num_elems)
    integer :: variables_at_node(num_nodes,0:2,2)
    real(dp) :: prq_solution(mesh_dof,2),inletbc,outletbc,density,grav_vect(3)
    logical:: ADD
    logical :: FIX(mesh_dof)
    character(len=60) ::bc_type,g_type
 
  ! Local variables
    integer :: nonode,np,noelem,ne,ny1,ny2,nj,np_in
    double precision :: grav
    character(len=60) :: sub_name

  sub_name = 'boundary_conditions'
  call enter_exit(sub_name,1)
  IF(.NOT.ADD)THEN
     ! Initial values
     FIX(1:mesh_dof)=.FALSE.
     prq_solution = 0
     ! Fixed boundary conditions  
     ! These are inlet BCs, apply to all inlet BCs (there should only be one)
     DO ne=1,num_elems
        !ne=elems(noelem)
        IF (elem_cnct(-1,0,ne) == 0) THEN !Entry element
           IF(BC_TYPE == 'pressure')THEN          
              np=elem_nodes(1,ne)
              ny1=variables_at_node(np,1,1) !for fixed pressure BC
              FIX(ny1)=.TRUE. !set fixed
              prq_solution(ny1,1)=inletbc !Putting BC value into solution array
              np_in=np !inlet node set here, gravity reference to this point
           ELSE IF(BC_TYPE == 'flow')THEN
              ny1=variables_at_elem(0,1,ne) !fixed
              FIX(ny1)=.TRUE. !set fixed
              prq_solution(ny1,1)=inletbc !Putting BC value into solution array
              np_in=elem_nodes(1,ne)
           ENDIF
        ENDIF
     ENDDO
  ELSE !Add terminal pressure BC for all terminal branches
     DO nonode=1,num_units
        np=elem_nodes(2,units(nonode)) !Second node in element = terminal node
        ny1=variables_at_node(np,1,1) !for fixed pressure BC
        FIX(ny1)=.TRUE. !set fixed
!! NB// Add gravitational factor in here
        IF(np_in.eq.0) then
           write(*,*),"Warning --> np_in is not set yet, setting to first node as default"
           np_in=1 !Setting to first node as default
        ENDIF
        IF(g_type.eq.'off') then
           grav=0.d0
        ELSE                    
           grav=0.d0
           DO nj=1,3
              grav=grav+density*grav_vect(nj)*9810.d0*(node_xyz(nj,np_in)-node_xyz(nj,np))
           ENDDO
           write(*,*),"gravity",ne,grav
        ENDIF
        prq_solution(ny1,1)=outletbc-grav !Putting BC value into solution array       
        !write(*,*),"BC----",ny1,outletbc,grav,prq_solution(ny1,1)
     ENDDO
  ENDIF
    call enter_exit(sub_name,2)
  end subroutine boundary_conditions
!
!###################################################################################
!
subroutine calculate_resistance(density,gamma,viscosity)
    use arrays,only: dp,num_elems,num_nodes,elems,elem_nodes,elem_cnct,node_xyz,units,&
        num_units,elem_field,elem_ordrs
    use other_consts
    use indices
    use diagnostics, only: enter_exit
    real(dp)::density,gamma,viscosity
!local variables
    integer :: ne,np1,np2,noelem
    real(dp) :: resistance,reynolds,zeta
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
       !write(*,*),"avrad",ne,elem_field(ne_radius,ne)
       elem_field(ne_radius,ne)=(elem_field(ne_radius_in,ne)+elem_field(ne_radius_out,ne))/2
       ! element Poiseuille (laminar) resistance in units of Pa.s.mm-3        
       resistance = 8.d0*viscosity*elem_field(ne_length,ne)/ &
            (PI*elem_field(ne_radius,ne)**4) !laminar resistance

       ! element turbulent resistance (flow in bifurcating tubes)
       !reynolds=DABS(elem_field(ne_flow,ne)*2.d0*density/ &
         !   (PI*elem_field(ne_radius,ne)*viscosity))
       zeta = MAX(1.d0,dsqrt(2.d0*elem_field(ne_radius,ne)* &
            reynolds/elem_field(ne_length,ne))*gamma)
       elem_field(ne_resist,ne) = resistance * zeta
       !write(*,*),"TESTING RESISTANCE",ne,elem_field(ne_resist,ne),elem_field(ne_radius,ne),elem_ordrs(2,ne)
    enddo 

    call enter_exit(sub_name,2)
  end subroutine calculate_resistance
!
!##################################################################
!
!*calc_variable_maps:* calculates the mapping between the nodes and elements and
!the problem variables that are needed for matrix setup and solution
  subroutine calc_variable_maps(mesh_from_variables,variables_at_elem,&
                variable_totals,variables_at_node,mesh_dof,num_vars)
    use arrays,only: dp,num_elems,num_nodes,elems,elem_nodes
    use diagnostics, only: enter_exit
    character(len=60) :: sub_name

     integer :: mesh_dof
     integer :: mesh_from_variables(0:2,mesh_dof,0:2)
     integer :: variables_at_elem(0:2,2,num_elems)
     integer :: variables_at_node(num_nodes,0:2,2)
     integer :: variable_totals(0:2,2)
     integer :: num_vars
!     Local Variables
    integer :: ny_start=0  
    integer :: nat=1 !Number of auxially basis functions, i.e. 1 element solution...
    integer :: i,nc,ne,nn,noelem,nonode,np,nrc,ny

    sub_name = 'calc_variable_maps'
    call enter_exit(sub_name,1)
     
     variable_totals = 0
     mesh_from_variables = 0
     variables_at_elem = 0
     variables_at_node = 0
     
!nrc = loops from 0,1,2

!  Set up mapping arrays for current region:
!  includes nested loop creating variables_at_node & variables_at_elem consecutively
     DO nrc=0,2 !row or column no
        ny=0 !variable tag
        DO nc=1,2 !no of dependent variable types
          ! IF(nrc.NE.0) ny=ny_start !--> resets to ny_start only for nrc=1,2. Maybe should also do for nrc=1??
           ny=ny_start
           DO ne=1,num_elems
              !ne=elems(noelem)
              DO nn=1,2 !Loop over number of nodes in each element: always=2 in 1D element
                 np=elem_nodes(nn,ne)
                 IF(variables_at_node(np,nrc,nc).EQ.0)THEN
                    !Checks if this node already been done
                    ny=ny+1
                     IF(ny.GT.mesh_dof) THEN
                       write(*,*),"1. Need to increase mesh_dof! ny=",ny
                       call exit(1)
                    ENDIF
                    IF(nrc.NE.0) THEN
                       IF(ny.GT.variable_totals(nrc,nc)) variable_totals(nrc,nc)=ny
                    ENDIF
                    variables_at_node(np,nrc,nc)=ny
                    IF(nrc.NE.1.OR.nc.EQ.1) THEN
                       ! Don't set mesh_from_variables for nrc=1(rows), nc<>1(LHS) cases
                       mesh_from_variables(0,ny,nrc)=1 !mesh dof is nodal based
                       mesh_from_variables(1,ny,nrc)=np
                       mesh_from_variables(2,ny,nrc)=nc
                    ENDIF !nrc.NE.1
                 ENDIF !variables_at_node.EQ.0
              ENDDO !nn (np)
              ny=ny+1
              IF(ny.GT.mesh_dof) THEN
                 write(*,*),"2. Need to increase mesh_dof! ny=",ny
                 call exit(1)
              ENDIF
              IF(nrc.NE.0) THEN
                 IF(ny.GT.variable_totals(nrc,nc)) variable_totals(nrc,nc)=ny
              ENDIF
              variables_at_elem(nrc,nc,ne)=ny
              IF(nrc.NE.1.OR.nc.EQ.1) THEN
                 !                 Don't set mesh_from_variables for nrc=1(rows), nc<>1(LHS) cases
                 mesh_from_variables(0,ny,nrc)=2 !mesh dof is element based
                 mesh_from_variables(1,ny,nrc)=nc
                 mesh_from_variables(2,ny,nrc)=ne
              ENDIF
           ENDDO !noelem (ne)          
        ENDDO !nc
     ENDDO !nrc
     num_vars=ny
     !write(*,*),"Max ny number=",ny


    call enter_exit(sub_name,2)
  end subroutine calc_variable_maps
!
!##################################################################
!
!*tree_resistance:* Calculates the total resistance of a tree

  subroutine tree_resistance(resistance)
    use indices
    use arrays,only: dp,num_elems,elem_cnct,elem_field,elems
    use diagnostics, only: enter_exit
    character(len=60) :: sub_name
!local variables
    real(dp), intent(out) :: resistance
    real(dp) :: invres,elem_res(num_elems)
    integer :: num,num2,ne,ne2

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
    variables_at_node,variables_at_elem,FIX)

    use indices
    use arrays,only: dp,num_elems,elem_ordrs,elems,elem_nodes,num_nodes
    use diagnostics, only: enter_exit
    integer, intent(in) :: mesh_dof
    integer,intent(in) :: variables_at_elem(0:2,2,num_elems)
    integer,intent(in) :: variables_at_node(num_nodes,0:2,2)
    real(dp), intent(in) :: pressure_in, pressure_out,cardiac_output
    real(dp) :: prq_solution(mesh_dof,2)
    logical, intent(in) :: FIX(mesh_dof)
!local variables
    integer :: nn,ne,np,noelem,n_variable
    character(len=60) :: sub_name
   sub_name = 'intialise_solution'
    call enter_exit(sub_name,1)
    do ne=1,num_elems
       !ne=elems(noelem)
       do nn=1,2 !Loop over number of nodes in each element: always=2 in 1D element
         np=elem_nodes(nn,ne) !Node number
         n_variable=variables_at_node(np,0,1) !--> This will be a pressure because it's the variable at a node
         if(.NOT.FIX(n_variable))then
           prq_solution(n_variable,1)=(pressure_in+pressure_out)/2.0
         endif
         n_variable=variables_at_elem(0,1,ne) !--> This will be a flow because it's the variable at the element
         if(.NOT.FIX(n_variable))then
           prq_solution(n_variable,1)=cardiac_output/(2**(elem_ordrs(1,ne)-1)) !Here you can use the generation number to split flow
         endif
       enddo
    enddo
    call enter_exit(sub_name,2)
  end subroutine initialise_solution
!
!##################################################################
!
!*calc_sparse_1d_tree:* Calculates sparsity structure for 1d tree problems

subroutine calc_sparse_1dtree(mesh_dof,variables_at_elem,variables_at_node,FIX,NonZeros,&
        MatrixSize,SparseCol,SparseRow,SparseVal,RHS,prq_solution,update_resistance_entries)
    use indices
    use arrays,only: dp,num_elems,elem_ordrs,elems,elem_nodes,num_nodes,elems_at_node,elem_field
    use diagnostics, only: enter_exit
    integer, intent(in) :: mesh_dof
    integer,intent(in) :: variables_at_elem(0:2,2,num_elems)
    integer,intent(in) :: variables_at_node(num_nodes,0:2,2)
    logical, intent(in) :: FIX(mesh_dof)
    integer :: NonZeros,MatrixSize
    integer :: SparseCol(NonZeros)
    integer :: SparseRow(MatrixSize+1)
    integer :: update_resistance_entries(num_elems)
    real(dp) :: SparseVal(NonZeros)
    real(dp) :: RHS(MatrixSize)
    real(dp) :: prq_solution(mesh_dof,2)
!local variables
    integer :: noelem,ne,np1,variable,variable1,NumDiag,nhs,np2,variable2,ost2,nzz,nzz_val,&
      nzz_row,ost1,ny2c,nn,np,ny,ny2,ne2,noelem2,nonode
    integer :: NPLIST(0:num_nodes)
    logical :: FLOW_BALANCED
    real(dp) :: grav,flow_term
    character(len=60) :: sub_name
   sub_name = 'calc_sparse_1dtree'
    call enter_exit(sub_name,1)

  SparseCol=0
  SparseRow=1
  SparseVal=0.0_dp
  RHS=0.0_dp
  NPLIST=0
  NumDiag=0!number of diagonal entries
  !NonZeros=0
  nzz=1 !position in SparseCol
  nzz_val=1 !position in SparseVal
  nzz_row=1 !position in SparseRow
  ost1=0!offset
  ost2=0!offset

  do ne=1,num_elems
     !ne=elems(noelem)!elem no
     np1=elem_nodes(2,ne) !second node in that element
     variable1=variables_at_node(np1,1,1) !which entry
     if((.NOT.FIX(variable1)).OR.(elems_at_node(np1,0).LE.1))then !if its not fixed, or it only has one element connected
       NumDiag=NumDiag+1!one more diagonal entry
        do nhs=1,2 !for each row entry
              np2=elem_nodes(nhs,ne) 
              variable2=variables_at_node(np2,1,1) 
              if(.NOT.FIX(variable2))then !if variable for node2 is not fixed
                 ost2=0!zero offset
                 do variable=variable2-1,1,-1 !loop over previous columns
                    if(FIX(variable))then !new line
                       ost2=ost2+1 !count # of fixed BC
                    endif!if variable us fixed
                 enddo!end looping over pervious cols
                 SparseCol(nzz)=(variable2-ost2) !store the col # offset by ost2 -correct
                 nzz=nzz+1
                 !if(GRAVITY.eq.'nogravity') then
                   grav=0.d0
                 !else
                 !   grav=0.d0
                 !   DO nj=1,3                                  
                 !      grav=grav+density*grav_vect(nj)*9810.d0*(node_xyz(nj,elem_nodes(1,ne))-node_xyz(nj,elem_nodes(2,ne)))
                 !   ENDDO
                 !endif
		 IF (nhs.EQ.1) THEN !row entry one
                    SparseVal(nzz_val)=1.0_dp
                    nzz_val=nzz_val+1
                    SparseVal(nzz_val)=-prq_solution(variable2,1)+grav !minus inlet pressure
                 ELSEIF (nhs.EQ.2) THEN !row entry two
                    SparseVal(nzz_val)=-1.0_dp
                    nzz_val=nzz_val+1
                 ENDIF
             else!notfix variable2
!                 if(GRAVITY.eq.'nogravity') then
                    grav=0.d0
 !                else
  !                  grav=0.d0
  !                  DO nj=1,3                                  
 !                      grav=grav+density*grav_vect(nj)*9810.d0*(node_xyz(nj,elem_nodes(1,ne))-node_xyz(nj,elem_nodes(2,ne)))
  !                  ENDDO
   !              endif
                 IF(nhs.EQ.1) THEN
                    RHS(nzz_row)=-prq_solution(variable2,1)+grav
                 ELSEIF(nhs.EQ.2) THEN
                    RHS(nzz_row)=prq_solution(variable2,1)+grav
                 ENDIF
             endif
        enddo!for each row entry,nhs
           !now looking at flow variables
           ny2c=variables_at_elem(0,1,ne)
           IF(.NOT.FIX(ny2c)) THEN !if not fixed
              ost2=0
              DO ny=ny2c-1,1,-1 !loop over previous columns
                 IF(FIX(ny))THEN !new line
                    ost2=ost2+1 !count # of fixed BC
                    ENDIF
              ENDDO !ny
              SparseCol(nzz)=(ny2c-ost2) !store column #
              nzz=nzz+1
              SparseVal(nzz_val)=-elem_field(ne_resist,ne) !resistance components of dP=QR
              update_resistance_entries(ne)=nzz_val! needed?
              nzz_val=nzz_val+1
           ELSE
              RHS(nzz_row)=prq_solution(ny2c,1)
           ENDIF !FIX

        nzz_row=nzz_row+1
        SparseRow(nzz_row)=nzz !stores row #
        !Now flow balance equation
        DO nn=1,2 !balances at each node of ne
           FLOW_BALANCED=.FALSE. !initialise
           np=elem_nodes(nn,ne)
           DO nonode=1,NPLIST(0) !junctions balanced at already
              IF(np.EQ.NPLIST(nonode)) THEN
                 FLOW_BALANCED=.TRUE. !already balanced at np
              ENDIF
           ENDDO !nonode
           IF((elems_at_node(np,0).GT.1).AND.(.NOT.FLOW_BALANCED))THEN !if there is more than one node at a junction an we arent flow balanced
              !at an unbalanced junction
              DO noelem2=1,elems_at_node(np,0) !for elems with
                 ne2=elems_at_node(np,noelem2) !subtended branches only
!                 IF(nzz.LT.NISC_GKM) THEN
                    ny2=variables_at_elem(1,1,ne2)
                    IF(.NOT.FIX(ny2))THEN
                       ost2=0
                       DO ny=ny2-1,1,-1 !loop over previous columns
                          IF(FIX(ny))THEN !new line
                             ost2=ost2+1 !count # of fixed BC
                          ENDIF
                       ENDDO !ny
                       SparseCol(nzz)=(ny2-ost2)
                       nzz=nzz+1

                       !! NEED TO TEST THIS PART!!
                       !IF(NORD(5,ne2).EQ.0.OR.ne.EQ.ne2)  THEN !capillary
                          flow_term=1.0_dp !Q1-Q2-Q3=0
                       !ELSEIF(NORD(5,ne2).EQ.1) THEN
                       !   !symmetric arteriole - flow decreases by 2: Q1-2*Q2=0
                       !   flow_term=2.d0
                       !ELSEIF(NORD(5,ne2).EQ.-1) THEN !venule
                       !   !symmetric venule - flow increases by 2: Q1-0.5*Q2=0 
                       !   flow_term=0.5d0
                       !ENDIF
                       IF(np.EQ.elem_nodes(2,ne2)) THEN !end node
                          SparseVal(nzz_val)=flow_term
                          nzz_val=nzz_val+1
                       ELSEIF(np.EQ.elem_nodes(1,ne2)) THEN !start node
                         SparseVal(nzz_val)=-flow_term
                          nzz_val=nzz_val+1
                       ENDIF
                    ENDIF
              ENDDO !noelem2
              NPLIST(0)=NPLIST(0)+1 !stores nodes where flow
              NPLIST(NPLIST(0))=np !balance been done
           ENDIF !elem_from_node(np,0).GT.1
           IF((.NOT.FLOW_BALANCED).AND.(elems_at_node(np,0).GT.1))THEN
              nzz_row=nzz_row+1
              SparseRow(nzz_row)=nzz !store the row #
              NumDiag=NumDiag+1 !# entries on diagonal
           ENDIF !.NOT.FLOW_BALANCED
        ENDDO !nn

     else
       ost1=ost1+1
     endif

  enddo

    call enter_exit(sub_name,2)
  end subroutine calc_sparse_1dtree


!
!##################################################################
!
!*calc_sparse_size:* Calculates sparsity sizes

subroutine calc_sparse_size(mesh_dof,variables_at_elem,variables_at_node,FIX,NonZeros,MatrixSize)
    use indices
    use arrays,only: dp,num_elems,elem_ordrs,elems,elem_nodes,num_nodes,elems_at_node,elem_field
    use diagnostics, only: enter_exit
    integer, intent(in) :: mesh_dof
    integer,intent(in) :: variables_at_elem(0:2,2,num_elems)
    integer,intent(in) :: variables_at_node(num_nodes,0:2,2)
    logical, intent(in) :: FIX(mesh_dof)
    integer :: NonZeros,MatrixSize
!local variables
    integer :: noelem,ne,np1,variable,variable1,NumDiag,nhs,np2,variable2,ost2,nzz,nzz_val,&
      nzz_row,ost1,ny2c,nn,np,ny,ny2,ne2,noelem2,nonode
    integer :: NPLIST(0:num_nodes)
    logical :: FLOW_BALANCED
    real(dp) :: grav,flow_term
    character(len=60) :: sub_name
   sub_name = 'calc_sparse_size'
    call enter_exit(sub_name,1)

 
  NPLIST=0
  NumDiag=0!number of diagonal entries
  !NonZeros=0
  nzz=1 !position in SparseCol
  nzz_val=1 !position in SparseVal
  nzz_row=1 !position in SparseRow
  ost1=0!offset
  ost2=0!offset

  do ne=1,num_elems
     !ne=elems(noelem)!elem no
     np1=elem_nodes(2,ne) !second node in that element
     variable1=variables_at_node(np1,1,1) !which entry
     if((.NOT.FIX(variable1)).OR.(elems_at_node(np1,0).LE.1))then !if its not fixed, or it only has one element connected
       NumDiag=NumDiag+1!one more diagonal entry
        do nhs=1,2 !for each row entry
              np2=elem_nodes(nhs,ne) 
              variable2=variables_at_node(np2,1,1) 
              if(.NOT.FIX(variable2))then !if variable for node2 is not fixed
                 nzz=nzz+1
		 IF (nhs.EQ.1) THEN !row entry one
                    nzz_val=nzz_val+1
                 ELSEIF (nhs.EQ.2) THEN !row entry two
                    nzz_val=nzz_val+1
                 ENDIF
             else!notfix variable2
             endif
        enddo!for each row entry,nhs
           !now looking at flow variables
           ny2c=variables_at_elem(0,1,ne)
           IF(.NOT.FIX(ny2c)) THEN !if not fixed
              nzz=nzz+1
              nzz_val=nzz_val+1
           ELSE
           ENDIF !FIX
        nzz_row=nzz_row+1
        !Now flow balance equation
        DO nn=1,2 !balances at each node of ne
           FLOW_BALANCED=.FALSE. !initialise
           np=elem_nodes(nn,ne)
           DO nonode=1,NPLIST(0) !junctions balanced at already
              IF(np.EQ.NPLIST(nonode)) THEN
                 FLOW_BALANCED=.TRUE. !already balanced at np
              ENDIF
           ENDDO !nonode
           IF((elems_at_node(np,0).GT.1).AND.(.NOT.FLOW_BALANCED))THEN !if there is more than one node at a junction an we arent flow balanced
              !at an unbalanced junction
              DO noelem2=1,elems_at_node(np,0) !for elems with
                 ne2=elems_at_node(np,noelem2) !subtended branches only
!                 IF(nzz.LT.NISC_GKM) THEN
                    ny2=variables_at_elem(1,1,ne2)
                    IF(.NOT.FIX(ny2))THEN
                       nzz=nzz+1
                       IF(np.EQ.elem_nodes(2,ne2)) THEN !end node
                          nzz_val=nzz_val+1
                       ELSEIF(np.EQ.elem_nodes(1,ne2)) THEN !start node
                          nzz_val=nzz_val+1
                       ENDIF
                    ENDIF
              ENDDO !noelem2
              NPLIST(0)=NPLIST(0)+1 !stores nodes where flow
              NPLIST(NPLIST(0))=np !balance been done
           ENDIF !elem_from_node(np,0).GT.1
           IF((.NOT.FLOW_BALANCED).AND.(elems_at_node(np,0).GT.1))THEN
              nzz_row=nzz_row+1
              NumDiag=NumDiag+1 !# entries on diagonal
           ENDIF !.NOT.FLOW_BALANCED
        ENDDO !nn
     else
     endif

  enddo 
    NonZeros=nzz-1
    MatrixSize=nzz_row-1

    call enter_exit(sub_name,2)
  end subroutine calc_sparse_size

!##################################################################
!
!*calc_press_area:* Calculates new radii based on pressure area relnships

subroutine calc_press_area(KOUNT,variables_at_node,prq_solution,mesh_dof,ptrans,beta,Go_artery,Ptm_max)
use indices
use arrays,only: dp,num_nodes,num_elems,elem_field,elem_nodes,elem_cnct!,num_elems,elem_ordrs,elems,elem_nodes,num_nodes,elems_at_node,elem_field
use diagnostics, only: enter_exit
integer,intent(in) :: KOUNT,mesh_dof
integer,intent(in) :: variables_at_node(num_nodes,0:2,2)
real(dp) ::  prq_solution(mesh_dof,2)
real(dp) :: ptrans,beta,Go_artery,Ptm_max

!local variables
integer :: np,de,ne,ny,numelem,nn
double precision :: HEIGHT(3),G_PLEURAL,Ptm,R0,Pblood,Ppl

character(len=60) :: sub_name
sub_name = 'calc_press_area'
call enter_exit(sub_name,1)

IF(KOUNT.EQ.1)THEN !store initial, unstressed radius values
  DO ne=1,num_elems
    elem_field(ne_radius_in0,ne)=elem_field(ne_radius_in,ne)
    elem_field(ne_radius_out0,ne)=elem_field(ne_radius_in,ne)
ENDDO !elems
ENDIF

DO ne=1,num_elems
  do nn=1,2
  np=elem_nodes(1,ne)
  ny=variables_at_node(np,0,1)
  Ppl=(ptrans*98.07d0/1000.0d0)!-G_PLEURAL !cmH2O->kPa
  Pblood=prq_solution(ny,1)/1000.0d0 ! Pa->kPa
  Ptm=Pblood+Ppl     ! kPa
  R0=elem_field(ne_radius_in0,ne)
!...ARC: giving a maximum distension
  IF(Ptm.LT.Ptm_max.and.Go_artery.gt.0.d0)THEN
     if(nn.eq.1) elem_field(ne_radius_in,ne)=R0*((Ptm/Go_artery)+1.d0)**(1.d0/beta)
     if(nn.eq.2) elem_field(ne_radius_out,ne)=R0*((Ptm/Go_artery)+1.d0)**(1.d0/beta)
  ELSEIF(Ptm.lt.0.or.Go_artery.eq.0.d0)THEN
    IF(Ptm.lt.0)write(*,*) 'Transmural pressure < zero',ne,Ptm,Pblood,Ppl
     if(nn.eq.1) elem_field(ne_radius_in,ne)=R0
     if(nn.eq.2) elem_field(ne_radius_out,ne)=R0
  ELSE!ptm>ptmmax
    if(nn.eq.1) elem_field(ne_radius_in,ne)=R0*((Ptm_max/Go_artery)+1.d0)**(1.d0/beta)
    if(nn.eq.2) elem_field(ne_radius_in,ne)=R0*((Ptm_max/Go_artery)+1.d0)**(1.d0/beta)
  ENDIF
  enddo!nn
enddo!ne
!ny=variables_at_node(np,0,1)

!...  KSB: Adding gravity-dependent pleural pressure
!G_PLEURAL=0.d0    !gravitational force
!DO nj=1,3
!HEIGHT(nj)=node_xyz(nj,np)-node_xyz(nj,np_in)
!G_PLEURAL=G_PLEURAL+PLEURAL_DENSITY*grav_vect(nj)*9810.d0*HEIGHT(nj) !kg
!ENDDO
!Ppl=(ptrans*98.07d0/1000.0d0)!-G_PLEURAL !cmH2O->kPa
!Pblood=solution_prq(ny,1)/1000.0d0 ! Pa->kPa
!Ptm=Pblood+Ppl     ! kPa


!R0=elem_field(nj_radius0,np)
!...ARC: giving a maximum distension
!IF(Ptm.LT.Ptm_max.and.Go_vein.gt.0.d0)THEN
!node_field(nj_radius,np)=R0*((Ptm/Go)+1.d0)**(1.d0/beta)
!ELSEIF(Ptm.lt.0.or.Go.eq.0.d0)THEN
!IF(Ptm.lt.0)write(*,*) 'Transmural pressure < zero',ne,Ptm,Pblood,Ppl
!node_field(nj_radius,np)=R0
!ELSE
!node_field(nj_radius,np)=R0*((Ptm_max/Go)+1.d0)**(1.d0/beta)
!write(*,*) 'Transmural pressure < zero',ne,Ptm,Pblood,Ppl,Go,Go_vein
!ENDIF
!Hypoxia removed
!ENDDO     !nonode



call enter_exit(sub_name,2)
end subroutine calc_press_area

end module pressure_resistance_flow
