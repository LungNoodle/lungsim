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
    use arrays,only: dp,num_elems,num_nodes
    use diagnostics, only: enter_exit
    !local variables
    integer :: mesh_dof,variable_types
    integer, allocatable :: node_to_variables(:,:,:)
    integer, allocatable :: variables_to_node(:,:,:)
    integer, allocatable :: variables_to_elem(:,:,:)
    integer, dimension(0:2,2) :: variable_totals
    integer :: num_vars
    integer :: AllocateStatus

    real(dp), allocatable :: prq_solution(:,:)
    real(dp) :: viscosity,density,inletbc,outletbc,grav_vect(3),gamma
    logical, allocatable :: FIX(:)
    logical :: ADD=.FALSE.
    character(len=60) :: sub_name,mesh_type,vessel_type,mechanics_type,bc_type,g_type
    integer :: g_dirn

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
vessel_type='rigid'
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

!!---------DESCRIPTION OF IMPORTANT PARAMETERS-----------
!viscosity: fluid viscosity
!density:fluid density
!

!!!set the default values for the parameters that control the prq simulation, these should be controlled by user input
call read_params_evaluate_prq(viscosity,density,gamma)


!! Allocate memory to solution arrays
mesh_dof=num_elems+num_nodes
variable_types=2 !pressure/flow
    allocate (node_to_variables(0:2,mesh_dof,0:2), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for node_to_variables array ***"
    allocate (variables_to_elem(0:2,variable_types,num_elems), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for variables_to_elem array ***"
    allocate (variables_to_node(num_nodes,0:2,variable_types), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for variables_to_node array ***"
    allocate (prq_solution(mesh_dof,2), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for prq_solution array ***"
    allocate (FIX(mesh_dof), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory for FIX array ***"


!! Setting up mappings between nodes, elements and solution variables
    call calc_variable_maps(node_to_variables,variables_to_elem,&
                variable_totals,variables_to_node,mesh_dof,num_vars)

!! Define boundary conditions
    !first call to define inlet boundary conditions
    call boundary_conditions(ADD,FIX,bc_type,g_type,grav_vect,density,inletbc,outletbc,&
      variables_to_node,variables_to_elem,prq_solution,mesh_dof)
    !second call if simple tree need to define pressure bcs at all terminal branches
    if(mesh_type.eq.'simple_tree')then
        ADD=.TRUE.
        call boundary_conditions(ADD,FIX,bc_type,g_type,grav_vect,density,inletbc,outletbc,&
            variables_to_node,variables_to_elem,prq_solution,mesh_dof)   
    endif
    print *,prq_solution(:,1)  
 
!! Calculate resistance of each element
   call calculate_resistance(density,gamma,viscosity)
        
!! Calculate sparsity structure for solution matrices
  ! call calc_sparse_1dtree(SparseCol,SparseRow,FIX,SparseVal,RHS,NDIAG,NPLIST,&
   !     variables_to_node,variables_to_elem,nz_esed,prq_solution)


    deallocate (node_to_variables, STAT = AllocateStatus)
    deallocate (variables_to_elem, STAT = AllocateStatus)
    deallocate (variables_to_node, STAT = AllocateStatus)
    deallocate (prq_solution, STAT = AllocateStatus)
    deallocate (FIX, STAT = AllocateStatus)
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
       variables_to_node,variables_to_elem,prq_solution,mesh_dof)
 use arrays,only: dp,num_elems,num_nodes,elems,elem_nodes,elem_cnct,node_xyz,units,&
        num_units
    use diagnostics, only: enter_exit

    integer :: mesh_dof
    integer :: variables_to_elem(0:2,2,num_elems)
    integer :: variables_to_node(num_nodes,0:2,2)
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
              ny1=variables_to_node(np,1,1) !for fixed pressure BC
              FIX(ny1)=.TRUE. !set fixed
              prq_solution(ny1,1)=inletbc !Putting BC value into solution array
              np_in=np !inlet node set here, gravity reference to this point
           ELSE IF(BC_TYPE == 'flow')THEN
              ny1=variables_to_elem(0,1,ne) !fixed
              FIX(ny1)=.TRUE. !set fixed
              prq_solution(ny1,1)=inletbc !Putting BC value into solution array
              np_in=elem_nodes(1,ne)
           ENDIF
        ENDIF
     ENDDO
  ELSE !Add terminal pressure BC for all terminal branches
    write(*,*) num_units,units
     DO nonode=1,num_units
        np=elem_nodes(2,units(nonode)) !Second node in element = terminal node
        ny1=variables_to_node(np,1,1) !for fixed pressure BC
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
        write(*,*),"BC----",ny1,outletbc,grav,prq_solution(ny1,1)
     ENDDO
  ENDIF
    call enter_exit(sub_name,2)
  end subroutine boundary_conditions
!
!###################################################################################
!
subroutine calculate_resistance(density,gamma,viscosity)
    use arrays,only: dp,num_elems,num_nodes,elems,elem_nodes,elem_cnct,node_xyz,units,&
        num_units,elem_field
    use other_consts
    use indices
    use diagnostics, only: enter_exit
    real(dp)::density,gamma,viscosity
!local variables
    integer :: ne,np1,np2
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
      
       reynolds=DABS(elem_field(ne_flow,ne)*2.d0*density/ &
            (PI*elem_field(ne_radius,ne)*viscosity))
       zeta = MAX(1.d0,dsqrt(2.d0*elem_field(ne_radius,ne)* &
            reynolds/elem_field(ne_length,ne))*gamma)
       elem_field(ne_resist,ne) = resistance * zeta
       write(*,*),"TESTING RESISTANCE",ne,elem_field(ne_resist,ne),elem_field(ne_radius,ne),zeta
    enddo 

    call enter_exit(sub_name,2)
  end subroutine calculate_resistance
!
!##################################################################
!
!*calc_variable_maps:* calculates the mapping between the nodes and elements and
!the problem variables that are needed for matrix setup and solution
  subroutine calc_variable_maps(node_to_variables,variables_to_elem,&
                variable_totals,variables_to_node,mesh_dof,num_vars)
    use arrays,only: dp,num_elems,num_nodes,elems,elem_nodes
    use diagnostics, only: enter_exit
    character(len=60) :: sub_name

     integer :: mesh_dof
     integer :: node_to_variables(0:2,mesh_dof,0:2)
     integer :: variables_to_elem(0:2,2,num_elems)
     integer :: variables_to_node(num_nodes,0:2,2)
     integer :: variable_totals(0:2,2)
     integer :: num_vars
!     Local Variables
    integer :: ny_start=0  
    integer :: nat=1 !Number of auxially basis functions, i.e. 1 element solution...
    integer :: i,nc,ne,nn,noelem,nonode,np,nrc,ny

    sub_name = 'calc_variable_maps'
    call enter_exit(sub_name,1)
     
     variable_totals = 0
     node_to_variables = 0
     variables_to_elem = 0
     variables_to_node = 0
     
!nrc = loops from 0,1,2

!  Set up mapping arrays for current region:
!  includes nested loop creating variables_to_node & variables_to_elem consecutively
     DO nrc=0,2 !row or column no
        ny=0 !variable tag
        DO nc=1,2 !no of dependent variable types
          ! IF(nrc.NE.0) ny=ny_start !--> resets to ny_start only for nrc=1,2. Maybe should also do for nrc=1??
           ny=ny_start
           DO ne=1,num_elems
              !ne=elems(noelem)
              DO nn=1,2 !Loop over number of nodes in each element: always=2 in 1D element
                 np=elem_nodes(nn,ne)
                 IF(variables_to_node(np,nrc,nc).EQ.0)THEN
                    !Checks if this node already been done
                    ny=ny+1
                     IF(ny.GT.mesh_dof) THEN
                       write(*,*),"1. Need to increase mesh_dof! ny=",ny
                       call exit(1)
                    ENDIF
                    IF(nrc.NE.0) THEN
                       IF(ny.GT.variable_totals(nrc,nc)) variable_totals(nrc,nc)=ny
                    ENDIF
                    variables_to_node(np,nrc,nc)=ny
                    IF(nrc.NE.1.OR.nc.EQ.1) THEN
                       ! Don't set node_to_variables for nrc=1(rows), nc<>1(LHS) cases
                       node_to_variables(0,ny,nrc)=1 !mesh dof is nodal based
                       node_to_variables(1,ny,nrc)=np
                       node_to_variables(2,ny,nrc)=nc
                    ENDIF !nrc.NE.1
                 ENDIF !variables_to_node.EQ.0
              ENDDO !nn (np)
              ny=ny+1
              IF(ny.GT.mesh_dof) THEN
                 write(*,*),"2. Need to increase mesh_dof! ny=",ny
                 call exit(1)
              ENDIF
              IF(nrc.NE.0) THEN
                 IF(ny.GT.variable_totals(nrc,nc)) variable_totals(nrc,nc)=ny
              ENDIF
              variables_to_elem(nrc,nc,ne)=ny
              IF(nrc.NE.1.OR.nc.EQ.1) THEN
                 !                 Don't set node_to_variables for nrc=1(rows), nc<>1(LHS) cases
                 node_to_variables(0,ny,nrc)=2 !mesh dof is element based
                 node_to_variables(1,ny,nrc)=nc
                 node_to_variables(2,ny,nrc)=ne
              ENDIF
           ENDDO !noelem (ne)          
        ENDDO !nc
     ENDDO !nrc
     num_vars=ny
     !write(*,*),"Max ny number=",ny


    call enter_exit(sub_name,2)
  end subroutine calc_variable_maps
!
end module pressure_resistance_flow
