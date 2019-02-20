!> \file
!> \author Merryn Tawhai
!> \brief This module contains code specific to running gas mixing problems
!>
!> \section LICENSE
!>
!>
!> Contributor(s):
!>
!>\Description
!> This module contains code specific to running gas mixing problems

module gasmix
  use arrays,only: dp
  implicit none
  private airway_mesh_deform

  integer,public :: inlet_node = 1
  integer,private :: NonZeros_unreduced
  real(dp),private :: ideal_mass,initial_mass
  real(dp),private :: total_volume_change

  integer,allocatable :: sparsity_col(:),reduced_col(:)
  integer,allocatable :: sparsity_row(:),reduced_row(:)
  real(dp),allocatable :: global_K(:),global_M(:),global_AA(:),global_BB(:)
  real(dp),allocatable :: global_R(:)

contains

!!!#########################################################################

  subroutine assemble_gasmix(diffusion_coeff,nonzeros_unreduced)

    use arrays,only: dp,elem_nodes,num_elems
    use geometry, only: volume_of_mesh
    use diagnostics, only: enter_exit
    implicit none

    integer,intent(in) :: nonzeros_unreduced
    real(dp),intent(in) :: diffusion_coeff

    integer :: i,j,ncol,ne,nentry,nrow
    real(dp) :: elem_K(2,2),elem_M(2,2),elem_R(2)
    logical :: found
    character(len=60) :: sub_name

!!!................................................................

    sub_name = 'assemble_gasmix'
    call enter_exit(sub_name,1)

    global_K(1:nonzeros_unreduced) = 0.0_dp
    global_M(1:nonzeros_unreduced) = 0.0_dp

    do ne=1,num_elems
       call element_gasmix(ne,elem_K,elem_M,elem_R,diffusion_coeff)
       do i=1,2
          nrow = elem_nodes(i,ne)
          do j=1,2
             ncol = elem_nodes(j,ne)
             found=.false.
             nentry = sparsity_row(nrow) ! start check at start of row
             do while (.not.found)
                if(ncol.eq.sparsity_col(nentry))then
                   found = .true.
                else
                   nentry = nentry+1
                endif
             enddo
             global_K(nentry) = global_K(nentry) + elem_K(i,j)
             global_M(nentry) = global_M(nentry) + elem_M(i,j)
          enddo !j
       enddo !i
    enddo !noelem

    call enter_exit(sub_name,2)

  end subroutine assemble_gasmix

!!!################################################################################

  subroutine calc_mass(nj,nu_field,gas_mass)
    use arrays,only: dp,elem_cnct,elem_nodes,elem_symmetry,&
         node_field,num_elems,num_nodes,num_units,elem_field,units,unit_field
    use indices,only: ne_vol,nu_vol
    use diagnostics, only: enter_exit
    implicit none

    integer,intent(in) :: nj,nu_field
    real(dp) :: gas_mass
    !     Local Variables
    integer :: ne,ne0,np1,np2,nunit
    real(dp) :: average_conc
    real(dp),allocatable :: tree_mass(:)
    character(len=60):: sub_name

    !...........................................................................

    sub_name = 'calc_mass'
    call enter_exit(sub_name,1)

    if(.not.allocated(tree_mass)) allocate(tree_mass(num_nodes))

    ! initialise to the mass in each element
    do ne=1,num_elems
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       average_conc = (node_field(nj,np1)+node_field(nj,np2))/2.0_dp
       tree_mass(ne) = average_conc*elem_field(ne_vol,ne)
    enddo

    ! add the mass in each elastic unit to terminal elements
    do nunit=1,num_units
       ne=units(nunit)
       tree_mass(ne) = tree_mass(ne) + &
            unit_field(nu_vol,nunit)*unit_field(nu_field,nunit)
    enddo

    ! sum mass recursively up the tree
    do ne=num_elems,2,-1 ! not for the stem branch; parent = 0
       ne0=elem_cnct(-1,1,ne)
       tree_mass(ne0) = tree_mass(ne0) + DBLE(elem_symmetry(ne))*tree_mass(ne)
    enddo !noelem

    gas_mass = tree_mass(1)

    deallocate(tree_mass)

    call enter_exit(sub_name,2)

  end subroutine calc_mass

!!!###################################################################

  subroutine element_gasmix(ne,elem_K,elem_M,elem_R,diffusion_coeff)
    use arrays,only: dp,elem_field,elem_symmetry
    use indices,only: ne_a_A,ne_length,ne_radius
    use other_consts
    implicit none

    integer,intent(in) :: ne
    real(dp) :: elem_K(2,2),elem_M(2,2),elem_R(2)
    real(dp),intent(in) :: diffusion_coeff

    !Local variables
    real(dp) :: a_A_ratio,inner_area,length,outer_area,radius

    radius = elem_field(ne_radius,ne)
    length = elem_field(ne_length,ne)
    a_A_ratio = elem_field(ne_a_A,ne)
    outer_area=PI*radius**2
    inner_area=outer_area*a_A_ratio

    elem_M(1,1) = outer_area*length/3.0_dp*DBLE(elem_symmetry(ne))
    elem_M(1,2) = outer_area*length/3.0_dp/2.0_dp*DBLE(elem_symmetry(ne))
    elem_M(2,1) = outer_area*length/3.0_dp/2.0_dp
    elem_M(2,2) = outer_area*length/3.0_dp

    elem_K(1,1) = (inner_area*diffusion_coeff/length)*DBLE(elem_symmetry(ne))
    elem_K(1,2) = (-inner_area*diffusion_coeff/length)*DBLE(elem_symmetry(ne))
    elem_K(2,1) = -inner_area*diffusion_coeff/length
    elem_K(2,2) = inner_area*diffusion_coeff/length

    elem_R = 0.0_dp

  end subroutine element_gasmix

!!!########################################################################

  subroutine initial_gasmix(initial_concentration,inlet_concentration)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_INITIAL_GASMIX" :: INITIAL_GASMIX

    use arrays,only: dp,node_field,num_nodes
    use indices,only: nj_conc1,nu_conc1
    use diagnostics, only: enter_exit
    implicit none

    real(dp),intent(in) :: initial_concentration,inlet_concentration

    character(len=60):: sub_name

    ! #########################################################################

    sub_name = 'initial_gasmix'
    call enter_exit(sub_name,1)

    node_field(nj_conc1,1:num_nodes) = initial_concentration

    ! initialise the 'ideal mass' to the mass of gas initially in model
    call calc_mass(nj_conc1,nu_conc1,ideal_mass)

    node_field(nj_conc1,1) = inlet_concentration
    total_volume_change = 0.0_dp ! records the volume change from FRC

! allocate the arrays for solving
    if(.not.allocated(sparsity_col)) allocate(sparsity_col(1+3*(num_nodes-1)))
    if(.not.allocated(reduced_col))  allocate(reduced_col(1+3*(num_nodes-1)))
    if(.not.allocated(sparsity_row)) allocate(sparsity_row(num_nodes+1))
    if(.not.allocated(reduced_row))  allocate(reduced_row(num_nodes+1))
    if(.not.allocated(global_K))     allocate(global_K(1+3*(num_nodes-1)))
    if(.not.allocated(global_M))     allocate(global_M(1+3*(num_nodes-1)))
    if(.not.allocated(global_AA))    allocate(global_AA(1+3*(num_nodes-1)))
    if(.not.allocated(global_BB))    allocate(global_BB(num_nodes))
    if(.not.allocated(global_R))     allocate(global_R(num_nodes))

    ! calculate the sparsity pattern for unreduced system
    call sparse_gasmix

    call enter_exit(sub_name,2)

  end subroutine initial_gasmix

!!!##########################################################################

  subroutine reduce_gasmix(MatrixSize,NonZeros,noffset_entry,noffset_row,&
       inspiration)
    use arrays,only: num_nodes
    implicit none

    integer :: MatrixSize,NonZeros,&
         noffset_entry,noffset_row
    logical,intent(in) :: inspiration

    integer :: i

    if(inspiration)then !remove first row and column (note: also for breath-hold)

       do i=1,num_nodes ! one more than # of rows
          reduced_row(i) = sparsity_row(i+1)-3
       enddo
       NonZeros = NonZeros_unreduced - 3
       do i=1,NonZeros
          reduced_col(i) = sparsity_col(i+3)-1
       enddo
       reduced_row(1)=1
       MatrixSize = num_nodes - 1
       noffset_entry = 3
       noffset_row = 1

    else !expiration

       do i=1,num_nodes+1
          reduced_row(i) = sparsity_row(i)
       enddo
       NonZeros = NonZeros_unreduced
       do i=1,NonZeros
          reduced_col(i) = sparsity_col(i)
       enddo
       reduced_row(1)=1
       MatrixSize = num_nodes
       noffset_entry = 0
       noffset_row = 0

    endif

  end subroutine reduce_gasmix

!!!#######################################################################

  subroutine solve_gasmix(fileid,inr_itr_max,out_itr_max,diffusion_coeff,&
       dt,initial_volume,inlet_concentration,inlet_flow,solve_tolerance,time_end,&
       time_start,inspiration)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SOLVE_GASMIX" :: SOLVE_GASMIX

!!! Assemble matrices for 1D inert gas mixing equation, and solve. The sparsity
!!! structure is first calculated. At each time step the mesh size is changed
!!! incrementally. The equations are solved using a Lagrange-Galerkin method which separates
!!! solution of the advective and diffusive components of the equation. The advective
!!! component is solved using method of characteristics (in this case a simple tracking
!!! of material point movement); and this is followed by solution of the diffusion equation
!!! with the advective solution as an initial condition. See Tawhai, M.H., PhD thesis for
!!! further explanation.

    use arrays,only: dp,node_field,num_nodes
    use exports,only: export_node_field
    use indices,only: nj_conc1,nu_conc1
    use other_consts
    use diagnostics, only: enter_exit
    use geometry, only: volume_of_mesh
    use solve,only: pmgmres_ilu_cr
    implicit none

    integer,intent(in) :: fileid,inr_itr_max,out_itr_max
    real(dp),intent(in) :: diffusion_coeff,dt,initial_volume,&
         inlet_concentration,inlet_flow,solve_tolerance,time_end,time_start
    logical,intent(in) :: inspiration

    ! Local variables
    real(dp),allocatable :: solution(:)

    integer :: MatrixSize,nonzeros,ncol,nentry, &
         noffset_entry,noffset_row,np,nrow,nrow_BB,SolverFlag
    real(dp) :: AA,BB,current_mass,current_volume, &
         mass_error,theta,time,volume_error,volume_tree,zero_tol=1.0e-8_dp,&
         mass0,mass1,mass_error_deform,mass_error_track,mass_error_solve
    logical :: carryon
    character(len=60) :: sub_name

    ! #############################################################################

    sub_name = 'solve_gasmix'
    call enter_exit(sub_name,1)

    ! allocatable array to store the current solution
    if(.not.allocated(solution)) allocate(solution(num_nodes))

    ! the weighting for matrices in the reduced system: A = M+K*dt*theta; B = -K*c^(n)*dt
    ! can change this such that fully implicit or fully explicit, however 2/3 generally most stable
    theta=2.0_dp/3.0_dp

    ! get the sparsity arrays for the reduced system. uses compressed row format.
    call reduce_gasmix(MatrixSize,nonzeros,noffset_entry,noffset_row,inspiration)

    time = time_start ! initialise the time

    carryon = .true. ! logical for whether solution continues

    ! main time-stepping loop:  time-stepping continues while 'carryon' is true
    do while (carryon) !

       time = time + dt ! increment time

       if(time_end + 0.5_dp*dt - time .gt. zero_tol)then !continue

          call calc_mass(nj_conc1,nu_conc1,mass0) ! calculate model mass, for error check
          call airway_mesh_deform(dt,initial_volume,inlet_flow) ! change model size by dV
          call calc_mass(nj_conc1,nu_conc1,mass1) ! calculate model mass, for error check
          mass_error_deform = 100 * (mass1-mass0)/mass0 ! error from deforming
          call update_unit_mass(dt,inlet_concentration,inlet_flow) ! track gas into/out units

          call track_back(dt,inlet_concentration,inlet_flow,inspiration) ! track gas along branches

          call calc_mass(nj_conc1,nu_conc1,mass0) ! calculate model mass, for error check
          mass1 = mass1 + inlet_flow*dt*node_field(nj_conc1,1) ! expected mass following tracking
          mass_error_track = 100 * (mass0-mass1)/mass1 ! error from tracking units and branches

          ! assemble the element matrices. Element matrix calculation can be done directly
          ! (based on assumption of interpolation functions) or using Gaussian interpolation.
          call assemble_gasmix(diffusion_coeff,nonzeros_unreduced)

          ! initialise the values in the solution matrices
          global_AA(1:nonzeros) = 0.0_dp ! equivalent to M in Tawhai thesis
          global_BB(1:num_nodes) = 0.0_dp ! equivalent to K in Tawhai thesis
          global_R(1:num_nodes) = 0.0_dp

          ! Assemble the reduced system of matrices
          do np=1,num_nodes ! Loop over rows of unreduced system
             nrow = np ! conveniently true for the way we set up our models
             ! different boundary conditions are applied during inspiration and
             ! expiration: Dirichlet at model entry during inspiration (concentration
             ! = inlet_concentration), and Neumann at model entry during expiration (dcdt = 0)
             if(.not.inspiration)then
                BB = global_R(nrow)         !get reduced R.H.S.vector
                do nentry = sparsity_row(nrow),sparsity_row(nrow+1)-1  !each row entry
                   ncol = sparsity_col(nentry)
                   BB = BB-global_K(nentry)*node_field(nj_conc1,ncol)*dt ! -K*c^(n)*dt
                   AA = global_M(nentry) + dt*theta*global_K(nentry) ! M+K*dt*theta
                   global_AA(nentry) = global_AA(nentry) + AA
                enddo
                global_BB(nrow) =  global_BB(nrow) + BB
             elseif(inspiration.and.np.ne.1)then !not first row
                BB = global_R(nrow)         !get reduced R.H.S.vector
                do nentry = sparsity_row(nrow),sparsity_row(nrow+1)-1  !each row entry
                   ncol = sparsity_col(nentry)
                   BB = BB-global_K(nentry)*node_field(nj_conc1,ncol)*dt
                   AA = global_M(nentry) + dt*theta*global_K(nentry) !M+K*dt*theta
                   if(ncol.ne.1)then ! not first column
                      global_AA(nentry-noffset_entry) = &
                           global_AA(nentry-noffset_entry) + AA
                   endif
                enddo
                global_BB(nrow-noffset_row) = &
                     global_BB(nrow-noffset_row) + BB
             endif
          enddo

          solution(1:num_nodes) = 0.0_dp

          ! Call a solver to solve the system of reduced equations.
          ! Here we use an iterative solver (GMRES == Generalised Minimal
          ! RESidual method). The solver requires the solution matrices to
          ! be represented in compressed row format.
          call pmgmres_ilu_cr (MatrixSize,NonZeros,reduced_row,&
               reduced_col,global_AA,solution,global_BB,&
               out_itr_max,inr_itr_max,solve_tolerance,solve_tolerance,SolverFlag)

          ! transfer the solver solution (in 'Solution') to the node field array
          do np = 1,num_nodes
             if(.not.inspiration.or.(inspiration.and.np.gt.1))then
                if(inspiration)then
                   nrow_BB=np-1
                else
                   nrow_BB=np
                endif
                node_field(nj_conc1,np) = node_field(nj_conc1,np) &
                     + Solution(nrow_BB) !c^(n+1)=c^(n)+dc
                node_field(nj_conc1,np) = MAX(0.d0,node_field(nj_conc1,np))
                node_field(nj_conc1,np) = MIN(1.d0,node_field(nj_conc1,np))
             endif
          enddo

!          if(.not.inspiration)then
!             call smooth_expiration(nj_conc1)
!          endif
          ! estimate the volume and mass errors
          call calc_mass(nj_conc1,nu_conc1,current_mass)
          mass_error_solve = 100 * (current_mass-mass0)/mass0

          call volume_of_mesh(current_volume,volume_tree)
          ideal_mass = ideal_mass + inlet_flow*dt*node_field(nj_conc1,1)
          volume_error = 1.0e+2_dp*(current_volume - (initial_volume +  &
               total_volume_change))/(initial_volume + total_volume_change)

          if(ideal_mass.gt.0.0_dp)then
             mass_error = 1.0e+2_dp*(current_mass - ideal_mass)/ideal_mass
          else
             mass_error = 0.0_dp
          endif
          ! output current solution to screen
          write(*,'(F10.3,7(F10.2),F12.3)') &
               time,inlet_flow/1.0e+6_dp,total_volume_change/1.0e+6_dp, &
               current_volume/1.0e+6_dp,volume_error,mass_error_deform,&
               mass_error_track,mass_error_solve,node_field(nj_conc1,1)
          if(.not.inspiration)then
             write(fileid,'(F10.3,3(D14.5),4(F10.2),D14.5)') &
                  time,inlet_flow,total_volume_change, &
                  current_volume,volume_error,mass_error_deform,&
                  mass_error_track,mass_error_solve,node_field(nj_conc1,1)
          endif
       else
          carryon = .false.
       endif
    enddo !carryon

    call enter_exit(sub_name,2)

  end subroutine solve_gasmix

!!!########################################################################

  subroutine sparse_gasmix
    use arrays,only: elem_cnct,elem_nodes,num_elems

    implicit none

    integer :: n_unreduced,i,ncol,ne,ne2,np1,np2,nrow

    sparsity_row(1) = 1
    n_unreduced = 1

    do ne=1,num_elems ! note using local numbering
       if(elem_cnct(-1,0,ne).eq.0)then !at the inlet
          np1=elem_nodes(1,ne) ! start node
          nrow=np1
          do i=1,2
             np2=elem_nodes(i,ne)
             ncol=np2
             sparsity_col(n_unreduced)=ncol
             n_unreduced=n_unreduced+1
          enddo
          sparsity_row(nrow+1)=n_unreduced
       endif

       np1=elem_nodes(2,ne) !end node
       nrow=np1
       do i=1,2
          np2=elem_nodes(i,ne)
          ncol=np2
          sparsity_col(n_unreduced)=ncol
          n_unreduced=n_unreduced+1
       enddo
       do i=1,elem_cnct(1,0,ne) ! for each child branch
          ne2=elem_cnct(1,i,ne)
          np2=elem_nodes(2,ne2)
          ncol=np2
          sparsity_col(n_unreduced)=ncol
          n_unreduced=n_unreduced+1
       enddo
       sparsity_row(nrow+1)=n_unreduced
    enddo !noelem

    NonZeros_unreduced = n_unreduced - 1

  end subroutine sparse_gasmix

!!!####################################################################

  subroutine update_unit_mass(dt,inlet_concentration,inlet_flow)
    use arrays,only: dp,elem_field,elem_nodes,num_units,units,unit_field
    use indices,only: ne_Vdot,nu_conc1,nu_vol
    use geometry, only: volume_of_mesh
    use diagnostics, only: enter_exit
    implicit none

    real(dp),intent(in) :: dt,inlet_concentration,inlet_flow

    ! Local parameters
    integer :: ne,np,nunit
    real(dp) :: concentration,mass_change,&
         unit_mass,volume_change,sum_mass_change
    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_unit_mass'
    call enter_exit(sub_name,1)

    ! calculate the change in volume of elastic units
    ! note that conducting airways are assumed rigid
    sum_mass_change=0.0_dp
    do nunit=1,num_units
       ne = units(nunit) ! local element number
       np = elem_nodes(2,ne) ! end node
       unit_mass = unit_field(nu_vol,nunit)*unit_field(nu_conc1,nunit)
       volume_change = elem_field(ne_Vdot,ne)*inlet_flow*dt
       if(inlet_flow.ge.0.0_dp)then !filling units
          call track_to_location(ne,concentration,mass_change,&
               dt,inlet_concentration,inlet_flow)
        else ! emptying units
          mass_change = unit_field(nu_conc1,nunit)*volume_change
       endif
       unit_mass = unit_mass + mass_change
       unit_field(nu_conc1,nunit) = unit_mass/unit_field(nu_vol,nunit)
       sum_mass_change = sum_mass_change+mass_change
    enddo

    call enter_exit(sub_name,2)

  end subroutine update_unit_mass

!!!####################################################################

  subroutine airway_mesh_deform(dt,initial_volume,inlet_flow)

!!! changes the size of elastic units and airways, based on pre-computed
!!! fields for flow into the units (ne_Vdot) and element volume change (ne_dvdt).
!!! The concentration in the units is adjusted to make sure that mass is conserved
!!! when the unit changes volume. If a unit changes in volume, then both the radius
!!! and length scale as the cube root of volume change. For alveolated airways (where
!!! a/A is < 1) the outer radius (indexed as nj_radius) is scaled and a/A unchanged.

    use arrays,only: dp,elem_field,elem_nodes,&
         num_elems,num_units,units,unit_field
    use indices,only: ne_dvdt,ne_Vdot,ne_length,ne_radius,&
         ne_vol,nu_conc1,nu_vol
    use geometry, only: volume_of_mesh
    use diagnostics, only: enter_exit
    implicit none

    real(dp),intent(in) :: dt,initial_volume,inlet_flow

    ! Local parameters
    integer :: ne,np,nunit
    real(dp) :: current_volume,ratio,tree_volume,unit_mass,volume_change
    character(len=60) :: sub_name

    ! #################################################################

    sub_name = 'airway_mesh_deform'
    call enter_exit(sub_name,1)

!!! calculate the change in volume of elastic units:
    ! mass(unit) = volume(unit) * concentration(unit)
    ! dv(unit) = flow(unit) * dt
    ! volume(new) = volume(old) + dv
    ! concentration(new) = mass/volume(new)
    do nunit=1,num_units !for each of the elastic units
       ne = units(nunit) ! local element number
       np = elem_nodes(2,ne) ! end node, attaches to unit
       unit_mass = unit_field(nu_vol,nunit)*unit_field(nu_conc1,nunit) ! m = v(unit)*c
       volume_change = elem_field(ne_Vdot,ne)*inlet_flow*dt ! dv = q(unit)*dt
       unit_field(nu_vol,nunit) = unit_field(nu_vol,nunit) + volume_change ! v(new)
       ! important: adjust the concentration to maintain mass conservation
       unit_field(nu_conc1,nunit) = unit_mass/unit_field(nu_vol,nunit)
    enddo

!!! calculate the inflation/deflation of multi-branching acini (or other branches)
    ! dv(element) = dvdt(element) * total_volume_change
    ! ratio_of_volumes = (volume(element)+dv(element))/volume(element)
    ! volume(new) = volume(old) * ratio_of_volumes
    ! length(new) = length(old) * cuberoot(ratio_of_volumes)
    ! radius(new) = radius(old) * cuberoot(ratio_of_volumes)
    do ne=1,num_elems
       volume_change = elem_field(ne_dvdt,ne)*inlet_flow*dt
       ratio = (volume_change+elem_field(ne_vol,ne))/elem_field(ne_vol,ne)
       elem_field(ne_vol,ne)=elem_field(ne_vol,ne)*ratio
       elem_field(ne_length,ne)=elem_field(ne_length,ne)*(ratio)**(1/3)
       elem_field(ne_radius,ne)=elem_field(ne_radius,ne)*(ratio)**(1/3)
    enddo

!!! calculate the total model volume
    call volume_of_mesh(current_volume,tree_volume)

    total_volume_change = current_volume-initial_volume

    call enter_exit(sub_name,2)

  end subroutine airway_mesh_deform

!!! ######################################################################

  subroutine smooth_expiration(nj_field)
    use arrays,only: elem_cnct,elem_ordrs,num_elems
    implicit none

    integer,intent(in) :: nj_field

    integer :: ne,ne_next,ngen,ngen_next,num_smooth
    integer,allocatable :: nsmooth_list(:)
    logical,allocatable :: smoothed(:)

    allocate(nsmooth_list(num_elems))
    allocate(smoothed(num_elems))
    smoothed(1:num_elems)=.FALSE.

    do ne=1,num_elems
       ngen=elem_ordrs(1,ne)
       if(ngen.gt.1)then
          if(elem_cnct(1,0,ne).eq.1)then !do only for discretised branch
             if(.not.smoothed(ne))then  !hasn't been done as part of branch
                ne_next=elem_cnct(1,1,ne) !element # of adjacent element
                ngen_next=elem_ordrs(1,ne_next) !generation of adjacent element
                if(ngen.eq.ngen_next)then !same generation, therefore same branch
                   num_smooth = 1
                   nsmooth_list(1) = ne
                   do while(ne_next.ne.0)
                      num_smooth = num_smooth+1
                      nsmooth_list(num_smooth)=ne_next
                      smoothed(ne_next)=.TRUE.
                      if(elem_cnct(1,0,ne_next).eq.1)then
                         ne_next=elem_cnct(1,1,ne_next)
                         ngen_next=elem_ordrs(1,ne_next)
                         if(ngen.ne.ngen_next)  ne_next=0
                      else
                         ne_next=0
                      endif
                   enddo     !WHILE
                   call smooth_expiration_linear(nj_field,nsmooth_list,num_smooth)
                endif          !gen
             endif             !SMOOTHED
          else                 !(NXI=0 or 2)
             smoothed(ne)=.TRUE.
          endif                !NXI
       endif
    enddo !noelem(ne)
    
    deallocate(nsmooth_list)
    deallocate(smoothed)

  end subroutine smooth_expiration

!!! ######################################################################

  subroutine smooth_expiration_linear(nj_field,nsmooth_list,num_smooth)
    use arrays,only: dp,elem_field,elem_nodes,node_field,num_elems
    use indices,only: ne_length
    implicit none

!!! Parameter List
    integer,intent(in) :: nj_field,num_smooth,nsmooth_list(*)

!!! Local Variables
    integer :: n,ne,np
    real(dp) :: conc_end,conc_start,intercept,slope
    real(dp),allocatable :: length(:)

    allocate(length(num_elems))

    !for first node
    ne=nsmooth_list(1)
    np=elem_nodes(1,ne)
    conc_start=node_field(nj_field,np)

    !for last node
    ne=nsmooth_list(num_smooth)
    np=elem_nodes(2,ne)
    conc_end=node_field(nj_field,np)

    do n=1,num_smooth
       ne=nsmooth_list(n)
       if(n.eq.1)then
          length(n)=elem_field(ne_length,ne)
       else
          length(n)=length(n-1)+elem_field(ne_length,ne)
       endif
    enddo
    !calculate linear equation
    intercept=conc_start
    slope=(conc_end-intercept)/length(num_smooth)

    !update solution values for the branch
    ne=nsmooth_list(1)
    np=elem_nodes(1,ne)

    intercept=node_field(nj_field,np)
    ne=nsmooth_list(num_smooth)
    np=elem_nodes(2,ne)

    slope=(node_field(nj_field,np)-intercept)/length(num_smooth)

    node_field(nj_field,np)=intercept
    do n=1,num_smooth
       ne=nsmooth_list(n)
       np=elem_nodes(2,ne)
       node_field(nj_field,np)=intercept+slope*length(n)
    enddo
    
    deallocate(length)

  end subroutine smooth_expiration_linear

!!! ##################################################################

  subroutine transfer_flow_vol_from_units()
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_TRANSFER_FLOW_VOL_FROM_UNITS" :: TRANSFER_FLOW_VOL_FROM_UNITS

    use arrays,only: dp,elem_cnct,elem_field,elem_symmetry,&
         elem_units_below,expansile,units,num_elems,num_units,&
         unit_field
    use indices,only: ne_dvdt,ne_Vdot,ne_length,ne_radius,ne_vol,nu_vol
    implicit none

!!! Parameters

!!! Local variables
    integer :: nchild,ne_child,ne_parent,ne_stem,&
         nparent,nunit,num_list,num_list_temp,num_list_total
    integer,allocatable :: elem_list(:),elem_list_temp(:),&
         elem_list_total(:)
    real(dp),allocatable :: vol_below(:)
    real(dp) :: ratio


    allocate(vol_below(num_elems))
    vol_below(1:num_elems) = elem_field(ne_vol,1:num_elems)
    elem_field(ne_dvdt,1:num_elems) = 0.0_dp !rate of change of element volume
    expansile(1:num_elems) = .false.

    allocate(elem_list(num_elems))
    allocate(elem_list_temp(num_elems))
    allocate(elem_list_total(num_elems))

    do nunit = 1,num_units
       ne_stem=units(nunit) ! the element number supplying a unit
       num_list=1
       elem_list(1)=ne_stem
       elem_list(2:elem_units_below(ne_stem)) = 0

       num_list_temp=0
       elem_list_temp(1:elem_units_below(ne_stem)) = 0

       num_list_total = 0
       elem_list_total = 0 !record all elements below the unit

       do while(num_list.gt.0)
          do nparent=1,num_list ! for each 'parent' element in the list
             ne_parent=elem_list(nparent) ! the 'parent' element number
             do nchild=1,elem_cnct(1,0,ne_parent) ! for each child branch of 'parent'
                ne_child = elem_cnct(1,nchild,ne_parent)
                num_list_temp=num_list_temp+1
                elem_list_temp(num_list_temp)=ne_child
                num_list_total=num_list_total+1
                elem_list_total(num_list_total)=ne_child
                expansile(ne_child)=.true.
             enddo !nchild
          enddo !nparent
          num_list = num_list_temp
          num_list_temp = 0
          elem_list(1:num_list) =elem_list_temp(1:num_list)
          elem_list_temp(1:elem_units_below(ne_stem)) = 0
       enddo !while

       ! calculate the volume below each branch
       do nchild = num_list_total,1,-1
          ne_child=elem_list_total(nchild)
          ne_parent=elem_cnct(-1,1,ne_child)
          vol_below(ne_parent) = vol_below(ne_parent) &
               + DBLE(elem_symmetry(ne_child))*vol_below(ne_child)
       enddo !nchild

!!! scale the elements so same total size as original unit size, and set the
!!! rate of volume change for an element to be its volume/unit volume*stem flow
       ratio = unit_field(nu_vol,nunit)/(vol_below(ne_stem)&
            -elem_field(ne_vol,ne_stem))
       do nchild=1,num_list_total
          ne_child=elem_list_total(nchild)
          elem_field(ne_vol,ne_child)=elem_field(ne_vol,ne_child)*ratio
          elem_field(ne_length,ne_child)=elem_field(ne_length,ne_child)*(ratio)**(1/3)
          elem_field(ne_radius,ne_child)=elem_field(ne_radius,ne_child)*(ratio)**(1/3)
          elem_field(ne_dvdt,ne_child)=elem_field(ne_Vdot,ne_stem)*&
               elem_field(ne_vol,ne_child)/&
               (unit_field(nu_vol,nunit)-elem_field(ne_vol,ne_stem))
       enddo
!!! calculate the flow field. flow(in) = flow(out) + DV
       do nparent=num_list_total,1,-1
          ne_parent=elem_list_total(nparent)
          elem_field(ne_Vdot,ne_parent)=elem_field(ne_dvdt,ne_parent)
          do nchild=1,elem_cnct(1,0,ne_parent)
             ne_child=elem_cnct(1,nchild,ne_parent)
               elem_field(ne_Vdot,ne_parent)=elem_field(ne_Vdot,ne_parent)+&
                  elem_field(ne_Vdot,ne_child)*elem_symmetry(ne_child)
          enddo !nchild
       enddo !nparent
    enddo !nunit

    num_units = 0

    deallocate(elem_list)
    deallocate(elem_list_temp)
    deallocate(elem_list_total)
    deallocate(vol_below)

  end subroutine transfer_flow_vol_from_units

!!! ##################################################################

  subroutine track_back(dt,inlet_concentration,inlet_flow,inspiration)
    use arrays,only: dp,elem_cnct,elem_field,elem_nodes,elem_symmetry,&
         elems_at_node,node_field,num_nodes,num_units,unit_field,units
    use indices,only: ne_dvdt,ne_Vdot,ne_vol,nj_conc1,nu_conc1
    implicit none

!!! Parameters
    real(dp),intent(in) :: dt,inlet_concentration,inlet_flow
    logical,intent(in) :: inspiration

!!! Local variables
    integer :: elems_below(1024),i,j,ne,ne_parent,nextra,np,np1,np2,nunit,num_to_check
    real(dp) :: branch_fraction(1024),branch_time(1024),concentration,cumulative_mass,&
         flow_fraction,flow_parent,local_xi,time_through_element,total_time
    real(dp),allocatable :: initial_conc(:)

    allocate(initial_conc(num_nodes))
    initial_conc = 0.0_dp

    if(inspiration)then
       do np=2,num_nodes
          ne=elems_at_node(np,1) ! first element; will be most proximal branch
          call track_to_location(ne,concentration,cumulative_mass,&
               dt,inlet_concentration,inlet_flow)
          initial_conc(np) = concentration
       enddo !np
       initial_conc(1) = inlet_concentration

    else !tracking for expiration
       ! set the concentration at terminal branch to the same as in elastic unit (if attached)
       do nunit=1,num_units
          ne=units(nunit)
          np=elem_nodes(2,ne)
          initial_conc(np) = unit_field(nu_conc1,nunit)
          node_field(nj_conc1,np) = unit_field(nu_conc1,nunit)
       enddo
       do np=1,num_nodes
          if(np.eq.1)then
             elems_below(1) = elems_at_node(np,1)
             branch_time(1) = 0.0_dp
             branch_fraction(1) = 1.0_dp
             num_to_check = 1
             ne_parent = elems_at_node(np,1)
          else
             if(elems_at_node(np,0).eq.1)then  !only one adjoining element, so a terminal node
                num_to_check = 0
                initial_conc(np) = node_field(nj_conc1,np)
             else
                ne_parent = elems_at_node(np,1)
                do j=2,elems_at_node(np,0)
                   elems_below(j-1) = elems_at_node(np,j) ! check each element
                   branch_time(j-1) = 0.0_dp
                   branch_fraction(j-1) = 1.0_dp
                enddo
                num_to_check = elems_at_node(np,0)-1
             endif
          endif
          do while(num_to_check.gt.0) ! while still some to check
             nextra=0            !records the addition of elements to check
             do j=1,num_to_check
                ne=elems_below(j)
                if(np.gt.1) ne_parent = elem_cnct(-1,1,ne)
                total_time=branch_time(j)
                time_through_element = dabs(elem_field(ne_vol,ne)/&
                     (inlet_flow*elem_field(ne_Vdot,ne)))
                if(total_time+time_through_element.ge.dt)then
                   !     A material point is within this element
                   local_xi = (dt-total_time)/time_through_element
                   np1 = elem_nodes(1,ne)
                   np2 = elem_nodes(2,ne)
                   concentration = node_field(nj_conc1,np1)*(1.0_dp-local_xi)&
                        +node_field(nj_conc1,np2)*local_xi
!!! this isn't quite right: needs to take into account the rate of change of parent branch size
!!! also should do this when estimating time to pass through a branch.
                   flow_parent = elem_field(ne_Vdot,ne_parent) - elem_field(ne_dvdt,ne_parent)
                   flow_fraction = branch_fraction(j)*(elem_field(ne_Vdot,ne)*&
                        elem_symmetry(ne))/flow_parent
                   initial_conc(np) = initial_conc(np) + concentration*flow_fraction
                else
                   !     need to check the next elements
                   if(elem_cnct(1,0,ne).eq.0)then !set to terminal node concentration
                      np2 = elem_nodes(2,ne)
                      concentration = node_field(nj_conc1,np2)
                      flow_parent = elem_field(ne_Vdot,ne_parent) - elem_field(ne_dvdt,ne_parent)
                      flow_fraction = branch_fraction(j)*(elem_field(ne_Vdot,ne)*&
                           elem_symmetry(ne))/flow_parent
                      initial_conc(np) = initial_conc(np) + &
                           concentration*flow_fraction
                   else
                      do i=1,elem_cnct(1,0,ne) ! for each child element
                         nextra=nextra+1
                         elems_below(num_to_check+nextra) = elem_cnct(1,i,ne)
                         branch_time(num_to_check+nextra) = &
                              total_time+time_through_element
                         flow_parent = elem_field(ne_Vdot,ne_parent) - elem_field(ne_dvdt,ne_parent)
                         branch_fraction(num_to_check+nextra) = branch_fraction(j)*&
                              (elem_field(ne_Vdot,ne)*elem_symmetry(ne))/flow_parent
                      enddo   ! i
                   endif      ! terminal or not
                endif         ! exceed the time criteria?
             enddo               ! j

             if(nextra.gt.64)then
                write(*,*) 'Need to decrease the time step'
             endif
             do j=1,nextra
                elems_below(j) = elems_below(num_to_check+j)
                branch_time(j) = branch_time(num_to_check+j)
                branch_fraction(j)= branch_fraction(num_to_check+j)
             enddo               !j, for nextra
             num_to_check=nextra
          enddo                  ! WHILE num_to_check.GT.0
       enddo                     !np
    endif !inspiration

    node_field(nj_conc1,1:num_nodes) = initial_conc(1:num_nodes)

    deallocate(initial_conc)

  end subroutine track_back

!!! ##################################################################

  subroutine track_to_location(ne,concentration,cumulative_mass,&
       dt,inlet_concentration,inlet_flow)

    use arrays,only: dp,elem_cnct,elem_field,elem_nodes,node_field
    use indices,only: ne_Vdot,ne_vol,nj_conc1
    implicit none

!!! Parameters
    integer,intent(in) :: ne
    real(dp),intent(in) :: dt,inlet_concentration,inlet_flow
    real(dp) :: concentration,cumulative_mass

!!! Local variables
    integer :: ne0,np1,np2
    real(dp) :: local_xi,mean_concentration,proportion_flow,&
         time_through_element,total_time
    logical :: continue

    cumulative_mass = 0.0_dp !used only if there is a 'unit' appended

    time_through_element = elem_field(ne_vol,ne)/(inlet_flow*elem_field(ne_Vdot,ne))
    total_time = time_through_element
    np1 = elem_nodes(1,ne)
    np2 = elem_nodes(2,ne)

    if(total_time.ge.dt)then !location is within this element
       local_xi = 1.0_dp-dt/time_through_element
       concentration = (1.0_dp-local_xi)*node_field(nj_conc1,np1) + &
            local_xi*node_field(nj_conc1,np2)

       mean_concentration = 0.5_dp*(concentration+node_field(nj_conc1,np2))
       cumulative_mass = cumulative_mass + mean_concentration* &
            elem_field(ne_vol,ne)*(1.0_dp-local_xi)
    else

       mean_concentration = 0.5_dp*(node_field(nj_conc1,np1)+node_field(nj_conc1,np2))
       cumulative_mass = cumulative_mass + mean_concentration*elem_field(ne_vol,ne)
       continue = .true.
       ne0=ne

       do while (continue)

          if(elem_cnct(-1,1,ne0).gt.0)then
             ne0 = elem_cnct(-1,1,ne0) !parent element
             np1 = elem_nodes(1,ne0)
             np2 = elem_nodes(2,ne0)

             time_through_element = elem_field(ne_vol,ne0)/&
                  (inlet_flow*elem_field(ne_Vdot,ne0))
             total_time = total_time + time_through_element

             proportion_flow = elem_field(ne_Vdot,ne)/elem_field(ne_Vdot,ne0)

             if(total_time.ge.dt)then !location is within this element
                local_xi = (total_time-dt)/time_through_element !logic correct
                concentration = (1.0_dp-local_xi)*node_field(nj_conc1,np1) + &
                     local_xi*node_field(nj_conc1,np2)

                mean_concentration = 0.5_dp*(concentration+node_field(nj_conc1,np2))
                cumulative_mass = cumulative_mass + mean_concentration* &
                     elem_field(ne_vol,ne0)*(1.0_dp-local_xi)*proportion_flow
                continue=.false.
             else
                concentration = 0.5_dp*(node_field(nj_conc1,np1)+node_field(nj_conc1,np2))
                cumulative_mass = cumulative_mass + concentration* &
                     elem_field(ne_vol,ne0)*proportion_flow
             endif
          else
             continue=.false.
             concentration = inlet_concentration

             proportion_flow = elem_field(ne_Vdot,ne)/elem_field(ne_Vdot,1)
             cumulative_mass = cumulative_mass + concentration* &
                  inlet_flow*(dt-total_time)*proportion_flow
          endif
       enddo
    endif

  end subroutine track_to_location

!!!##############################################################################

  subroutine normalised_slope(num_breaths,dt,t_fit_0,t_fit_1,&
       time_expiration,vol_expired,vol_frc,resultsfile)

    use arrays,only: dp

    implicit none

!!! Parameters
    integer,intent(in) :: num_breaths
    real(dp),intent(in) :: dt,t_fit_0,t_fit_1,time_expiration,vol_expired,vol_frc
    character(len=*),intent(in):: resultsfile

!!! Local variables
    integer :: fh,i,ibeg,iend,ios=0, i_ss_end,nbreath,nlabel,num_data
    real(dp) :: exp_vol,exp_conc,exp_time,exp_vol_last,flow,&
         norm_slope,Sacin,Scond,slope,Sn0,Sn1,Sn2,sum_gas_expired,sum_n2_expired,&
         sum_x,sum_xsq,sum_xy,sum_y,turnovers,&
         turnover_0,turnover_1,turnover_2,vol_end,vol_start
    character(len=100) :: buffer
    logical :: first,found_1,found_2

    vol_start = t_fit_0 * vol_expired ! the expired volume where fit starts from
    vol_end = t_fit_1 * vol_expired !the expired volume where fit goes to

    fh=20
    open(fh, file=resultsfile, status='old') !open existing file

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.
    ios=0
    first = .true.
    found_1 = .false.
    found_2 = .false.

    do nbreath = 1,num_breaths
       exp_time = 0.0_dp
       exp_vol_last = 0.0_dp
       sum_gas_expired = 0.0_dp
       sum_x = 0.0_dp
       sum_y = 0.0_dp
       sum_xsq = 0.0_dp
       sum_xy = 0.0_dp
       i = 0
       do while (exp_time < time_expiration)
          read(fh, '(A)', iostat=ios) buffer
          if(ios .ne. 0)then
             write(*,'('' Premature end of file, at breath'',I3,'' time'',F8.3)') nbreath,exp_time
             stop
          endif
          ! line contains: time, flow, DV, vol, 4 x %errors, expired conc
          i_ss_end = len(buffer)
          do nlabel = 1,9
             ibeg = index(buffer," ") + 1 !get location of first integer beyond ws in string
             buffer = adjustl(buffer(ibeg:i_ss_end)) ! get info beyond ws, remove leading ws
             iend = index(buffer," ") !get location of next ws in sub-string
             select case (nlabel)
             case(2)
                read (buffer(1:iend-1), *, iostat=ios) flow
             case(9)
                read (buffer(1:iend-1), *, iostat=ios) exp_conc
             end select
          enddo
          exp_vol = exp_vol_last + dabs(flow * dt)
          if(exp_vol.ge.vol_start.and.exp_vol.le.vol_end)then
             i=i+1
             sum_x = sum_x + exp_vol
             sum_y = sum_y + exp_conc
             sum_xsq = sum_xsq + exp_vol**2
             sum_xy = sum_xy + exp_vol * exp_conc
          endif
          sum_gas_expired = sum_gas_expired + dabs(flow) * dt * exp_conc
          exp_time = exp_time + dt
          exp_vol_last = exp_vol
       enddo ! while time
       num_data = i

       slope = -(num_data*sum_xy-sum_x*sum_y)/(num_data*sum_xsq-sum_x**2)*1.0e+6_dp !fractional conc/L
       sum_n2_expired = 1.0_dp - sum_gas_expired/vol_expired !fractional
       norm_slope=slope/sum_n2_expired

       turnovers = nbreath*vol_expired/vol_frc

       write(*,'(I3,F7.2,F12.4,D13.4,F10.4)') nbreath,turnovers,norm_slope,slope,sum_n2_expired

       if(first)then
          first = .false.
          turnover_0 = turnovers
          Sn0 = norm_slope
       endif
       if(.not.found_1.and.turnovers.ge.1.5_dp)then
          found_1 = .true.
          turnover_1 = turnovers
          Sn1 = norm_slope
       else if(.not.found_2.and.turnovers.ge.6.0_dp)then
          found_2 = .true.
          turnover_2 = turnovers
          Sn2 = norm_slope
       endif

    enddo ! num_breaths

    if(found_1.and.found_2)then
       Scond=(Sn2-Sn1)/(turnover_2-turnover_1)
       Sacin=Sn0-Scond*turnover_0
       write(*,'('' Sacin = '',F10.4,'' L.-1, Scond = '',F10.4)') Sacin,Scond
       write(*,'('' Calculated using TO = '',F10.4,'' and '',F10.4)') turnover_1,turnover_2
    else
       write(*,'('' Turnovers is'',F10.4,''; too small for calculating Scond and Sacin'')') turnovers
    endif

    close(fh)

  end subroutine normalised_slope

end module gasmix
