module wave_transmission

!*Brief Description:* Simulating wave propagation in a 1D tree structure
!
!*LICENSE:*
!
!
!
!*Full Description:*
!Simulating wave propagation in a 1D tree structure
!
  use arrays
  use capillaryflow
  use diagnostics
  use indices
  use other_consts
  use math_utilities
  use pressure_resistance_flow

  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private
  public evaluate_wave_transmission


contains
!
!##############################################################################
!
subroutine evaluate_wave_transmission(grav_dirn,grav_factor,n_time,heartrate,a0,no_freq,a,b,n_adparams,&
  admittance_param,n_model,model_definition,cap_model,remodeling_grade,bc_type,lobe_imped)

  integer, intent(in) :: n_time
  real(dp), intent(in) :: heartrate
  real(dp), intent(in) :: a0
  integer, intent(in):: no_freq
  real(dp), intent(in) :: a(no_freq)
  real(dp), intent(in) :: b(no_freq)
  integer, intent(in) :: n_adparams
  real(dp), intent(in) :: admittance_param(n_adparams)
  integer, intent(in) :: n_model
  real(dp), intent(in) :: model_definition(n_model)
  integer, intent(in) :: grav_dirn
  integer, intent(in) :: cap_model
  real(dp), intent(in) :: remodeling_grade
  character(len=60) :: bc_type
  character(len=60) :: lobe_imped

  type(all_admit_param) :: admit_param
  type(fluid_properties) :: fluid
  type(elasticity_param) :: elast_param

  character(len=60) :: mesh_type
  real(dp) :: viscosity
  real(dp) :: density
  real(dp) :: harmonic_scale
  real(dp) :: steady_flow
  complex(dp), allocatable :: eff_admit(:,:)
  complex(dp), allocatable :: char_admit(:,:)
  complex(dp), allocatable :: reflect(:,:)
  complex(dp), allocatable :: prop_const(:,:)
  complex(dp), allocatable :: p_factor(:,:)
  complex(dp), allocatable :: q_factor(:,:)
  real(dp), allocatable :: forward_pressure(:)
  real(dp), allocatable :: reflected_pressure(:)
  real(dp), allocatable :: forward_pressure_previous(:)
  real(dp), allocatable :: reflected_pressure_previous(:)
  real(dp), allocatable :: forward_flow(:)
  real(dp), allocatable :: reflected_flow(:)
  real(dp), allocatable :: terminals_radius(:),WSS(:),elem_rad(:),alpha(:)
  real(dp), allocatable :: p_terminal(:),p_previous(:),terminal_flow(:),elem_flow(:)
  integer :: min_art,max_art,min_ven,max_ven,min_cap,max_cap,ne,nu,nt,nf,np,np_previous,ne_previous
  character(len=30) :: tree_direction,mechanics_type
  real(dp) start_time,end_time,dt,time,omega,delta_p,Ppl
  real(dp) grav_vect(3),grav_factor,mechanics_parameters(2)
  integer :: AllocateStatus,fid=10,fid2=20,fid3=30,fid4=40,fid5=50,fid6=60,fid7=70,fid8=80,fid9=90
  integer :: fid10=100,fid11=110,fid12=120,fid13=130
  character(len=60) :: sub_name
  logical :: vein_found=.False.
  integer :: vein_elem=0
  integer, parameter :: num_freq = 101, num_vessels = 17
  integer :: i, j
  real :: freq(num_freq)
  character(len=5) :: vessel_names(num_vessels) = ["LUL_A", "LUL_V", "LLL_A", "LLL_V", "RUL_A", "RUL_V", "RLL_A", "RLL_V",&
  "RML_A", "RML_V", "MPA_A", "LPA_A", "RPA_A", "RBS_A", "RBS_A", "LBS_A", "LBS_V"]
  real, dimension(num_vessels,num_freq) :: impedance, phase
  ! character(len=10) :: units(num_units) = ["dyne.s/cm5", "radians"]
  real(dp), dimension(num_freq) :: imped

  sub_name = 'evalulate_wave_transmission'
  call enter_exit(sub_name,1)
  !!MODEL TYPE AND FLUID PROPERTIES
  !mesh_type: can be simple_tree, full_plus_ladder, full_sheet, full_tube The first can be airways, arteries, veins but no special features at the terminal level, the last one has arteries and veins connected by capillary units of some type (lung ladder acinus, lung sheet capillary bed, capillaries are just tubes represented by an element)
  if(model_definition(1).eq.1.0_dp)then
    mesh_type='simple_tree'
  elseif(model_definition(1).eq.2.0_dp)then
    mesh_type='full_plus_ladder'
  !elseif(model_definition(1).eq.3.0_dp)then
  !  full_sheet
  !elseif(model_definition(1).eq.4.0_dp)then
  ! full_tube
  else
    print *, 'ERROR: Your geometry choice has not yet been implemented'
    call exit(0)
  endif
  !viscosity and density of fluid
  if(model_definition(2).eq.1.0_dp)then !BLOOD
    viscosity=fluid%blood_viscosity
    density=fluid%blood_density
  elseif(model_definition(2).eq.2.0_dp)then !AIR
    viscosity=fluid%air_viscosity
    density=fluid%air_density
  else
    viscosity=model_definition(3)
    density=model_definition(4)
  endif

  !!SET UP ADMITTANCE MODEL
  if(admittance_param(1).eq.1.0_dp)then
    admit_param%admittance_type='lachase_standard'
    elast_param%vessel_type='elastic_hooke'
    elast_param%elasticity_parameters(1)=admittance_param(2)!Pa
    elast_param%elasticity_parameters(2)=admittance_param(3)!Unitless
    elast_param%elasticity_parameters(3)=admittance_param(4)!dummy
  elseif(admittance_param(1).eq.2.0_dp)then
    admit_param%admittance_type='lachase_modified'
    elast_param%vessel_type='elastic_hooke'
    elast_param%elasticity_parameters(1)=admittance_param(2)!Pa
    elast_param%elasticity_parameters(2)=admittance_param(3)!Unitless
    elast_param%elasticity_parameters(3)=admittance_param(4)!dummy
  elseif(admittance_param(1).eq.3.0_dp)then
    admit_param%admittance_type='zhu_chesler'
    elast_param%vessel_type='elastic_hooke'
    elast_param%elasticity_parameters(1)=admittance_param(2)!Pa
    elast_param%elasticity_parameters(2)=admittance_param(3)!Unitless
    elast_param%elasticity_parameters(3)=admittance_param(4)!dummy
  elseif(admittance_param(1).eq.4.0_dp)then
    admit_param%admittance_type='duan_zamir'
    elast_param%vessel_type='elastic_alpha'
    elast_param%elasticity_parameters(1)=admittance_param(2)!/Pa
    elast_param%elasticity_parameters(2)=admittance_param(3)!Unitless
    elast_param%elasticity_parameters(3)=admittance_param(4)!dummy
  else
    print *, 'ERROR: Your admittance model choice has not yet been implemented'
    call exit(0)
  endif

  !! SET UP PARAMETERS DEFINING OUTLET BOUNDARY CONDITIONS
  if(admittance_param(5).eq.1.0_dp)then !note we need to check that the right number of parameters have been input
    admit_param%bc_type='two_unit_wk'
    admit_param%two_parameter%admit_P1=admittance_param(6)
    admit_param%two_parameter%admit_P2=admittance_param(7)
  elseif(admittance_param(5).eq.2.0_dp)then
    admit_param%bc_type='three_unit_wk'
    admit_param%three_parameter%admit_P1=admittance_param(6)
    admit_param%three_parameter%admit_P2=admittance_param(7)
    admit_param%three_parameter%admit_P3=admittance_param(8)
  elseif(admittance_param(5).eq.4.0_dp)then
    admit_param%bc_type='two_wk_plus'
    admit_param%four_parameter%admit_P1=admittance_param(6)
    admit_param%four_parameter%admit_P2=admittance_param(7)
    admit_param%four_parameter%admit_P3=admittance_param(8)
    admit_param%four_parameter%admit_P4=admittance_param(9)
  elseif(admittance_param(5).eq.5.0_dp)then
    admit_param%bc_type='zero_reflection'
  else
    print *, 'ERROR: Your boundary condition choice has not yet been implemented'
    call exit(0)
  endif

  if(remodeling_grade.ne.0.0_dp) then
    write(*,*) 'Solving remodeling case, grade',remodeling_grade,' - make sure you are using elastic_alpha vessel type'
  endif
  mechanics_type='linear'
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

  !!Determine steady component of flow
  if(a0.eq.0.0_dp)then !Using steady flow solution at inlet as a0
    steady_flow=elem_field(ne_Qdot,1)!ASSUMING FIRST ELEMENT
  else !otherwise input a0 is used
    steady_flow=a0
  endif
  !! SET UP PARAMETERS DEFINING COMPLIANCE MODEL
  harmonic_scale=heartrate/60.0_dp !frequency of first harmonic (Hz)
  !!ALLOCATE MEMORY
  allocate (eff_admit(1:no_freq,num_elems), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for eff_admit array ***"
  allocate (char_admit(1:no_freq,num_elems), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for char_admit array ***"
  allocate (reflect(1:no_freq,num_elems), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for reflect array ***"
  allocate (prop_const(1:no_freq,num_elems), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for prop_const array ***"
  allocate (p_factor(1:no_freq,num_elems), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for p_factor array ***"
  allocate (q_factor(1:no_freq,num_elems), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for q_factor array ***"
  allocate (forward_pressure(n_time), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for forward_p array ***"
  allocate (reflected_pressure(n_time), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for reflected_p array ***"
  allocate (forward_flow(n_time), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for forward_q array ***"
  allocate (reflected_flow(n_time), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for reflected_q array ***"
  allocate (forward_pressure_previous(n_time))
  if (AllocateStatus /= 0) STOP "*** Not enough memory for forward_p_p array ***"
  allocate (reflected_pressure_previous(n_time))
  if (AllocateStatus /= 0) STOP "*** Not enough memory for reflected_p_p array ***"
  allocate (terminals_radius(n_time))
  if (AllocateStatus /= 0) STOP "*** Not enough memory for terminals_radius array ***"
  allocate (WSS(n_time))
  if (AllocateStatus /= 0) STOP "*** Not enough memory for WSS array ***"
  allocate (elem_rad(n_time))
  if (AllocateStatus /= 0) STOP "*** Not enough memory for elem_rad array ***"
  allocate (elem_flow(n_time))
  if (AllocateStatus /= 0) STOP "*** Not enough memory for elem_flow array ***"
  allocate (alpha(num_elems))
  if (AllocateStatus /= 0) STOP "*** Not enough memory for elasticity array ***"

  !initialise admittance
  char_admit=0.0_dp
  eff_admit=0.0_dp
  alpha=0.0_dp
  !calculate characteristic admittance of each branch
  call characteristic_admittance(no_freq,char_admit,prop_const,alpha,harmonic_scale, &
    density,viscosity,admit_param,elast_param,mechanics_parameters,grav_vect,remodeling_grade)

  !Apply boundary conditions to terminal units
  call boundary_admittance(no_freq,eff_admit,char_admit,admit_param,harmonic_scale,&
    density,viscosity,elast_param,mesh_type)


    ! calculate effective admittance through the tree
    if(mesh_type.eq.'full_plus_ladder')then
        min_art=1
        ne=1
        do while(elem_field(ne_group,ne).eq.0.0_dp)
            max_art=ne
            ne=ne+1
        enddo
        min_ven=ne
        do while(elem_field(ne_group,ne).eq.2.0_dp)
            max_ven=ne
            ne=ne+1
        enddo
        min_cap=ne
        max_cap=num_elems

        !vein admittance
        tree_direction='converging'
        call tree_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale,&
            min_ven,max_ven,tree_direction)
!        !cap admittance
        call capillary_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale,&
            min_cap,max_cap,elast_param,mechanics_parameters,grav_vect,cap_model)!
        !art admittance
        tree_direction='diverging'
        call tree_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale,&
            min_art,max_art,tree_direction)
    else !Assume simple tree
        tree_direction='diverging'
        min_art=1
        max_art=num_elems
        call tree_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale,&
            min_art,max_art,tree_direction)
    endif

    !calculate pressure drop through arterial tree (note to do veins too need to implement this concept thro' whole ladder model)
    !Also need to implement in reverse for veins
    call pressure_flow_factor(no_freq,p_factor,q_factor,reflect,prop_const,char_admit,eff_admit,&
    harmonic_scale,min_art,max_art,bc_type)
    if(lobe_imped.eq.'ON') then ! export lobe imped
      ! Open output file
      open(unit=10, file='lobe_imped.json', status='replace')
      ! Write header and frequency values
      freq(1) = 0
      do j = 2, no_freq+1
        freq(j) = (j-1)*harmonic_scale
      enddo
      write(10, *) "{"
      write(10, *) " ""frequency"": [", freq(1), ",", (freq(i),",",i=2,no_freq), freq(no_freq+1), "],"
      write(10, *) '"vessel_names": ["LUL_A","LUL_V","LLL_A","LLL_V","RUL_A","RUL_V","RLL_A","RLL_V","RML_A","RML_V","MPA_A_1",&
      "MPA_A_2", "MPA_A_3", "MPA_A_4", "LPA_A_1","LPA_A_2", "LPA_A_3", "RPA_A_1", "RPA_A_2", "RPA_A_3","RBS_A","RBS_V","LBS_A"&
      ,"LBS_V"],'
      ! Write impedance and phase matrices
      write(10, *) " ""impedance"":{"

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculating LUL_A impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,11))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven+10)))/elem_field(ne_Qdot,11)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,11))
      enddo
      write(10, *) " ""LUL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LUL_V impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,min_ven+10))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven)))/elem_field(ne_Qdot,11)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,min_ven+10))
      enddo
      write(10, *) " ""LUL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LLL_A impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,20))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven+19)))/elem_field(ne_Qdot,20)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,20))
      enddo
      write(10, *) " ""LLL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LLL_V impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,min_ven+19))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven)))/elem_field(ne_Qdot,20)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,min_ven+19))
      enddo
      write(10, *) " ""LLL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RUL_A impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,15))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven+14)))/elem_field(ne_Qdot,15)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,15))
      enddo
      write(10, *) " ""RUL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RUL_V impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,min_ven+14))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven)))/elem_field(ne_Qdot,15)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,min_ven+14))
      enddo
      write(10, *) " ""RUL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RLL_A impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,23))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven+22)))/elem_field(ne_Qdot,23)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,23))
      enddo
      write(10, *) " ""RLL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RLL_V impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,min_ven+22))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven)))/elem_field(ne_Qdot,23)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,min_ven+22))
      enddo
      write(10, *) " ""RLL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RML_A impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,24))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven+23)))/elem_field(ne_Qdot,24)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,24))
      enddo
      write(10, *) " ""RML_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RML_V impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,min_ven+23))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven)))/elem_field(ne_Qdot,24)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,min_ven+23))
      enddo
      write(10, *) " ""RML_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,1))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven)))/elem_field(ne_Qdot,1)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,1))
      enddo
      write(10, *) " ""MPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,2))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven)))/elem_field(ne_Qdot,1)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,2))
      enddo
      write(10, *) " ""MPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,3))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven)))/elem_field(ne_Qdot,1)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,3))
      enddo
      write(10, *) " ""MPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,4))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven)))/elem_field(ne_Qdot,1)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,4))
      enddo
      write(10, *) " ""MPA_A_4"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,5))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven+4)))/elem_field(ne_Qdot,5)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,5))
      enddo
      write(10, *) " ""LPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,6))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven+5)))/elem_field(ne_Qdot,6)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,6))
      enddo
      write(10, *) " ""LPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,7))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven+6)))/elem_field(ne_Qdot,7)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,7))
      enddo
      write(10, *) " ""LPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,8))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven+7)))/elem_field(ne_Qdot,8)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,8))
      enddo
      write(10, *) " ""RPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,9))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven+8)))/elem_field(ne_Qdot,9)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,9))
      enddo
      write(10, *) " ""RPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,10))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven+9)))/elem_field(ne_Qdot,10)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,10))
      enddo
      write(10, *) " ""RPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RBS_A impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,16))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven+15)))/elem_field(ne_Qdot,16)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,16))
      enddo
      write(10, *) " ""RBS_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RBS_V impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,min_ven+15))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven)))/elem_field(ne_Qdot,16)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,min_ven+15))
      enddo
      write(10, *) " ""RBS_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LBS_A impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,12))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven+11)))/elem_field(ne_Qdot,12)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,12))
      enddo
      write(10, *) " ""LBS_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LBS_V impedances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 10000*(node_field(nj_bv_press,elem_nodes(1,min_ven+11))-node_field&
      (nj_bv_press,elem_nodes(2,min_ven)))/elem_field(ne_Qdot,12)
      do i = 1, no_freq
        imped(i+1) = 10000.0/abs(eff_admit(i,min_ven+11))
      enddo
      write(10, *) " ""LBS_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"

      write(10, *) " ""unit"": ""dyne.s/cm5"""
      write(10, *) " },"
      write(10, *) " ""phase"":{"

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculating LUL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,11)),real(eff_admit(i,11), 8)) ! -1 is to make it impedance phase
      enddo
      write(10, *) " ""LUL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LUL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,min_ven+10)),real(eff_admit(i,min_ven+10), 8))
      enddo
      write(10, *) " ""LUL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LLL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,20)),real(eff_admit(i,20), 8)) ! -1 is to make it impedance phase
      enddo
      write(10, *) " ""LLL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LLL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,min_ven+19)),real(eff_admit(i,min_ven+19), 8))
      enddo
      write(10, *) " ""LLL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RUL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,15)),real(eff_admit(i,15), 8))
      enddo
      write(10, *) " ""RUL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RUL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,min_ven+14)),real(eff_admit(i,min_ven+14), 8))
      enddo
      write(10, *) " ""RUL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RLL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,23)),real(eff_admit(i,23), 8))
      enddo
      write(10, *) " ""RLL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RLL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,min_ven+22)),real(eff_admit(i,min_ven+22), 8))
      enddo
      write(10, *) " ""RLL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RML_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,24)),real(eff_admit(i,24), 8))
      enddo
      write(10, *) " ""RML_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RML_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,min_ven+23)),real(eff_admit(i,min_ven+23), 8))
      enddo
      write(10, *) " ""RML_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,1)),real(eff_admit(i,1), 8))
      enddo
      write(10, *) " ""MPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,2)),real(eff_admit(i,2), 8))
      enddo
      write(10, *) " ""MPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,3)),real(eff_admit(i,3), 8))
      enddo
      write(10, *) " ""MPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,4)),real(eff_admit(i,4), 8))
      enddo
      write(10, *) " ""MPA_A_4"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,5)),real(eff_admit(i,5), 8))
      enddo
      write(10, *) " ""LPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,6)),real(eff_admit(i,6), 8))
      enddo
      write(10, *) " ""LPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,7)),real(eff_admit(i,7), 8))
      enddo
      write(10, *) " ""LPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,8)),real(eff_admit(i,8), 8))
      enddo
      write(10, *) " ""RPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,9)),real(eff_admit(i,9), 8))
      enddo
      write(10, *) " ""RPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,10)),real(eff_admit(i,10), 8))
      enddo
      write(10, *) " ""RPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RBS_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,16)),real(eff_admit(i,16), 8))
      enddo
      write(10, *) " ""RBS_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RBS_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,min_ven+15)),real(eff_admit(i,min_ven+15), 8))
      enddo
      write(10, *) " ""RBS_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LBS_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,12)),real(eff_admit(i,12), 8))
      enddo
      write(10, *) " ""LBS_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LBS_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = -1*atan2(dimag(eff_admit(i,min_ven+11)),real(eff_admit(i,min_ven+11), 8))
      enddo
      write(10, *) " ""LBS_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      write(10, *) " ""unit"": ""radians"""
      write(10, *) " },"
      ! write(10, *) " ""Flow phase"":{"
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ! Calculating LUL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,11)),real(q_factor(i,11), 8))+b(i)
      ! enddo
      ! write(10, *) " ""LUL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating LUL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,min_ven+10)),real(q_factor(i,min_ven+10), 8))+b(i)
      ! enddo
      ! write(10, *) " ""LUL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating LLL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,20)),real(q_factor(i,20), 8))+b(i)
      ! enddo
      ! write(10, *) " ""LLL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating LLL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,min_ven+19)),real(q_factor(i,min_ven+19), 8))+b(i)
      ! enddo
      ! write(10, *) " ""LLL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RUL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,15)),real(q_factor(i,15), 8))+b(i)
      ! enddo
      ! write(10, *) " ""RUL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RUL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,min_ven+14)),real(q_factor(i,min_ven+14), 8))+b(i)
      ! enddo
      ! write(10, *) " ""RUL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RLL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,23)),real(q_factor(i,23), 8))+b(i)
      ! enddo
      ! write(10, *) " ""RLL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RLL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,min_ven+22)),real(q_factor(i,min_ven+22), 8))+b(i)
      ! enddo
      ! write(10, *) " ""RLL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RML_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,24)),real(q_factor(i,24), 8))+b(i)
      ! enddo
      ! write(10, *) " ""RML_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RML_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,min_ven+23)),real(q_factor(i,min_ven+23), 8))+b(i)
      ! enddo
      ! write(10, *) " ""RML_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,1)),real(q_factor(i,1), 8))+b(i)!-atan2(dimag(char_admit(i,1))&
      !   !,real(char_admit(i,1), 8))
      ! enddo
      ! write(10, *) " ""MPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,2)),real(q_factor(i,2), 8))+b(i)
      ! enddo
      ! write(10, *) " ""MPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,3)),real(q_factor(i,3), 8))+b(i)
      ! enddo
      ! write(10, *) " ""MPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,4)),real(q_factor(i,4), 8))+b(i)
      ! enddo
      ! write(10, *) " ""MPA_A_4"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating LPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,5)),real(q_factor(i,5), 8))+b(i)
      ! enddo
      ! write(10, *) " ""LPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating LPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,6)),real(q_factor(i,6), 8))+b(i)
      ! enddo
      ! write(10, *) " ""LPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating LPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,7)),real(q_factor(i,7), 8))+b(i)
      ! enddo
      ! write(10, *) " ""LPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,8)),real(q_factor(i,8), 8))+b(i)
      ! enddo
      ! write(10, *) " ""RPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,9)),real(q_factor(i,9), 8))+b(i)
      ! enddo
      ! write(10, *) " ""RPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,10)),real(q_factor(i,10), 8))+b(i)
      ! enddo
      ! write(10, *) " ""RPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RBS_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,16)),real(q_factor(i,16), 8))+b(i)
      ! enddo
      ! write(10, *) " ""RBS_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RBS_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,min_ven+15)),real(q_factor(i,min_ven+15), 8))+b(i)
      ! enddo
      ! write(10, *) " ""RBS_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating LBS_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,12)),real(q_factor(i,12), 8))+b(i)
      ! enddo
      ! write(10, *) " ""LBS_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating LBS_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = 0
      ! do i = 1, no_freq
      !   imped(i+1) = atan2(dimag(q_factor(i,min_ven+11)),real(q_factor(i,min_ven+11), 8))+b(i)
      ! enddo
      ! write(10, *) " ""LBS_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! write(10, *) " ""unit"": ""radians"""
      ! write(10, *) " },"
      write(10, *) " ""Pressure phase"":{"

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculating LUL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,11)),real(p_factor(i,11), 8))+b(i) ! -1 is to make it impedance phase
      enddo
      write(10, *) " ""LUL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LUL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,min_ven+10)),real(p_factor(i,min_ven+10), 8))+b(i)
      enddo
      write(10, *) " ""LUL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LLL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,20)),real(p_factor(i,20), 8))+b(i) ! -1 is to make it impedance phase
      enddo
      write(10, *) " ""LLL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LLL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,min_ven+19)),real(p_factor(i,min_ven+19), 8))+b(i)
      enddo
      write(10, *) " ""LLL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RUL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,15)),real(p_factor(i,15), 8))+b(i)
      enddo
      write(10, *) " ""RUL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RUL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,min_ven+14)),real(p_factor(i,min_ven+14), 8))+b(i)
      enddo
      write(10, *) " ""RUL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RLL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,23)),real(p_factor(i,23), 8))+b(i)
      enddo
      write(10, *) " ""RLL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RLL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,min_ven+22)),real(p_factor(i,min_ven+22), 8))+b(i)
      enddo
      write(10, *) " ""RLL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RML_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,24)),real(p_factor(i,24), 8))+b(i)
      enddo
      write(10, *) " ""RML_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RML_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,min_ven+23)),real(p_factor(i,min_ven+23), 8))+b(i)
      enddo
      write(10, *) " ""RML_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,1)),real(p_factor(i,1), 8))+b(i)
      enddo
      write(10, *) " ""MPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,2)),real(p_factor(i,2), 8))+b(i)
      enddo
      write(10, *) " ""MPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,3)),real(p_factor(i,3), 8))+b(i)
      enddo
      write(10, *) " ""MPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,4)),real(p_factor(i,4), 8))+b(i)
      enddo
      write(10, *) " ""MPA_A_4"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,5)),real(p_factor(i,5), 8))+b(i)
      enddo
      write(10, *) " ""LPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,6)),real(p_factor(i,6), 8))+b(i)
      enddo
      write(10, *) " ""LPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,7)),real(p_factor(i,7), 8))+b(i)
      enddo
      write(10, *) " ""LPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,8)),real(p_factor(i,8), 8))+b(i)
      enddo
      write(10, *) " ""RPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,9)),real(p_factor(i,9), 8))+b(i)
      enddo
      write(10, *) " ""RPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,10)),real(p_factor(i,10), 8))+b(i)
      enddo
      write(10, *) " ""RPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RBS_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,16)),real(p_factor(i,16), 8))+b(i)
      enddo
      write(10, *) " ""RBS_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RBS_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,min_ven+15)),real(p_factor(i,min_ven+15), 8))+b(i)
      enddo
      write(10, *) " ""RBS_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LBS_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,12)),real(p_factor(i,12), 8))+b(i)
      enddo
      write(10, *) " ""LBS_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LBS_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = 0
      do i = 1, no_freq
        imped(i+1) = atan2(dimag(p_factor(i,min_ven+11)),real(p_factor(i,min_ven+11), 8))+b(i)
      enddo
      write(10, *) " ""LBS_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      write(10, *) " ""unit"": ""radians"""
      write(10, *) " },"
      ! write(10, *) " ""Flow amplitude"":{"
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ! Calculating LUL_A Flow ampl !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,11)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,11))*a(i)!*abs(char_admit(i,11))*a(i)
      ! enddo
      ! write(10, *) " ""LUL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating LUL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,min_ven+10)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,min_ven+10))*a(i)!*abs(char_admit(i,min_ven+10))*a(i)
      ! enddo
      ! write(10, *) " ""LUL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating LLL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,20)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,20))*a(i)!*abs(char_admit(i,20))*a(i)
      ! enddo
      ! write(10, *) " ""LLL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating LLL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,min_ven+19)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,min_ven+19))*a(i)!*abs(char_admit(i,min_ven+19))*a(i)
      ! enddo
      ! write(10, *) " ""LLL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RUL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,15)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,15))*a(i)!*abs(char_admit(i,15))*a(i)
      ! enddo
      ! write(10, *) " ""RUL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RUL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,min_ven+14)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,min_ven+14))*a(i)!*abs(char_admit(i,min_ven+14))*a(i)
      ! enddo
      ! write(10, *) " ""RUL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RLL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,23)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,23))*a(i)!*abs(char_admit(i,23))*a(i)
      ! enddo
      ! write(10, *) " ""RLL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RLL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,min_ven+22)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,min_ven+22))*a(i)!*abs(char_admit(i,min_ven+22))*a(i)
      ! enddo
      ! write(10, *) " ""RLL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RML_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,24)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,24))*a(i)!*abs(char_admit(i,24))*a(i)
      ! enddo
      ! write(10, *) " ""RML_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RML_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,min_ven+23)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,min_ven+23))*a(i)!*abs(char_admit(i,min_ven+23))*a(i)
      ! enddo
      ! write(10, *) " ""RML_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,1)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i, 1))*a(i)!*abs(char_admit(i,1))*a(i)
      ! enddo
      ! write(10, *) " ""MPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,2)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,2))*a(i)!*abs(char_admit(i,2))*a(i)
      ! enddo
      ! write(10, *) " ""MPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,3)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,3))*a(i)!*abs(char_admit(i,3))*a(i)
      ! enddo
      ! write(10, *) " ""MPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,4)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,4))*a(i)!*abs(char_admit(i,4))*a(i)
      ! enddo
      ! write(10, *) " ""MPA_A_4"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating LPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,5)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,5))*a(i)!*abs(char_admit(i,5))*a(i)
      ! enddo
      ! write(10, *) " ""LPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating LPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,6)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,6))*a(i)!*abs(char_admit(i,6))*a(i)
      ! enddo
      ! write(10, *) " ""LPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating LPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,7)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,7))*a(i)!*abs(char_admit(i,7))*a(i)
      ! enddo
      ! write(10, *) " ""LPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,8)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,8))*a(i)!*abs(char_admit(i,8))*a(i)
      ! enddo
      ! write(10, *) " ""RPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,9)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,9))*a(i)!*abs(char_admit(i,9))*a(i)
      ! enddo
      ! write(10, *) " ""RPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,10)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,10))*a(i)!*abs(char_admit(i,10))*a(i)
      ! enddo
      ! write(10, *) " ""RPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RBS_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,16)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,16))*a(i)!*abs(char_admit(i,16))*a(i)
      ! enddo
      ! write(10, *) " ""RBS_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating RBS_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,min_ven+15)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,min_ven+15))*a(i)!*abs(char_admit(i,min_ven+15))*a(i)
      ! enddo
      ! write(10, *) " ""RBS_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating LBS_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,12)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,12))*a(i)!*abs(char_admit(i,12))*a(i)
      ! enddo
      ! write(10, *) " ""LBS_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !Calculating LBS_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! imped(1) = elem_field(ne_Qdot,min_ven+11)
      ! do i = 1, no_freq
      !   imped(i+1) = abs(q_factor(i,min_ven+11))*a(i)!*abs(char_admit(i,min_ven+11))*a(i)
      ! enddo
      ! write(10, *) " ""LBS_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      ! write(10, *) " ""unit"": ""mm3/s"""
      ! write(10, *) " },"
      write(10, *) " ""Pressure amplitude"":{"

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculating LUL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,11))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,11))*a(i) ! -1 is to make it impedance phase
      enddo
      write(10, *) " ""LUL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LUL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,min_ven+10))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,min_ven+10))*a(i)
      enddo
      write(10, *) " ""LUL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LLL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,20))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,20))*a(i)
      enddo
      write(10, *) " ""LLL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LLL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,min_ven+19))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,min_ven+19))*a(i)
      enddo
      write(10, *) " ""LLL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RUL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,15))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,15))*a(i)
      enddo
      write(10, *) " ""RUL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RUL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,min_ven+14))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,min_ven+14))*a(i)
      enddo
      write(10, *) " ""RUL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RLL_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,23))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,23))*a(i)
      enddo
      write(10, *) " ""RLL_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RLL_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,min_ven+22))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,min_ven+22))*a(i)
      enddo
      write(10, *) " ""RLL_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RML_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,24))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,24))*a(i)
      enddo
      write(10, *) " ""RML_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RML_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,min_ven+23))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,min_ven+23))*a(i)
      enddo
      write(10, *) " ""RML_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,1))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,1))*a(i)
      enddo
      write(10, *) " ""MPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,2))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,2))*a(i)
      enddo
      write(10, *) " ""MPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,3))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,3))*a(i)
      enddo
      write(10, *) " ""MPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,4))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,4))*a(i)
      enddo
      write(10, *) " ""MPA_A_4"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,5))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,5))*a(i)
      enddo
      write(10, *) " ""LPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,6))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,6))*a(i)
      enddo
      write(10, *) " ""LPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,7))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,7))*a(i)
      enddo
      write(10, *) " ""LPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,8))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,8))*a(i)
      enddo
      write(10, *) " ""RPA_A_1"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,9))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,9))*a(i)
      enddo
      write(10, *) " ""RPA_A_2"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,10))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,10))*a(i)
      enddo
      write(10, *) " ""RPA_A_3"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RBS_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,16))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,16))*a(i)
      enddo
      write(10, *) " ""RBS_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RBS_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,min_ven+15))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,min_ven+15))*a(i)
      enddo
      write(10, *) " ""RBS_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LBS_A phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,12))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,12))*a(i)
      enddo
      write(10, *) " ""LBS_A"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LBS_V phase !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      imped(1) = node_field(nj_bv_press,elem_nodes(1,min_ven+11))
      do i = 1, no_freq
        imped(i+1) = abs(p_factor(i,min_ven+11))*a(i)
      enddo
      write(10, *) " ""LBS_V"": [", imped(1), ",", (imped(i),",",i=2,no_freq), imped(no_freq+1), "],"
      write(10, *) " ""unit"": ""Pa"""
      write(10, *) " },"
      write(10, *) " ""radius"":{"

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculating LUL_A radius  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LUL_A"": [", elem_field(ne_radius_out,11), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LUL_V radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LUL_V"": [", elem_field(ne_radius_out,min_ven+10), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LLL_A radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LLL_A"": [", elem_field(ne_radius_out,20), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LLL_V radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LLL_V"": [", elem_field(ne_radius_out,min_ven+19), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RUL_A radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RUL_A"": [", elem_field(ne_radius_out,15), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RUL_V radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RUL_V"": [", elem_field(ne_radius_out,min_ven+14), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RLL_A radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RLL_A"": [", elem_field(ne_radius_out,23), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RLL_V radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RLL_V"": [", elem_field(ne_radius_out,min_ven+22), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RML_A radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RML_A"": [", elem_field(ne_radius_out,24), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RML_V radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RML_V"": [", elem_field(ne_radius_out,min_ven+23), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""MPA_A_1"": [", elem_field(ne_radius_out,1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""MPA_A_2"": [", elem_field(ne_radius_out,2), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""MPA_A_3"": [", elem_field(ne_radius_out,3), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""MPA_A_4"": [", elem_field(ne_radius_out,4), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LPA_A_1"": [", elem_field(ne_radius_out,5), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LPA_A_2"": [", elem_field(ne_radius_out,6), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LPA_A_3"": [", elem_field(ne_radius_out,7), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RPA_A_1"": [", elem_field(ne_radius_out,8), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RPA_A_2"": [", elem_field(ne_radius_out,9), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RPA_A_3"": [", elem_field(ne_radius_out,10), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RBS_A radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RBS_A"": [", elem_field(ne_radius_out,16), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RBS_V radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RBS_V"": [", elem_field(ne_radius_out,min_ven+15), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LBS_A radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LBS_A"": [", elem_field(ne_radius_out,12), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LBS_V radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LBS_V"": [", elem_field(ne_radius_out,min_ven+11), "],"

      write(10, *) " ""unit"": ""mm"""
      write(10, *) " },"
      write(10, *) " ""unstrained radius"":{"

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculating LUL_A unstrained radius  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LUL_A"": [", elem_field(ne_radius_out0,11), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LUL_V unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LUL_V"": [", elem_field(ne_radius_out0,min_ven+10), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LLL_A unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LLL_A"": [", elem_field(ne_radius_out0,20), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LLL_V unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LLL_V"": [", elem_field(ne_radius_out0,min_ven+19), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RUL_A unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RUL_A"": [", elem_field(ne_radius_out0,15), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RUL_V unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RUL_V"": [", elem_field(ne_radius_out0,min_ven+14), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RLL_A unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RLL_A"": [", elem_field(ne_radius_out0,23), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RLL_V unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RLL_V"": [", elem_field(ne_radius_out0,min_ven+22), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RML_A unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RML_A"": [", elem_field(ne_radius_out0,24), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RML_V unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RML_V"": [", elem_field(ne_radius_out0,min_ven+23), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""MPA_A_1"": [", elem_field(ne_radius_out0,1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""MPA_A_2"": [", elem_field(ne_radius_out0,2), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""MPA_A_3"": [", elem_field(ne_radius_out0,3), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""MPA_A_4"": [", elem_field(ne_radius_out0,4), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LPA_A_1"": [", elem_field(ne_radius_out0,5), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LPA_A_2"": [", elem_field(ne_radius_out0,6), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LPA_A_3"": [", elem_field(ne_radius_out0,7), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RPA_A_1"": [", elem_field(ne_radius_out0,8), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RPA_A_2"": [", elem_field(ne_radius_out0,9), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RPA_A_3"": [", elem_field(ne_radius_out0,10), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RBS_A unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RBS_A"": [", elem_field(ne_radius_out0,16), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RBS_V unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RBS_V"": [", elem_field(ne_radius_out0,min_ven+15), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LBS_A unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LBS_A"": [", elem_field(ne_radius_out0,12), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LBS_V unstrained radius !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LBS_V"": [", elem_field(ne_radius_out0,min_ven+11), "],"

      write(10, *) " ""unit"": ""mm"""
      write(10, *) " },"
      write(10, *) " ""mean flow"":{"

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculating LUL_A mean flow  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LUL_A"": [", elem_field(ne_Qdot,11), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LUL_V mean flow  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LUL_V"": [", elem_field(ne_Qdot,min_ven+10), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LLL_A mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LLL_A"": [", elem_field(ne_Qdot,20), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LLL_V mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LLL_V"": [", elem_field(ne_Qdot,min_ven+19), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RUL_A mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RUL_A"": [", elem_field(ne_Qdot,15), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RUL_V mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RUL_V"": [", elem_field(ne_Qdot,min_ven+14), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RLL_A mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RLL_A"": [", elem_field(ne_Qdot,23), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RLL_V mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RLL_V"": [", elem_field(ne_Qdot,min_ven+22), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RML_A mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RML_A"": [", elem_field(ne_Qdot,24), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RML_V mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RML_V"": [", elem_field(ne_Qdot,min_ven+23), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""MPA_A_1"": [", elem_field(ne_Qdot,1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""MPA_A_2"": [", elem_field(ne_Qdot,2), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""MPA_A_3"": [", elem_field(ne_Qdot,3), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""MPA_A_4"": [", elem_field(ne_Qdot,4), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LPA_A_1"": [", elem_field(ne_Qdot,5), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LPA_A_2"": [", elem_field(ne_Qdot,6), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LPA_A_3"": [", elem_field(ne_Qdot,7), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RPA_A_1"": [", elem_field(ne_Qdot,8), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RPA_A_2"": [", elem_field(ne_Qdot,9), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RPA_A_3"": [", elem_field(ne_Qdot,10), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RBS_A mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RBS_A"": [", elem_field(ne_Qdot,16), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RBS_V mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RBS_V"": [", elem_field(ne_Qdot,min_ven+15), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LBS_A mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LBS_A"": [", elem_field(ne_Qdot,12), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LBS_V mean flow !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LBS_V"": [", elem_field(ne_Qdot,min_ven+11), "],"

      write(10, *) " ""unit"": ""mm^3/s"""
      write(10, *) " },"
      write(10, *) " ""Length"":{"

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculating LUL_A Length  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LUL_A"": [", elem_field(ne_length,11), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LUL_V Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LUL_V"": [", elem_field(ne_length,min_ven+10), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LLL_A Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LLL_A"": [", elem_field(ne_length,20), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LLL_V Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LLL_V"": [", elem_field(ne_length,min_ven+19), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RUL_A Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RUL_A"": [", elem_field(ne_length,15), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RUL_V Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RUL_V"": [", elem_field(ne_length,min_ven+14), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RLL_A Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RLL_A"": [", elem_field(ne_length,23), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RLL_V Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RLL_V"": [", elem_field(ne_length,min_ven+22), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RML_A Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RML_A"": [", elem_field(ne_length,24), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RML_V Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RML_V"": [", elem_field(ne_length,min_ven+23), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""MPA_A_1"": [", elem_field(ne_length,1), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""MPA_A_2"": [", elem_field(ne_length,2), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""MPA_A_3"": [", elem_field(ne_length,3), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating MPA Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""MPA_A_4"": [", elem_field(ne_length,4), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LPA_A_1"": [", elem_field(ne_length,5), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LPA_A_2"": [", elem_field(ne_length,6), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LPA Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LPA_A_3"": [", elem_field(ne_length,7), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RPA_A_1"": [", elem_field(ne_length,8), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RPA_A_2"": [", elem_field(ne_length,9), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RPA Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RPA_A_3"": [", elem_field(ne_length,10), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RBS_A Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RBS_A"": [", elem_field(ne_length,16), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating RBS_V Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""RBS_V"": [", elem_field(ne_length,min_ven+15), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LBS_A Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LBS_A"": [", elem_field(ne_length,12)+elem_field(ne_length,13)+elem_field(ne_length,14), "],"
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Calculating LBS_V Length !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(10, *) " ""LBS_V"": [", elem_field(ne_length,min_ven+11)+elem_field(ne_length,min_ven+12)&
      +elem_field(ne_length,min_ven+13), "],"

      write(10, *) " ""unit"": ""mm"""
      write(10, *) " }"
      write(10, *) "}"
      close(10)
    endif ! export lobe admittance

    open(fid5, file = 'inputadmittance.txt',action='write')
    write(fid5,fmt=*) 'input admittance:'
    do nf=1,no_freq
        omega=nf*harmonic_scale
        write(fid5,fmt=*) omega,abs(eff_admit(nf,1)),&
            atan2(dimag(eff_admit(nf,1)),real(eff_admit(nf,1), 8))
    enddo
    close(fid5)

    start_time=0.0_dp
    end_time=60.0_dp/heartrate
    dt=(end_time-start_time)/n_time
    time=start_time
    !consider first pressure and flow into the vessel (at x=0)
    open(fid, file = 'incident_pressure.txt', action='write')
    open(fid2, file = 'incident_flow.txt',action='write')
    open(fid3, file = 'terminal_pressure.txt',action='write')
    open(fid4, file = 'terminal_flow.txt',action='write')
    open(fid6, file = 'terminal_radii.txt', action='write')
    open(fid7, file = 'terminal_WSS.txt',action='write')
    open(200, file = 'WSS_arterial.txt', action='write')
    open(210, file = 'radii.txt', action='write')
    open(220, file='flows.txt', status='replace', action='write')
    open(230, file='pressures.txt', status='replace', action='write')
    call write_time_dependent_perfusion_header(220,'flows','mm^3/s', &
      start_time,dt,n_time)
    call write_time_dependent_perfusion_header(230,'pressures','Pa', &
      start_time,dt,n_time)
    do nu =1,num_units
        ne=units(nu)  ! terminal elements
        ne_previous=elem_cnct(-1,1,ne)
        ! initialisation
        forward_pressure=0.0_dp
        reflected_pressure=0.0_dp
        forward_pressure_previous=0.0_dp
        reflected_pressure_previous=0.0_dp
        forward_flow=0.0_dp
        reflected_flow=0.0_dp
        terminals_radius=0.0_dp
        WSS=0.0_dp
        do nt=1,n_time
            do nf=1,no_freq
                omega=2*pi*nf*harmonic_scale
                forward_pressure(nt)=forward_pressure(nt)+abs(p_factor(nf,ne))*a(nf)*cos(omega*time+b(nf)+&
                    atan2(dimag(p_factor(nf,ne)),real(p_factor(nf,ne), 8)))
                forward_pressure_previous(nt)=forward_pressure_previous(nt)+abs(p_factor(nf,ne_previous))*&
                a(nf)*cos(omega*time+b(nf)+atan2(dimag(p_factor(nf,ne_previous)),real(p_factor(nf,ne_previous), 8)))

                reflected_pressure(nt)=reflected_pressure(nt)+abs(p_factor(nf,ne))*a(nf)*&
                    abs(reflect(nf,ne))*exp((-2*elem_field(ne_length,ne))*(real(prop_const(nf,ne), 8)))*&
                    cos(omega*time+b(nf)+&
                    atan2(dimag(p_factor(nf,ne)),real(p_factor(nf,ne), 8))+&
                    (-2*elem_field(ne_length,ne))*(dimag(prop_const(nf,ne)))+&
                    atan2(dimag(reflect(nf,ne)),real(reflect(nf,ne), 8)))
                reflected_pressure_previous(nt)=reflected_pressure_previous(nt)+abs(p_factor(nf,ne_previous))*&
                a(nf)*abs(reflect(nf,ne_previous))*exp((-2*elem_field(ne_length,ne_previous))*&
                (real(prop_const(nf,ne_previous), 8)))*cos(omega*time+b(nf)+&
                    atan2(dimag(p_factor(nf,ne_previous)),real(p_factor(nf,ne_previous), 8))+&
                    (-2*elem_field(ne_length,ne_previous))*(dimag(prop_const(nf,ne_previous)))+&
                    atan2(dimag(reflect(nf,ne_previous)),real(reflect(nf,ne_previous), 8)))

                forward_flow(nt)=forward_flow(nt)+abs(char_admit(nf,ne))*abs(p_factor(nf,ne))*a(nf)*&
                    cos(omega*time+b(nf)+&
                    atan2(dimag(p_factor(nf,ne)),real(p_factor(nf,ne), 8))+&
                    atan2(dimag(char_admit(nf,ne)),real(char_admit(nf,ne), 8)))

                reflected_flow(nt)=reflected_flow(nt)+abs(char_admit(nf,ne))*abs(p_factor(nf,ne))*a(nf)*&
                    abs(reflect(nf,ne))*exp((-2*elem_field(ne_length,ne))*(real(prop_const(nf,ne), 8)))*&
                    cos(omega*time+b(nf)+&
                    atan2(dimag(p_factor(nf,ne)),real(p_factor(nf,ne), 8))+&
                    (-2*elem_field(ne_length,ne))*(dimag(prop_const(nf,ne)))+&
                    atan2(dimag(reflect(nf,ne)),real(reflect(nf,ne), 8))+&
                    atan2(dimag(char_admit(nf,ne)),real(char_admit(nf,ne), 8)))

            enddo
            time=time+dt
        enddo

        np=elem_nodes(2,ne) ! terminals
        np_previous=elem_nodes(1,ne) ! upstream node to terminal nodes
        if(.not.allocated(p_terminal)) allocate (p_terminal(n_time))
        if(.not.allocated(p_previous)) allocate (p_previous(n_time))
        if(.not.allocated(terminal_flow)) allocate (terminal_flow(n_time))
        terminal_flow=0.0_dp
        p_terminal = forward_pressure+reflected_pressure + node_field(nj_bv_press,np)
        p_previous = forward_pressure_previous+reflected_pressure_previous +&
        node_field(nj_bv_press,np_previous)
        terminal_flow = forward_flow - reflected_flow + elem_field(ne_Qdot,ne)
        call calculate_ppl(np,grav_vect,mechanics_parameters,Ppl)
        do nt=1,n_time
          delta_p = ((p_terminal(nt) + p_previous(nt))/2) - Ppl
          terminals_radius(nt) = elem_field(ne_radius_out0,ne)*(1+alpha(ne)*delta_p)
          WSS(nt) = (4.0_dp * viscosity * terminal_flow(nt))/(pi * terminals_radius(nt)**3) !wall shear stress at terminals
        enddo
        write(fid,fmt=*) ne, forward_pressure+node_field(nj_bv_press,np) ! incident pressure
        write(fid2,fmt=*) ne, forward_flow+elem_field(ne_Qdot,ne) ! incident flow
        write(fid3,fmt=*) ne, forward_pressure+reflected_pressure + node_field(nj_bv_press,np) !terminal total pressure
        write(fid4,fmt=*) ne, forward_flow-reflected_flow + elem_field(ne_Qdot,ne) !terminal total flow
        write(fid6,fmt=*) ne, terminals_radius ! terminal radii
        write(fid7,fmt=*) ne, WSS ! terminal elements wall shear stress

    enddo

    !!! Doing the same for all frequencies for MPA
    !!! Export MPA flow components (flow/pressure)
    open(fid8, file = 'Inlet_pressure.txt', action='write')
    open(fid9, file = 'Inlet_flow.txt',action='write')
    open(fid10, file = 'Inlet_forward_flow.txt',action='write')
    open(fid11, file = 'Inlet_reflected_flow.txt',action='write')
    open(fid12, file = 'Inlet_forward_pressure.txt',action='write')
    open(fid13, file = 'Inlet_reflected_pressure.txt',action='write')
    open(140, file = 'Outlet_pressure.txt',action='write')
    open(150, file = 'Outlet_flow.txt', action='write')
    open(160, file = 'Outlet_forward_pressure.txt',action='write')
    open(170, file = 'Outlet_reflected_pressure.txt',action='write')
    open(180, file = 'Outlet_forward_flow.txt',action='write')
    open(190, file = 'Outlet_reflected_flow.txt',action='write')
    ne = 1 ! MPA inlet element
    forward_pressure=0.0_dp ! reseting the for MPA
    reflected_pressure=0.0_dp ! reseting the for MPA
    forward_flow=0.0_dp ! reseting the for MPA
    reflected_flow=0.0_dp ! reseting the for MPA
    elem_rad=0.0_dp
    WSS=0.0_dp
    do nt=1,n_time
        do nf=1,no_freq
            omega=2*pi*nf*harmonic_scale
            forward_pressure(nt)=forward_pressure(nt)+abs(p_factor(nf,ne))*a(nf)*cos(omega*time+b(nf)+&
                atan2(dimag(p_factor(nf,ne)),real(p_factor(nf,ne), 8)))

            reflected_pressure(nt)=reflected_pressure(nt)+abs(p_factor(nf,ne))*a(nf)*&
                abs(reflect(nf,ne))*exp((-2*elem_field(ne_length,ne))*(real(prop_const(nf,ne), 8)))*&
                cos(omega*time+b(nf)+&
                atan2(dimag(p_factor(nf,ne)),real(p_factor(nf,ne), 8))+&
                (-2*elem_field(ne_length,ne))*(dimag(prop_const(nf,ne)))+&
                atan2(dimag(reflect(nf,ne)),real(reflect(nf,ne), 8)))

            forward_flow(nt)=forward_flow(nt)+abs(char_admit(nf,ne))*abs(p_factor(nf,ne))*a(nf)*&
                cos(omega*time+b(nf)+&
                atan2(dimag(p_factor(nf,ne)),real(p_factor(nf,ne), 8))+&
                atan2(dimag(char_admit(nf,ne)),real(char_admit(nf,ne), 8)))

            reflected_flow(nt)=reflected_flow(nt)+abs(char_admit(nf,ne))*abs(p_factor(nf,ne))*a(nf)*&
                abs(reflect(nf,ne))*exp((-2*elem_field(ne_length,ne))*(real(prop_const(nf,ne), 8)))*&
                cos(omega*time+b(nf)+&
                atan2(dimag(p_factor(nf,ne)),real(p_factor(nf,ne), 8))+&
                (-2*elem_field(ne_length,ne))*(dimag(prop_const(nf,ne)))+&
                atan2(dimag(reflect(nf,ne)),real(reflect(nf,ne), 8))+&
                atan2(dimag(char_admit(nf,ne)),real(char_admit(nf,ne), 8)))

        enddo
        time=time+dt
    enddo

    np=elem_nodes(2,ne) ! end node of an elem
    np_previous=elem_nodes(1,ne) ! upstream node of the element
    if(.not.allocated(p_terminal)) allocate (p_terminal(n_time))
    if(.not.allocated(p_previous)) allocate (p_previous(n_time))
    if(.not.allocated(elem_flow)) allocate (elem_flow(n_time))
    if(.not.allocated(elem_rad)) allocate (elem_rad(n_time))
    p_terminal = forward_pressure + reflected_pressure + node_field(nj_bv_press,np)
    p_previous = node_field(nj_bv_press,np_previous)
    elem_flow = forward_flow - reflected_flow + elem_field(ne_Qdot,ne)
    call calculate_ppl(np,grav_vect,mechanics_parameters,Ppl)

    do nt=1,n_time
      delta_p = ((p_terminal(nt) + p_previous(nt))/2)-Ppl
      elem_rad(nt) = elem_field(ne_radius_out0,ne)*(1+alpha(ne)*delta_p)
      WSS(nt) = (4.0_dp * viscosity * elem_flow(nt))/(pi * elem_rad(nt)**3) !wall shear stress at terminals
    enddo

    write(fid8,fmt=*) forward_pressure+reflected_pressure + node_field(nj_bv_press,np_previous) !Inlet total pressure
    write(fid9,fmt=*) forward_flow - reflected_flow + elem_field(ne_Qdot,ne) !Inlet MPA flow
    write(fid10,fmt=*) forward_flow !Inlet forward flow
    write(fid11,fmt=*) reflected_flow !Inlet reflected flow
    write(fid12,fmt=*) forward_pressure !Inlet forward pressure
    write(fid13,fmt=*) reflected_pressure !Inlet reflected pressure
    write(200,fmt=*) ne, WSS ! WSS for element 1 to be added to WSS_arterial.txt
    write(210,fmt=*) ne, elem_rad ! radius for element 1 to be added to radii.txt
    write(220,fmt=*) ne, forward_flow - reflected_flow + elem_field(ne_Qdot,ne) ! Flow for element 1 to be added to flows.txt
    write(230,fmt=*) ne, p_terminal ! Pressure for element 1 to be added to pressures.txt


    do ne =2,min_ven-1  ! reporting time-dependent radius and WSS for all the arteries
      ne_previous=elem_cnct(-1,1,ne)
      forward_pressure=0.0_dp
      reflected_pressure=0.0_dp
      forward_pressure_previous=0.0_dp
      reflected_pressure_previous=0.0_dp
      forward_flow=0.0_dp
      reflected_flow=0.0_dp
      WSS=0.0_dp
      elem_rad=0.0_dp
      elem_flow=0.0_dp
      do nt=1,n_time
          do nf=1,no_freq
              omega=2*pi*nf*harmonic_scale
              forward_pressure(nt)=forward_pressure(nt)+abs(p_factor(nf,ne))*a(nf)*cos(omega*time+b(nf)+&
                  atan2(dimag(p_factor(nf,ne)),real(p_factor(nf,ne), 8)))
              forward_pressure_previous(nt)=forward_pressure_previous(nt)+abs(p_factor(nf,ne_previous))*&
              a(nf)*cos(omega*time+b(nf)+atan2(dimag(p_factor(nf,ne_previous)),real(p_factor(nf,ne_previous), 8)))

              reflected_pressure(nt)=reflected_pressure(nt)+abs(p_factor(nf,ne))*a(nf)*&
                  abs(reflect(nf,ne))*exp((-2*elem_field(ne_length,ne))*(real(prop_const(nf,ne), 8)))*&
                  cos(omega*time+b(nf)+&
                  atan2(dimag(p_factor(nf,ne)),real(p_factor(nf,ne), 8))+&
                  (-2*elem_field(ne_length,ne))*(dimag(prop_const(nf,ne)))+&
                  atan2(dimag(reflect(nf,ne)),real(reflect(nf,ne), 8)))
              reflected_pressure_previous(nt)=reflected_pressure_previous(nt)+abs(p_factor(nf,ne_previous))*&
              a(nf)*abs(reflect(nf,ne_previous))*exp((-2*elem_field(ne_length,ne_previous))*&
              (real(prop_const(nf,ne_previous), 8)))*cos(omega*time+b(nf)+&
                  atan2(dimag(p_factor(nf,ne_previous)),real(p_factor(nf,ne_previous), 8))+&
                  (-2*elem_field(ne_length,ne_previous))*(dimag(prop_const(nf,ne_previous)))+&
                  atan2(dimag(reflect(nf,ne_previous)),real(reflect(nf,ne_previous), 8)))

              forward_flow(nt)=forward_flow(nt)+abs(char_admit(nf,ne))*abs(p_factor(nf,ne))*a(nf)*&
                  cos(omega*time+b(nf)+&
                  atan2(dimag(p_factor(nf,ne)),real(p_factor(nf,ne), 8))+&
                  atan2(dimag(char_admit(nf,ne)),real(char_admit(nf,ne), 8)))

              reflected_flow(nt)=reflected_flow(nt)+abs(char_admit(nf,ne))*abs(p_factor(nf,ne))*a(nf)*&
                  abs(reflect(nf,ne))*exp((-2*elem_field(ne_length,ne))*(real(prop_const(nf,ne), 8)))*&
                  cos(omega*time+b(nf)+&
                  atan2(dimag(p_factor(nf,ne)),real(p_factor(nf,ne), 8))+&
                  (-2*elem_field(ne_length,ne))*(dimag(prop_const(nf,ne)))+&
                  atan2(dimag(reflect(nf,ne)),real(reflect(nf,ne), 8))+&
                  atan2(dimag(char_admit(nf,ne)),real(char_admit(nf,ne), 8)))

          enddo
          time=time+dt
      enddo

      np=elem_nodes(2,ne) ! end node of an elem
      np_previous=elem_nodes(1,ne) ! upstream node of the element
      if(.not.allocated(p_terminal)) allocate (p_terminal(n_time))
      if(.not.allocated(p_previous)) allocate (p_previous(n_time))
      if(.not.allocated(elem_flow)) allocate (elem_flow(n_time))
      if(.not.allocated(elem_rad)) allocate (elem_rad(n_time))
      p_terminal = forward_pressure + reflected_pressure + node_field(nj_bv_press,np)
      p_previous = forward_pressure_previous + reflected_pressure_previous +&
      node_field(nj_bv_press,np_previous)
      elem_flow = forward_flow - reflected_flow + elem_field(ne_Qdot,ne) ! forwards wave - reflected + S.S flow = total elem flow
      call calculate_ppl(np,grav_vect,mechanics_parameters,Ppl)
      do nt=1,n_time
        delta_p = ((p_terminal(nt) + p_previous(nt))/2) - Ppl
        elem_rad(nt) = elem_field(ne_radius_out0,ne)*(1+alpha(ne)*delta_p)
        WSS(nt) = (4.0_dp * viscosity * elem_flow(nt))/(pi * elem_rad(nt)**3) !wall shear stress at terminals
      enddo
      write(200,fmt=*) ne, WSS ! elements wall shear stress in cardiac cycle
      write(210,fmt=*) ne, elem_rad ! elements radii in cardiac cycle
      write(220,fmt=*) ne, forward_flow-reflected_flow+elem_field(ne_Qdot,ne) ! elem flows in cardiac cycle
      write(230,fmt=*) ne, p_terminal ! element pressures in cardiac cycle
    enddo

    do ne = 1, num_elems
      if((elem_field(ne_group,ne).eq.2.0_dp).and.(vein_found.eqv..False.))  then !only applying on veins
        vein_found = .True.
        vein_elem = ne
      endif
    end do

!! Extracting the waveforms at the venous outlet !!!!!!!!!!!!!!!!!!!!

    ne = vein_elem  !!!
    forward_pressure=0.0_dp ! resetting the array for outlet
    reflected_pressure=0.0_dp ! resetting the array for outlet
    forward_flow=0.0_dp ! resetting the array for outlet
    reflected_flow=0.0_dp ! resetting the array for outlet
    if (bc_type.eq.'pressure')then
      do nt=1,n_time
          do nf=1,no_freq
              omega=2*pi*nf*harmonic_scale
              forward_pressure(nt)=forward_pressure(nt)+abs(p_factor(nf,ne))*a(nf)*cos(omega*time+b(nf)+&
                  atan2(dimag(p_factor(nf,ne)),real(p_factor(nf,ne), 8)))

              reflected_pressure(nt)=reflected_pressure(nt)+abs(p_factor(nf,ne))*a(nf)*&
                  abs(reflect(nf,ne))*exp((-2*elem_field(ne_length,ne))*(real(prop_const(nf,ne), 8)))*&
                  cos(omega*time+b(nf)+&
                  atan2(dimag(p_factor(nf,ne)),real(p_factor(nf,ne), 8))+&
                  (-2*elem_field(ne_length,ne))*(dimag(prop_const(nf,ne)))+&
                  atan2(dimag(reflect(nf,ne)),real(reflect(nf,ne), 8)))

              forward_flow(nt)=forward_flow(nt)+abs(char_admit(nf,ne))*abs(p_factor(nf,ne))*a(nf)*&
                  cos(omega*time+b(nf)+&
                  atan2(dimag(p_factor(nf,ne)),real(p_factor(nf,ne), 8))+&
                  atan2(dimag(char_admit(nf,ne)),real(char_admit(nf,ne), 8)))

              reflected_flow(nt)=reflected_flow(nt)+abs(char_admit(nf,ne))*abs(p_factor(nf,ne))*a(nf)*&
                  abs(reflect(nf,ne))*exp((-2*elem_field(ne_length,ne))*(real(prop_const(nf,ne), 8)))*&
                  cos(omega*time+b(nf)+&
                  atan2(dimag(p_factor(nf,ne)),real(p_factor(nf,ne), 8))+&
                  (-2*elem_field(ne_length,ne))*(dimag(prop_const(nf,ne)))+&
                  atan2(dimag(reflect(nf,ne)),real(reflect(nf,ne), 8))+&
                  atan2(dimag(char_admit(nf,ne)),real(char_admit(nf,ne), 8)))

          enddo
          time=time+dt
      enddo
    elseif(bc_type.eq.'flow')then
      do nt=1,n_time
          do nf=1,no_freq
              omega=2*pi*nf*harmonic_scale

              forward_pressure(nt)=forward_pressure(nt)+(abs(q_factor(nf,ne))/abs(char_admit(nf,ne)))*a(nf)*&
              cos(omega*time+b(nf)+atan2(dimag(q_factor(nf,ne)),real(q_factor(nf,ne), 8))-&
              atan2(dimag(char_admit(nf,ne)),real(char_admit(nf,ne), 8)))

              reflected_pressure(nt)=reflected_pressure(nt)+(abs(q_factor(nf,ne))/abs(char_admit(nf,ne)))*a(nf)*&
                  abs(reflect(nf,ne))*exp((-2*elem_field(ne_length,ne))*(real(prop_const(nf,ne), 8)))*&
                  cos(omega*time+b(nf)+&
                  atan2(dimag(q_factor(nf,ne)),real(q_factor(nf,ne), 8))+&
                  (-2*elem_field(ne_length,ne))*(dimag(prop_const(nf,ne)))+&
                  atan2(dimag(reflect(nf,ne)),real(reflect(nf,ne), 8))-&
                  atan2(dimag(char_admit(nf,ne)),real(char_admit(nf,ne), 8)))

              forward_flow(nt)=forward_flow(nt)+abs(q_factor(nf,ne))*a(nf)*&
                  cos(omega*time+b(nf)+&
                  atan2(dimag(q_factor(nf,ne)),real(q_factor(nf,ne), 8)))

              reflected_flow(nt)=reflected_flow(nt)+abs(q_factor(nf,ne))*a(nf)*&
                  abs(reflect(nf,ne))*exp((-2*elem_field(ne_length,ne))*(real(prop_const(nf,ne), 8)))*&
                  cos(omega*time+b(nf)+&
                  atan2(dimag(q_factor(nf,ne)),real(q_factor(nf,ne), 8))+&
                  (-2*elem_field(ne_length,ne))*(dimag(prop_const(nf,ne)))+&
                  atan2(dimag(reflect(nf,ne)),real(reflect(nf,ne), 8)))

          enddo
          time=time+dt
      enddo
    else
      print *, "ERROR: Boundary condition type not recognised. Choose flow or pressure type."
      call exit(0)
    endif
    np=elem_nodes(2,ne) ! outlet node
    write(140,fmt=*) forward_pressure+reflected_pressure + node_field(nj_bv_press,np) !Outlet total pressure
    write(150,fmt=*) forward_flow - reflected_flow + elem_field(ne_Qdot,ne) !Outlet total flow waveform
    write(160,fmt=*) forward_pressure !Outlet forward flow
    write(170,fmt=*) reflected_pressure !Outlet reflected flow
    write(180,fmt=*) forward_flow !Outlet forward pressure
    write(190,fmt=*) reflected_flow !Outlet reflected pressure

    close(fid)
    close(fid2)
    close(fid3)
    close(fid4)
    close(fid6)
    close(fid7)
    close(fid8)
    close(fid9)
    close(140)
    close(150)
    close(160)
    close(170)
    close(180)
    close(190)
    close(200)
    close(210)
    close(220)
    close(230)


  !!DEALLOCATE MEMORY
  deallocate (p_terminal)
  deallocate (p_previous)
  deallocate (terminal_flow)
  deallocate (eff_admit, STAT = AllocateStatus)
  deallocate (char_admit, STAT = AllocateStatus)
  deallocate (reflect, STAT = AllocateStatus)
  deallocate (prop_const, STAT=AllocateStatus)
  deallocate (p_factor, STAT=AllocateStatus)
  deallocate (q_factor, STAT=AllocateStatus)
  deallocate (forward_pressure, STAT=AllocateStatus)
  deallocate (reflected_pressure, STAT=AllocateStatus)
  deallocate (forward_pressure_previous, STAT=AllocateStatus)
  deallocate (reflected_pressure_previous, STAT=AllocateStatus)
  deallocate (forward_flow, STAT=AllocateStatus)
  deallocate (reflected_flow, STAT=AllocateStatus)
  deallocate(terminals_radius, STAT=AllocateStatus)
  deallocate(WSS, STAT=AllocateStatus)
  call enter_exit(sub_name,2)
end subroutine evaluate_wave_transmission
!
!##############################################################################
!
subroutine write_time_dependent_perfusion_header(output_unit,field_name,field_units, &
  start_time,dt,n_time)

  integer, intent(in) :: output_unit,n_time
  character(len=*), intent(in) :: field_name,field_units
  real(dp), intent(in) :: start_time,dt

  integer :: nt

  write(output_unit,'(A,A,A)') '# Time-dependent perfusion ',trim(field_name), &
    ' for one cardiac cycle'
  write(output_unit,'(A)') '# Time units: s'
  write(output_unit,'(A,A)') '# Value units: ',trim(field_units)
  write(output_unit,'(A)',advance='no') '# Element'
  do nt=1,n_time
    write(output_unit,'(1X,ES16.8)',advance='no') start_time+(nt-1)*dt
  enddo
  write(output_unit,*)

end subroutine write_time_dependent_perfusion_header
!
!##############################################################################
!
!*boundary_admittance* applies chosen admittance boundary conditions at the terminal units
subroutine boundary_admittance(no_freq,eff_admit,char_admit,admit_param,harmonic_scale,&
  density,viscosity,elast_param,mesh_type)

  integer, intent(in) :: no_freq
  complex(dp), intent(inout) :: eff_admit(1:no_freq,num_elems)
  complex(dp), intent(inout) :: char_admit(1:no_freq,num_elems)
  type(all_admit_param) :: admit_param
  real(dp), intent(in) :: harmonic_scale
  real(dp), intent(in) :: density
  real(dp), intent(in) :: viscosity
  character(len=60), intent(in) :: mesh_type

  type(elasticity_param) :: elast_param

  !local variables
  integer :: nf,ne,nunit
  real(dp) :: omega,R1,R2,C,length,radius,C_term,E,h_bar
  real(dp) ::  h,L_term,R_term,vein_res
  complex(dp) :: term_admit

  character(len=60) :: sub_name
  sub_name = 'boundary_admittance'
  call enter_exit(sub_name,1)
    if(admit_param%bc_type.eq.'two_unit_wk')then
      R1=admit_param%two_parameter%admit_P1
      C=admit_param%two_parameter%admit_P2
      do nf=1,no_freq !step through frequencies
        omega=nf*2*PI*harmonic_scale
        if(mesh_type.eq.'simple_tree')then
          do nunit=1,num_units
            ne=units(nunit)
           !temporarily store in eff_admit, to be added to the char admit
            eff_admit(nf,ne)=(1.0_dp+cmplx(0.0_dp,1.0_dp)*omega*R1*C)/R1
          enddo
        elseif(mesh_type.eq.'full_plus_ladder')then
          do ne=1,num_elems
            if(elem_cnct(1,0,ne).eq.0)then
              !temporarily store in eff_admit, to be added to the char admit
              eff_admit(nf,ne)=(1.0_dp+cmplx(0.0_dp,1.0_dp)*omega*R1*C)/R1
            endif
          enddo
        endif
      enddo
    elseif(admit_param%bc_type.eq.'three_unit_wk')then
      R1=admit_param%three_parameter%admit_P1
      R2=admit_param%three_parameter%admit_P2
      C=admit_param%three_parameter%admit_P3
      do nf=1,no_freq !step through frequencies
        omega=nf*2*PI*harmonic_scale
        if(mesh_type.eq.'simple_tree')then
          do nunit=1,num_units
            ne=units(nunit)
            !temporarily store in eff_admit, to be added to the char admit
            eff_admit(nf,ne)=(1+cmplx(0.0_dp,1.0_dp)*omega*R2*C)/(R1+R2+cmplx(0,1)*omega*R1*R2*C)
          enddo
        elseif(mesh_type.eq.'full_plus_ladder')then
          do ne=1,num_elems
            if(elem_cnct(1,0,ne).eq.0)then
              !temporarily store in eff_admit, to be added to the char admit
              eff_admit(nf,ne)=(1+cmplx(0,1)*omega*R2*C)/(R1+R2+cmplx(0,1)*omega*R1*R2*C)
            endif
          enddo
        endif
      enddo
    elseif(admit_param%bc_type.eq.'two_wk_plus')then
      !special case for uterine arteries which are in parallel with shunts
      E=elast_param%elasticity_parameters(1) !Pa
      h_bar=elast_param%elasticity_parameters(1)!this is a fraction of the radius so is unitless
      vein_res=0.45_dp
      R1=admit_param%four_parameter%admit_P1
      C=admit_param%four_parameter%admit_P2
      length=admit_param%four_parameter%admit_P3
      radius=admit_param%four_parameter%admit_P4
      do nf=1,no_freq !step through frequencies
        omega=nf*2*PI*harmonic_scale
        do nunit=1,num_units
          ne=units(nunit)
         !temporarily store in eff_admit, to be added to the char admit
         !ADMITTANCE DUE TO THE TERMINAL LOAD
          eff_admit(nf,ne)=(1.0_dp+cmplx(0,1.0_dp)*omega*R1*C)/R1
          ! A SECOND ADMITTANCE IN PARALLEL REPRESENTING SHUNTS
          h=h_bar*radius
          C_term=3.0_dp*PI*radius**3/(2.0_dp*h*E)!
          L_term=9.0_dp*density&
             /(4.0_dp*PI*radius**2)!per unit length
          R_term=81.0_dp*viscosity/ &
             (8.0_dp*PI*radius**4) !laminar resistance per unit length
         !G=0.0_dp
          term_admit=sqrt(cmplx(0.0_dp,1.0_dp,8)*omega*C_term)&
            /sqrt(R_term+cmplx(0.0_dp,1.0_dp,8)*omega*L_term)*50.0_dp*1.0_dp
                        term_admit=term_admit/(1+term_admit*vein_res)
          eff_admit(nf,ne)=term_admit+eff_admit(nf,ne)
        enddo
      enddo
    elseif(admit_param%bc_type.eq.'zero_reflection')then
     do nf=1,no_freq
      if(mesh_type.eq.'simple_tree')then
        do nunit=1,num_units
          ne=units(nunit)
          eff_admit(nf,ne)=char_admit(nf,ne) !in effective admittance subroutine this will become a 'dummy' daughter admittance and zero the reflection coefficient
        enddo
       elseif(mesh_type.eq.'full_plus_ladder')then
          do ne=1,num_elems
            if(elem_cnct(1,0,ne).eq.0)then
              !temporarily store in eff_admit
              eff_admit(nf,ne)=char_admit(nf,ne)
            endif
          enddo
        endif
      enddo
    endif
      call enter_exit(sub_name,2)
end subroutine boundary_admittance

!
!##############################################################################
!
!*characteristic_admittance* calculates the characteristic admittance of each
subroutine characteristic_admittance(no_freq,char_admit,prop_const,alpha,harmonic_scale,&
  density,viscosity,admit_param,elast_param,mechanics_parameters,grav_vect,remodeling_grade)

  integer, intent(in) :: no_freq
  complex(dp), intent(inout) :: char_admit(1:no_freq,num_elems)
  complex(dp), intent(inout) :: prop_const(1:no_freq,num_elems)
  real(dp), intent(inout) :: alpha(num_elems)
  real(dp), intent(in) :: harmonic_scale
  real(dp), intent(in) :: density
  real(dp), intent(in) :: viscosity
  real(dp),intent(in) :: mechanics_parameters(2),grav_vect(3)
  real(dp), intent(in) :: remodeling_grade
  type(elasticity_param) :: elast_param
  type(all_admit_param) :: admit_param

  !local variables
  real(dp) :: L,C,R, G,omega,gen_factor
  real(dp) :: E,h_bar,h,wavespeed,wolmer !should be global - maybe express as alpha (i.e. pre multiply)
  complex(dp) :: f10,bessel0,bessel1
  integer :: ne,nf,nn,np
  integer :: exit_status=0
  real(dp) :: R0,Ppl,Ptm,Rg_in,Rg_out
  real(dp) :: alt_hyp,alt_fib,prox_fib,narrow_rad_one,narrow_rad_two,narrow_factor,prune_rad,prune_fraction ! Remodeling parameters
  character(len=60) :: sub_name
  sub_name = 'characteristic_admittance'
  call enter_exit(sub_name,1)
  alpha=elast_param%elasticity_parameters(1)
  if(remodeling_grade.eq.0.0_dp) then  ! Solving for Healthy
    do ne=1,num_elems
      do nn=1,2
        if(nn.eq.1) np=elem_nodes(1,ne)
        if(nn.eq.2) np=elem_nodes(2,ne)
        call calculate_ppl(np,grav_vect,mechanics_parameters,Ppl)
        Ptm=Ppl     ! Pa
        if(nn.eq.1)R0=elem_field(ne_radius_in0,ne)
        if(nn.eq.2)R0=elem_field(ne_radius_out0,ne)
        if(admit_param%admittance_type.eq.'duan_zamir')then!alpha controls elasticity
         if(nn.eq.1)Rg_in=R0*(Ptm*elast_param%elasticity_parameters(1)+1.d0)
         if(nn.eq.2)Rg_out=R0*(Ptm*elast_param%elasticity_parameters(1)+1.d0)
        else !Hooke type elasticity
           h=elast_param%elasticity_parameters(2)*R0
          if(nn.eq.1) Rg_in=R0+3.0_dp*R0**2*Ptm/(4.0_dp*elast_param%elasticity_parameters(1)*h)
          if(nn.eq.2) Rg_out=R0+3.0_dp*R0**2*Ptm/(4.0_dp*elast_param%elasticity_parameters(1)*h)
        endif
      enddo
      elem_field(ne_radius_out,ne)=(Rg_in+Rg_out)/2.0_dp
    enddo
    E=elast_param%elasticity_parameters(1) !Pa
    h_bar=elast_param%elasticity_parameters(2)!this is a fraction of the radius so is unitless
    do ne=1,num_elems
      if(admit_param%admittance_type.eq.'lachase_standard')then
        h=h_bar*elem_field(ne_radius_out,ne)
        C=3.0_dp*PI*elem_field(ne_radius_out,ne)**3*elem_field(ne_length,ne)/(2.0_dp*h*E)
        L=density*elem_field(ne_length,ne)/(4*PI*elem_field(ne_radius_out,ne)**2)
        R=8.0_dp*viscosity*elem_field(ne_length,ne)/ &
            (PI*elem_field(ne_radius_out,ne)**4) !laminar resistance
        G=0.0_dp
      elseif(admit_param%admittance_type.eq.'lachase_modified')then
        h=h_bar*elem_field(ne_radius_out,ne)
        C=3.0_dp*PI*elem_field(ne_radius_out,ne)**3/(2.0_dp*h*E)!
        L=9.0_dp*density&
           /(4.0_dp*PI*elem_field(ne_radius_out,ne)**2)!per unit length
        R=81.0_dp*viscosity/ &
             (8.0_dp*PI*elem_field(ne_radius_out,ne)**4) !laminar resistance per unit length
        G=0.0_dp
      elseif(admit_param%admittance_type.eq.'zhu_chesler')then
        h=h_bar*elem_field(ne_radius_out,ne)
        C=3.0_dp*PI*elem_field(ne_radius_out,ne)**3*elem_field(ne_length,ne)/(2.0_dp*h*E)
        L=9.0_dp*density*elem_field(ne_length,ne)/(4.0_dp*PI*elem_field(ne_radius_out,ne)**2)
        R=8.0_dp*viscosity*elem_field(ne_length,ne)/ &
              (PI*elem_field(ne_radius_out,ne)**4) !laminar resistance
        G=0.0_dp
      elseif(admit_param%admittance_type.eq.'duan_zamir')then
        do nf=1,no_freq !radius needs to  be multipled by 1000 to go to mm (units of rest of model)
         omega=nf*2*PI*harmonic_scale!q/s
         wolmer=(elem_field(ne_radius_out,ne))*sqrt(omega*density/viscosity)
         call bessel_complex(wolmer*cmplx(0.0_dp,1.0_dp,8)**(3.0_dp/2.0_dp),bessel0,bessel1)
         f10=2*bessel1/(wolmer*cmplx(0.0_dp,1.0_dp,8)**(3.0_dp/2.0_dp)*bessel0)!no units
         wavespeed=sqrt(1.0_dp/(2*density*elast_param%elasticity_parameters(1)))*sqrt(1-f10)! !mm/s
         char_admit(nf,ne)=PI*(elem_field(ne_radius_out,ne))**2/(density*wavespeed/(1-f10))*sqrt(1-f10)!mm3/Pa
         prop_const(nf,ne)=cmplx(0.0_dp,1.0_dp,8)*omega/(wavespeed)!1/mm
        enddo
      else !Unrecognised admittance model
        print *, "EXITING"
        print *, "Unrecognised admittance model, please check inputs"
        call exit(exit_status)
      endif
      if(admit_param%admittance_type.eq.'duan_zamir')then
      else
        do nf=1,no_freq
          omega=nf*2*PI*harmonic_scale
          char_admit(nf,ne)=sqrt(G+cmplx(0.0_dp,1.0_dp,8)*omega*C)/sqrt(R+cmplx(0.0_dp,1.0_dp,8)*omega*L)!mm3/Pa.s
          prop_const(nf,ne)=sqrt((G+cmplx(0.0_dp,1.0_dp,8)*omega*C)*(R+cmplx(0.0_dp,1.0_dp,8)*omega*L))!1/mm
        enddo!nf
      endif
    enddo!ne
  else ! Solving for remodeling case - only implemented for elastic_alpha

    ! hypertophy effect
    if(remodeling_grade.ge.60)then
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
      prune_fraction = (1.0_dp/16) * remodeling_grade - 1.0_dp/8
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

    do ne=1,num_elems
      do nn=1,2
        if(nn.eq.1) np=elem_nodes(1,ne)
        if(nn.eq.2) np=elem_nodes(2,ne)
        call calculate_ppl(np,grav_vect,mechanics_parameters,Ppl)
        Ptm=Ppl     ! Pa
        if(nn.eq.1)R0=elem_field(ne_radius_in0,ne)
        if(nn.eq.2)R0=elem_field(ne_radius_out0,ne)
        if(admit_param%admittance_type.eq.'duan_zamir')then!alpha controls elasticity
           if(elem_field(ne_group,ne).eq.0.0_dp)then !applying remodeling factors on arteries only
             if(nn.eq.1) then
               if((R0.gt.narrow_rad_one).and.(R0.lt.0.5)) then !Hypertophy+narrowing effect
                 if(R0.lt.0.05) then !only narrowing factor
                   Rg_in=narrow_factor*R0*(Ptm*elast_param%elasticity_parameters(1)+1.d0)
                 elseif(R0.gt.narrow_rad_two) then ! Hypertrophy only
                   Rg_in=R0*(Ptm*alt_hyp*elast_param%elasticity_parameters(1)+1.d0)
                   alpha(ne)=alt_hyp*elast_param%elasticity_parameters(1)
                 else ! both hypertophy and narrowing
                   Rg_in=narrow_factor*R0*(Ptm*alt_hyp*alt_fib*elast_param%elasticity_parameters(1)+1.d0)
                   alpha(ne)=alt_hyp*alt_fib*elast_param%elasticity_parameters(1)
                 endif
               else ! out of range of target vessels,hence, No remodeling
                 Rg_in=R0*(Ptm*elast_param%elasticity_parameters(1)+1.d0)
               endif ! radius condition
             endif ! nn=1
             if(nn.eq.2) then
               if((R0.gt.narrow_rad_one).and.(R0.lt.0.5)) then !Hypertophy+narrowing effect
                if(R0.lt.0.05) then !only narrowing factor
                  Rg_out=narrow_factor*R0*(Ptm*elast_param%elasticity_parameters(1)+1.d0)
                elseif(R0.gt.narrow_rad_two) then ! only hypertophy
                  Rg_out=R0*(Ptm*alt_hyp*elast_param%elasticity_parameters(1)+1.d0)
                  alpha(ne)=alt_hyp*elast_param%elasticity_parameters(1)
                else ! both hypertophy and narrowing
                  Rg_out=narrow_factor*R0*(Ptm*alt_hyp*alt_fib*elast_param%elasticity_parameters(1)+1.d0)
                  alpha(ne)=alt_hyp*alt_fib*elast_param%elasticity_parameters(1)
                endif
               else ! out of range of target vessels,hence, No remodeling
                 Rg_out=R0*(Ptm*elast_param%elasticity_parameters(1)+1.d0)
               endif ! radius condition
             endif ! nn=2
           else !everything except arteries is treated as normal
             if(nn.eq.1)Rg_in=R0*(Ptm*elast_param%elasticity_parameters(1)+1.d0)
             if(nn.eq.2)Rg_out=R0*(Ptm*elast_param%elasticity_parameters(1)+1.d0)
           endif
        else ! Hooke type elasticity
           h=elast_param%elasticity_parameters(2)*R0
          if(nn.eq.1) Rg_in=R0+3.0_dp*R0**2*Ptm/(4.0_dp*elast_param%elasticity_parameters(1)*h)
          if(nn.eq.2) Rg_out=R0+3.0_dp*R0**2*Ptm/(4.0_dp*elast_param%elasticity_parameters(1)*h)
        endif
      enddo
      elem_field(ne_radius_out,ne)=(Rg_in+Rg_out)/2.0_dp
    enddo
    E=elast_param%elasticity_parameters(1) !Pa
    h_bar=elast_param%elasticity_parameters(2)!this is a fraction of the radius so is unitless
    do ne=1,num_elems
      if(admit_param%admittance_type.eq.'lachase_standard')then
        h=h_bar*elem_field(ne_radius_out,ne)
        C=3.0_dp*PI*elem_field(ne_radius_out,ne)**3*elem_field(ne_length,ne)/(2.0_dp*h*E)
        L=density*elem_field(ne_length,ne)/(4*PI*elem_field(ne_radius_out,ne)**2)
        R=8.0_dp*viscosity*elem_field(ne_length,ne)/ &
            (PI*elem_field(ne_radius_out,ne)**4) !laminar resistance
        G=0.0_dp
      elseif(admit_param%admittance_type.eq.'lachase_modified')then
        h=h_bar*elem_field(ne_radius_out,ne)
        C=3.0_dp*PI*elem_field(ne_radius_out,ne)**3/(2.0_dp*h*E)!
        L=9.0_dp*density&
           /(4.0_dp*PI*elem_field(ne_radius_out,ne)**2)!per unit length
        R=81.0_dp*viscosity/ &
             (8.0_dp*PI*elem_field(ne_radius_out,ne)**4) !laminar resistance per unit length
        G=0.0_dp
      elseif(admit_param%admittance_type.eq.'zhu_chesler')then
        h=h_bar*elem_field(ne_radius_out,ne)
        C=3.0_dp*PI*elem_field(ne_radius_out,ne)**3*elem_field(ne_length,ne)/(2.0_dp*h*E)
        L=9.0_dp*density*elem_field(ne_length,ne)/(4.0_dp*PI*elem_field(ne_radius_out,ne)**2)
        R=8.0_dp*viscosity*elem_field(ne_length,ne)/ &
              (PI*elem_field(ne_radius_out,ne)**4) !laminar resistance
        G=0.0_dp
      elseif(admit_param%admittance_type.eq.'duan_zamir')then
       do nf=1,no_freq !radius needs to  be multipled by 1000 to go to mm (units of rest of model)
         omega=nf*2*PI*harmonic_scale!q/s
         wolmer=(elem_field(ne_radius_out,ne))*sqrt(omega*density/viscosity) !radii is already affected by a factor
         call bessel_complex(wolmer*cmplx(0.0_dp,1.0_dp,8)**(3.0_dp/2.0_dp),bessel0,bessel1)
         f10=2*bessel1/(wolmer*cmplx(0.0_dp,1.0_dp,8)**(3.0_dp/2.0_dp)*bessel0)!no units
         if(elem_field(ne_group,ne).eq.0.0_dp)then !applying elasticity factor on wavespeed for arteries
           if((R0.gt.narrow_rad_one).and.(R0.lt.0.5))then
             if(R0.lt.0.05) then ! only narrowing
               wavespeed=sqrt(1.0_dp/(2*density*elast_param%elasticity_parameters(1)))*sqrt(1-f10)! !mm/s
             elseif(R0.gt.narrow_rad_two) then ! only Hypertrophy
               wavespeed=sqrt(1.0_dp/(2*density*alt_hyp*elast_param%elasticity_parameters(1)))*sqrt(1-f10)! !mm/s
             else ! narrowing and hypertrophy (also fibrosis depending on grade)
               wavespeed=sqrt(1.0_dp/(2*density*alt_hyp*alt_fib*elast_param%elasticity_parameters(1)))*sqrt(1-f10)! !mm/s
             endif
           else ! not in range of remodeling target radii, hence, no remodeling
             wavespeed=sqrt(1.0_dp/(2*density*elast_param%elasticity_parameters(1)))*sqrt(1-f10)! !mm/s
           endif
         else !apply normal elasticity on everything except arteries
           wavespeed=sqrt(1.0_dp/(2*density*elast_param%elasticity_parameters(1)))*sqrt(1-f10)! !mm/s
         endif
         char_admit(nf,ne)=PI*(elem_field(ne_radius_out,ne))**2/(density*wavespeed/(1-f10))*sqrt(1-f10)!mm3/Pa
         prop_const(nf,ne)=cmplx(0.0_dp,1.0_dp,8)*omega/(wavespeed)!1/mm
       enddo
      else !Unrecognised admittance model
        print *, "EXITING"
        print *, "Unrecognised admittance model, please check inputs"
        call exit(exit_status)
      endif
      if(admit_param%admittance_type.eq.'duan_zamir')then
      else
        do nf=1,no_freq
          omega=nf*2*PI*harmonic_scale
          char_admit(nf,ne)=sqrt(G+cmplx(0.0_dp,1.0_dp,8)*omega*C)/sqrt(R+cmplx(0.0_dp,1.0_dp,8)*omega*L)!mm3/Pa.s
          prop_const(nf,ne)=sqrt((G+cmplx(0.0_dp,1.0_dp,8)*omega*C)*(R+cmplx(0.0_dp,1.0_dp,8)*omega*L))!1/mm
        enddo!nf
      endif
    enddo!ne
  endif
    call enter_exit(sub_name,2)
  end subroutine characteristic_admittance

!##################################################################
!
!*tree_admittance:* Calculates the total admittance of a tree
subroutine tree_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale,&
  min_elem,max_elem,tree_direction)

  integer, intent(in) :: no_freq
  complex(dp), intent(inout) :: eff_admit(1:no_freq,num_elems)
  complex(dp), intent(in) :: char_admit(1:no_freq,num_elems)
  complex(dp), intent(inout) :: reflect(1:no_freq,num_elems)
  complex(dp), intent(in) :: prop_const(1:no_freq,num_elems)
  real(dp), intent(in) :: harmonic_scale
  integer, intent(in) :: min_elem,max_elem
  character(len=30), intent(in) :: tree_direction

  character(len=60) :: sub_name
!local variables
  real(dp) :: invres,elem_res(num_elems),omega
  integer :: num2,ne,ne2,nf,num3,ne3,ne_sist
  complex(dp) :: daughter_admit,sister_admit,sister_current,a,b,m1,m2

    sub_name = 'tree_admittance'
    call enter_exit(sub_name,1)
    reflect(:,:)=cmplx(0.0_dp,0.0_dp,8)

    if(tree_direction.eq.'diverging')then
        do nf=1,no_freq
            omega=nf*2*PI*harmonic_scale
            do ne=max_elem,min_elem,-1!step backward through elements
                daughter_admit=cmplx(0.0_dp,0.0_dp,8)!

                do num2=1,elem_cnct(1,0,ne)!will only do stuff to non-terminals will add one daughter if no branching
                    ne2=elem_cnct(1,num2,ne)!for each downstream element
                    daughter_admit=daughter_admit+eff_admit(nf,ne2)!sum of two child elements
                enddo
                if(elem_cnct(1,0,ne).gt.0)then !not a terminal as it has a downstream
                    reflect(nf,ne)=(char_admit(nf,ne)-daughter_admit)/&
                        (char_admit(nf,ne)+daughter_admit)!double checked
                    eff_admit(nf,ne)=char_admit(nf,ne)*(1&
                        -reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))/&
                        (1+reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))
                else !a terminal
                    daughter_admit=eff_admit(nf,ne) !a boundary condition is applied here
                    reflect(nf,ne)=(char_admit(nf,ne)-daughter_admit)/&
                        (char_admit(nf,ne)+daughter_admit)
                     !now we overwrite the effective admittance of the terminal to include reflection from the daughter.
                    eff_admit(nf,ne)=char_admit(nf,ne)*(1&
                        -reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))/&
                        (1&
                        +reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))
                endif
            enddo!ne
        enddo!nf
    elseif(tree_direction.eq.'converging')then
        do nf=1,no_freq
            omega=nf*2*PI*harmonic_scale
            do ne=min_elem,max_elem!step forward through elements
                daughter_admit=cmplx(0.0_dp,0.0_dp,8)!
                sister_admit=cmplx(0.0_dp,0.0_dp,8)!
                do num2=1,elem_cnct(1,0,ne)!will only do stuff to non-terminals
                    ne2=elem_cnct(1,num2,ne)!for each downstream element
                    daughter_admit=daughter_admit+eff_admit(nf,ne2)!sum of two child elements
                    do num3=1,elem_cnct(-1,0,ne2)!sisters
                        ne3=elem_cnct(-1,num3,ne2)!for each upstream element of the daughter
                        if(ne3.ne.ne)then
                            ne_sist=ne3
                            sister_admit=sister_admit+char_admit(nf,ne3)
                            sister_current=exp(-1.0_dp*prop_const(nf,ne3)*elem_field(ne_length,ne3))/&
                                exp(-1.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne))
                            a=exp(-1.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne))
                            b=exp(-1.0_dp*prop_const(nf,ne3)*elem_field(ne_length,ne3))
                        endif
                    enddo
                enddo
                m1=1.0_dp+2.0_dp*(b/a)*((2*a*b-a**2-1)*char_admit(nf,ne)+&
                    (a**2-1)*sister_admit+(a**2-1)*daughter_admit)/((1-b**2)*char_admit(nf,ne)+&
                    (1+b**2-2*a*b)*sister_admit+(1-b**2)*daughter_admit)
                if(elem_cnct(1,0,ne).gt.0)then !not a terminal
                    reflect(nf,ne)=(char_admit(nf,ne)-m1*sister_admit-daughter_admit)/&
                        (char_admit(nf,ne)+daughter_admit+sister_admit)!ARC- to check
                    eff_admit(nf,ne)=char_admit(nf,ne)*(1&
                        -reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))/&
                        (1+reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))
                      else !a terminal
                    daughter_admit=eff_admit(nf,ne) !a boundary condition is applied here
                    reflect(nf,ne)=(char_admit(nf,ne)-daughter_admit)/&
                        (char_admit(nf,ne)+daughter_admit)
                    eff_admit(nf,ne)=char_admit(nf,ne)*(1&
                        -reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))/&
                        (1+reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))
                endif
            enddo!ne
        enddo!nf
    endif

  call enter_exit(sub_name,2)
end subroutine tree_admittance
!##################################################################
!
!*capillaryadmittance:* Calculates the total admittance of a tree
subroutine capillary_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale,&
  min_elem,max_elem,elast_param,mechanics_parameters,grav_vect,cap_model)

  integer, intent(in) :: no_freq
  complex(dp), intent(inout) :: eff_admit(1:no_freq,num_elems)
  complex(dp), intent(inout) :: char_admit(1:no_freq,num_elems)
  complex(dp), intent(inout) :: reflect(1:no_freq,num_elems)
  complex(dp), intent(in) :: prop_const(1:no_freq,num_elems)
  real(dp), intent(in) :: harmonic_scale
  integer, intent(in) :: min_elem,max_elem
  real(dp),intent(in) :: mechanics_parameters(2),grav_vect(3)
  integer, intent(in) :: cap_model

  type(capillary_bf_parameters) :: cap_param
  type(elasticity_param) :: elast_param

  character(len=60) :: sub_name
!local variables
  integer :: ne, ne2,num2,nf,ne0,ne1
  real(dp) :: daughter_admit,omega,length,Hart,Hven,Gamma_sheet,alpha_c,Ptp
  integer :: grav_dirn,i
  real (dp) :: Ppl,P1,P2
  real(dp) :: Q01,Rin,Rout,Lin,Lout,x_cap,y_cap,z_cap
  complex(dp) :: eff_admit_downstream(no_freq)

  sub_name = 'capillary_admittance'
  call enter_exit(sub_name,1)
  reflect(:,:)=cmplx(0.0_dp,0.0_dp,8)

  do ne=min_elem,max_elem
    ne0=elem_cnct(-1,1,ne)!upstream element number
    ne1=elem_cnct(1,1,ne) !downstream element number
    P1=node_field(nj_bv_press,elem_nodes(2,ne0)) !pressure at start node of capillary element
    P2=node_field(nj_bv_press,elem_nodes(1,ne1))
    Q01=elem_field(ne_Qdot,ne0) !flow in element upstream of capillary element !mm^3/s
    Rin=elem_field(ne_radius_out0,ne0)!radius of upstream element
    Rout=elem_field(ne_radius_out0,ne1) !radius of downstream element
    x_cap=node_xyz(1,elem_nodes(1,ne))
    y_cap=node_xyz(2,elem_nodes(1,ne))
    z_cap=node_xyz(3,elem_nodes(1,ne))
    call calculate_ppl(elem_nodes(1,ne),grav_vect,mechanics_parameters,Ppl)
    Lin=elem_field(ne_length,ne0)
    Lout=elem_field(ne_length,ne1)
    Ptp=(cap_param%Palv-(-Ppl))/98.06d0 !Pa -> cmH2O
    do i=1,no_freq
      eff_admit_downstream(i)=eff_admit(i,ne1)
    enddo
    call cap_flow_admit(ne,eff_admit(:,ne),eff_admit_downstream,Lin,Lout,P1,P2,&
                        Ppl,Q01,Rin,Rout,x_cap,y_cap,z_cap,no_freq,harmonic_scale,&
                        elast_param,cap_model)
   enddo!ne

end subroutine capillary_admittance
!
!################################################
!
!*pressure_factor:* Calculates change in pressure and flow through tree
subroutine pressure_flow_factor(no_freq,p_factor,q_factor,reflect,prop_const,char_admit,eff_admit&
  ,harmonic_scale,ne_min,ne_max,bc_type)

  integer, intent(in) :: no_freq
  complex(dp), intent(inout) :: p_factor(1:no_freq,num_elems)
  complex(dp), intent(inout) :: q_factor(1:no_freq,num_elems)
  complex(dp), intent(inout) :: reflect(1:no_freq,num_elems)
  complex(dp), intent(in) :: prop_const(1:no_freq,num_elems)
  complex(dp), intent(in) :: char_admit(1:no_freq,num_elems)
  complex(dp), intent(in) :: eff_admit(1:no_freq,num_elems)
  real(dp), intent(in) :: harmonic_scale
  integer, intent(in) :: ne_min,ne_max
  character(len=60) :: bc_type



  character(len=60) :: sub_name
!local variables
  integer :: ne, nf,ne_up
  real(dp) :: omega

  sub_name = 'pressure_flow_factor'
  call enter_exit(sub_name,1)

  p_factor=1.0_dp
  q_factor=1.0_dp
    do nf=1,no_freq
      omega=nf*2*PI*harmonic_scale
      do ne=ne_min,ne_max
        !look for upstream element
        if(elem_cnct(-1,0,ne).eq.0)then !no upstream elements, inlet, ignore
        ne_up=ne_min
          p_factor(nf,ne)=(1.0_dp)!* &!assumes input admittance is the same as characteristic admittance for this vessel
        else
          ne_up=elem_cnct(-1,1,ne)
          p_factor(nf,ne)=p_factor(nf,ne_up)*(1+reflect(nf,ne_up))* &
            exp(-1.0_dp*elem_field(ne_length,ne_up)*prop_const(nf,ne_up))/&
            (1+reflect(nf,ne)*exp(-2.0_dp*elem_field(ne_length,ne)*prop_const(nf,ne)))
          q_factor(nf,ne)=p_factor(nf,ne)*char_admit(nf,ne)
        endif!neup
      enddo!ne
    enddo!nf

  call enter_exit(sub_name,2)
end subroutine pressure_flow_factor

end module wave_transmission
