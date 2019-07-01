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
  use arrays, only: dp
  use other_consts, only: PI
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
subroutine evaluate_wave_transmission(grav_dirn,grav_factor,&
    n_time,heartrate,a0,no_freq,a,b,n_adparams,admittance_param,n_model,model_definition,cap_model)
!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_WAVE_TRANSMISSION: EVALUATE_WAVE_PROPAGATION
  use indices
  use arrays, only: dp,all_admit_param,num_elems,elem_field,fluid_properties,elasticity_param,num_units,&
    units,node_xyz,elem_cnct,elem_nodes,node_field
  use diagnostics, only: enter_exit



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
  integer, intent(in) :: cap_model ! This determines the capillary model (1: ha=hv, 2:non-constant sheet height)

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
  real(dp), allocatable :: forward_pressure(:)
  real(dp), allocatable :: reflected_pressure(:)
  real(dp), allocatable :: forward_flow(:)
  real(dp), allocatable :: reflected_flow(:)
  integer :: min_art,max_art,min_ven,max_ven,min_cap,max_cap,ne,nu,nt,nf,np
  character(len=30) :: tree_direction,mechanics_type
  real(dp) start_time,end_time,dt,time,omega
  real(dp) grav_vect(3),grav_factor,mechanics_parameters(2)
  integer :: AllocateStatus,fid=10,fid2=20,fid3=30,fid4=40,fid5=50
  character(len=60) :: sub_name
  logical :: constant_sheet_height

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
  
  !! Determining the capillary model (Constant sheet height / Non-constant sheet height)
  if (cap_model.eq.1) then
     constant_sheet_height = .True.
  else
     constant_sheet_height = .False.
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
  else!otherwise input a0 is used
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
  allocate (forward_pressure(n_time), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for forward_p array ***"
  allocate (reflected_pressure(n_time), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for reflected_p array ***"
    allocate (forward_flow(n_time), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for forward_q array ***"
  allocate (reflected_flow(n_time), STAT = AllocateStatus)
  if (AllocateStatus /= 0) STOP "*** Not enough memory for reflected_q array ***"

  !initialise admittance
  char_admit=0.0_dp
  eff_admit=0.0_dp
  !calculate characteristic admittance of each branch
  call characteristic_admittance(no_freq,char_admit,prop_const,harmonic_scale, &
    density,viscosity,admit_param,elast_param,mechanics_parameters,grav_vect)

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
            min_cap,max_cap,elast_param,mechanics_parameters,grav_vect,constant_sheet_height)
        !art admittance
        tree_direction='diverging'
        call tree_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale,&
            min_art,max_art,tree_direction)
    else!Assume simple tree
        tree_direction='diverging'
        min_art=1
        max_art=num_elems
        call tree_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale,&
            min_art,max_art,tree_direction)
    endif

!
!    !calculate pressure drop through arterial tree (note to do veins too need to implement this concept thro' whole ladder model)
!    !Also need to implement in reverse for veins
    call pressure_factor(no_freq,p_factor,reflect,prop_const,harmonic_scale,min_art,max_art)
    open(fid5, file = 'inputimpedance.txt',action='write')
    write(fid5,fmt=*) 'input impedance:'
    do nf=1,no_freq
        omega=nf*harmonic_scale
        write(fid5,fmt=*) omega,abs(eff_admit(nf,1)),&
            atan2(imagpart(eff_admit(nf,1)),realpart(eff_admit(nf,1)))
    enddo
    close(fid5)

    start_time=0.0_dp
    end_time=60.0_dp/heartrate
    dt=(end_time-start_time)/n_time
    time=start_time
    !consider first pressure and flow into the vessel (at x=0)
    open(fid, file = 'incident_pressure.txt', action='write')
    open(fid2, file = 'incident_flow.txt',action='write')
    open(fid3, file = 'total_pressure.txt',action='write')
    open(fid4, file = 'total_flow.txt',action='write')
    do nu =1,num_units
        ne=units(nu)
        forward_pressure=0.0_dp
        reflected_pressure=0.0_dp
        forward_flow=0.0_dp
        reflected_flow=0.0_dp
        do nt=1,n_time
            do nf=1,no_freq
                omega=2*pi*nf*harmonic_scale
                forward_pressure(nt)=forward_pressure(nt)+abs(p_factor(nf,ne))*a(nf)*cos(omega*time+b(nf)+&
                    atan2(imagpart(p_factor(nf,ne)),realpart(p_factor(nf,ne))))

                reflected_pressure(nt)=reflected_pressure(nt)+abs(p_factor(nf,ne))*a(nf)*&
                    abs(reflect(nf,ne))*exp((-2*elem_field(ne_length,ne))*realpart(prop_const(nf,ne)))*&
                    cos(omega*time+b(nf)+&
                    atan2(imagpart(p_factor(nf,ne)),realpart(p_factor(nf,ne)))+&
                    (-2*elem_field(ne_length,ne))*imagpart(prop_const(nf,ne))+&
                    atan2(imagpart(reflect(nf,ne)),realpart(reflect(nf,ne))))

                forward_flow(nt)=forward_flow(nt)+abs(char_admit(nf,ne))*abs(p_factor(nf,ne))*a(nf)*&
                    cos(omega*time+b(nf)+&
                    atan2(imagpart(p_factor(nf,ne)),realpart(p_factor(nf,ne)))+&
                    atan2(imagpart(char_admit(nf,ne)),realpart(char_admit(nf,ne))))

                reflected_flow(nt)=reflected_flow(nt)+abs(char_admit(nf,ne))*abs(p_factor(nf,ne))*a(nf)*&
                    abs(reflect(nf,ne))*exp((-2*elem_field(ne_length,ne))*realpart(prop_const(nf,ne)))*&
                    cos(omega*time+b(nf)+&
                    atan2(imagpart(p_factor(nf,ne)),realpart(p_factor(nf,ne)))+&
                    (-2*elem_field(ne_length,ne))*imagpart(prop_const(nf,ne))+&
                    atan2(imagpart(reflect(nf,ne)),realpart(reflect(nf,ne)))+&
                    atan2(imagpart(char_admit(nf,ne)),realpart(char_admit(nf,ne))))

            enddo
            time=time+dt
        enddo
        np=elem_nodes(2,ne)
        write(fid,fmt=*) ne, forward_pressure+node_field(nj_bv_press,np)
        write(fid2,fmt=*) ne, forward_flow+elem_field(ne_Qdot,ne)

        write(fid3,fmt=*) ne, forward_pressure+reflected_pressure + node_field(nj_bv_press,np)
        write(fid4,fmt=*) ne, forward_flow-reflected_flow + elem_field(ne_Qdot,ne)


    enddo
    close(fid)
    close(fid2)
    close(fid3)
    close(fid4)


  !!DEALLOCATE MEMORY
  deallocate (eff_admit, STAT = AllocateStatus)
  deallocate (char_admit, STAT = AllocateStatus)
  deallocate (reflect, STAT = AllocateStatus)
  deallocate (prop_const, STAT=AllocateStatus)
  deallocate (p_factor, STAT=AllocateStatus)
  deallocate (forward_pressure, STAT=AllocateStatus)
  deallocate (reflected_pressure, STAT=AllocateStatus)
  deallocate (forward_flow, STAT=AllocateStatus)
  deallocate (reflected_flow, STAT=AllocateStatus)
  call enter_exit(sub_name,2)
end subroutine evaluate_wave_transmission
!
!##############################################################################
!
!*boundary_admittance* applies chosen admittance boundary conditions at the terminal units
subroutine boundary_admittance(no_freq,eff_admit,char_admit,admit_param,harmonic_scale,&
  density,viscosity,elast_param,mesh_type)
!DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_boundary_admittance: boundary_admittance
  use arrays,only: num_elems,all_admit_param,units,num_units,elasticity_param,elem_cnct
  use diagnostics, only: enter_exit

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
subroutine characteristic_admittance(no_freq,char_admit,prop_const,harmonic_scale,&
  density,viscosity,admit_param,elast_param,mechanics_parameters,grav_vect)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAD:"SO_characteristic_admittance: characteristic_admittance
  use other_consts, only: MAX_STRING_LEN
  use indices
  use arrays, only: num_elems,elem_field,elasticity_param,all_admit_param,elem_nodes
  use pressure_resistance_flow, only: calculate_ppl
  use math_utilities, only: bessel_complex
  use diagnostics, only: enter_exit

  integer, intent(in) :: no_freq
  complex(dp), intent(inout) :: char_admit(1:no_freq,num_elems)
  complex(dp), intent(inout) :: prop_const(1:no_freq,num_elems)
  real(dp), intent(in) :: harmonic_scale
  real(dp), intent(in) :: density
  real(dp), intent(in) :: viscosity
  real(dp),intent(in) :: mechanics_parameters(2),grav_vect(3)

  type(elasticity_param) :: elast_param
  type(all_admit_param) :: admit_param

  !local variables
  real(dp) :: L,C,R, G,omega,gen_factor
  real(dp) :: E,h_bar,h,wavespeed,wolmer !should be global - maybe express as alpha (i.e. pre multiply)
  complex(dp) :: f10,bessel0,bessel1
  integer :: ne,nf,nn,np
  integer :: exit_status=0
  real(dp) :: R0,Ppl,Ptm,Rg_in,Rg_out
  character(len=60) :: sub_name
  sub_name = 'characteristic_admittance'
  call enter_exit(sub_name,1)

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
      else!Hooke type elasticity
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

  call enter_exit(sub_name,2)
end subroutine characteristic_admittance

!##################################################################
!
!*tree_admittance:* Calculates the total admittance of a tree
subroutine tree_admittance(no_freq,eff_admit,char_admit,reflect,prop_const,harmonic_scale,&
  min_elem,max_elem,tree_direction)
  use indices
  use arrays,only: dp,num_elems,elem_cnct,elem_field
  use diagnostics, only: enter_exit
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
                else!a terminal
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
                else!a terminal
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
  min_elem,max_elem,elast_param,mechanics_parameters,grav_vect,constant_sheet_height)
  use indices
  use arrays,only: dp,num_elems,elem_cnct,elem_field,capillary_bf_parameters,elem_nodes,&
    node_field,node_xyz,elasticity_param
  use pressure_resistance_flow,only: calculate_ppl
  use capillaryflow, only:cap_flow_admit
  use diagnostics, only: enter_exit

  integer, intent(in) :: no_freq
  complex(dp), intent(inout) :: eff_admit(1:no_freq,num_elems)
  complex(dp), intent(inout) :: char_admit(1:no_freq,num_elems)
  complex(dp), intent(inout) :: reflect(1:no_freq,num_elems)
  complex(dp), intent(in) :: prop_const(1:no_freq,num_elems)
  real(dp), intent(in) :: harmonic_scale
  integer, intent(in) :: min_elem,max_elem
  real(dp),intent(in) :: mechanics_parameters(2),grav_vect(3)
  logical, intent(in) :: constant_sheet_height

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
Ppl,Q01,Rin,Rout,x_cap,y_cap,z_cap,no_freq,harmonic_scale,elast_param,constant_sheet_height)
   enddo!ne

end subroutine capillary_admittance
!
!################################################
!
!*pressure_factor:* Calculates change in pressure through tree
  subroutine pressure_factor(no_freq,p_factor,reflect,prop_const,harmonic_scale,ne_min,ne_max)
    use indices
    use arrays,only: dp,num_elems,elem_cnct,elem_field
    use diagnostics, only: enter_exit
    integer, intent(in) :: no_freq
    complex(dp), intent(inout) :: p_factor(1:no_freq,num_elems)
    complex(dp), intent(inout) :: reflect(1:no_freq,num_elems)
    complex(dp), intent(in) :: prop_const(1:no_freq,num_elems)
    real(dp), intent(in) :: harmonic_scale
    integer, intent(in) :: ne_min,ne_max



    character(len=60) :: sub_name
!local variables
    integer :: ne, nf,ne_up
    real(dp) :: omega

    sub_name = 'pressure_factor'
    call enter_exit(sub_name,1)

    p_factor=1.0_dp
    do nf=1,no_freq
      omega=nf*2*PI*harmonic_scale
      do ne=ne_min,ne_max
        !look for upstream element
        if(elem_cnct(-1,0,ne).eq.0)then !no upstream elements, inlet, ignore
        ne_up=ne_min
          p_factor(nf,ne)=(1.0_dp)!* &!assumes input admittance is the same as characteristic admittance for this vessel
            !exp(-1.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne))!/&
            !(1+reflect(nf,ne)*exp(-2.0_dp*prop_const(nf,ne)*elem_field(ne_length,ne)))
        else
          ne_up=elem_cnct(-1,1,ne)
          p_factor(nf,ne)=p_factor(nf,ne_up)*(1+reflect(nf,ne_up))* &
            exp(-1.0_dp*elem_field(ne_length,ne_up)*prop_const(nf,ne_up))/&
            (1+reflect(nf,ne)*exp(-2.0_dp*elem_field(ne_length,ne)*prop_const(nf,ne)))
        endif!neup
      enddo!ne
    enddo!nf

    call enter_exit(sub_name,2)
end subroutine pressure_factor
end module wave_transmission
