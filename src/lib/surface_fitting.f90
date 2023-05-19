module surface_fitting
!!! This is a non-generalised adaptation of the general geometric fitting
!!! methods from the Auckland Bioengineering Institute's legacy 'CMISS' code.

  use arrays
  use diagnostics
  use geometry
  use other_consts
  use mesh_utilities
  use precision
  use solve

  implicit none

  public &
       define_data_fit_group, &
       fit_surface_geometry, &
       initialise_fit_mesh, &
       pxi, &
       reset_fitting

  integer,parameter :: nmax_data_elem = 4000     ! max # data points on an element
  integer,parameter :: nmax_versn = 6            ! max # versions of node (derivative)
  integer,parameter :: num_coords = 3            ! x,y,z coordinates
  integer,parameter :: num_deriv = 4             ! coordinate + 3 derivatives
  integer,parameter :: num_deriv_elem = 16       ! nodes * num_deriv = 16
  integer,parameter :: num_elem_nodes = 4        ! code is unashamedly for 4-noded surface elements!
  integer,parameter :: num_fit = 3               ! number of parameters in fit (x,y,z)
  integer,parameter :: num_gauss = 9             ! number of gauss points

  integer,parameter :: nsize_gkk = 10000000        ! size of A matrix, until allocation sorted out

  !real(dp),parameter :: loose_tol = 1.0e-6_dp

contains

!!! ##########################################################################      

  subroutine fit_surface_geometry(niterations,fitting_file)
    !*fit_surface_geometry:* completes 'niterations' of geometry fitting to a  
    ! surface, via minimising the least squares distance between a list of 
    ! data points (3D RC coordinates) and a surface mesh (assumed bi-cubic 
    ! Hermite only). 'fitting_file' lists the nodes/derivatives that are fixed,
    !  and any mapping of nodes and/or derivatives

    integer,intent(in) :: niterations                   ! user-specified number of fitting iterations
    character(len=255),intent(in) :: fitting_file ! file that lists versions/mapping/BCs
    ! Local variables
    integer  :: nfit,nk,NOT_1,NOT_2,np,num_depvar,nv,ny_max
    logical :: first = .true., writefile = .false.
    ! local allocatable arrays
    integer,allocatable :: data_elem(:)              
    integer,allocatable :: data_on_elem(:,:)     ! list of data closest to elements
    integer,allocatable :: elem_list(:)          ! list of elements in fit (all)
    integer,allocatable :: ndata_on_elem(:)      ! # of data points closest to element
    integer,allocatable :: npny(:,:)             ! maps deriv, vers, node, etc for a dep. variable
    integer,allocatable :: np_list_redist(:,:)   ! lists of nodes to be uniformly redistributed
    integer,allocatable :: nynp(:,:,:,:)         ! dep. variable # for node, deriv, version etc.
    integer,allocatable :: nynr(:)               ! list of all dep. variables
    integer,allocatable :: nyny(:,:)             ! maps dep. variable to another dep. variable
    real(dp),allocatable :: cyny(:,:)            ! weighting for mapped dep. variables
    real(dp),allocatable :: data_xi(:,:)         ! xi location of data points
    real(dp),allocatable :: fit_soln(:,:,:,:)    ! current fit geometry solution
    real(dp),allocatable :: sobelov_wts(:,:)     ! Scaling factor for, and Sobelov weights
    logical,allocatable :: fix_bcs(:)            ! logical for boundary conditions

    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'fit_surface_geometry'
    call enter_exit(sub_name,1)

!!! allocate element-sized arrays
    allocate(elem_list(num_elems_2d))
    allocate(sobelov_wts(0:6,num_elems_2d))
    allocate(np_list_redist(num_elems_2d*5,20)) ! sizing is arbitrary

!!! allocate data arrays
    allocate(data_elem(num_data))
    allocate(data_on_elem(num_elems_2d,nmax_data_elem))
    allocate(ndata_on_elem(num_elems_2d))
    allocate(data_xi(2,num_data))
    data_elem = 0

!!! allocate dependent variable arrays
    ny_max = sum(node_versn_2d(:))*num_nodes_2d*num_fit*num_deriv  ! nodes * coordinates * #derivatives+1
    allocate(nynr(0:ny_max))
    allocate(npny(1:6,ny_max))
    allocate(nynp(num_deriv,nmax_versn,num_fit,num_nodes_2d))
    allocate(fix_bcs(ny_max))

    write(*,'('' Define boundary conditions and mapping '')')
    call define_geometry_fit(elem_list,np_list_redist,npny,num_depvar,nynp,nynr,nyny,&
         cyny,sobelov_wts,fit_soln,fitting_file,fix_bcs)
    call set_linear_derivatives
    scale_factors_2d = 1.0_dp

!!! find the closest surface to each data point, and calculate the Xi
!!! coordinates of the data point to the surface
    write(*,'('' Calculating normal projections (slow first time) for '',i8,'' data points'')') &
         num_data
    call define_xi_closest(data_elem,data_on_elem,elem_list,ndata_on_elem,data_xi,first)
    first = .false.
    
!!! list the total error between the data points and the surface
    write(*,'(/'' TOTAL RMS ERROR PRIOR TO FITTING:'')') 
    call list_data_error(data_on_elem,ndata_on_elem,data_xi)
!!! read 'fitting_file' to define the fitting constraints. set up mapping
!!! arrays, dependent variable arrays, Sobelov smoothing weights (hardcoded)
    do nfit = 1,niterations ! user-defined number of iterations
       write(*,'(/'' FITTING ITERATION'',I3)') nfit
!!!    solve for new nodal coordinates and derivatives
       write(*,'('' Solve fitting problem '')')
       call solve_geometry_fit(data_on_elem,ndata_on_elem,num_depvar,&
            elem_list,not_1,not_2,npny,nynp,&
            nyny,data_xi,cyny,sobelov_wts,fit_soln,fix_bcs)
       call update_versions(nynp,fit_soln,fix_bcs)
       call calc_arclengths
       write(*,'('' Update pseudo-landmarks locations '')')
       ! update the node locations on base, fissures, anterior and posterior lines
       call distribute_surface_node_fit(np_list_redist,nynp,fit_soln,fix_bcs) ! lateral-base

!!!    update the scale factors for new geometry if NOT unit scale factors
!       write(*,'('' Update scale factors '')')
!       call update_scale_factor_norm
!!!    update the data point projections and their Xi coordinates
       write(*,'('' Calculating normal projections '')')
       call define_xi_closest(data_elem,data_on_elem,elem_list,ndata_on_elem,data_xi,first)
!!!    calculated the updated error between data and surface
       write(*,'(/'' CURRENT PROJECTION ERROR FOR ALL DATA:'')')
       if(num_groups <= 1)then
          call list_data_error(data_on_elem,ndata_on_elem,data_xi)
       else
          if(nfit.eq.niterations)then
             writefile = .true.
          endif
          call list_data_error_by_group(data_on_elem,ndata_on_elem,data_xi,writefile)
       endif
    enddo

    call calc_data_field_distance(data_elem,data_xi)

    deallocate(elem_list)
    deallocate(sobelov_wts)
    deallocate(data_elem)
    deallocate(data_on_elem)
    deallocate(ndata_on_elem)
    deallocate(data_xi)
    deallocate(cyny)
    deallocate(nynr)
    deallocate(npny)
    deallocate(np_list_redist)
    deallocate(nynp)
    deallocate(nyny)
    deallocate(fix_bcs)
    deallocate(fit_soln)

    call enter_exit(sub_name,2)
    
  end subroutine fit_surface_geometry

!!! ##########################################################################      
  
  subroutine reset_fitting()

    if(allocated(nelem_groups)) deallocate(nelem_groups)
    if(allocated(data_xyz)) deallocate(data_xyz)
    if(allocated(data_field)) deallocate(data_field)
    if(allocated(data_weight)) deallocate(data_weight)

  end subroutine reset_fitting
  
!!! ##########################################################################      
  
  subroutine define_geometry_fit(elem_list,np_list_redist,npny,num_depvar,nynp,nynr,nyny,&
       cyny,sobelov_wts,fit_soln,fitting_file,fix_bcs)

!!! read information from 'fitting_file' to determine the boundary conditions
!!! and mapping of dependent variables for geometric fitting. set up the 
!!! dependent variable-to-mapping arrays  
  
!!! dummy arguments
    integer :: elem_list(:),np_list_redist(:,:),npny(:,:),num_depvar,nynp(:,:,:,:),nynr(0:)
    integer,allocatable :: nyny(:,:)
    real(dp) :: sobelov_wts(0:,:)
    real(dp),allocatable :: cyny(:,:),fit_soln(:,:,:,:)
    character(len=*),intent(in) :: fitting_file
    logical :: fix_bcs(:)
!!! local variables
    integer :: i,ibeg,idx,iend,ierror,IPFILE=10,i_ss_end,L,ne,nelem,nh,nj,nk, &
         node,np,np_global,n_group,number_of_maps,number_of_fixed,num_redists, &
         nv,nv_fix,ny,stat
    integer,allocatable :: nmap_info(:,:)
    real(dp) :: temp_weights(7)
    character(len=300) :: readfile,string,sub_string,temp
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'define_geometry_fit'
    call enter_exit(sub_name,1)

    if(.not.allocated(fit_soln)) allocate(fit_soln(4,10,16,num_nodes_2d))
    allocate(nmap_info(7,num_nodes_2d*num_deriv*nmax_versn))
    allocate(nelem_groups(20,num_elems_2d))
    nelem_groups = 0
    
    ! linear fitting for 3 geometric variables. solution stored in fields 1,2,3
    ! includes Sobelov smoothing on the geometry field
    
    !***Set up dependent variable interpolation information
    fit_soln = node_xyz_2d    
 
    elem_list(1:num_elems_2d) = elems_2d(1:num_elems_2d) ! global element numbers
    !elem_list = 0

    !*** Calculate ny maps
    call calculate_ny_maps(npny,num_depvar,nynp)
    
    fix_bcs = .false. !initialise, default
   
    ! *** Specify smoothing constraints on each element
    ! set some default values in case smoothing not specified
    sobelov_wts(0,:) = 1.0_dp 
    sobelov_wts(1,:) = 1.0_dp !the scaling factor for the Sobolev weights
    !  The 5 weights on derivs wrt Xi_1/_11/_2/_22/'_12 are:
    sobelov_wts(2,:) = 1.0e-2_dp !weight for deriv wrt Xi_1
    sobelov_wts(3,:) = 0.4_dp
    sobelov_wts(4,:) = 1.0e-2_dp
    sobelov_wts(5,:) = 0.4_dp
    sobelov_wts(6,:) = 0.8_dp

    if(index(fitting_file, ".ipmap")> 0) then !full filename is given
       readfile = fitting_file
    else ! need to append the correct filename extension
       readfile = trim(fitting_file)//'.ipmap'
    endif
    open(IPFILE, file = readfile, status='old')

!!! read the element list for fitting, the fixed node locations and/or derivatives,
!!! and the nodal derivative mapping for versions of nodes. Node locations for
!!! multiple versions are assumed to be mapped to version 1.
    read(unit=IPFILE, fmt="(a)", iostat=ierror) string
    if(index(string, "Element group")> 0) then
       n_group = 0
       read_element_groups : do
          n_group = n_group + 1
          ibeg = index(string,"(")+1 
          iend = index(string,")")-1 
          sub_string = adjustl(string(ibeg:iend)) ! get the items between brackets
          read(sub_string, fmt=*, iostat=ierror) temp
          elem_group_names(n_group) = trim(temp)
          ibeg = index(string,"{")+1 
          iend = index(string,"}")-1 
          sub_string = string(ibeg:iend) ! get the items between brackets
          nelem = 0
          do while (index(trim(sub_string), ' ').ne.0)
             idx = index(trim(sub_string), ' ')
             nelem = nelem + 1
             read(sub_string(1:idx-1),*,iostat=stat) nelem_groups(n_group,nelem)
             sub_string = trim(sub_string(idx+1:iend))
          enddo
          nelem = nelem + 1
          read(sub_string(1:iend),*,iostat=stat) nelem_groups(n_group,nelem)
          read(unit=IPFILE, fmt="(a)", iostat=ierror) string
          if(index(string, "Element group") == 0) exit ! move to fixed nodes
       enddo read_element_groups
    endif

    if(index(string, "Fixed node")> 0) then
       read_fixed_nodes : do
          ibeg = index(string,"(")+1 
          iend = index(string,")")-1 
          sub_string = adjustl(string(ibeg:iend)) ! get the characters between brackets
          read(sub_string, fmt=*, iostat=ierror) np_global
          ibeg = index(string,"{")+1 
          iend = index(string,"}")-1 
          sub_string = adjustl(string(ibeg:iend)) ! get the characters between brackets
          read(sub_string, fmt=*, iostat=ierror) nv_fix, nk
          np = get_local_node_f(2,np_global)
          nk = nk+1 !read in 0 for coordinate, 1 for 1st deriv, 2 for 2nd deriv
          if(nv_fix.eq.0)then ! do for all versions
             do nv = 1,node_versn_2d(np)
                do nh = 1,3
                   ny = nynp(nk,nv,nh,np)
                   fix_bcs(ny) = .true.
                enddo !nh
             enddo !nv
          else
             nv = nv_fix
             do nh = 1,3
                ny = nynp(nk,nv,nh,np)
                fix_bcs(ny) = .true.
                !fit_soln(nk,nv,nh,np) = 0.0_dp ! shouldn't be zero
             enddo !nh
          endif
          read(unit=IPFILE, fmt="(a)", iostat=ierror) string
          if(index(string, "Fixed node") == 0) exit ! move to fixed nodes
       enddo read_fixed_nodes
    endif
       
    number_of_maps = 0
    nmap_info = 0
    if(index(string, "Mapped node")> 0) then
       read_mapped_nodes : do
          number_of_maps = number_of_maps + 1
          ibeg = index(string,"(")+1 
          iend = index(string,")")-1 
          sub_string = adjustl(string(ibeg:iend)) ! get the characters between brackets
          read(sub_string, fmt=*, iostat=ierror) nmap_info(1,number_of_maps)
          ibeg = index(string,"{")+1 
          iend = index(string,"}")-1 
          sub_string = adjustl(string(ibeg:iend)) ! get the characters between brackets
          read(sub_string, fmt=*, iostat=ierror) nmap_info(2:7,number_of_maps)
          read(unit=IPFILE, fmt="(a)", iostat=ierror) string
          if(index(string, "Mapped node") == 0) exit  ! go to the enxt section
       enddo read_mapped_nodes
    endif

    np_list_redist = 0
    num_redists = 0
    if(index(string, "Redistribute node")> 0) then
       read_redistribute_nodes : do
          ibeg = index(string,"{")+1 
          iend = index(string,"}")-1 
          sub_string = adjustl(string(ibeg:iend)) ! get the characters between brackets
          num_redists = num_redists + 1
          read(sub_string, fmt=*, iostat=ierror) np_list_redist(num_redists,:)
          read(unit=IPFILE, fmt="(a)", iostat=ierror) string
          if(index(string, "Redistribute node") == 0) exit  ! go to next section
       enddo read_redistribute_nodes
    endif

    if(index(string, "Sobelov weights")> 0) then
       read_smoothing : do
          ibeg = index(string,"(")+1 
          iend = index(string,")")-1 
          sub_string = adjustl(string(ibeg:iend)) ! get the characters between brackets
          read(sub_string, fmt=*, iostat=ierror) ne
          ibeg = index(string,"{")+1 
          iend = index(string,"}")-1 
          sub_string = adjustl(string(ibeg:iend)) ! get the characters between brackets
          read(sub_string, fmt=*, iostat=ierror) temp_weights(1:7)
          if(ne.eq.0)then
             forall(i = 0:6) sobelov_wts(i,1:num_elems_2d) = temp_weights(i+1)
          else
             ne = get_local_elem_2d(ne)
             forall(i = 0:6) sobelov_wts(i,ne) = temp_weights(i+1)
          endif
          read(unit=IPFILE, fmt="(a)", iostat=ierror) string
          if(index(string, "Sobelov weights") == 0) exit  ! end of file
       enddo read_smoothing
    endif

    close(IPFILE)
    
    call map_versions(nmap_info,number_of_maps,num_depvar,nynp,nyny,cyny,fit_soln,fix_bcs)
    
    ! fix ALL of the cross derivatives, and set to zero
    nk = 4
    do np = 1,num_nodes_2d
       do nv = 1,node_versn_2d(np)
          do nj = 1,num_coords
             ny = nynp(4,nv,nj,np)
             fix_bcs(ny) = .TRUE.
             node_xyz_2d(nk,nv,nj,np) = 0.0_dp
          enddo !nj
       enddo !nv
    enddo !np

    deallocate(nmap_info)
    
    call enter_exit(sub_name,2)

  end subroutine define_geometry_fit
  
!!! ##########################################################################      
  
  subroutine initialise_fit_mesh()
    !*initialise_fit_mesh:* scale and translate the mesh to align with a data
    ! cloud. uses the centre of mass and the range of data coordinates.

    ! Local variables
    integer :: i
    real(dp) :: datacofm(3),meshcofm(3),datarange(3),meshrange(3), &
         movemesh(3),scalemesh(3)
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'initialise_fit_mesh'
    call enter_exit(sub_name,1)

    do i = 1,3
       datacofm(i) = sum(data_xyz(i,:))/real(num_data)
       datarange(i) = maxval(data_xyz(i,:)) - minval(data_xyz(i,:))
       meshcofm(i) = sum(node_xyz_2d(1,1,i,:))/real(num_nodes_2d)
       meshrange(i) = maxval(node_xyz_2d(1,1,i,:)) - minval(node_xyz_2d(1,1,i,:))
    enddo

    scalemesh = datarange/meshrange
    movemesh = datacofm - meshcofm * scalemesh

    forall (i=1:3) node_xyz_2d(1,:,i,:) = node_xyz_2d(1,:,i,:) * &
         scalemesh(i) + movemesh(i)

    call enter_exit(sub_name,2)
    
  end subroutine initialise_fit_mesh
  
!!! ##########################################################################      
  
  subroutine gauss1(PG)
    
!!! dummy arguments
    real(dp) :: PG(:,:,:)
!!! local variables
    integer :: I,J,ng,nk,nn,ns,nu
    real(dp) :: D(3),XI(3),XIGG(3,3,2)
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'gauss1'
    call enter_exit(sub_name,1)
    
    D = [-0.3872983346207410_dp, 0.0_dp, 0.3872983346207410_dp]
    
    do j=1,3
       do i=1,3
          XIGG(i,j,1) = 0.50_dp+D(i)
          XIGG(i,j,2) = 0.50_dp+D(j)
          ng = I+(J-1)*3
          XI(1:2) = XIGG(i,j,1:2)
          ns=0
          do nn=1,num_elem_nodes !number of local nodes
             do nk=1,num_deriv !number of derivatives at node
                ns=ns+1 
                do nu=1,6 !number of xi coord derivs at node
                   PG(ns,nu,ng) = PSI1(nu,nk,nn,XI)
                enddo !nu
             enddo !nk
          enddo !nn
       enddo !i
    enddo !j

    call enter_exit(sub_name,2)

  end subroutine gauss1
  
!!! ##########################################################################      
  
  subroutine globalf(nony,not_1,not_2,npny,nyno,nynp,nyny,cony,cyno,cyny,fix_bcs)

!!! calculates the mapping arrays nyno/nony/cyno/cony

!!! dummy arguments
    integer :: nony(0:,:),not_1,not_2,npny(:,:),nyno(0:,:,:),nynp(:,:,:,:),nyny(0:,:)
    real(dp) :: cony(:),cyno(0:,:),cyny(0:,:)
    logical :: fix_bcs(:)
!!! local variables
    integer :: nh,nv,nk,no,no_tot(2),np,nrc,ny,nyo,nyr,nyr2,nyy2(2),ny2
    real(dp) :: RATIO
    logical :: done
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'globalf'
    call enter_exit(sub_name,1)
    
!!!***  Initialise mapping arrays
    nony = 0
    cony = 0.0_dp
    nyno = 0
    cyno = 0.0_dp
    no_tot = 0
    
!!!*** Calculate mapping arrays
    do np = 1,num_nodes_2d
       do nh = 1,num_fit
          do nv = 1,node_versn_2d(np)
             do nk = 1,num_deriv
                ny = nynp(nk,nv,nh,np)
                if(.not.fix_bcs(ny)) then 
                   ! variable needs to be solved for
                   done = .false.
                   ny2 = ny !the default
                   if(nyny(0,ny).ne.0) then ! a special mapping
                      ny2 = nyny(1,ny)
                      RATIO = cyny(1,ny) ! the weighting of the mapping
                   endif
                   if(ny2.ne.ny) then ! for mapping, where ny exists
                      done = .true.
                      if(ny2.eq.0) then !dof not used
                         fix_bcs(ny) = .TRUE.
                      else if(fix_bcs(ny2)) then
                         fix_bcs(ny) = .TRUE.
                      else            ! no mapping
                         do nrc = 1,1 ! nrc = 1,2 local row and local column
                            nyr = ny
                            nony(0,ny) = 1
                            no = nony(1,ny2)
                            nony(1,ny) = no
                            cony(ny) = RATIO*cony(ny2)
                            nyo = nyno(0,no,nrc)+1
                            nyno(0,no,nrc) = nyo
                            nyno(nyo,no,nrc) = ny
                            cyno(nyo,no) = RATIO*cony(ny2)
                         enddo !nrc
                      endif !ny2=0/fix_bcs
                   endif !ny.NE.ny2
                   
                   if(.not.done) then
                      do nrc=1,1 !rows and columns
                         no_tot(nrc) = no_tot(nrc)+1
                         nony(0,ny) = 1
                         nony(1,ny) = no_tot(nrc)
                         cony(ny) = 1.0_dp
                         nyno(0,no_tot(nrc),nrc) = 1
                         nyno(1,no_tot(nrc),nrc) = ny
                         cyno(0,no_tot(nrc)) = 0.0_dp
                         cyno(1,no_tot(nrc)) = 1.0_dp
                      enddo !nrc
                   endif !not done
                endif !fix
             enddo !nk
          enddo !nv
       enddo !nh
    enddo !np
    
    NOT_1 = no_tot(1)
    NOT_2 = no_tot(1)
    nyno(:,:,2) = nyno(:,:,1)

    call enter_exit(sub_name,2)

  end subroutine globalf

!!! ##########################################################################      

  subroutine list_data_error(data_on_elem,ndata_on_elem,data_xi)

!!! calculate and write out the RMS error for distance between data points
!!! and 2d mesh surface

!!! dummy arguments
    integer :: data_on_elem(:,:),ndata_on_elem(:) 
    real(dp) :: data_xi(:,:)
!!! local variables
    integer elem,nd,nde,num_data_infit,ne,nj
    real(dp) :: data_xi_local(2),EDD,elem_xyz(num_deriv_elem,num_coords), &
         SAED,SMED,SUM,SQED,X(6)
         
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'list_data_error'
    call enter_exit(sub_name,1)
    
    SMED=0.0_dp
    SAED=0.0_dp
    SQED=0.0_dp
    num_data_infit=0
    
    do ne=1,num_elems_2d
       call node_to_local_elem(ne,elem_xyz)
       elem=ne
       do nde=1,ndata_on_elem(elem) !for each data point on element
          nd=data_on_elem(elem,nde) !the data point number
          data_xi_local(1:2) = data_xi(1:2,nd)
          do nj=1,num_coords 
             X(nj)=PXI(1,data_xi_local,elem_xyz(1,nj))
          enddo
          SUM=0.0_dp
          do nj=1,num_coords
             SUM=SUM+(X(nj)-data_xyz(nj,nd))**2
          enddo !nj
          EDD = sqrt(SUM)  ! distance of the point from the surface
          SMED=SMED+EDD
          SAED=SAED+abs(EDD)
          SQED=SQED+EDD**2
          num_data_infit=num_data_infit+1
       enddo !nde
    enddo !list of elements
    
    if(num_data_infit.GT.1) then
       write(*,'('' Number of data points in fit ='',I8)') num_data_infit
       write(*,'('' Average absolute error  : '',D12.6,'' +/- '',D12.6)') &
            SAED/real(num_data_infit),sqrt((SQED-SAED**2/real(num_data_infit))/ &
            real(num_data_infit-1))
       write(*,'('' Root mean squared error : '',D12.6)') &
            sqrt(SQED/DBLE(num_data_infit))
    else
       WRITE(*,'('' No data points in any elements'')')
       stop
    endif !ndtot>1
    
    call enter_exit(sub_name,2)

  end subroutine list_data_error
  
!!! ##########################################################################      

  subroutine list_data_error_by_group(data_on_elem,ndata_on_elem,data_xi,writefile)

!!! calculate and write out the RMS error for distance between data points
!!! and 2d mesh surface

!!! dummy arguments
    integer :: data_on_elem(:,:),ndata_on_elem(:) 
    real(dp) :: data_xi(:,:)
    logical :: writefile
!!! local variables
    integer elem,i,nd,nde,ngroup,num_data_infit,ne,nj
    real(dp) :: data_xi_local(2),EDD,elem_xyz(num_deriv_elem,num_coords), &
         SAED,SMED,SUM,SQED,X(6)
         
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'list_data_error'
    call enter_exit(sub_name,1)

    if(writefile)then
       open(17, file = 'fitting.err', status='replace')
    endif
    do ngroup = 1,num_groups
       SMED=0.0_dp
       SAED=0.0_dp
       SQED=0.0_dp
       num_data_infit=0
    
       do i = 1,count(nelem_groups(ngroup,:)/=0)
          ne = get_local_elem_2d(nelem_groups(ngroup,i))
          call node_to_local_elem(ne,elem_xyz)
          elem = ne
          do nde = 1,ndata_on_elem(elem) !for each data point on element
             nd = data_on_elem(elem,nde) !the data point number
             data_xi_local(1:2) = data_xi(1:2,nd)
             do nj=1,num_coords 
                X(nj)=PXI(1,data_xi_local,elem_xyz(1,nj))
             enddo
             SUM=0.0_dp
             do nj=1,num_coords
                SUM=SUM+(X(nj)-data_xyz(nj,nd))**2
             enddo !nj
             EDD = sqrt(SUM)  ! distance of the point from the surface
             SMED=SMED+EDD
             SAED=SAED+abs(EDD)
             SQED=SQED+EDD**2
             num_data_infit = num_data_infit+1
          enddo !nde
       enddo !list of elements
       if(num_data_infit.GT.1) then
          if(writefile)then
             write(17,'('' Group'',a13,'' (n='',i5,''): err_av='',f6.2,'' +/-'',f6.2,'' mm;''&
                  &'' RMS='',f6.2,'' mm'')') trim(elem_group_names(ngroup)), &
                  num_data_infit, SAED/real(num_data_infit), &
                  sqrt((SQED-SAED**2/real(num_data_infit))/ real(num_data_infit-1)), &
                  sqrt(SQED/DBLE(num_data_infit))
          endif
          write(*,'('' Group'',a13,'' (n='',i5,''): err_av='',f6.2,'' +/-'',f6.2,'' mm;''&
               &'' RMS='',f6.2,'' mm'')') trim(elem_group_names(ngroup)), &
               num_data_infit, SAED/real(num_data_infit), &
               sqrt((SQED-SAED**2/real(num_data_infit))/ real(num_data_infit-1)), &
               sqrt(SQED/DBLE(num_data_infit))
       else
          if(writefile)then
             write(17,'('' Group'',a13,'' (n='',i5,''): no data points in any elements'')') &
                  trim(elem_group_names(ngroup)), num_data_infit
          endif
          write(*,'('' Group'',a13,'' (n='',i5,''): no data points in any elements'')') &
               trim(elem_group_names(ngroup)), num_data_infit
          !WRITE(*,'('' No data points in any elements'')')
          !stop
       endif !ndtot>1
    enddo !ngroup

    if(writefile) close(17)
    
    call enter_exit(sub_name,2)

  end subroutine list_data_error_by_group
  
!!! ##########################################################################      
  
  subroutine map_versions(nmap_info,number_of_maps,num_depvar,nynp,nyny, &
       cyny,fit_soln,fix_bcs)

!!! dummy arguments
    integer, intent(in) :: nmap_info(:,:),num_depvar
    integer :: nynp(:,:,:,:)
    integer,allocatable :: nyny(:,:)
    real(dp),allocatable :: cyny(:,:)
    real(dp) :: fit_soln(:,:,:,:)
    logical :: fix_bcs(:)
!!! local variables
    integer :: i,ibeg,iend,ierror,i_ss_end,nj,node,np,number_of_maps,nv, &
         NV_MAX,ny,nk_t,nv_t,nj_t,np_t,nk_m,nv_m,nj_m,np_m, &
         ny_t
    real(dp) :: r_map_coef
    character(len=132) :: string
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'map_versions'
    call enter_exit(sub_name,1)

!!! fix the boundary conditions for coordinates for nodes with versions, such that
!!! versions higher than 1 map to version 1     
    do np = 1,num_nodes_2d
       if(node_versn_2d(np).gt.1)then
          do nv = 2,node_versn_2d(np)
             do nj = 1,num_coords
                node_xyz_2d(1,nv,nj,np) = node_xyz_2d(1,1,nj,np)
                fit_soln(1,nv,nj,np) = fit_soln(1,1,nj,np)
                ny = nynp(1,nv,nj,np)
                fix_bcs(ny) = .TRUE.
             enddo !nj
          enddo !nv
       endif
    enddo !node
    
!!! read in the following for mapping:
!!!    node, version, derivative >> node, version, derivative, mapping coefficient

    ! allocate memory for dependent variable mapping arrays
    if(.not.allocated(cyny)) allocate(cyny(0:number_of_maps,num_depvar)) 
    if(.not.allocated(nyny)) allocate(nyny(0:number_of_maps,num_depvar))
    nyny = 0       ! initialise depvar to depvar mapping
    cyny = 0.0_dp  ! initialise weighting for mappings
    
!!! note that the global node numbers are used in the mapping file, whereas we need to use
!!! local numbering for the computation. Read in as global and then map to local below.
    do node = 1,number_of_maps ! for the number of nodes with mappings
       do nj = 1,num_coords
          nk_m = nmap_info(3,node)+1 !derivative
          nv_m = nmap_info(2,node) !version
          nj_m = nj !coordinate
          np_m = get_local_node_f(2,nmap_info(1,node)) !global node mapped to local node
          ny = nynp(nk_m,nv_m,nj_m,np_m)
          nk_t = nmap_info(6,node)+1 !derivative
          nv_t = nmap_info(5,node) !version
          nj_t = nj !coordinate
          np_t = get_local_node_f(2,nmap_info(4,node)) !global node mapped to local node
          ny_t = nynp(nk_t,nv_t,nj_t,np_t)
          r_map_coef = REAL(nmap_info(7,node)) !mapping coefficient, +1 or -1
          
          if(ny.gt.0) then
             nyny(0,ny) = nyny(0,ny)+1 ! increment array size
             nyny(nyny(0,ny),ny) = ny_t
             cyny(0,ny) = 0.0_dp
             cyny(nyny(0,ny),ny) = r_map_coef
             node_xyz_2d(nk_m,nv_m,nj_m,np_m) = node_xyz_2d(nk_t,nv_t,nj_t,np_t)*r_map_coef
             fit_soln(nk_m,nv_m,nj_m,np_m) = node_xyz_2d(nk_t,nv_t,nj_t,np_t)*r_map_coef
          endif ! ny.GT.0
       enddo !nj
    enddo
    
!!! make all coordinates for different versions be the same
    do np=1,num_nodes_2d
       do nj=1,3
          ny_t=nynp(1,1,nj,np)
          do nv=2,node_versn_2d(np) !for each version
             ny=nynp(1,nv,nj,np)
             if(ny > 0) then
                nyny(0,ny)=nyny(0,ny)+1 ! increment array size
                nyny(nyny(0,ny),ny)=ny_t
                cyny(0,ny)=0.0_dp
                cyny(nyny(0,ny),ny)=1.0_dp
                node_xyz_2d(1,nv,nj,np)=node_xyz_2d(1,1,nj,np)
                fit_soln(1,nv,nj,np)=node_xyz_2d(1,1,nj,np)
             endif ! ny.GT.0
          enddo !nv
       enddo !nj
    enddo

    call enter_exit(sub_name,2)

  end subroutine map_versions
  
!!! ##########################################################################

  subroutine update_versions(nynp,fit_soln,fix_bcs)

    integer :: nynp(:,:,:,:)  
    real(dp) :: fit_soln(:,:,:,:)
    logical :: fix_bcs(:)
    ! Local variables
    integer :: nj,nk,np,nv,ny
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'update_versions'
    call enter_exit(sub_name,1)
    
!!! fix the boundary conditions for coordinates for nodes with versions, such that
!!! versions higher than 1 map to version 1     
    do np = 1,num_nodes_2d
       if(node_versn_2d(np) > 1)then
          do nv = 2,node_versn_2d(np)
             do nj = 1,3
                node_xyz_2d(1,nv,nj,np) = node_xyz_2d(1,1,nj,np)
                fit_soln(1,nv,nj,np) = fit_soln(1,1,nj,np)
                ny = nynp(1,nv,nj,np)
                fix_bcs(ny) = .TRUE.
             enddo !nj
          enddo !nv
       endif
    enddo !node

!!! fix ALL of the cross derivatives, and set to zero
    do np = 1,num_nodes_2d
       nk = 4 !index for 1-2 cross-derivative
       do nv = 1,node_versn_2d(np)
          do nj = 1,3
             ny = nynp(nk,nv,nj,np)
             fix_bcs(ny) = .TRUE.
             node_xyz_2d(nk,nv,nj,np) = 0.0_dp
          enddo !nj
       enddo !nv
    enddo !np

    call enter_exit(sub_name,2)
    
  end subroutine update_versions
  
!!! ##########################################################################      
  
  subroutine local_dof(n_dof,ne,ny_local,nynp)

    integer :: ny_local(:),n_dof,ne,nynp(:,:,:,:)
!!! local variables
    integer nh,nk,nn,np,nv
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'local_dof'
    call enter_exit(sub_name,1)
    
    n_dof = 0
    do nh = 1,num_fit
       do nn = 1,num_elem_nodes !nodal variables
          np = elem_nodes_2d(nn,ne)
          nv = elem_versn_2d(nn,ne)
          do nk = 1,num_deriv
             n_dof = n_dof + 1
             ny_local(n_dof) = nynp(nk,nv,nh,np)
          enddo !nk
       enddo !nn
    enddo !nhj
    
    call enter_exit(sub_name,2)

  end subroutine local_dof

!!! ##########################################################################      

  subroutine distribute_surface_node_fit(np_list,nynp,fit_soln,fix_bcs)

    integer,intent(in) :: np_list(:,:),nynp(:,:,:,:)
    real(dp) :: fit_soln(:,:,:,:)
    logical,intent(in) :: fix_bcs(:)
    ! Local variables
    integer :: i,in_line,iredist,j,k,line_numbers(20),ne,nline,nlist,nn(4),node_2,np, &
         np1,np2,np_between(20),np_end,np_start,num_lines,nv,nv1,nv2,n_xi_dctn,ny
    real(dp) :: new_length,new_xyz(20,3),segment_length,sum_length,xi

    iredist = 1
    distribute_line : do
       if(np_list(iredist,1).eq.0) exit  distribute_line ! no more lines to distribute
       np_start = get_local_node_f(2,np_list(iredist,1))

       ! get the start and end nodes, and between nodes
       np_between = 0
       i = 1
       get_line_nodes : do
          i = i+1
          if(np_list(iredist,i).eq.0) exit get_line_nodes
          np_between(i-1) = get_local_node_f(2,np_list(iredist,i))
       enddo get_line_nodes
       nlist = i-2
       np_end = np_between(i-2)
       
       ! find the first element that has np_start and np_between(1)
       find_element : do i = 1,elems_at_node_2d(np_start,0)
          ne = elems_at_node_2d(np_start,i)
          if(inlist(np_between(1),elem_nodes_2d(1:4,ne))) exit find_element
       enddo find_element

       ! get the list of line segments
       num_lines = 0
       segment_length = 0.0_dp

       do i = 1,nlist
          np = np_between(i) ! the next 'between' node
          do j = 1,num_lines_2d
             nline = lines_2d(j)
             if((np_start.eq.nodes_in_line(2,1,nline).and.np.eq.nodes_in_line(3,1,nline)).or.&
                  (np_start.eq.nodes_in_line(3,1,nline).and.np.eq.nodes_in_line(2,1,nline)))then
                num_lines = num_lines + 1
                line_numbers(num_lines) = nline
                segment_length = segment_length + arclength(nline)
             endif
          enddo
          np_start = np
       enddo ! i
       segment_length = segment_length/real(num_lines)

       ! redistribute along the line segments
       do i = 1,nlist-1
          node_2 = np_between(i)
          nline = line_numbers(i)
          new_length = segment_length*real(i)
          sum_length = 0.0_dp
          check_lines: do j = 1,num_lines ! check each line segment
             if(new_length.gt.sum_length.and.new_length.le. &
                  sum_length+arclength(line_numbers(j)))then
                ! in this segment
                in_line = line_numbers(j)
                np1 = nodes_in_line(2,1,in_line)
                np2 = nodes_in_line(3,1,in_line)
                nv1 = line_versn_2d(1,1,in_line)
                nv2 = line_versn_2d(2,1,in_line)
                n_xi_dctn = nodes_in_line(1,0,nline)
                if(i.gt.1.and.j.gt.1)then
                   if(np1.ne.nodes_in_line(3,1,line_numbers(j-1))) then
                      xi = 1.0_dp - (new_length-sum_length)/arclength(line_numbers(j))
                   else
                      xi = (new_length-sum_length)/arclength(line_numbers(j))
                   endif
                endif
                exit check_lines
             else
                sum_length = sum_length+arclength(line_numbers(j))
             endif
          enddo check_lines

          ! x(xi) = phi_10*x1 + phi_20*x2 + phi_11*x1' + phi_21*x2'
          new_xyz(i,:) = hermite(1,1,1,xi)*node_xyz_2d(1,nv1,:,np1) + &
               hermite(2,1,1,xi)*node_xyz_2d(1,nv2,:,np2) + &
               hermite(1,2,1,xi)*node_xyz_2d(n_xi_dctn+1,nv1,:,np1) + &
               hermite(2,2,1,xi)*node_xyz_2d(n_xi_dctn+1,nv2,:,np2)
       enddo ! i

       ! update the node coordinates
       do i = 1,num_lines-1
          node_2 = np_between(i)
          node_xyz_2d(1,1,:,node_2) = new_xyz(i,:)
          do nv = 1,node_versn_2d(np)
             do j = 1,3
                node_xyz_2d(1,nv,j,node_2) = node_xyz_2d(1,1,j,node_2)
                ny = nynp(1,nv,j,node_2)
                fit_soln(1,nv,j,node_2) = node_xyz_2d(1,nv,j,node_2)
             enddo ! j
          enddo ! nv
       enddo

       iredist = iredist + 1
       
    enddo distribute_line
    call calc_arclengths
    
  end subroutine distribute_surface_node_fit
  
!!! ##########################################################################      

  subroutine distribute_nodes_between(gnode_1,gnode_2,n_xi_dctn)
    !*distribute_nodes_between*: update the location of all nodes on the line
    ! in the n_xi_dctn direction between node_1 and node_2, so that they are
    ! uniformly spread out wrt arclength.
    
    integer,intent(in) :: gnode_1,gnode_2,n_xi_dctn
    ! Local variables
    integer :: count_checks,i,in_line,j,ne_check,nline,nn,node_1,node_2,np1,np2, &
         nthline,num_lines,nv1,nv2
    integer,allocatable :: line_numbers(:)
    real(dp) :: new_length,new_xyz(100,3),segment_length,sum_length, &
         total_length,xi
    logical :: carry_on
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'distribute_nodes_between'
    call enter_exit(sub_name,1)

    allocate(line_numbers(num_lines_2d))
    line_numbers = 0
    num_lines = 0
    count_checks = 0
    carry_on = .true.

    node_1 = get_local_node_f(2,gnode_1)
    node_2 = get_local_node_f(2,gnode_2)

    if(elems_at_node_2d(node_1,0).eq.0.or.elems_at_node_2d(node_2,0).eq.0) &
         ! one or both of the nodes is not in an element
         carry_on = .false.

    ! get a list of all nodes that are on the line between node_1 and node_2
    if(carry_on)then
       ne_check = elems_at_node_2d(node_1,1)
       do i = 1,4
          if(node_1.eq.elem_nodes_2d(i,ne_check)) nn = i
       enddo
       if(nn.eq.1.and.n_xi_dctn.eq.1)then
          nthline = 1
          num_lines = num_lines + 1
          line_numbers(num_lines) = elem_lines_2d(nthline,ne_check)
       else if(nn.eq.1.and.n_xi_dctn.eq.2)then
          nthline = 3
          num_lines = num_lines + 1
          line_numbers(num_lines) = elem_lines_2d(nthline,ne_check)
       else if(nn.eq.2.and.n_xi_dctn.eq.1)then
          nthline = 1
       else if(nn.eq.2.and.n_xi_dctn.eq.2)then
          nthline = 4
          num_lines = num_lines + 1
          line_numbers(num_lines) = elem_lines_2d(nthline,ne_check)
       else if(nn.eq.3.and.n_xi_dctn.eq.1)then
          nthline = 2
          num_lines = num_lines + 1
          line_numbers(num_lines) = elem_lines_2d(nthline,ne_check)
       else if(nn.eq.3.and.n_xi_dctn.eq.2)then
          nthline = 3
       else if(nn.eq.4.and.n_xi_dctn.eq.1)then
          nthline = 2
       else if(nn.eq.4.and.n_xi_dctn.eq.2)then
          nthline = 4
       endif

       if(num_lines.eq.1)then ! check that we don't just have one line!
          nline = line_numbers(num_lines)
          if(nodes_in_line(3,1,nline).eq.node_2) carry_on = .false.
       endif
       
       ne_check = elem_cnct_2d(n_xi_dctn,1,ne_check)
       
       do while(carry_on)
          nline = elem_lines_2d(nthline,ne_check)
          if(nodes_in_line(3,1,nline).eq.node_2)then
             num_lines = num_lines + 1
             line_numbers(num_lines) = elem_lines_2d(nthline,ne_check)
             carry_on = .false.
          else
             ! for collapsed elements, go to the next one
             do while (nodes_in_line(2,1,nline).eq.nodes_in_line(3,1,nline))
                ne_check = elem_cnct_2d(n_xi_dctn,1,ne_check)
                nline = elem_lines_2d(nthline,ne_check)
             enddo
             num_lines = num_lines + 1
             line_numbers(num_lines) = elem_lines_2d(nthline,ne_check)
             count_checks = count_checks + 1
             if(count_checks.gt.num_elems_2d) carry_on = .false.
             ne_check = elem_cnct_2d(n_xi_dctn,1,ne_check)
          endif
       enddo ! while
    endif

    total_length = 0.0_dp
    do i = 1,num_lines
       nline = line_numbers(i)
       total_length = total_length + arclength(nline)
    enddo ! nline
    segment_length = total_length/real(num_lines)

    do i = 1,num_lines-1
       nline = line_numbers(i)
       node_2 = nodes_in_line(3,1,nline)
       new_length = segment_length*real(i)
       sum_length = 0.0_dp
       check_lines: do j = 1,num_lines ! check each line segment
          if(new_length.gt.sum_length.and.new_length.le. &
               sum_length+arclength(line_numbers(j)))then
             ! in this segment
             in_line = line_numbers(j)
             xi = (new_length-sum_length)/arclength(line_numbers(j))
             np1 = nodes_in_line(2,1,in_line)
             np2 = nodes_in_line(3,1,in_line)
             nv1 = line_versn_2d(1,1,in_line)
             nv2 = line_versn_2d(2,1,in_line)
             exit check_lines
          else
             sum_length = sum_length+arclength(line_numbers(j))
          endif
       enddo check_lines

       ! x(xi) = phi_10*x1 + phi_20*x2 + phi_11*x1' + phi_21*x2'
       new_xyz(i,:) = hermite(1,1,1,xi)*node_xyz_2d(1,nv1,:,np1) + &
            hermite(2,1,1,xi)*node_xyz_2d(1,nv2,:,np2) + &
            hermite(1,2,1,xi)*node_xyz_2d(n_xi_dctn+1,nv1,:,np1) + &
            hermite(2,2,1,xi)*node_xyz_2d(n_xi_dctn+1,nv2,:,np2)
          
    enddo

    do i = 1,num_lines-1
       nline = line_numbers(i)
       node_2 = nodes_in_line(3,1,nline)
       node_xyz_2d(1,1,:,node_2) = new_xyz(i,:)
       forall (j=1:6) node_xyz_2d(1,j,:,node_2) = node_xyz_2d(1,1,:,node_2)
    enddo
    
    deallocate(line_numbers)
    
    call enter_exit(sub_name,2)
    
  end subroutine distribute_nodes_between

!!! ##########################################################################      

  subroutine centre_a_node(update_node,n_xi_dctn)
    !*centre_a_node*: updates the coordinates of the given node so that it sits
    ! at xi between the two adjacent nodes in the n_xi_dctn direction

    integer,intent(in) :: n_xi_dctn,update_node
    ! Local variables
    integer :: i,iline_1,iline_2,nj,nline,np,np1,np2,nv,nv1,nv2
    real(dp) :: line_length,line_xyz(2,3,2),location = 0.5_dp, new_xyz(3),xi_on_line
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'centre_a_node'
    call enter_exit(sub_name,1)

    np = get_local_node_f(2,update_node)
    line_length = 0.0_dp
    
    do i = 1,num_lines_2d
       nline = lines_2d(i)
       if(nodes_in_line(1,0,nline).eq.n_xi_dctn)then ! only check lines in the right Xi direction
          if(nodes_in_line(3,1,nline).eq.np.or. &
               nodes_in_line(2,1,nline).eq.np)then
             line_length = line_length + arclength(i) ! calculated the total arclength
             if(nodes_in_line(3,1,nline).eq.np) iline_1 = i  ! np is the second node
             if(nodes_in_line(2,1,nline).eq.np) iline_2 = i
          endif
       endif
    enddo

    line_length = line_length * location ! get the arclength for the adjusted node
    if(line_length.le.arclength(iline_1))then
       xi_on_line = line_length/arclength(iline_1)
       nline = iline_1
    else
       xi_on_line = (arclength(iline_2)-line_length)/arclength(iline_2)
       nline = iline_2
    endif

    np1 = nodes_in_line(2,1,iline_1)
    np2 = nodes_in_line(3,1,iline_1)
    nv1 = line_versn_2d(1,1,iline_1)
    nv2 = line_versn_2d(2,1,iline_1)
    line_xyz(1,:,1) = node_xyz_2d(1,1,:,np1)
    line_xyz(1,:,2) = node_xyz_2d(1,1,:,np2)
    line_xyz(2,:,1) = node_xyz_2d(n_xi_dctn+1,nv1,:,np1)
    line_xyz(2,:,2) = node_xyz_2d(n_xi_dctn+1,nv2,:,np2)
    
    ! x(xi) = phi_10*x1 + phi_20*x2 + phi_11*x1' + phi_21*x2'
    new_xyz(:) = hermite(1,1,1,xi_on_line)*line_xyz(1,:,1) + &
         hermite(2,1,1,xi_on_line)*line_xyz(1,:,2) + &
         hermite(1,2,1,xi_on_line)*line_xyz(2,:,1) + &
         hermite(2,2,1,xi_on_line)*line_xyz(2,:,2)

    do nv = 1,node_versn_2d(np)
       node_xyz_2d(1,nv,:,np) = new_xyz(:)
    enddo
    write(*,'(''Node adjusted = '',3(f9.2))') node_xyz_2d(1,1,:,np)

    call enter_exit(sub_name,2)

  end subroutine centre_a_node
  
!!! ##########################################################################      
  
  subroutine set_linear_derivatives

    ! Local variables
    integer :: ne,nk,nl,nline,np1,np2,nv1,nv2

    do ne = 1,num_elems_2d
       np1 = elem_nodes_2d(1,ne)
       nv1 = elem_versn_2d(1,ne)
       np2 = elem_nodes_2d(2,ne)
       nv2 = elem_versn_2d(2,ne)
       node_xyz_2d(2,nv1,:,np1) = node_xyz_2d(1,1,:,np2) &
            - node_xyz_2d(1,1,:,np1)

       np2 = elem_nodes_2d(3,ne)
       nv2 = elem_versn_2d(3,ne)
       node_xyz_2d(3,nv1,:,np1) = node_xyz_2d(1,1,:,np2) &
            - node_xyz_2d(1,1,:,np1)

       np1 = elem_nodes_2d(2,ne)
       nv1 = elem_versn_2d(2,ne)
       np2 = elem_nodes_2d(4,ne)
       nv2 = elem_versn_2d(4,ne)
       node_xyz_2d(3,nv1,:,np1) = node_xyz_2d(1,1,:,np2) &
            - node_xyz_2d(1,1,:,np1)

       np1 = elem_nodes_2d(3,ne)
       nv1 = elem_versn_2d(3,ne)
       node_xyz_2d(2,nv1,:,np1) = node_xyz_2d(1,1,:,np2) &
            - node_xyz_2d(1,1,:,np1)
    enddo !ne
    
  end subroutine set_linear_derivatives
  
!!! ##########################################################################      

  function psi1(nu,nk,nn,XI)
    
!!! dummy arguments
    integer nk,nn,nu
    real(dp) :: XI(:)
!!! local variables
    integer :: ido(num_deriv,2),inp(4,2),ipu(6,2),ni
    real(dp) :: psi1
  
    ipu = reshape([1,2,3,1,1,2,1,1,1,2,3,2],shape(ipu))
    inp = reshape([1,2,1,2,1,1,2,2],shape(inp)) ! the node index positions
    ido = reshape([1,2,1,2,1,1,2,2],shape(ido))
    
    psi1 = 1.0_dp
    do ni=1,2
       psi1 = psi1*hermite(inp(nn,ni),ido(nk,ni),ipu(nu,ni),xi(ni))
    enddo

  end function psi1

!!! ###############################################################

  function pxi(nu,xi,x)
    
!!! dummy arguments
    integer :: nu
    real(dp) :: xi(2),x(num_deriv_elem)
!!! local variables
    integer :: nk,nn,ns
    real(dp) :: pxi
    
    pxi = 0.0_dp
    ns = 0
    do nn=1,num_elem_nodes
       do nk=1,num_deriv
          ns = ns+1
          pxi = pxi + psi1(nu,nk,nn,xi)*x(ns)
       enddo
    enddo

  end function pxi
  
!!! ##########################################################################      

  subroutine update_scale_factor_norm
  
!!! local variables
    integer :: nj,nk,nk1,np,nv
    real(dp) :: SCALE,XD(3),ZERO_TOL=1.0e-12_dp
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'update_scale_factor_norm'
    call enter_exit(sub_name,1)
    
    do np=1,num_nodes_2d
       do nv=1,node_versn_2d(np)
          do nk1=1,2 !loop over 1st derivs
             nk = nk1+1
             XD(1:3) = node_xyz_2d(nk,nv,1:3,np)
             SCALE = XD(1)**2+XD(2)**2+XD(3)**2
             if(SCALE.gt.ZERO_TOL) then
                ! Normalise the nodal derivatives
                SCALE = 1.0_dp/DSQRT(SCALE)
                do nj=1,3
                   node_xyz_2d(nk,nv,nj,np) = node_xyz_2d(nk,nv,nj,np)*SCALE
                   ! adjust cross-derivatives
                   node_xyz_2d(4,nv,nj,np) = node_xyz_2d(4,nv,nj,np)*SCALE
                enddo ! nj
             endif !>0
          enddo ! nk1
       enddo !nv
    enddo !no_np
    
    call calc_scale_factors_2d('arcl')
    
    call enter_exit(sub_name,2)

  end subroutine update_scale_factor_norm
  
!!! ##########################################################################      

  subroutine node_to_local_elem(ne,elem_xyz)

!!! copies geometry information from nodes into a local element array
    
!!! dummy arguments
    integer,intent(in) :: ne
    real(dp) :: elem_xyz(:,:)
!!! local variablesK     Local Variables
    integer :: nj,nk,nn,np,ns,nv
    
    ! --------------------------------------------------------------------------

    do nj=1,3
       ns=0
       do nn=1,num_elem_nodes
          np=elem_nodes_2d(nn,ne)
          nv=elem_versn_2d(nn,ne)
          do nk=1,num_deriv
             ns=ns+1
             elem_xyz(ns,nj)=node_xyz_2d(nk,nv,nj,np)*scale_factors_2d(ns,ne)
          enddo
       enddo
    enddo !nj
    
  end subroutine node_to_local_elem

!!! ##########################################################################      

  subroutine make_element_matrices(ne,fit_soln_local,fit_soln, &
       data_on_elem,ndata_on_elem,data_xi,ER,ES,sobelov_wts)

!!!    Evaluates element rhs, ER(ns), in calculation of least squares
!!!    fit of linear field variables, defined by nodal values
!!!    node_xyz_2d(nk,nv,nj,np), to the set of data values data_xyz(nj,nd) with
!!!    weights data_weight(nj,nd) at local coordinate values data_xi(ni,nd).

!!!    ZDES evaluates element stiffness matrix ES(ms,ns) in calculation
!!!    of least squares fit of linear field variables, defined by nodal
!!!    values node_xyz_2d(nk,nv,nj,np), to the set of data values XD(nj,nd) with
!!!    weights data_weight(nj,nd) at local coordinate values data_xi(ni,nd), where
!!!    nj=NJO.

!!! dummy arguments
    integer :: data_on_elem(:,:),ndata_on_elem(:),ne
    real(dp) :: data_xi(:,:),ER(:),ES(:,:),sobelov_wts(0:,:),&
         fit_soln(:,:,:,:),fit_soln_local(:,:)
!!! local variables
    integer nd,nde,ng,nh,nh1,nh2,nhj1,nhj2,nhs1,nhs1_for_nhj1,nhs2,nk,nk1,nk2, &
         nn,nn1,nn2,np,ns,ns1,ns2,nu,nv
    real(dp) :: SUM1,SUM2,SUM3,SUM4,X,ZDL(3,nmax_data_elem)
    real(dp) :: PD(num_deriv_elem),PG(16,6,9),WG(9)
    real(dp),dimension(3,nmax_data_elem) :: WDL
    real(dp),dimension(2,nmax_data_elem) :: XIDL
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'make_element_matrices'
    call enter_exit(sub_name,1)
    
    WG = [7.7160493827160628e-2_dp, 0.12345679012345677_dp, 7.7160493827160628e-2_dp,&
         0.12345679012345677_dp, 0.19753086419753044_dp, 0.12345679012345677_dp,&
         7.7160493827160628e-2_dp, 0.12345679012345677_dp, 7.7160493827160628e-2_dp]

    ER = 0.0_dp
    ES = 0.0_dp

    call gauss1(PG)

    do nh = 1,num_fit
       ns = 0
       do nn = 1,num_elem_nodes
          np = elem_nodes_2d(nn,ne)
          nv = elem_versn_2d(nn,ne)
          do nk = 1,num_deriv
             ns = ns+1
             fit_soln_local(ns,nh) = fit_soln(nk,nv,nh,np)*scale_factors_2d(ns,ne)
          enddo !nk
       enddo !nn
    enddo !nhx
    
!!! evaluate the element matrix ER
    
    do nde = 1,ndata_on_elem(ne) ! for each data point on the element
       nd = data_on_elem(ne,nde) ! the data point number
       XIDL(1:2,nde) = data_xi(1:2,nd)
       ZDL(1:3,nde) = data_xyz(1:3,nd)
       WDL(1:3,nde) = data_weight(1:3,nd)
    enddo !nde
    
    nhs1=0
    do nh = 1,num_fit
       do nde = 1,ndata_on_elem(ne)
          X = PXI(1,XIDL(1:2,nde),fit_soln_local(1:num_deriv_elem,nh))
          ZDL(nh,nde) = ZDL(nh,nde)-X
       enddo !nde
       ns1 = 0
       do nn1 = 1,num_elem_nodes
          do nk1 = 1,num_deriv
             nhs1 = nhs1+1
             ns1 = ns1+1
             SUM1 = 0.0_dp
             do nde = 1,ndata_on_elem(ne)
                SUM1 = SUM1 + PSI1(1,nk1,nn1,XIDL(1:2,nde))*ZDL(nh,nde)*WDL(nh,nde)
             enddo !nde
             SUM2 = 0.0_dp
             do ng = 1,num_gauss
                SUM3 = 0.0_dp
                do nu = 2,6 !for 2d elements
                   SUM4 = 0.0_dp
                   do ns2 = 1,num_deriv_elem
                      SUM4 = SUM4+fit_soln_local(ns2,nh)*PG(ns2,nu,ng)
                   enddo !ns2
                   SUM3 = SUM3+SUM4*PG(ns1,nu,ng)*sobelov_wts(nu,ne) !*sobelov_wts(1,ne)
                enddo !nu
                SUM2 = SUM2-SUM3*WG(ng) !*RG(ng)
             enddo !ng
             ER(nhs1) = ER(nhs1)+(SUM1+SUM2*sobelov_wts(0,ne))*scale_factors_2d(ns1,ne)
          enddo !nk1
       enddo !nn1
    enddo !nhj1
    
!!! evaluate the element matrix ES
    
    ES = 0.0_dp
    nhs1 = 0
    ! for each of the 3 dependent variables to be fitted 
    do nhj1 = 1,num_fit !nhj are vars for the fit problem njj
       nh1 = nhj1
       nhs1_for_nhj1 = nhs1
       do nde = 1,ndata_on_elem(ne)
          nhs1 = nhs1_for_nhj1
          ns1 = 0
          do nn1 = 1,num_elem_nodes
             do nk1 = 1,num_deriv
                nhs1 = nhs1+1
                ns1 = ns1+1
                PD(ns1) = PSI1(1,nk1,nn1,XIDL(1:2,nde))
             enddo !nk1
          enddo !nn1
          nhs1 = nhs1_for_nhj1
          do ns1 = 1,num_deriv_elem
             nhs1 = nhs1+1
             nhs2 = 0
             do nhj2 = 1,num_fit !columns
                nh2 = nhj2
                do ns2 = 1,num_deriv_elem
                   nhs2 = nhs2+1
                   if(nhj2.EQ.nhj1) then !to avoid coupling for now
                      ES(nhs1,nhs2) = ES(nhs1,nhs2)+PD(ns1)*PD(ns2) &
                           *WDL(nh1,nde)*scale_factors_2d(ns1,ne)*scale_factors_2d(ns2,ne)
                   endif !nhj2=nhj1
                enddo !ns2
             enddo !nhj2
          enddo !ns1
       enddo !nde
       
       ns1 = 0
       nhs1 = nhs1_for_nhj1
       do nn1 = 1,num_elem_nodes
          do nk1 = 1,num_deriv
             nhs1 = nhs1+1
             ns1 = ns1+1
             nhs2 = 0
             do nhj2 = 1,num_fit !columns
                ns2 = 0
                do nn2 = 1,num_elem_nodes
                   do nk2 = 1,num_deriv
                      nhs2 = nhs2+1
                      ns2 = ns2+1
                      if(nhj2.EQ.nhj1) then !to avoid coupling for now
                         SUM2 = 0.0_dp
                         do ng = 1,num_gauss
                            SUM3 = 0.0_dp
                            do nu = 2,6
                               SUM3 = SUM3+ &
                                    PG(ns1,nu,ng)*PG(ns2,nu,ng)*sobelov_wts(nu,ne)
                            enddo !nu
                            SUM2 = SUM2+SUM3*WG(ng)
                         enddo !ng
                         ES(nhs1,nhs2) = ES(nhs1,nhs2)+(SUM2*sobelov_wts(0,ne))* &
                              scale_factors_2d(ns1,ne)*scale_factors_2d(ns2,ne)
                      endif !nhj2=nhj1
                   enddo !nk2
                enddo !nn2
             enddo !nhj2
          enddo !nk1
       enddo !nn1
    enddo !nhj1
    
    call enter_exit(sub_name,2)

  end subroutine make_element_matrices
  
!!! ##########################################################################      
  
  subroutine calculate_ny_maps(npny,num_depvar,nynp)
  
!!! dummy arguments
    integer :: npny(:,:),num_depvar,nynp(:,:,:,:)
!!! local variables
    integer nh,nk,np,nv,ny
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'calculate_ny_maps'
    call enter_exit(sub_name,1)
    
    !***  Initialise mapping arrays 
    nynp = 0
    npny = 0
    
    !***  Set up mapping arrays
    ny = 0
    do nh = 1,num_fit
       do np = 1,num_nodes_2d
          do nv = 1,node_versn_2d(np)
             do nk = 1,num_deriv
                ny = ny+1
                nynp(nk,nv,nh,np) = ny
                npny(1,ny) = nk
                npny(2,ny) = nv
                npny(3,ny) = nh
                npny(4,ny) = np
             enddo !nk
          enddo !nv
       enddo !np
    enddo !njj
    num_depvar = ny

    call enter_exit(sub_name,2)

  end subroutine calculate_ny_maps

!!!#############################################################################

  subroutine define_data_fit_group(datafile, groupname)
    !*define_data_fit_group:* reads data points from a file and associates with a named group

    character(len=*) :: datafile
    character(len=*) :: groupname
    ! Local variables
    integer :: iend,ierror,length_string,ncount,n_add_data,nj,itemp
    real(dp),allocatable :: temp_data_field(:,:),temp_data_xyz(:,:),temp_data_weight(:,:)
    character(len=132) :: buffer,readfile
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'define_data_fit_group'
    call enter_exit(sub_name,1)
    
    if(index(datafile, ".ipdata")> 0) then !full filename is given
       readfile = datafile
    else ! need to append the correct filename extension
       readfile = trim(datafile)//'.ipdata'
    endif

    open(10, file=readfile, status='old')
    read(unit=10, fmt="(a)", iostat=ierror) buffer

!!! first run through to count the number of data points
    !set the counted number of data points to zero
    ncount = 0
    read_line_to_count : do
       read(unit=10, fmt="(a)", iostat=ierror) buffer
       if(ierror<0) exit !ierror<0 means end of file
       ncount = ncount + 1
    end do read_line_to_count
    n_add_data = ncount
    close (10)
    if(.not.allocated(data_xyz))then ! first time reading, so allocate and initialise num_data
       allocate(data_field(4,n_add_data))
       allocate(data_xyz(3,n_add_data))
       allocate(data_weight(3,n_add_data))
       num_data = 0
       num_groups = 0
    endif

    num_groups = num_groups + 1
    data_group_names(num_groups) = trim(groupname)

    write(*,'('' Read'',I7,'' data points from file'')') n_add_data !num_data
!!! allocate arrays now that we know the size required
    if(num_data.gt.0)then
       allocate(temp_data_field(4,num_data))
       allocate(temp_data_xyz(3,num_data))
       allocate(temp_data_weight(3,num_data))
       temp_data_field(:,1:num_data) = data_field(:,1:num_data)
       temp_data_xyz(:,1:num_data) = data_xyz(:,1:num_data)
       temp_data_weight(:,1:num_data) = data_weight(:,1:num_data)
       deallocate(data_field)
       deallocate(data_xyz)
       deallocate(data_weight)
       allocate(data_field(4,num_data+n_add_data))
       allocate(data_xyz(3,num_data+n_add_data))
       allocate(data_weight(3,num_data+n_add_data))
       data_field(:,1:num_data) = temp_data_field(:,1:num_data)
       data_xyz(:,1:num_data) = temp_data_xyz(:,1:num_data)
       data_weight(:,1:num_data) = temp_data_weight(:,1:num_data)
       deallocate(temp_data_field)
       deallocate(temp_data_xyz)
       deallocate(temp_data_weight)
    endif

    ! record the start and end data numbers for this group
    ndata_groups(num_groups,1) = num_data + 1
    ndata_groups(num_groups,2) = num_data + n_add_data
    !set the counted number of data points to previous total
    ncount = num_data
    num_data = num_data+n_add_data

!!! read the data point information
    open(10, file=readfile, status='old')
    read(unit=10, fmt="(a)", iostat=ierror) buffer

    read_line_of_data : do
       
       ! read the data #; z; y; z; wd1; wd2; wd3 for each data point
       read(unit=10, fmt="(a)", iostat=ierror) buffer
       if(ierror<0) exit !ierror<0 means end of file
       length_string = len_trim(buffer) !length of buffer, and removed trailing blanks
       
       ! read data number
       buffer=adjustl(buffer) !remove leading blanks
       iend=index(buffer," ",.false.)-1 !index returns location of first blank
       if(length_string == 0) exit
       ncount=ncount+1
       read (buffer(1:iend), '(i6)') itemp
       
       do nj=1,3
          ! read x,y,z coordinates
          buffer = adjustl(buffer(iend+1:length_string)) !remove data number from string
          buffer = adjustl(buffer) !remove leading blanks
          length_string = len(buffer) !new length of buffer
          iend=index(buffer," ",.false.)-1 !index returns location of first blank
          read (buffer(1:iend), '(D25.17)') data_xyz(nj,ncount)
       enddo !nj
       data_weight(1:3,ncount)=1.0_dp
    enddo read_line_of_data
    
    close(10)
    
    call enter_exit(sub_name,2)

  end subroutine define_data_fit_group
  
!!! ##########################################################################

  subroutine calc_data_field_distance(data_elem,data_xi)
    implicit none

    integer :: data_elem(:)
    real(dp) :: data_xi(:,:)
!!! local variables
    integer :: nd,ne,nj
    real(dp) :: dz(3),elem_xyz(num_deriv_elem,num_coords),sq,xi(2),z(3)

    data_field = 0.0_dp
    do nd = 1,num_data
       ne = data_elem(nd)
       if(ne.ne.0)then
          call node_to_local_elem(ne,elem_xyz)
          xi(1:2) = data_xi(1:2,nd)
          sq = 0.0_dp
          do nj = 1,3
             z(nj) = pxi(1,xi,elem_xyz(1,nj))
             dz(nj) = z(nj) - data_xyz(nj,nd)
             sq = sq + dz(nj)**2.0_dp
          enddo
          if(abs(sq) < loose_tol)then
             forall (nj = 1:4) data_field(nj,nd) = 0.0_dp
          else
             sq = sqrt(sq)
             forall (nj = 1:3) data_field(nj,nd) = dz(nj)/sq
             data_field(4,nd) = sq
          endif
       endif
    enddo !nd

  end subroutine calc_data_field_distance

!!! ##########################################################################      

  subroutine define_xi_closest(data_elem,data_on_elem,elem_list,ndata_on_elem,data_xi,first)
!!! find the closest xi location on a 2d mesh surface to each data point
    implicit none

    integer,intent(in) :: elem_list(:)
    integer :: data_elem(:),data_on_elem(:,:),ndata_on_elem(:)
    real(dp) :: data_xi(:,:)
    logical,intent(in) :: first
!!! local variables
    integer :: egroup,i,n_check,ne_checklist(5),ngroup,IT,ITMAX=20, &
         nd,ne,neadj,nelast,neold,ni,nj
    real(dp) :: sqmax,sqnd,temp,elem_xyz(num_deriv_elem,num_coords),xi(3)
    real(dp),allocatable :: sq(:)
    logical :: found,not_converged,not_converged_at_all

    integer :: n_data,np
    character(len=200) :: exfile
    character(len=1) :: string_ne1
    character(len=2) :: string_ne2
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'define_xi_closest'
    call enter_exit(sub_name,1)
    
    allocate(sq(num_data))    

    sqmax = 1.0e4_dp*1.0e4_dp
    
    !  initialise
    sq = 0.0_dp
    xi = 0.5_dp

    if(num_groups <= 0)then ! fitting to all elements at once
       if(first)then ! check every element for every data point
          do nd = 1,num_data
             sqmax = 1.0e4_dp*1.0e4_dp
             not_converged_at_all = .false.
             do i = 1,count(elem_list(:)/=0)
                ne = get_local_elem_2d(elem_list(i))
                xi = 0.5_dp
                call node_to_local_elem(ne,elem_xyz)
                found = .false.
                call project_orthogonal(nd,SQND,elem_xyz,xi,found,not_converged)
                if(.not.not_converged)then
                   not_converged_at_all = .true.
                   if(sqnd.lt.sqmax)then
                      sqmax = sqnd
                      data_xi(1:2,nd) = xi(1:2)
                      data_elem(nd) = ne
                      SQ(nd) = SQND
                   endif
                endif !FOUND
             enddo !i
          enddo ! nd
       else

          do nd = 1,num_data
             not_converged_at_all = .false.
             ne = data_elem(nd)
             if(ne.ne.0)then ! i.e. only for data points that have projected to an element
                ne_checklist(1) = ne
                n_check = 1
                if(elem_cnct_2d(-1,0,ne).ne.0)then
                   n_check = n_check + 1
                   ne_checklist(n_check) = elem_cnct_2d(-1,1,ne)
                endif
                if(elem_cnct_2d(1,0,ne).ne.0)then
                   n_check = n_check + 1
                   ne_checklist(n_check) = elem_cnct_2d(1,1,ne)
                endif
                if(elem_cnct_2d(-2,0,ne).ne.0)then
                   n_check = n_check + 1
                   ne_checklist(n_check) = elem_cnct_2d(-2,1,ne)
                endif
                if(elem_cnct_2d(2,0,ne).ne.0)then
                   n_check = n_check + 1
                   ne_checklist(n_check) = elem_cnct_2d(2,1,ne)
                endif
                sqmax = 1.0e4_dp*1.0e4_dp
                do i = 1,n_check
                   ne = ne_checklist(i)
                   if(i.eq.1)then
                      xi(1:2) = data_xi(1:2,nd)
                   else
                      xi = 0.5_dp
                   endif
                   call node_to_local_elem(ne,elem_xyz)
                   found = .false.
                   call project_orthogonal(nd,sqnd,elem_xyz,xi,found,not_converged)
                   if(.not.not_converged)then
                      not_converged_at_all = .true.
                      if(sqnd.lt.sqmax)then
                         sqmax = sqnd
                         data_xi(1:2,nd) = xi(1:2)
                         data_elem(nd)=ne
                         sq(nd) = sqnd
                      endif
                   endif
                enddo !i
             endif ! ne.ne.0
          enddo ! nd
       endif
    else
       data_elem = 0
       do ngroup = 1,num_groups
          do i = 1,num_groups
             if(data_group_names(ngroup) == elem_group_names(i)) egroup = i
          enddo
          do nd = ndata_groups(ngroup,1),ndata_groups(ngroup,2)
             sqmax = 1.0e4_dp*1.0e4_dp
             not_converged_at_all = .false.
             do i = 1,count(nelem_groups(egroup,:)/=0)
                ne = get_local_elem_2d(nelem_groups(egroup,i))
                xi = 0.5_dp
                call node_to_local_elem(ne,elem_xyz)
                found = .false.
                call project_orthogonal(nd,SQND,elem_xyz,xi,found,not_converged)
                if(abs(xi(1)).gt.zero_tol.and.abs(xi(1)).lt.1.0_dp-zero_tol.and. &
                     abs(xi(2)).gt.zero_tol.and.abs(xi(2)).lt.1.0_dp-zero_tol) then
                   if(sqnd.lt.sqmax)then
                      sqmax = sqnd
                      data_xi(1:2,nd) = xi(1:2)
                      data_elem(nd) = ne
                      SQ(nd) = SQND
                   endif
                endif !FOUND
             enddo !i

             ! check the distance: if too far from the surface then remove from fit
             if (sqrt(sq(nd)) > 150.0_dp)then
                data_xi(1:2,nd) = 0.0_dp
                data_elem(nd) = 0
                SQ(nd) = 0.0_dp
             endif
          enddo ! nd
       enddo ! ngroup
    endif
    ndata_on_elem=0
    do nd = 1,num_data
       ne = data_elem(nd)
       if(ne > 0) then
          ndata_on_elem(ne) = ndata_on_elem(ne)+1
          if(ndata_on_elem(ne) > nmax_data_elem)then
             write(*,'(''Number of data points on element'',i6,'' exceeds'',i6)') ne,nmax_data_elem
             write(*,'(''Reduce the data point density: no advantage in using this many'')')
             stop
          endif
          data_on_elem(ne,ndata_on_elem(ne)) = nd
       endif
    enddo
    deallocate(sq)
    
    call enter_exit(sub_name,2)

  end subroutine define_xi_closest

!!! ##########################################################################      

  subroutine solve_geometry_fit(data_on_elem,ndata_on_elem,num_depvar,&
       elem_list,not_1,not_2,npny,nynp,nyny,data_xi,cyny,sobelov_wts,&
       fit_soln,fix_bcs)

!!! dummy arguments
    integer :: data_on_elem(:,:),ndata_on_elem(:),not_1,not_2,num_depvar,&
         elem_list(:),npny(:,:),nynp(:,:,:,:),nyny(0:,:)
    real(dp) :: data_xi(:,:),cyny(0:,:),sobelov_wts(0:,:),fit_soln(:,:,:,:)
    logical :: fix_bcs(:)
!!! local variables
    integer :: l,ny_local(3*16),n_dof,ne,nh,nh1,nh2,nhs1,nhs2, &
         nk,no1,no2,no_nynr1,no_nynr2,noy1,noy2,np,nv,ny1,ny2,ny3,nyo1,nz,nzz
    integer,allocatable :: nony(:,:)
    integer,allocatable :: nyno(:,:,:)
    real(dp) :: co1,co2,ER(num_fit*num_deriv_elem),ES(3*16,3*16),&
         fit_soln_local(16,3)
    real(dp),allocatable :: cony(:)
    real(dp),allocatable :: cyno(:,:)
    real(dp),allocatable :: GR(:)      ! right-hand-side vector
    real(dp),allocatable :: GRR(:)     ! reduced right-hand-side vector
    real(dp),allocatable :: incr_soln(:)         ! current solution returned from solver
    logical :: FIRST_A

    ! make all of these allocatable!
    real(dp),dimension(nsize_gkk) :: GKK
    real(dp),dimension(nsize_gkk) :: GK
! doesn't like allocating these! gives different answer for errors
!    real(dp),allocatable :: GKK(:)
!    real(dp),allocatable :: GK(:)
    integer :: np_temp,i

    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'solve_geometry_fit'
    call enter_exit(sub_name,1)
    
    allocate(incr_soln(num_depvar))
    allocate(nony(0:1,num_depvar))
    allocate(nyno(0:5,num_depvar,2))
    allocate(cony(num_depvar))
    allocate(cyno(0:5,num_depvar))
    allocate(GR(num_depvar))
    allocate(GRR(num_depvar))

    GK = 0.0_dp

    !*** Calculate solution mapping arrays for the current fit variable
    call globalf(nony,not_1,not_2,npny,nyno,nynp,nyny,cony,cyno,cyny,fix_bcs)
    
    if(NOT_2.EQ.0) then
       write(*,'('' >>The number of unknowns is zero'')')
       stop
    endif

    FIRST_A = .TRUE.
    
    GR = 0.0_dp
    
    do l = 1,count(elem_list/=0) !loop over elements in the fit
       ne = get_local_elem_2d(elem_list(l)) ! elem_list stores global elements
       
       call make_element_matrices(ne,fit_soln_local,fit_soln, &
            data_on_elem,ndata_on_elem,data_xi,ER,ES,sobelov_wts)
       
       call local_dof(n_dof,ne,ny_local,nynp)

       ! Assemble element  matrices into global matrices
       do nh1 = 1,n_dof 
          ny1 = ny_local(nh1)
          if(ny1.eq.0)then
             write(*,'('' No dependent variable for node in element'',i6,'': are &
                  &you sure you have set up versions correctly?'')') ne
             stop
          endif
          GR(ny1) = GR(ny1) + ER(nh1)
          do nh2 = 1,n_dof 
             ny2 = ny_local(nh2)
             nz = ny1+(ny2-1)*num_depvar
             GK(nz) = GK(nz) + ES(nh1,nh2)
          enddo !nh2
       enddo !nh1
    enddo !l (ne)
    
    !----------------------- generate reduced system -----------------------
    
    GKK = 0.0_dp
    GRR = 0.0_dp
    
    !*** generate the reduced system of equations
    do ny1 = 1,num_depvar !loop global rows of GK
       do noy1 = 1,nony(0,ny1) !loop over #no's attached to ny1
          no1 = nony(noy1,ny1) !no# attached to row ny1
          co1 = cony(ny1) !coupling coeff for row mapping
          !                     ie row_no1=a*row_ny1+b*row_ny2
          GRR(no1) = GRR(no1)+GR(ny1)*co1 !get reduced R.H.S.vector
          do ny2 = 1,num_depvar !loop over #cols of GK
             !local GK var #
             nz = ny1+(ny2-1)*num_depvar
             if(nz.ne.0) then
                do noy2 = 1,nony(0,ny2) !loop over #no's for ny2
                   no2 = nony(noy2,ny2) !no# attached to ny2
                   co2 = cony(ny2) !coup coeff col mapping
                   !                     i.e. var_no1=a*var_ny1+b*var_ny2
                   nzz = no1+(no2-1)*NOT_1
                   if(nzz.ne.0) GKK(nzz) = GKK(nzz)+GK(nz)*co1*co2
                enddo !noy2
             endif
          enddo !no_nynr2
       enddo !noy1
    enddo !no_nynr1
    
    !-------------- solve reduced system of linear equations ---------------
    call solve_fit_system(NOT_1,NOT_2,num_depvar,GKK,GRR,incr_soln)

    do no1 = 1,NOT_2 ! for each unknown
       do nyo1 = 1,nyno(0,no1,1)
          ny1 = nyno(nyo1,no1,1) ! the dependent variable number
          co1 = cyno(nyo1,no1)   ! the weighting for mapped variables
          nk = npny(1,ny1)       ! derivative number
          nv = npny(2,ny1)       ! version number
          nh = npny(3,ny1)       ! dependent variable number
          np = npny(4,ny1)       ! node number
          !current fit solution =        previous       +     increment
          fit_soln(nk,nv,nh,np) = fit_soln(nk,nv,nh,np) + incr_soln(no1)*co1
       enddo !nyo1
    enddo !no1

    !*** Copy the fitted solution back into node_xyz_2d
    node_xyz_2d = fit_soln

    deallocate(incr_soln)
    deallocate(nony)
    deallocate(nyno)
    deallocate(cony)
    deallocate(cyno)
    deallocate(GR)
    deallocate(GRR)

    call enter_exit(sub_name,2)
  
  end subroutine solve_geometry_fit
  
!!! ##########################################################################        
  
  subroutine project_orthogonal(nd,SQ,elem_xyz,xi,inelem,not_converged)

!!! dummy arguments        
    integer :: nd
    real(dp) :: sq,elem_xyz(num_deriv_elem,num_coords),xi(:)
    logical :: inelem,not_converged
!!! local variables
    integer :: IT,ITMAX,BOUND(2),it2,ni,nifix,nj
    real(dp) :: LOOSE_TOL=1.0e-6_dp
    real(dp) :: DELTA,DET,D2SQV2,D2SQVW2,D2SQXI(2,2),D2ZXI(3,2,2),DSQXI(2), &
         DSQXI1,DSQXI2,DSQV,DSQVW,DZ(3),DZXI(3,2),EVMIN,EVMAX,H(2), &
         MU,SQLIN,SQDIFF,SQDPRED,TEMP,TEMP1,TEMP2,TOL, &
         TOL2,V(2),V1,V2,VMAX=1.0_dp,W,XILIN(2),Z(3)
    logical :: CONVERGED,ENFORCE(2),FREE,NEWTON
    
    ! --------------------------------------------------------------------------

    not_converged = .false.
    
    ITMAX = 10 ! max # iterations to use
    DELTA = VMAX/4.0_dp
    TOL = 5.0_dp*LOOSE_TOL !must be > sqrt(eps) or SQLIN<=SQ check may not work
    TOL2 = TOL**2
    SQ = 0.0_dp
    do nj=1,num_coords
       Z(nj) = PXI(1,XI,elem_xyz(1,nj))
       DZ(nj) = Z(nj)-data_xyz(nj,nd)
       SQ = SQ+DZ(nj)**2
    enddo !nj
    IT=0
    CONVERGED=.FALSE.
    do WHILE(.NOT.CONVERGED.AND.IT.LT.ITMAX)
       DSQXI = 0.0_dp
       do nj=1,num_coords
          DZXI(nj,1) = PXI(2,XI,elem_xyz(1,nj))
          DZXI(nj,2) = PXI(4,XI,elem_xyz(1,nj))
          DSQXI(1) = DSQXI(1)+DZXI(nj,1)*DZ(nj)
          DSQXI(2) = DSQXI(2)+DZXI(nj,2)*DZ(nj)
       enddo !nj
       do ni=1,2
          if(dabs(xi(ni)) < zero_tol) then
             BOUND(ni)=1
             ENFORCE(ni)=DSQXI(ni).GE.0.0_dp
          else if(dabs(xi(ni)-1.0_dp) < zero_tol) then
             BOUND(ni)=-1
             ENFORCE(ni)=DSQXI(ni).LE.0.0_dp
          else
             BOUND(ni)=0
             ENFORCE(ni)=.FALSE.
          endif
       enddo !ni
       if(ENFORCE(1).AND.ENFORCE(2)) exit
       
       D2SQXI = 0.0_dp
       do nj=1,num_coords
          D2ZXI(nj,1,1) = PXI(3,XI,elem_xyz(1,nj))
          D2ZXI(nj,1,2) = PXI(5,XI,elem_xyz(1,nj))
          D2ZXI(nj,2,2) = PXI(5,XI,elem_xyz(1,nj))
          D2SQXI(1,1) = D2SQXI(1,1)+DZXI(nj,1)*DZXI(nj,1)+D2ZXI(nj,1,1)*DZ(nj)
          D2SQXI(1,2) = D2SQXI(1,2)+DZXI(nj,1)*DZXI(nj,2)+D2ZXI(nj,1,2)*DZ(nj)
          D2SQXI(2,2) = D2SQXI(2,2)+DZXI(nj,2)*DZXI(nj,2)+D2ZXI(nj,2,2)*DZ(nj)
       enddo !nj
       
       !       A Newton step is taken if the condition of the Hessian
       !       guarantees that the step will be within the trust region.
       !       Otherwise the Hessian is shifted towrds a diagonal matrix to
       !       shift the step towards steepest descent.  Usually it is a much
       !       better direction than steepest descent.  I think it is close to
       !       the best direction in the trust region.
       
       !***    Find the smallest eigen value of the Hessian.
       DSQXI2 = DSQXI(1)**2+DSQXI(2)**2
       DSQXI1 = DSQRT(DSQXI2)
       TEMP1 = (D2SQXI(1,1)+D2SQXI(2,2))/2.0_dp
       TEMP2 = DSQRT(((D2SQXI(1,1)-D2SQXI(2,2))/2.0_dp)**2+D2SQXI(1,2)**2)
       EVMIN = TEMP1-TEMP2
       EVMAX = TEMP1+TEMP2
       if(DSQXI1.LT.TOL2) exit
       do it2=1,ITMAX
          TEMP = DSQXI1/DELTA
          NEWTON = EVMIN.GE.TEMP
          if(NEWTON) then !Newton is safe
             H(1) = D2SQXI(1,1)
             H(2) = D2SQXI(2,2)
             DET = EVMIN*EVMAX
          else
             !***        Shift eigenvalues to restrict step
             MU = TEMP-EVMIN
             H(1) = D2SQXI(1,1)+MU
             H(2) = D2SQXI(2,2)+MU
             DET = TEMP*(EVMAX+MU)
          endif
          V(1) = -(H(2)*DSQXI(1)-D2SQXI(1,2)*DSQXI(2))/DET
          V(2) = (D2SQXI(1,2)*DSQXI(1)-H(1)*DSQXI(2))/DET
          V2 = V(1)**2+V(2)**2
          DSQV = DSQXI(1)*V(1)+DSQXI(2)*V(2)
          !         This checks that numerical errors have not
          !         prevented the step direction being a descent direction.
          
          if(DSQV**2.LT.DSQXI2*V2*TOL2) then !try a smaller trust region
             DELTA = DELTA/10.0_dp
          else !step is good
             !***        Check feasible and limit step size
             FREE = .TRUE.
             do ni=1,2
                if(BOUND(ni)/=0 )then
                   if(BOUND(ni)>0.EQV.V(ni)<0.0_dp) then
                      FREE=.FALSE.
                      nifix=ni
                   endif
                endif
             enddo
             W=1.0_dp
             if(FREE) then
                V1=DSQRT(V2) !currently < DELTA
                D2SQV2=V(1)*(V(1)*D2SQXI(1,1)+2.0_dp*V(2)*D2SQXI(1,2)) &
                     +V(2)**2*D2SQXI(2,2)
                if(.NOT.NEWTON) then
                   !               Try to step to estimate of minimum along line
                   if(V1.GT.0.0_dp) then
                      W=DELTA/V1
                      if(D2SQV2.GT.0.0_dp) then !minimum exists
                         W=DMIN1(W,-DSQV/D2SQV2) !minimum if within trust region
                      endif
                   endif
                endif !newton
             else
                if(ENFORCE(2)) then !gradient suggests must use ni=1
                   nifix=2
                else if(ENFORCE(1)) then !gradient suggests must use ni=2
                   nifix=1
                   !else Gradient points into element.
                   !               Fix the direction that prevented the step.
                   !               P.d. Hessian guarantees there is only one of these
                   !               for this type of gradient.
                endif
                ni=3-nifix
                if(.NOT.INELEM) then
                   !***            If stepping predominantly out of element then exit
                   nifix=3-ni
                   if(DABS(V(nifix)).GT.DABS(DSQXI(ni)/H(ni))) then
                      XI(nifix)=XI(nifix)+V(nifix)
                      exit
                   endif
                endif
                V(nifix)=0.0_dp
                if(D2SQXI(ni,ni).GT.0.0_dp) then !minimum exists
                   V(ni)=-DSQXI(ni)/D2SQXI(ni,ni)
                   V1=DABS(V(ni))
                   NEWTON=V1.LE.DELTA
                endif
                if(.NOT.NEWTON) then
                   V(ni)=-DSIGN(DELTA,DSQXI(ni))
                   V1=DELTA
                endif
                V2=V1*V1
                DSQV=DSQXI(ni)*V(ni)
                D2SQV2=V2*D2SQXI(ni,ni)
             endif !free
             !***        First half of convergence test.
             !           Should be before boundary colllision check
             CONVERGED=V1*W.LT.TOL
             !***        Try the step.  (Name: XILIN is historical)
             XILIN(1)=XI(1)+V(1)*W
             XILIN(2)=XI(2)+V(2)*W
             !***        Test for boundary collision
             do ni=1,2
                if(XILIN(ni).LT.0.0_dp) then
                   XILIN(ni)=0.0_dp
                   W=XI(ni)/(-V(ni))
                   XILIN(3-ni)=XI(3-ni)+V(3-ni)*W
                else if(XILIN(ni).GT.1.0_dp) then
                   XILIN(ni)=1.0_dp
                   W=(1.0_dp-XI(ni))/V(ni)
                   XILIN(3-ni)=XI(3-ni)+V(3-ni)*W
                endif
             enddo !ni
             !***        Calculate new distance
             SQLIN=0.0_dp
             do nj=1,num_coords
                Z(nj)=PXI(1,XILIN,elem_xyz(1,nj))
                DZ(nj)=Z(nj)-data_xyz(nj,nd)
                SQLIN=SQLIN+DZ(nj)**2
             enddo
             !***        Second half of convergence test.
             CONVERGED=CONVERGED.AND.DABS(SQ-SQLIN)/(1.0_dp+SQ).LE.TOL
             if(CONVERGED) GO TO 5
             DSQVW=DSQV*W !<0
             D2SQVW2=0.5_dp*D2SQV2*W*W !1/2 for computational efficiency
             SQDIFF=0.5_dp*(SQLIN-SQ) !1/2 because derivs are for SQ/2
             SQDPRED=SQDIFF-DSQVW-D2SQVW2
             !***        Exit loop if decrease is satisfactory
             if(SQDIFF.LE.0.25_dp*DSQVW) then
                if(NEWTON) then
                   DELTA=V1 !next step smaller unless this is increased
                else if(W.GE.1.0_dp) then
                   !               If the quadratic model is good increase trust region size
                   if(SQDPRED.LT.-0.1_dp*SQDIFF) then
                      DELTA=DMIN1(VMAX,DELTA*2.0_dp)
                   endif
                endif
                GO TO 5
             endif
             !***        Calculate new trust region size from an estimate of the
             !***        minimum along the step direction using a cubic approximation.
             TEMP=-3.0_dp*SQDPRED !<0
             DELTA=W*V1*(D2SQVW2-DSQRT(D2SQVW2**2+TEMP*DSQVW))/TEMP !>0
!!! note: was getting delta=0 because d2sqvw2=0. the following avoids, but is it correct?
             if(abs(delta).le.zero_tol)then
                delta = 0.01_dp
             endif
          endif !DSQV**2.LT.DSQXI2*V2*TOL2
       enddo !it2
       
5      SQ=SQLIN
       XI(1)=XILIN(1)
       XI(2)=XILIN(2)
       IT=IT+1
    enddo
    
    if(.not.converged) then
       not_converged = .true.
    endif

    if(.NOT.inelem.AND.XI(1)>=0.0_dp.AND.XI(1)<=1.0_dp.AND.XI(2)>=0.0_dp.AND.XI(2)<=1.0_dp) then
       inelem=.TRUE.
    endif
    
  end subroutine project_orthogonal

!!! ##########################################################################      
  
  subroutine solve_fit_system(M,N,num_depvar,A,B,X)

    integer :: M,N,num_depvar
    real(dp) :: A(:),B(:),X(:)
!!! local variables
    integer,allocatable :: pivots(:)
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'solve_fit_system'
    call enter_exit(sub_name,1)
    
    allocate(pivots(num_depvar*num_depvar))

    ! LU Decomposition
    call lu_decomp(M,N,pivots,A)
    ! Solve the problem
    X(1:m) = B(1:m)
    call lu_solve(M,N,pivots,A,X)

    deallocate(pivots)
    
    call enter_exit(sub_name,2)

  end subroutine solve_fit_system
  
!!! ##########################################################################      

  subroutine lu_decomp(M,N,pivots,A)

    integer :: M,N,pivots(:)
    real(dp) :: A(M, * )
    integer :: j,j_pivot
    real(dp),allocatable :: temp_vec(:)
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'lu_decomp'
    call enter_exit(sub_name,1)

    allocate(temp_vec(M))

    do j = 1, min(M,N)
       j_pivot = j-1 + indexmax( M-j+1, A(j,j))
       pivots(j) = j_pivot
       if(abs(A(j_pivot,j)).gt.zero_tol )then
          if(j_pivot.ne.j)then
             ! swap rows a(j,:) and a(j_pivot,:)
             temp_vec(1:n) = a(j,1:n)
             a(j,1:n) = a(j_pivot,1:n)
             a(j_pivot,1:n) = temp_vec(1:n)
          endif
          if(j.lt.M) &
               ! scale vector by a constant (a(j+1,j) by 1/(a(j,j))
               call scale_row_lu(M-j,1.0_dp/A(j,j),A(j+1,j))
       endif
       if(j.lt.MIN(M,N))then
          ! A = A + x*y'
          call add_to_row_lu(M,M-j,N-j,A(j+1,j),A(j,j+1),A(j+1,j+1))
       endif
    enddo

    deallocate(temp_vec)
    
    call enter_exit(sub_name,2)

  end subroutine lu_decomp

!!! ##########################################################################      

  subroutine scale_row_lu(n,da,dx)

      real(dp) :: da,dx(*)
      integer :: i,n

      do i = 1,n
        dx(i) = da*dx(i)
     enddo
     
   end subroutine scale_row_lu

!!! ##########################################################################      

   subroutine add_to_row_lu(M,maxrow,N,X,Y,A)
     integer :: M,maxrow,N
     real(dp) :: A(M,*),X(*),Y(*)
     integer :: j, JY
     
     JY = 1
     do j = 1, N !n-j
        if(abs(Y(JY)).gt.zero_tol)then
           A(1:maxrow,j) = A(1:maxrow,j) - X(1:maxrow)*Y(jy)
        endif
        JY = JY + M
     enddo
     
   end subroutine add_to_row_lu

!!! ##########################################################################      

  subroutine lu_solve(M,N,pivots,A,B)

    integer :: M,N,pivots(:)
    real(dp) :: A(M,*),B(M, * )
    integer :: i,k,pivot_row
    real(dp) :: temp

    ! swap row i and pivot_row for each of rows 1..N
    do i = 1,N
       pivot_row = pivots(i)
       if(pivot_row.ne.i)then
          temp = B(i,1)
          B(i,1) = B(pivot_row,1)
          B(pivot_row,1) = temp
       endif
    enddo
    
    ! Solve L*X = B, overwriting B with X.
    do k = 1,N
       if(abs(B(k,1)).gt.zero_tol) &
            B(k+1:N,1) = B(k+1:N,1) - B(k,1)*A(k+1:N,k)
    enddo

    ! Solve U*X = B, overwriting B with X.
    do k = N,1,-1
       if(abs(B(k,1)).gt.zero_tol ) then
          B(k,1) = B(k,1)/A(k,k)
          B(1:k-1,1) = B(1:k-1,1) - B(k,1)*A(1:k-1,k)
       endif
    enddo

  end subroutine lu_solve

  function indexmax(n,dx)
    ! return the index of entry with largest abs value
    real(dp) :: dx(*),dmax
    integer :: i,n
    integer :: indexmax
     
    indexmax = 1
    if(n.gt.1)then
       dmax = abs(dx(1))
       do i = 2,n
          if(abs(dx(i)).gt.dmax)then
             indexmax = i
             dmax = abs(dx(i))
          endif
       enddo
    endif
     
   end function indexmax

!!! ##########################################################################      

end module surface_fitting
