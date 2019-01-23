module surface_fitting

  use arrays
  use other_consts
  use mesh_functions
  !use precision
  use solve

  implicit none

  public fit_surface_geometry,pxi

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

!!! completes 'niterations' of geometry fitting to a surface, via minimising 
!!! the least squares distance between a list of data points (3D RC coordinates)
!!! and a surface mesh (assumed bi-cubic Hermite only). 'fitting_file' lists
!!! the nodes/derivatives that are fixed, and any mapping of nodes and/or
!!! derivatives

    use geometry,only: get_local_node_f
!!! dummy arguments
    integer,intent(in) :: niterations             ! user-specified number of fitting iterations
    character(len=255),intent(in) :: fitting_file ! file that lists versions/mapping/BCs
!!! local variables
    integer  :: nfit,nk,NOT_1,NOT_2,np,num_depvar,nv,ny_max
    logical :: first = .true.
!!! local allocatable arrays
    integer,allocatable :: data_elem(:)              
    integer,allocatable :: data_on_elem(:,:)     ! list of data closest to elements
    integer,allocatable :: elem_list(:)          ! list of elements in fit (all)
    integer,allocatable :: ndata_on_elem(:)      ! # of data points closest to element
    integer,allocatable :: npny(:,:)             ! maps deriv, vers, node, etc for a dep. variable
    integer,allocatable :: nynp(:,:,:,:)         ! dep. variable # for node, deriv, version etc.
    integer,allocatable :: nynr(:)               ! list of all dep. variables
    integer,allocatable :: nyny(:,:)             ! maps dep. variable to another dep. variable
    real(dp),allocatable :: cyny(:,:)            ! weighting for mapped dep. variables
    real(dp),allocatable :: data_xi(:,:)         ! xi location of data points
    real(dp),allocatable :: fit_soln(:,:,:,:)    ! current fit geometry solution
    real(dp),allocatable :: sobelov_wts(:,:)     ! Scaling factor for, and Sobelov weights
    logical,allocatable :: fix_bcs(:)            ! logical for boundary conditions


!!! allocate element-sized arrays
    allocate(elem_list(0:num_elems_2d))
    allocate(sobelov_wts(0:6,num_elems_2d))

!!! allocate data arrays
    allocate(data_elem(num_data))
    allocate(data_on_elem(num_elems_2d,nmax_data_elem))
    allocate(ndata_on_elem(num_elems_2d))
    allocate(data_xi(2,num_data))
    data_elem = 0

!!! allocate dependent variable arrays
    ny_max = 0
    do np = 1,num_nodes_2d
       do nv = 1,node_versn_2d(np)
          ny_max = ny_max+1
       enddo
    enddo
    ny_max = ny_max*num_nodes_2d*num_fit*num_deriv  ! nodes * coordinates * #derivatives+1
    allocate(nynr(0:ny_max))
    allocate(npny(0:6,ny_max))
    allocate(nynp(num_deriv,nmax_versn,num_fit,num_nodes_2d))
    allocate(fix_bcs(ny_max))

!!! find the closest surface to each data point, and calculate the Xi
!!! coordinates of the data point to the surface
    write(*,'('' Calculating normal projections: slow first time '')')
    call define_xi_closest(data_elem,data_on_elem,ndata_on_elem,data_xi,first)
    first = .false.
    
!!! list the total error between the data points and the surface
    write(*,'(/'' TOTAL RMS ERROR PRIOR TO FITTING:'')') 
    call list_data_error(data_on_elem,ndata_on_elem,data_xi)

!!! read 'fitting_file' to define the fitting constraints. set up mapping
!!! arrays, dependent variable arrays, Sobelov smoothing weights (hardcoded)
    write(*,'('' Define fitting problem '')')
!    call define_geometry_fit(elem_list,npny,num_depvar,nynp,nynr,nyny,&
!         cyny,sobelov_wts,fit_soln,fitting_file,fix_bcs)

    do nfit = 1,niterations ! user-defined number of iterations
       call define_geometry_fit(elem_list,npny,num_depvar,nynp,nynr,nyny,&
            cyny,sobelov_wts,fit_soln,fitting_file,fix_bcs)
       write(*,'(/'' FITTING ITERATION'',I3)') nfit
!!!    solve for new nodal coordinates and derivatives
       write(*,'('' Solve fitting problem '')')
       call solve_geometry_fit(data_on_elem,ndata_on_elem,num_depvar,&
            elem_list,not_1,not_2,npny,nynp,nynr,&
            nyny,data_xi,cyny,sobelov_wts,fit_soln,fix_bcs)
!!!    update the scale factors for new geometry if NOT unit scale factors
!       write(*,'('' Update scale factors '')')
!       call update_scale_factor_norm
!!!    update the data point projections and their Xi coordinates
       write(*,'('' Calculating normal projections '')')
       call define_xi_closest(data_elem,data_on_elem,ndata_on_elem,data_xi,first)
!!!    calculated the updated error between data and surface
       write(*,'(/'' CURRENT RMS ERROR FOR ALL DATA:'')') 
       call list_data_error(data_on_elem,ndata_on_elem,data_xi)
    enddo

    deallocate(elem_list)
    deallocate(sobelov_wts)
    deallocate(data_on_elem)
    deallocate(ndata_on_elem)
    deallocate(data_xi)
    deallocate(cyny)
    deallocate(nynr)
    deallocate(npny)
    deallocate(nynp)
    deallocate(nyny)
    deallocate(fix_bcs)

  end subroutine fit_surface_geometry

!!! ##########################################################################      
  
  subroutine define_geometry_fit(elem_list,npny,num_depvar,nynp,nynr,nyny,&
       cyny,sobelov_wts,fit_soln,fitting_file,fix_bcs)

!!! read information from 'fitting_file' to determine the boundary conditions
!!! and mapping of dependent variables for geometric fitting. set up the 
!!! dependent variable-to-mapping arrays  
  
    use arrays,only: elems_2d,node_versn_2d,node_xyz_2d,num_elems_2d,num_nodes_2d
    use geometry,only: get_final_integer,get_local_node_f
!!! dummy arguments
    integer :: elem_list(0:),npny(0:,:),num_depvar,nynp(:,:,:,:),nynr(0:)
    integer,allocatable :: nyny(:,:)
    real(dp) :: sobelov_wts(0:,:)
    real(dp),allocatable :: cyny(:,:),fit_soln(:,:,:,:)
    character(len=*),intent(in) :: fitting_file
    logical :: fix_bcs(:)
!!! local variables
    integer :: i,ibeg,iend,ierror,IPFILE=10,i_ss_end,L,ne,nh,nj,nk,node,np,&
         np_global,number_of_fixed,nv,nv_fix,ny
    character(len=132) :: readfile,string
    

    if(.not.allocated(fit_soln)) allocate(fit_soln(4,10,16,num_nodes_2d))
    
    ! linear fitting for 3 geometric variables. solution stored in fields 1,2,3
    ! includes Sobelov smoothing on the geometry field
    
    !***Set up dependent variable interpolation information
    fit_soln = node_xyz_2d    
 
!!! the following not correct because it refers to the global element #s 
!!! use elem_list because we might want to fit only some of the elements  
!    elem_list(1:num_elems_2d) = elems_2d(1:num_elems_2d)
    forall (i=1:num_elems_2d) elem_list(i) = i
    elem_list(0) = num_elems_2d

    ! *** Specify smoothing constraints on each element
    do L=1,elem_list(0)
       ne=elem_list(L)
       sobelov_wts(0,ne) = 1.0_dp 
       sobelov_wts(1,ne) = 1.0_dp !the scaling factor for the Sobolev weights
       !  The 5 weights on derivs wrt Xi_1/_11/_2/_22/'_12 are:
       sobelov_wts(2,ne) = 1.0e-4_dp !weight for deriv wrt Xi_1
       sobelov_wts(3,ne) = 2.0e-3_dp
       sobelov_wts(4,ne) = 1.0e-4_dp
       sobelov_wts(5,ne) = 2.0e-3_dp
       sobelov_wts(6,ne) = 5.0e-3_dp
    enddo !L
    
    !*** Calculate ny maps
    call calculate_ny_maps(npny,num_depvar,nynp,nynr)
    
    fix_bcs = .false. !initialise, default
   
    readfile = trim(fitting_file)//'.ipmap'
    open(IPFILE, file = readfile, status='old')
    read_number_of_fixed : do
       read(unit=IPFILE, fmt="(a)", iostat=ierror) string
       if(index(string, "fixed")> 0) then
          call get_final_integer(string,number_of_fixed)
          exit read_number_of_fixed
       endif
    end do read_number_of_fixed
    
    do node = 1,number_of_fixed
       read(unit=IPFILE, fmt="(a)", iostat=ierror) string
       ibeg = 1
       i_ss_end = len(string) !get the end location of the sub-string
       iend=index(string," ") !get location of next blank in sub-string
       read (string(ibeg:iend-1), '(i6)' ) np_global
       np = get_local_node_f(2,np_global)

       string = adjustl(string(iend:i_ss_end)) ! get chars beyond " " and remove the leading blanks
       iend=index(string," ") !get location of next blank in sub-string
       read (string(ibeg:iend-1), '(i6)' ) nv_fix

       string = adjustl(string(iend:i_ss_end)) ! get chars beyond " " and remove the leading blanks
       read (string(ibeg:i_ss_end), '(i6)' ) nk

!       string = adjustl(string(iend:i_ss_end)) ! get chars beyond " " and remove the leading blanks
!       read (string(ibeg:i_ss_end), '(i6)' ) nk

       nk=nk+1 !read in 0 for coordinate, 1 for 1st deriv, 2 for 2nd deriv
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
             fit_soln(nk,nv,nh,np) = 0.0_dp
          enddo !nh
       endif
    enddo !node
    
    call map_versions(IPFILE,num_depvar,nynp,nyny,cyny,fit_soln,fix_bcs)
    close(IPFILE)
    
    ! fix ALL of the cross derivatives, and set to zero
    do np = 1,num_nodes_2d
       nk = 4 !index for 1-2 cross-derivative
       do nv = 1,node_versn_2d(np)
          do nj = 1,num_coords
             ny = nynp(nk,nv,nj,np)
             fix_bcs(ny) = .TRUE.
             node_xyz_2d(nk,nv,nj,np) = 0.0_dp
          enddo !nj
       enddo !nv
    enddo !np

  end subroutine define_geometry_fit
  
!!! ##########################################################################      
  
  subroutine gauss1(PG)
    
!!! dummy arguments
    real(dp) :: PG(:,:,:)
!!! local variables
    integer :: I,J,ng,nk,nn,ns,nu
    real(dp) :: D(3),XI(3),XIGG(3,3,2)
    
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

  end subroutine gauss1
  
!!! ##########################################################################      

  function getnyr(npny,ny,nynp)
    
!!! returns the dependent variable number
!!! dummy arguments
    integer :: npny(0:,:),ny,nynp(:,:,:,:)
!!! local variables
    integer :: nh,nk,np,nv
    integer :: getnyr
    
    getnyr = 0
    nk = npny(1,ny)
    nv = npny(2,ny)
    nh = npny(3,ny)
    np = npny(4,ny)
    getnyr = nynp(nk,nv,nh,np)
    
  end function getnyr
  
!!! ##########################################################################      
  
  subroutine globalf(nony,not_1,not_2,npny,nyno,nynp,nyny,cony,cyno,cyny,fix_bcs)

!!! calculates the mapping arrays nyno/nony/cyno/cony

    use arrays,only: node_versn_2d,num_nodes_2d

!!! dummy arguments
    integer :: nony(0:,:,:),not_1,not_2,npny(0:,:),nyno(0:,:,:),nynp(:,:,:,:),nyny(0:,:)
    real(dp) :: cony(0:,:,:),cyno(0:,:,:),cyny(0:,:)
    logical :: fix_bcs(:)
!!! local variables
    integer :: nh,nv,nk,no,no_tot(2),np,nrc,ny,nyy(2),nyo,nyr,nyr2,nyy2(2),ny2
    real(dp) :: COY,RATIO
    logical :: done
    
!!!***  Initialise mapping arrays
    nony = 0
    cony = 0.0_dp
    nyno = 0
    cyno = 0.0_dp
    no_tot = 0
    
!!!*** Calculate mapping arrays
    do np=1,num_nodes_2d
       do nh=1,num_fit
          do nv=1,node_versn_2d(np)
             do nk=1,num_deriv
                ny=nynp(nk,nv,nh,np)
                if(.not.fix_bcs(ny)) then 
                   ! variable needs to be solved for
                   done = .false.
                   nyy(1) = nynp(nk,nv,nh,np) !global row #
                   nyy(2) = ny !global variable #
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
                         nyy2(1) = getnyr(npny,ny2,nynp) !row#
                         nyy2(2) = ny2 !global col#
                         do nrc=1,2 !nrc=1,2 local row and local column
                            nyr = nyy(nrc)
                            nyr2 = nyy2(nrc)
                            nony(0,nyr,nrc) = 1
                            no = nony(1,nyr2,nrc)
                            nony(1,nyr,nrc) = no
                            COY = RATIO*cony(1,nyr2,nrc)
                            cony(1,nyr,nrc) = COY
!                            write(*,*) np,nh,nv,nk,ny,ny2,nrc
                            nyo = nyno(0,no,nrc)+1
                            nyno(0,no,nrc) = nyo
                            nyno(nyo,no,nrc) = nyr
                            cyno(nyo,no,nrc) = COY
                         enddo !nrc
                      endif !ny2=0/fix_bcs
                   endif !ny.NE.ny2
                   
                   if(.not.done) then
                      do nrc=1,2 !rows and columns
                         no_tot(nrc) = no_tot(nrc)+1
                         nony(0,nyy(nrc),nrc) = 1
                         nony(1,nyy(nrc),nrc) = no_tot(nrc)
                         cony(0,nyy(nrc),nrc) = 0.0_dp
                         cony(1,nyy(nrc),nrc) = 1.0_dp
                         nyno(0,no_tot(nrc),nrc) = 1
                         nyno(1,no_tot(nrc),nrc) = nyy(nrc)
                         cyno(0,no_tot(nrc),nrc) = 0.0_dp
                         cyno(1,no_tot(nrc),nrc) = 1.0_dp
                      enddo !nrc
                   endif !not done
                endif !fix
             enddo !nk
          enddo !nv
       enddo !nh
    enddo !np
    
    NOT_1 = no_tot(1)
    NOT_2 = no_tot(2)
    
  end subroutine globalf

!!! ##########################################################################      

  subroutine line_segments_for_2d_mesh
    
!!! sets up the line segment arrays for a 2d mesh
    
    use arrays,only: arclength,elem_cnct_2d,elem_nodes_2d,elem_versn_2d,elem_lines_2d,&
         lines_in_elem,line_versn_2d,&
         lines_2d,nodes_in_line,num_elems_2d,num_lines_2d,scale_factors_2d

!!! local variables
    integer :: ne,ne_adjacent,ni1,nj,npn(2)
    logical :: MAKE
    
    num_lines_2d=0
    ! estimate number of lines, for allocating memory to arrays
    do ne=1,num_elems_2d
       if(elem_cnct_2d(-1,0,ne) == 0) num_lines_2d=num_lines_2d+1
       if(elem_cnct_2d(-2,0,ne) == 0) num_lines_2d=num_lines_2d+1
       num_lines_2d=num_lines_2d+2 ! the minimum # of new lines for each element
    enddo

    if(.not.allocated(lines_2d)) allocate (lines_2d(0:num_lines_2d))
    if(.not.allocated(line_versn_2d)) allocate(line_versn_2d(2,3,num_lines_2d))
    if(.not.allocated(elem_lines_2d)) allocate (elem_lines_2d(4,num_elems_2d))
    if(.not.allocated(lines_in_elem)) allocate (lines_in_elem(0:4,num_lines_2d))
    if(.not.allocated(nodes_in_line)) allocate (nodes_in_line(3,0:3,num_lines_2d))
    if(.not.allocated(scale_factors_2d)) allocate(scale_factors_2d(16,num_elems_2d))
    if(.not.allocated(arclength)) allocate(arclength(3,num_lines_2d)) 

    lines_in_elem=0
    lines_2d=0
    elem_lines_2d=0
    nodes_in_line=0
    line_versn_2d=0
    num_lines_2d=0

    do ne=1,num_elems_2d
       !check whether to make a line
       MAKE=.FALSE.
       if(elem_cnct_2d(-1,0,ne) == 0) MAKE=.TRUE. !exterior, make line
       ne_adjacent=elem_cnct_2d(-1,1,ne)
       if(ne_adjacent.gt.0)then
          if(elem_lines_2d(4,ne_adjacent) == 0) MAKE=.TRUE.
       endif
       
       if(MAKE)then
          num_lines_2d=num_lines_2d+1
          lines_2d(num_lines_2d)=num_lines_2d !record a new line number
          lines_in_elem(0,num_lines_2d)=lines_in_elem(0,num_lines_2d)+1
          lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d)=ne !line num_lines_2d is in element ne
          elem_lines_2d(3,ne)=num_lines_2d !num_lines_2d is global line # corresponding to local line 3 of ne
          npn(1)=1
          npn(2)=3
          nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(1,ne) !records 1st node in line
          nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(3,ne) !records 2nd node in line
          nodes_in_line(1,0,num_lines_2d)=2 !Xi-direction of line segment num_lines_2d
          do nj=1,3
             nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
             do ni1=1,2
                line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
             enddo !n
          enddo !nj
       else !get adjacent element line number
          !WARNING:: this only works if all Xi directions are consistent!!!!
          ne_adjacent=elem_cnct_2d(-1,1,ne)
          elem_lines_2d(3,ne)=elem_lines_2d(4,ne_adjacent)
       endif
       
       !check whether to make a line
       MAKE=.FALSE.
       if(elem_cnct_2d(-2,0,ne) == 0) MAKE=.TRUE. !exterior, make line
       ne_adjacent=elem_cnct_2d(-2,1,ne)
       if(ne_adjacent.gt.0)then
          if(elem_lines_2d(2,ne_adjacent) == 0) MAKE=.TRUE.
       endif

       if(MAKE)then
          num_lines_2d=num_lines_2d+1
          lines_2d(num_lines_2d)=num_lines_2d !record a new line number
          lines_in_elem(0,num_lines_2d)=lines_in_elem(0,num_lines_2d)+1
          lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d)=ne !line num_lines_2d is in element ne
          elem_lines_2d(1,ne)=num_lines_2d !num_lines_2d is global line # corresponding to local line 1 of ne
          npn(1)=1
          npn(2)=2
          nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(1,ne) !records 1st node in line
          nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(2,ne) !records 2nd node in line
          !        write(*,*) 'line in -2 for ne',ne,num_lines_2d,' nodes',npne(1,ne),npne(2,ne)
          nodes_in_line(1,0,num_lines_2d)=1 !Xi-direction of line segment num_lines_2d
          do nj=1,3
             nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
             do ni1=1,2
                line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
             enddo !n
          enddo !nj
       else !get adjacent element line number
          !WARNING:: this only works if all Xi directions are consistent!!!!
          ne_adjacent=elem_cnct_2d(-2,1,ne)
          elem_lines_2d(1,ne)=elem_lines_2d(2,ne_adjacent)
          !        write(*,*) 'adjacent in -2',ne,ne_adjacent,elem_lines_2d(1,ne)
       endif
       
       num_lines_2d=num_lines_2d+1
       lines_2d(num_lines_2d)=num_lines_2d !record a new line number
       lines_in_elem(0,num_lines_2d)=lines_in_elem(0,num_lines_2d)+1
       lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d)=ne !line num_lines_2d is in element ne
       elem_lines_2d(4,ne)=num_lines_2d !num_lines_2d is global line # corresponding to local line 4 of ne
       npn(1)=2
       npn(2)=4
       nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(2,ne) !records 1st node in line
       nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(4,ne) !records 2nd node in line
       !     write(*,*) 'line in +2 for ne',ne,num_lines_2d,' nodes',npne(2,ne),npne(4,ne)
       nodes_in_line(1,0,num_lines_2d)=2 !Xi-direction of line segment num_lines_2d
       do nj=1,3
          nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
          do ni1=1,2
             line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
          enddo !n
       enddo !nj
       
       num_lines_2d = num_lines_2d+1
       lines_2d(num_lines_2d) = num_lines_2d !record a new line number
       lines_in_elem(0,num_lines_2d) = lines_in_elem(0,num_lines_2d)+1
       lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d) = ne !line num_lines_2d is in element ne
       elem_lines_2d(2,ne)=num_lines_2d !num_lines_2d is global line # corresponding to local line 2 of ne
       npn(1) = 3
       npn(2) = 4
       nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(3,ne) !records 1st node in line
       nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(4,ne) !records 2nd node in line
       nodes_in_line(1,0,num_lines_2d)=1 !Xi-direction of line segment num_lines_2d
       do nj=1,3
          nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
          do ni1=1,2
             line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
          enddo !n
       enddo !nj
       
    enddo !ne
    
    call calc_scale_factors_2d('arcl')
    
  end subroutine line_segments_for_2d_mesh

!!! ##########################################################################      

  subroutine list_data_error(data_on_elem,ndata_on_elem,data_xi)

!!! calculate and write out the RMS error for distance between data points
!!! and 2d mesh surface

    use arrays,only: data_xyz,num_elems_2d
!!! dummy arguments
    integer :: data_on_elem(:,:),ndata_on_elem(:) 
    real(dp) :: data_xi(:,:)
!!! local variables
    integer elem,nd,nde,num_data_infit,ne,nj
    real(dp) :: data_xi_local(2),EDD,SAED,SMED,SUM,SQED,X(6),&
         XE(num_deriv_elem,num_coords)
    

    SMED=0.0_dp
    SAED=0.0_dp
    SQED=0.0_dp
    num_data_infit=0
    
    do ne=1,num_elems_2d
       call xpxe(ne,xe)
       elem=ne
       do nde=1,ndata_on_elem(elem) !for each data point on element
          nd=data_on_elem(elem,nde) !the data point number
          data_xi_local(1:2) = data_xi(1:2,nd)
          do nj=1,num_coords 
             X(nj)=PXI(1,data_xi_local,XE(1,nj))
          enddo
          SUM=0.0_dp
          do nj=1,num_coords
             SUM=SUM+(X(nj)-data_xyz(nj,nd))**2
          enddo !nj
          EDD=DSQRT(SUM)
          SMED=SMED+EDD
          SAED=SAED+DABS(EDD)
          SQED=SQED+EDD**2
          num_data_infit=num_data_infit+1
       enddo !nde
    enddo !list of elements
    
    if(num_data_infit.GT.1) then
       write(*,'('' Number of data points in fit ='',I8)') num_data_infit
       !     write(*,'('' Average error           : '',D12.6,'' +/- '',D12.6)') &
       !          SMED/DBLE(num_data_infit), &
       !          DSQRT((SQED-SMED**2/DBLE(num_data_infit))/DBLE(num_data_infit-1))
       
       write(*,'('' Average absolute error  : '',D12.6,'' +/- '',D12.6)') &
            SAED/DBLE(num_data_infit),DSQRT((SQED-SAED**2/DBLE(num_data_infit))/ &
            DBLE(num_data_infit-1))
       write(*,'('' Root mean squared error : '',D12.6)') &
            DSQRT(SQED/DBLE(num_data_infit))
    else
       WRITE(*,'('' No data points in any elements'')')
       stop
    endif !ndtot>1
    
  end subroutine list_data_error
  
!!! ##########################################################################      
  
  subroutine map_versions(IPFILE,num_depvar,nynp,nyny,cyny,fit_soln,fix_bcs)

    use arrays,only: node_versn_2d,node_xyz_2d,num_nodes_2d
    use geometry,only: get_final_integer,get_local_node_f
!!! dummy arguments
    integer, intent(in) :: IPFILE,num_depvar
    integer :: nynp(:,:,:,:)
    integer,allocatable :: nyny(:,:)
    real(dp),allocatable :: cyny(:,:)
    real(dp) :: fit_soln(:,:,:,:)
    logical :: fix_bcs(:)
!!! local variables
    integer :: i,ibeg,iend,ierror,i_ss_end,nj,node,np,number_of_maps,nv, &
         NV_MAX,ny,nk_t,nv_t,nj_t,np_t,nk_m,nv_m,nj_m,np_m, &
         nmap_info(100,7),ny_t
    real(dp) :: r_map_coef
    character(len=132) :: string
    
!!! fix the boundary conditions for coordinates for nodes with versions, such that
!!! versions higher than 1 map to version 1     
    do np=1,num_nodes_2d
       NV_MAX=node_versn_2d(np)
       if(NV_MAX>1)then
          do nv=2,NV_MAX
             do nj=1,num_coords
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
!!!    default is that all versions are independent
    
    read_number_of_mappings : do
       read(unit=IPFILE, fmt="(a)", iostat=ierror) string
       ! read line containing "Number of mappings"
       if(index(string, "mappings")> 0) then
           call get_final_integer(string,number_of_maps)
          exit read_number_of_mappings
       endif
    end do read_number_of_mappings
    
    ! allocate memory for dependent variable mapping arrays
    write(*,*) 'Number of dependent variables =',num_depvar,'; squared =',num_depvar**2
    if(.not.allocated(cyny)) allocate(cyny(0:number_of_maps,num_depvar)) 
    if(.not.allocated(nyny)) allocate(nyny(0:number_of_maps,num_depvar))
    nyny = 0       ! initialise depvar to depvar mapping
    cyny = 0.0_dp  ! initialise weighting for mappings
    
!!! note that the global node numbers are used in the mapping file, whereas we need to use
!!! local numbering for the computation. Read in as global and then map to local below.
    do node=1,number_of_maps ! for the number of nodes with mappings
       read(unit=IPFILE, fmt="(a)", iostat=ierror) string
       ibeg=1
       i_ss_end=len(string) !get the end location of the sub-string
       do i=1,6
          iend=index(string," ") !get location of next blank in sub-string
          read (string(ibeg:iend-1), '(i6)' ) nmap_info(node,i)
          string = adjustl(string(iend:i_ss_end)) ! get the characters beyond " " and remove the leading blanks
       enddo !i
       read (string(ibeg:i_ss_end), '(i6)' ) nmap_info(node,7)
    enddo !node
    
    do node = 1,number_of_maps !for each mapping
       do nj = 1,num_coords
          nk_m = nmap_info(node,3)+1 !derivative
          nv_m = nmap_info(node,2) !version
          nj_m = nj !coordinate
          np_m = get_local_node_f(2,nmap_info(node,1)) !global node mapped to local node
          ny = nynp(nk_m,nv_m,nj_m,np_m)
          nk_t = nmap_info(node,6)+1 !derivative
          nv_t = nmap_info(node,5) !version
          nj_t = nj !coordinate
          np_t = get_local_node_f(2,nmap_info(node,4)) !global node mapped to local node
          ny_t = nynp(nk_t,nv_t,nj_t,np_t)
          r_map_coef = REAL(nmap_info(node,7)) !mapping coefficient, +1 or -1
          
          if(ny > 0) then
             nyny(0,ny) = nyny(0,ny)+1 ! increment array size
             nyny(nyny(0,ny),ny) = ny_t
             cyny(0,ny) = 0.0_dp
             cyny(nyny(0,ny),ny) = r_map_coef
             node_xyz_2d(nk_m,nv_m,nj_m,np_m) = node_xyz_2d(nk_t,nv_t,nj_t,np_t)*r_map_coef
             fit_soln(nk_m,nv_m,nj_m,np_m) = node_xyz_2d(nk_t,nv_t,nj_t,np_t)*r_map_coef
!             write(*,*) 'mapping ny',ny_t,' to',ny,' with',r_map_coef
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
    
  end subroutine map_versions
  
!!! ##########################################################################      

  subroutine melgef(LGE2,ne,NHST,nynp)

!!! calculates the row numbers (LGE(*,1)) and column numbers
!!! (LGE(*,2)) in the matrix for fitting for element variables nhs
!!! and fit variable njj in region nr.  It also returns the total
!!! number of element variables NHST(nrc).
    
    use arrays,only: elem_nodes_2d,elem_versn_2d
!!! dummy arguments    
    integer :: LGE2(num_fit*num_deriv_elem,2),ne,NHST(2),nynp(:,:,:,:)
!!! local variables
    integer nh,nk,nn,np,nrc,nv
    
    do nrc=1,2
       NHST(nrc)=0
       do nh=1,num_fit
          do nn=1,num_elem_nodes !nodal variables
             np=elem_nodes_2d(nn,ne)
             nv=elem_versn_2d(nn,ne)
             do nk=1,num_deriv
                NHST(nrc)=NHST(nrc)+1
                LGE2(NHST(nrc),nrc)=nynp(nk,nv,nh,np)
             enddo !nk
          enddo !nn
       enddo !nhj
    enddo !nrc
    
  end subroutine melgef

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
       psi1 = psi1*ph3(inp(nn,ni),ido(nk,ni),ipu(nu,ni),xi(ni))
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
  
    use arrays,only: node_versn_2d,node_xyz_2d,num_nodes_2d
    
!!! local variables
    integer :: nj,nk,nk1,np,nv
    real(dp) :: SCALE,XD(3),ZERO_TOL=1.0e-12_dp
    
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
    
  end subroutine update_scale_factor_norm
  
!!! ##########################################################################      

  subroutine xpxe(ne,xe)

!!! copies geometry information from nodes into a local element array
    
    use arrays,only: elem_nodes_2d,elem_versn_2d,node_xyz_2d,scale_factors_2d
!!! dummy arguments
    integer,intent(in) :: ne
    real(dp) :: xe(:,:)
!!! local variablesK     Local Variables
    integer :: nj,nk,nn,np,ns,nv
    
    do nj=1,3
       ns=0
       do nn=1,num_elem_nodes
          np=elem_nodes_2d(nn,ne)
          nv=elem_versn_2d(nn,ne)
          do nk=1,num_deriv
             ns=ns+1
             xe(ns,nj)=node_xyz_2d(nk,nv,nj,np)*scale_factors_2d(ns,ne)
          enddo
       enddo
    enddo !nj
    
  end subroutine xpxe

!!! ##########################################################################      
  
  subroutine zder(data_on_elem,ndata_on_elem,ne,data_xi,ER,PG,WDL,WG,sobelov_wts,&
       XIDL,fit_soln_local)
    
!!!    Evaluates element rhs, ER(ns), in calculation of least squares
!!!    fit of linear field variables, defined by nodal values
!!!    node_xyz_2d(nk,nv,nj,np), to the set of data values data_xyz(nj,nd) with
!!!    weights data_weight(nj,nd) at local coordinate values data_xi(ni,nd).

    use arrays,only: data_weight,data_xyz,scale_factors_2d
!!! dummy arguments
    integer :: data_on_elem(:,:),ndata_on_elem(:),ne
    real(dp) :: data_xi(:,:),ER(:),PG(:,:,:),WDL(:,:),WG(:),sobelov_wts(0:,:),&
         XIDL(:,:),fit_soln_local(:,:)
!!! local variables
    integer nd,nde,ng,nh,nhs1,nk1,nn1,ns1,ns2,nu
    real(dp) :: SUM1,SUM2,SUM3,SUM4,X,ZDL(3,nmax_data_elem)
    
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
                SUM1 = SUM1+PSI1(1,nk1,nn1,XIDL(1:2,nde))*ZDL(nh,nde)*WDL(nh,nde)
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
    
  end subroutine zder
  
!!! ##########################################################################      
  
  subroutine zdes(ndata_on_elem,ne,ES,PG,WDL,WG,sobelov_wts,XIDL)

!!!    ZDES evaluates element stiffness matrix ES(ms,ns) in calculation
!!!    of least squares fit of linear field variables, defined by nodal
!!!    values node_xyz_2d(nk,nv,nj,np), to the set of data values XD(nj,nd) with
!!!    weights data_weight(nj,nd) at local coordinate values data_xi(ni,nd), where
!!!    nj=NJO.

    use arrays,only: scale_factors_2d
!!! dummy arguments    
    integer :: ndata_on_elem(:),ne
    real(dp) :: ES(:,:),PG(:,:,:),WDL(:,:),WG(:),sobelov_wts(0:,:),XIDL(:,:)
!!! local variables
    integer nde,ng,nh1,nh2,nhj1,nhj2,nhs1,nhs1_for_nhj1,nhs2, &
         nk1,nk2,nn1,nn2,ns1,ns2,nu
    real(dp) :: PD(num_deriv_elem),SUM2,SUM3
    
    ES = 0.0_dp
    nhs1 = 0
    ! for each of the 3 dependent variables to be fitted (num_fit(1)=3)
    do nhj1=1,num_fit !nhj are vars for the fit problem njj
       nh1=nhj1
       nhs1_for_nhj1 = nhs1
       do nde=1,ndata_on_elem(ne)
          nhs1 = nhs1_for_nhj1
          ns1=0
          do nn1=1,num_elem_nodes
             do nk1=1,num_deriv
                nhs1=nhs1+1
                ns1=ns1+1
                PD(ns1)=PSI1(1,nk1,nn1,XIDL(1:2,nde))
             enddo !nk1
          enddo !nn1
          nhs1 = nhs1_for_nhj1
          do ns1=1,num_deriv_elem
             nhs1=nhs1+1
             nhs2=0
             do nhj2=1,num_fit !columns
                nh2=nhj2
                do ns2=1,num_deriv_elem
                   nhs2=nhs2+1
                   if(nhj2.EQ.nhj1) then !to avoid coupling for now
                      ES(nhs1,nhs2)=ES(nhs1,nhs2)+PD(ns1)*PD(ns2) &
                           *WDL(nh1,nde)*scale_factors_2d(ns1,ne)*scale_factors_2d(ns2,ne)
                   endif !nhj2=nhj1
                enddo !ns2
             enddo !nhj2
          enddo !ns1
       enddo !nde
       
       ns1=0
       nhs1 = nhs1_for_nhj1
       do nn1=1,num_elem_nodes
          do nk1=1,num_deriv
             nhs1=nhs1+1
             ns1=ns1+1
             nhs2=0
             do nhj2=1,num_fit !columns
                ns2=0
                do nn2=1,num_elem_nodes
                   do nk2=1,num_deriv
                      nhs2=nhs2+1
                      ns2=ns2+1
                      if(nhj2.EQ.nhj1) then !to avoid coupling for now
                         SUM2=0.0_dp
                         do ng=1,num_gauss
                            SUM3=0.0_dp
                            do nu=2,6
                               SUM3=SUM3+ &
                                    PG(ns1,nu,ng)*PG(ns2,nu,ng)*sobelov_wts(nu,ne)
                            enddo !nu
                            SUM2=SUM2+SUM3*WG(ng)
                         enddo !ng
                         ES(nhs1,nhs2)=ES(nhs1,nhs2)+(SUM2*sobelov_wts(0,ne))* &
                              scale_factors_2d(ns1,ne)*scale_factors_2d(ns2,ne)
                      endif !nhj2=nhj1
                   enddo !nk2
                enddo !nn2
             enddo !nhj2
          enddo !nk1
       enddo !nn1
    enddo !nhj1
    
  end subroutine zdes

!!! ##########################################################################      
    
  subroutine zpze_fit(ne,fit_soln_local,fit_soln)

    use arrays,only: elem_nodes_2d,elem_versn_2d,&
         scale_factors_2d
!!! dummy arguments
    integer,intent(in) :: ne
    real(dp) :: fit_soln_local(:,:)
    real(dp) :: fit_soln(:,:,:,:)
!!! local variables
    integer :: nh,nk,nn,np,ns,nv
    
    do nh=1,num_fit
       ns=0
       do nn=1,num_elem_nodes
          np=elem_nodes_2d(nn,ne)
          nv=elem_versn_2d(nn,ne)
          do nk=1,num_deriv
             ns=ns+1
             fit_soln_local(ns,nh)=fit_soln(nk,nv,nh,np)*scale_factors_2d(ns,ne)
          enddo !nk
       enddo !nn
    enddo !nhx
    
  end subroutine zpze_fit
  
!!! ##########################################################################      
  
  subroutine calculate_ny_maps(npny,num_depvar,nynp,nynr)
  
    use arrays,only: node_versn_2d,num_nodes_2d
!!! dummy arguments
    integer :: npny(0:,:),num_depvar,nynp(:,:,:,:),nynr(0:)
!!! local variables
    integer nh,nk,np,nv,ny
    
    !***  Initialise mapping arrays 
    nynp=0
    npny=0
    nynr=0
    
    !***  Set up mapping arrays
    ny = 0
    nynr(0)=0
    do nh=1,num_fit
       do np=1,num_nodes_2d
          do nv=1,node_versn_2d(np)
             do nk=1,num_deriv
                ny=ny+1
                nynr(0)=nynr(0)+1
                nynr(nynr(0))=ny
                nynp(nk,nv,nh,np) = ny
                npny(0,ny)=1 !mesh dof is node based
                npny(1,ny)=nk
                npny(2,ny)=nv
                npny(3,ny)=nh
                npny(4,ny)=np
                npny(5,ny)=1
             enddo !nk
          enddo !nv
       enddo !np
    enddo !njj
    num_depvar = ny

  end subroutine calculate_ny_maps

!!! ##########################################################################      
    
  subroutine define_2d_elements(ELEMFILE)

    use arrays,only: elems_2d,elem_nodes_2d,elem_versn_2d,node_versn_2d,num_elems_2d
    use geometry,only: get_final_integer,get_four_nodes,element_connectivity_2d
    character(len=*) :: ELEMFILE
    
    !     Local Variables
    integer :: ierror,ne,nn,noelem,np,number_of_elements
    character(len=132) :: ctemp1
    
    
    open(10, file=ELEMFILE, status='old')
    
    read_number_of_elements : do
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "elements")> 0) then
          call get_final_integer(ctemp1,number_of_elements)
          exit read_number_of_elements
       endif
    end do read_number_of_elements

    num_elems_2d=number_of_elements
    if(.not.allocated(elems_2d)) allocate(elems_2d(num_elems_2d))
    if(.not.allocated(elem_nodes_2d)) allocate(elem_nodes_2d(4,num_elems_2d))
    if(.not.allocated(elem_versn_2d)) allocate(elem_versn_2d(4,num_elems_2d))
    
    noelem=1
    
    read_an_element : do 
       !.......read element number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Element")> 0) then
          call get_final_integer(ctemp1,ne) !get element number
          elems_2d(noelem)=ne
          noelem=noelem+1
          
          read_element_nodes : do 
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "global")> 0) then !found the correct line
                call get_four_nodes(ne,ctemp1) !number of versions for node np
                ! note that only the ne'th data of elem_nodes_2d is passed to 'get_four_nodes'
                do nn=1,4
                   np=elem_nodes_2d(nn,ne)
                   if(node_versn_2d(np).gt.1)then
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !contains version# for njj=1
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !contains version# for njj=1
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !contains version# for njj=1
                      call get_final_integer(ctemp1,elem_versn_2d(nn,ne)) !get version#
                   else
                      elem_versn_2d(nn,ne)= 1
                   endif !nversions
                enddo !nn
                exit read_element_nodes
             endif !index
          end do read_element_nodes
          
          if(noelem.gt.number_of_elements) exit read_an_element
       endif
       
    end do read_an_element
    
    close(10)
    
    call element_connectivity_2d
    call line_segments_for_2d_mesh
    
  end subroutine define_2d_elements
  
!!! ##########################################################################      

  subroutine define_xi_closest(data_elem,data_on_elem,ndata_on_elem,data_xi,first)
!!! find the closest xi location on a 2d mesh surface to each data point
    implicit none

!!! dummy arguments    
    integer :: data_elem(:),data_on_elem(:,:),ndata_on_elem(:)
    real(dp) :: data_xi(:,:)
    logical,intent(in) :: first
!!! local variables
    integer :: i,n_check,ne_checklist(5),IT,ITMAX=20,nd,ne,neadj,nelast,neold,ni,nj
    real(dp) :: sqmax,sqnd,temp,xe(num_deriv_elem,num_coords),xi(3)
    real(dp),allocatable :: sq(:)
    logical :: found

    integer :: n_data
    character(len=200) :: exfile
    character(len=1) :: string_ne1
    character(len=2) :: string_ne2
    
    allocate(sq(num_data))    

    sqmax = 1.0e4_dp*1.0e4_dp
    
    !  initialise
    sq = 0.0_dp
    xi = 0.5_dp
    
!!! start by finding the closest centre of an element to each data point
!    do ne = 1,num_elems_2d
!       call xpxe(ne,xe)
!       do nd = 1,num_data
!          sqnd = 0.0_dp
!          do nj=1,num_coords
!             temp = pxi(1,xi,xe(1,nj))-data_xyz(nj,nd)
!             sqnd = sqnd+temp**2
!          enddo !nj
!          if(data_elem(nd).eq.0.or.sqnd.lt.sq(nd)) then
!             data_xi(1:2,nd) = xi(1:2)
!             if(nd.eq.3976)then
!                write(*,*) 'initial closest=',ne
!             endif
!             data_elem(nd) = ne
!             sq(nd) = sqnd
!          endif
!       enddo ! nd
!    enddo ! ne

    if(first)then ! check every element for every data point
       do nd = 1,num_data
          sqmax = 1.0e4_dp*1.0e4_dp
          do ne = 1,num_elems_2d
             xi = 0.5_dp
             call xpxe(ne,xe)
             found = .false.
             call project_orthogonal(nd,SQND,xe,xi,found)
             if(abs(xi(1)).ge.-zero_tol.and.abs(xi(1)).lt.1.0_dp+zero_tol.and. &
                  abs(xi(2)).ge.zero_tol.and.abs(xi(2)).lt.1.0_dp+zero_tol) then
                if(sqnd.lt.sqmax)then
                   sqmax = sqnd
                   data_xi(1:2,nd) = xi(1:2)
                   data_elem(nd) = ne
                   SQ(nd) = SQND
                endif
             endif !FOUND
          enddo
       enddo ! nd
    else

       do nd = 1,num_data
          ne = data_elem(nd)
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
             call xpxe(ne,xe)
             found = .true. !find nearest point in element
             call project_orthogonal(nd,sqnd,xe,xi,found)
             if(abs(xi(1)).ge.-zero_tol.and.abs(xi(1)).lt.1.0_dp+zero_tol.and. &
                  abs(xi(2)).ge.zero_tol.and.abs(xi(2)).lt.1.0_dp+zero_tol) then
                if(sqnd.lt.sqmax)then
                   sqmax = sqnd
                   data_xi(1:2,nd) = xi(1:2)
                   data_elem(nd)=ne
                   sq(nd) = sqnd
                endif
             endif
          enddo !i
       enddo ! nd
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

    exfile = 'temp.exdata'
    open(10, file = exfile, status = 'replace')

    do ne = 1,num_elems_2d
       !**   write the group name
       if(ne.lt.10)then
          write(string_ne1,'(i1)') ne
          write(10,'( '' Group name: '',A)') 'datapoints_'//string_ne1
       else
          write(string_ne2,'(i2)') ne
          write(10,'( '' Group name: '',A)') 'datapoints_'//string_ne2
       endif
       write(10,'(1X,''#Fields=1'')')
       write(10,'(1X,''1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
       write(10,'(1X,''  x.  Value index= 1, #Derivatives=0'')')
       write(10,'(1X,''  y.  Value index= 2, #Derivatives=0'')')
       write(10,'(1X,''  z.  Value index= 3, #Derivatives=0'')')
    
       do n_data = 1,ndata_on_elem(ne)
          nd = data_on_elem(ne,n_data)
          write(10,'(1X,''Node: '',I9)') nd
          write(10,'(1X,3E13.5)')  (data_xyz(nj,nd),nj=1,num_coords)
       enddo !n_data
    enddo !ne
    close(10)
    
  end subroutine define_xi_closest
  

!!! ##########################################################################      

  subroutine solve_geometry_fit(data_on_elem,ndata_on_elem,num_depvar,&
       elem_list,not_1,not_2,npny,nynp,nynr,nyny,data_xi,cyny,sobelov_wts,&
       fit_soln,fix_bcs)

    use arrays,only: node_xyz_2d
!!! dummy arguments
    integer :: data_on_elem(:,:),ndata_on_elem(:),not_1,not_2,num_depvar,&
         elem_list(0:),npny(0:,:),nynp(:,:,:,:),nynr(0:),nyny(0:,:)
    real(dp) :: data_xi(:,:),cyny(0:,:),sobelov_wts(0:,:),fit_soln(:,:,:,:)
    logical :: fix_bcs(:)
!!! local variables
    integer :: l,LGE2(3*16,2),ne,nh,nhs1,nhs2,NHST(2), &
         nk,no1,no2,no_nynr1,no_nynr2,noy1,noy2,np,nv,ny1,ny2,ny3,nyo1,nz,nzz
    integer,allocatable :: nony(:,:,:)
    integer,allocatable :: nyno(:,:,:)
    real(dp) :: co1,co2,ER(num_fit*num_deriv_elem),ES(3*16,3*16),&
         fit_soln_local(16,3),PG(16,6,9),WG(9)
    real(dp),allocatable :: cony(:,:,:)
    real(dp),allocatable :: cyno(:,:,:)
    real(dp),allocatable :: GR(:)      ! right-hand-side vector
    real(dp),allocatable :: GRR(:)     ! reduced right-hand-side vector
    real(dp),allocatable :: incr_soln(:)         ! current solution returned from solver
    logical :: FIRST_A,UPDATE_MATRIX

    ! make all of these allocatable!
    real(dp),dimension(nsize_gkk) :: GKK
    real(dp),dimension(nsize_gkk) :: GK
! doesn't like allocating these! gives different answer for errors
!    real(dp),allocatable :: GKK(:)
!    real(dp),allocatable :: GK(:)
    real(dp),dimension(3,nmax_data_elem) :: WDL
    real(dp),dimension(2,nmax_data_elem) :: XIDL
  
    
    WG = [7.7160493827160628e-2_dp, 0.12345679012345677_dp, 7.7160493827160628e-2_dp,&
         0.12345679012345677_dp, 0.19753086419753044_dp, 0.12345679012345677_dp,&
         7.7160493827160628e-2_dp, 0.12345679012345677_dp, 7.7160493827160628e-2_dp]

    allocate(incr_soln(num_depvar))
    allocate(nony(0:1,num_depvar,2))
    allocate(nyno(0:5,num_depvar,2))
    allocate(cony(0:1,num_depvar,2))
    allocate(cyno(0:5,num_depvar,2))
    allocate(GR(num_depvar))
    allocate(GRR(num_depvar))
!    allocate(GK(num_depvar*num_depvar))
!    allocate(GKK(num_depvar*num_depvar))

    call gauss1(PG)

    UPDATE_MATRIX=.TRUE.
    FIRST_A=.TRUE.
    
    GR=0.0_dp
    
    do l=1,elem_list(0) !loop over elements in the fit
       ne=elem_list(l)
       call melgef(LGE2,ne,NHST,nynp)
       ER=0.0_dp
       ES=0.0_dp
       call zpze_fit(ne,fit_soln_local,fit_soln) !gets fit_soln_local for element ne
       call zder(data_on_elem,ndata_on_elem,ne,data_xi,ER,PG,WDL,WG,&
            sobelov_wts,XIDL,fit_soln_local) 
       call zdes(ndata_on_elem,ne,ES,PG,WDL,WG,sobelov_wts,XIDL)
       
       !*** Assemble element stiffness matrix into global system.
       do nhs1=1,NHST(1) !3 dependent variables
          ny1=IABS(LGE2(nhs1,1))
          if(ny1.eq.0)then
             write(*,'('' No dependent variable for node in element'',i6,'': are &
                  &you sure you have set up versions correctly?'')') ne
             stop
          endif
          GR(ny1)=GR(ny1)+ER(nhs1)
          do nhs2=1,NHST(2) !3 dependent variables
             ny2=IABS(LGE2(nhs2,2))
             nz=ny1+(ny2-1)*num_depvar
             GK(nz)=GK(nz)+ES(nhs1,nhs2)
          enddo !nhs2
       enddo !nhs1
    enddo !l (ne)
    
    !*** Calculate solution mapping arrays for the current fit variable
    call globalf(nony,not_1,not_2,npny,nyno,nynp,nyny,cony,cyno,cyny,fix_bcs)
    
    if(NOT_2.EQ.0) then
       write(*,'('' >>The number of unknowns is zero'')')
       stop
    endif
    
    !----------------------- generate reduced system -----------------------
    
    GKK=0.0_dp
    GRR=0.0_dp
    
    !*** generate the reduced system of equations
    do no_nynr1=1,nynr(0) !loop global rows of GK
       ny1=nynr(no_nynr1) !is row #
       do noy1=1,nony(0,ny1,1) !loop over #no's attached to ny1
          no1=nony(noy1,ny1,1) !no# attached to row ny1
          co1=cony(noy1,ny1,1) !coupling coeff for row mapping
          !                     ie row_no1=a*row_ny1+b*row_ny2
          GRR(no1)=GRR(no1)+GR(ny1)*co1 !get reduced R.H.S.vector
          do no_nynr2=1,nynr(0) !loop over #cols of GK
             ny2=nynr(no_nynr2) !is global variable #
             ny3=getnyr(npny,ny2,nynp)
             !local GK var #
             nz=ny1+(ny3-1)*num_depvar
             if(nz.NE.0) then
                do noy2=1,nony(0,ny2,2) !loop over #no's for ny2
                   no2=nony(noy2,ny2,2) !no# attached to ny2
                   co2=cony(noy2,ny2,2) !coup coeff col mapping
                   !                     i.e. var_no1=a*var_ny1+b*var_ny2
                   nzz=no1+(no2-1)*NOT_1
                   write(*,*) 'nzz',nzz
                   if(nzz.NE.0) GKK(nzz)=GKK(nzz)+GK(nz)*co1*co2
                enddo !noy2
             endif
          enddo !no_nynr2
       enddo !noy1
    enddo !no_nynr1
    
    !-------------- solve reduced system of linear equations ---------------
    !Commented out since subroutines called further are temporarily unavailable
    write(*,*) NOT_1,NOT_2, num_depvar, size(GKK), GKK(1:10), size(GRR), GRR(1:10)
    write(*,*) size(incr_soln), incr_soln(1:10)
   !pause
    !call direct_solver(NOT_1,NOT_1,NOT_2,num_depvar,GKK,GRR,incr_soln,FIRST_A)
    
    do no1=1,NOT_2 ! for each unknown
       do nyo1=1,nyno(0,no1,2)
          ny1=nyno(nyo1,no1,2) ! the dependent variable number
          co1=cyno(nyo1,no1,2) ! the weighting for mapped variables
          nk=npny(1,ny1)     ! derivative number
          nv=npny(2,ny1)     ! version number
          nh=npny(3,ny1)     ! dependent variable number
          np=npny(4,ny1)     ! node number
          fit_soln(nk,nv,nh,np) = fit_soln(nk,nv,nh,np) + incr_soln(no1)*co1 
          ! current fit solution = previous + increment
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
!    deallocate(GK)
!    deallocate(GKK)

  end subroutine solve_geometry_fit
  
!!! ##########################################################################        
  
  subroutine project_orthogonal(nd,SQ,xe,xi,inelem)

    use arrays,only: data_xyz
!!! dummy arguments        
    integer :: nd
    real(dp) :: sq,xe(num_deriv_elem,num_coords),xi(:)
    logical :: inelem
!!! local variables
    integer :: IT,ITMAX,BOUND(2),it2,ni,nifix,nj
    real(dp) :: LOOSE_TOL=1.0e-6_dp
    real(dp) :: DELTA,DET,D2SQV2,D2SQVW2,D2SQXI(2,2),D2ZXI(3,2,2),DSQXI(2), &
         DSQXI1,DSQXI2,DSQV,DSQVW,DZ(3),DZXI(3,2),EVMIN,EVMAX,H(2), &
         MU,SQLIN,SQDIFF,SQDPRED,TEMP,TEMP1,TEMP2,TOL, &
         TOL2,V(2),V1,V2,VMAX=1.0_dp,W,XILIN(2),Z(3)
    logical :: CONVERGED,ENFORCE(2),FREE,NEWTON
    
    ITMAX = 10 ! max # iterations to use
    DELTA = VMAX/4.0_dp
    TOL = 5.0_dp*LOOSE_TOL !must be > sqrt(eps) or SQLIN<=SQ check may not work
    TOL2 = TOL**2
    SQ = 0.0_dp
    do nj=1,num_coords
       Z(nj) = PXI(1,XI,XE(1,nj))
       DZ(nj) = Z(nj)-data_xyz(nj,nd)
       SQ = SQ+DZ(nj)**2
    enddo !nj
    IT=0
    CONVERGED=.FALSE.
    do WHILE(.NOT.CONVERGED.AND.IT.LT.ITMAX)
       DSQXI = 0.0_dp
       do nj=1,num_coords
          DZXI(nj,1) = PXI(2,XI,XE(1,nj))
          DZXI(nj,2) = PXI(4,XI,XE(1,nj))
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
          D2ZXI(nj,1,1) = PXI(3,XI,XE(1,nj))
          D2ZXI(nj,1,2) = PXI(5,XI,XE(1,nj))
          D2ZXI(nj,2,2) = PXI(5,XI,XE(1,nj))
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
                Z(nj)=PXI(1,XILIN,XE(1,nj))
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
    
    if(.NOT.inelem.AND.XI(1)>=0.0_dp.AND.XI(1)<=1.0_dp.AND.XI(2)>=0.0_dp.AND.XI(2)<=1.0_dp) then
       inelem=.TRUE.
    endif
    
  end subroutine project_orthogonal

!!! ##########################################################################      
  


end module surface_fitting
