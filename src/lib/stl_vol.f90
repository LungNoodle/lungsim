module geometry
  !*Brief Description:* This module handles all geometry read/write/generation.
  !
  !*LICENSE:*
  !
  !
  !
  !*Full Description:*
  !
  !This module handles all geometry read/write/generation.
  use arrays
  use diagnostics
  use indices
  use mesh_utilities
  use other_consts ! currently has pi
  use precision ! sets dp for precision
  
  implicit none
  
  !Module parameters
  
  !Module types
  
  !Module variables
  
  !Interfaces
  private
  public add_mesh
  public add_matching_mesh
  public append_units
  public define_1d_elements
  public define_elem_geometry_2d
  public define_mesh_geometry_test
  public define_node_geometry
  public define_node_geometry_2d
  public define_data_geometry
  public group_elem_parent_term
  public define_rad_from_file
  public define_rad_from_geom
  public element_connectivity_1d
  public element_connectivity_2d
  public inlist
  public split_datacloud
  public split_datacloudfromSTL
  public evaluate_ordering
  public get_final_real
  public splitwords
  public get_local_node_f
  public make_data_grid
  public make_2d_vessel_from_1d
  public reallocate_node_elem_arrays
  public set_initial_volume
  public triangles_from_surface
  public volume_of_mesh
  public write_geo_file
  public get_final_integer
  public get_four_nodes
  public write_elem_geometry_2d
  public write_node_geometry_2d
  
contains

!!!#############################################################################
  
  subroutine allocate_node_arrays(num_nodes)
    !*allocate_node_arrays:* allocate memory for arrays associated with 1D trees
    
    integer,intent(in) :: num_nodes
    ! Local variables
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'allocate_node_arrays'
    call enter_exit(sub_name,1)
    
    if(.not.allocated(nodes)) allocate (nodes(num_nodes))
    if(.not.allocated(node_xyz)) allocate (node_xyz(3,num_nodes))
    if(.not.allocated(node_field)) allocate (node_field(num_nj,num_nodes))
    if(.not.allocated(elems_at_node)) allocate(elems_at_node(num_nodes,0:3))
    nodes = 0 !initialise node index values
    node_xyz = 0.0_dp !initialise
    node_field = 0.0_dp !initialise
    
    call enter_exit(sub_name,2)
    
  end subroutine allocate_node_arrays
  
!!!#############################################################################

  subroutine add_mesh(AIRWAY_MESHFILE)
    !*add_mesh:* Reads in an ipmesh file and adds this mesh to the terminal
    ! branches of an existing tree geometry
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ADD_MESH" :: ADD_MESH

    character(len=MAX_FILENAME_LEN), intent(in) :: AIRWAY_MESHFILE
    ! Local parameters
    character(len=100) :: buffer
    integer, parameter :: fh = 15
    integer :: ios
    integer :: line
    integer :: i,ibeg,iend,i_ss_end,j,ne,ne0,ne_global,ne_parent,ne_start, &
         ngen_parent,np,np0,np_global,&
         num_elems_new,num_elems_to_add,num_nodes_new,nunit,nlabel
    integer,dimension(1000) :: element_temp,generation, &
         parent_element,symmetry_temp
    real(dp),dimension(1000) :: length,radius,a_A
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'add_mesh'
    call enter_exit(sub_name,1)

    ios = 0
    line = 0
    open(fh, file=AIRWAY_MESHFILE)
    
    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.
    ios=0
    line=0
    i=0 ! count the number of elements read in
    
    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       ! line contains: element, parent element, generation,
       !                symmetry, length, outer radius, a/A ratio
       ! note that a/A ratio is always 1 for the conducting airways
       if (ios == 0) then
          line = line + 1
          i=i+1
          i_ss_end = len(buffer)
          
          do nlabel = 1,7
             ibeg = index(buffer," ") + 1 !get location of first integer beyond ws in string
             buffer = adjustl(buffer(ibeg:i_ss_end)) ! get info beyond ws, remove leading ws
             iend = index(buffer," ") !get location of next ws in sub-string
             select case (nlabel)
             case (1)
                read (buffer(1:iend-1), *, iostat=ios) element_temp(i) ! not used??
             case(2)
                read (buffer(1:iend-1), *, iostat=ios) parent_element(i)
             case(3)
                read (buffer(1:iend-1), *, iostat=ios) generation(i)
             case(4)
                read (buffer(1:iend-1), *, iostat=ios) symmetry_temp(i)
             case(5)
                read (buffer(1:iend-1), *, iostat=ios) length(i)
             case(6)
                read (buffer(1:iend-1), *, iostat=ios) radius(i)
             case(7)
                read (buffer(1:iend-1), *, iostat=ios) a_A(i)
             end select
          enddo
       endif
    enddo
    close(fh)
    
    num_elems_to_add = i
    
!!! increase the size of node and element arrays to accommodate the additional elements
    ! the number of nodes after adding mesh will be:
    num_nodes_new = num_nodes + num_units*num_elems_to_add
    ! the number of elems after adding mesh will be:
    num_elems_new = num_elems + num_units*num_elems_to_add
    call reallocate_node_elem_arrays(num_elems_new,num_nodes_new)
    
    ne = num_elems ! the starting local element number
    ne_global = elems(ne) ! assumes this is the highest element number (!!!)
    np = num_nodes ! the starting local node number
    np_global = nodes(np) ! assumes this is the highest node number (!!!)
    
    do nunit = 1,num_units ! for all terminal branches, append the mesh
       
       ne_parent = units(nunit) ! local element number of terminal, to append to
       ngen_parent = elem_ordrs(1,ne_parent)
       ne_start = ne !starting element number for the unit
       
       do i=1,num_elems_to_add
          
          if(parent_element(i).eq.0)then
             ne_parent = units(nunit)
          else
             ne_parent = ne_start+parent_element(i)
          endif
          
          ne0 = ne_parent
          np0 = elem_nodes(2,ne0)
          
          ne_global = ne_global + 1 ! new global element number
          ne = ne + 1 ! new local element number
          np_global = np_global + 1 !new global node number
          np = np + 1 ! new local node number
          
          nodes(np) = np_global
          elems(ne) = ne_global
          
          elem_nodes(1,ne) = np0
          elem_nodes(2,ne) = np
          
          elem_ordrs(1,ne) = ngen_parent + generation(i)
          elem_ordrs(no_type,ne) = 1   ! ntype ! 0 for respiratory, 1 for conducting
          elem_symmetry(ne) = symmetry_temp(i)+1 ! uses 0/1 in file; 1/2 in code
          
          ! record the element connectivity
          elem_cnct(-1,0,ne) = 1 ! one parent branch
          elem_cnct(-1,1,ne) = ne0 ! store parent element
          elem_cnct(1,0,ne0) = elem_cnct(1,0,ne0) + 1
          elem_cnct(1,elem_cnct(1,0,ne0),ne0) = ne
          
          ! record the direction and location of the branch
          do j=1,3
             elem_direction(j,ne) = elem_direction(j,ne0)
             node_xyz(j,np) = node_xyz(j,np0) + &
                  elem_direction(j,ne)*length(i)
          enddo !j
          
          elem_field(ne_length,ne) = length(i)
          elem_field(ne_radius,ne) = radius(i)
          elem_field(ne_a_A,ne) = a_A(i)
          elem_field(ne_vol,ne) = PI*radius(i)**2*length(i)
          
       enddo !i
    enddo !nunit
    
    num_nodes = np
    num_elems = ne
    
    call element_connectivity_1d
    call evaluate_ordering ! calculate new ordering of tree
    
    call enter_exit(sub_name,2)

  end subroutine add_mesh

!!!#############################################################################

  subroutine add_matching_mesh()
    !*add_matching_mesh:* 
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ADD_MATCHING_MESH" :: ADD_MATCHING_MESH

    !Parameters to become inputs
    real(dp) :: offset(3)
    logical :: REVERSE=.TRUE.
    character(len=60) :: mesh_type='terminal'
    !local variables
    integer :: num_nodes_new,num_elems_new,ne,ne_global,np,np_global,np0, &
         nonode,np_m
    integer :: nj,ne_m,noelem,ne0,n,nindex,ne1,noelem0,nu,cap_conns, &
         cap_term,np1,np2
    integer, allocatable :: np_map(:)
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'add_matching_mesh'
    call enter_exit(sub_name,1)
    !Ultimately offset should be an input argument
    offset(1)=0.0_dp
    offset(2)=1.0e-6_dp
    offset(3)=0.0_dp

    allocate(np_map(num_nodes))
!!! increase the size of node and element arrays to accommodate the additional elements
    ! the number of nodes after adding mesh will be:
    num_nodes_new = 2*num_nodes
    ! the number of elems after adding mesh will be:
    if(mesh_type.eq.'basic')then
       num_elems_new = 2*num_elems
    elseif(mesh_type.eq.'terminal')then
       num_elems_new = 2*num_elems + num_units
    elseif(mesh_type.eq.'ladder')then
       print *, "NOT YET IMPLEMENTED"
       call exit(0)
    endif
    call reallocate_node_elem_arrays(num_elems_new,num_nodes_new)
    noelem0=0
    ne0 = num_elems ! the starting local element number
    ne_global = elems(ne0) ! assumes this is the highest element number (!!!)
    np0 = num_nodes ! the starting local node number
    np_global = nodes(np0) ! assumes this is the highest node number (!!!)
    
    do nonode=1,num_nodes
       np=np_global+nonode
       np_m=nodes(nonode)
       np_map(np_m)=np !maps new to old node numbering
       nodes(np0+nonode)=np
       do nj=1,3
          node_xyz(nj,np)=node_xyz(nj,np_m)+offset(nj)
       enddo
       elems_at_node(np,0)=0 !initialise
       !Doesnt map versions, would be added here
    enddo
    
    do noelem=1,num_elems
       ne=ne_global+noelem
       elem_field(ne_group,ne)=2.0_dp!VEIN
       ne_m=elems(noelem)
       elem_field(ne_group,ne_m)=0.0_dp!ARTERY
       elems(ne0+noelem)=ne
       if(.NOT.REVERSE)then
          elem_nodes(1,ne)=np_map(elem_nodes(1,ne_m))
          elem_nodes(2,ne)=np_map(elem_nodes(2,ne_m))
          elem_cnct(1,0,ne)=elem_cnct(1,0,ne_m)!The numberdownstream are the number downstream
          elem_cnct(-1,0,ne)=elem_cnct(-1,0,ne_m)
          do n=1,elem_cnct(1,0,ne)
             elem_cnct(1,n,ne)=elem_cnct(1,n,ne_m)+ne0
          enddo
          do n=1,elem_cnct(-1,0,ne)
             elem_cnct(-1,n,ne)=elem_cnct(-1,n,ne_m)+ne0
          enddo
       else
          elem_nodes(1,ne)=np_map(elem_nodes(2,ne_m))
          elem_nodes(2,ne)=np_map(elem_nodes(1,ne_m))
          elem_cnct(-1,0,ne)=elem_cnct(1,0,ne_m) !The number upstream are the number downstream
          elem_cnct(1,0,ne)=elem_cnct(-1,0,ne_m)!The number downstream are the number upstream
          do n=1,elem_cnct(1,0,ne)
             elem_cnct(1,n,ne)=elem_cnct(-1,n,ne_m)+ne0
          enddo
          do n=1,elem_cnct(-1,0,ne)
             elem_cnct(-1,n,ne)=elem_cnct(1,n,ne_m)+ne0
          enddo
       endif
       !if worrying about regions and versions do it here
       elems_at_node(elem_nodes(1,ne),0)=elems_at_node(elem_nodes(1,ne),0)+1
       elems_at_node(elem_nodes(1,ne),elems_at_node(elem_nodes(1,ne),0))=ne
       elems_at_node(elem_nodes(2,ne),0)=elems_at_node(elem_nodes(2,ne),0)+1
       elems_at_node(elem_nodes(2,ne),elems_at_node(elem_nodes(2,ne),0))=ne
       nindex=no_gen
       elem_ordrs(nindex,ne)=elem_ordrs(nindex,ne_m)
       nindex=no_sord
       elem_ordrs(nindex,ne)=elem_ordrs(nindex,ne_m)
       nindex=no_hord
       elem_ordrs(nindex,ne)=elem_ordrs(nindex,ne_m)
    enddo
    
    !update current no of nodes and elements to determine connectivity
    np0=np !current highest node
    ne1=ne !current highest element
    noelem0=num_elems+noelem0
    if(mesh_type.eq.'ladder')then
       !To be implemented
    elseif(mesh_type.eq.'terminal')then
       cap_conns=0
       cap_term=0
       do nu=1,num_units
          ne=units(nu)
          cap_term=cap_term+1
          np1=elem_nodes(2,ne)
          np2=np_map(np1)
          noelem0=noelem0+1
          ne1=ne1+1
          elems(noelem0)=ne1
          elem_nodes(1,ne1)=np1
          elem_nodes(2,ne1)=np2
          elems_at_node(np1,0)=elems_at_node(np1,0)+1
          elems_at_node(np1,elems_at_node(np1,0))=ne1
          elems_at_node(np2,0)=elems_at_node(np2,0)+1
          elems_at_node(np2,elems_at_node(np2,0))=ne1
          elem_cnct(1,elem_cnct(1,0,ne)+1,ne)=ne1
          elem_cnct(1,0,ne)=elem_cnct(1,0,ne)+1
          elem_cnct(-1,elem_cnct(-1,0,ne+ne_global)+1,ne+ne_global)=ne1
          elem_cnct(-1,0,ne+ne_global)=elem_cnct(-1,0,ne+ne_global)+1
          elem_cnct(-1,0,ne1)=1
          elem_cnct(1,0,ne1)=1
          elem_cnct(-1,1,ne1)=ne
          elem_cnct(1,1,ne1)=ne+ne0
          nindex=no_gen
          elem_ordrs(nindex,ne1)=elem_ordrs(nindex,ne_m)
          nindex=no_sord
          elem_ordrs(nindex,ne1)=elem_ordrs(nindex,ne_m)
          nindex=no_hord
          elem_ordrs(nindex,ne1)=elem_ordrs(nindex,ne_m)
          elem_field(ne_group,ne1)=1.0_dp!connection between meshes
       enddo
       print *, 'Number of connections', cap_term
    endif
    num_nodes=num_nodes_new
    num_elems=num_elems_new
    deallocate(np_map)

    call enter_exit(sub_name,2)
    
  end subroutine add_matching_mesh

!!!#############################################################################

  subroutine append_units()
    !*append_units:* Appends terminal units at the end of a tree structure
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_APPEND_UNITS" :: APPEND_UNITS

    ! Local parameters
    integer :: ne,ne0,nu
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'append_units'
    call enter_exit(sub_name,1)

    num_units = 0
    do ne=1,num_elems
       if(elem_cnct(1,0,ne).eq.0)THEN
          num_units=num_units+1
       endif
    enddo
    
    if(allocated(units))then !increasing the array size; just overwrite
       deallocate(units)
       deallocate(unit_field)
    endif
    allocate(units(num_units))
    allocate(unit_field(num_nu,num_units))
    
    unit_field=0.0_dp
    units=0
    elem_units_below(1:num_elems) = 0 !initialise the number of terminal units below a branch
    
    nu=0
    do ne=1,num_elems
       if(elem_cnct(1,0,ne).eq.0)THEN
          nu=nu+1
          units(nu)=ne     !Set up units array containing terminals
          elem_units_below(ne)=1
       endif
    enddo
    
    ! count the effective number of elements below each branch
    do ne=num_elems,2,-1
       ne0=elem_cnct(-1,1,ne)
       elem_units_below(ne0) = elem_units_below(ne0) &
            + elem_units_below(ne)*elem_symmetry(ne)
    enddo !ne
    
    call enter_exit(sub_name,2)

  end subroutine append_units

!!!#############################################################################

  subroutine define_1d_elements(ELEMFILE)
    !*define_1d_elements:* Reads in an 1D element ipelem file to define a geometry
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_1D_ELEMENTS" :: DEFINE_1D_ELEMENTS
    
    character(len=MAX_FILENAME_LEN), intent(in) :: ELEMFILE
    !     Local Variables
    integer :: ibeg,iend,ierror,i_ss_end,j,ne,ne_global,&
         nn,np,np1,np2,np_global
    character(LEN=132) :: ctemp1
    character(len=300) :: readfile
    character(LEN=40) :: sub_string
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'define_1d_elements'
    call enter_exit(sub_name,1)
    
    if(index(ELEMFILE, ".ipelem")> 0) then !full filename is given
       readfile = ELEMFILE
    else ! need to append the correct filename extension
       readfile = trim(ELEMFILE)//'.ipelem'
    endif
    
    open(10, file=readfile, status='old')
    
    read_number_of_elements : do
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "elements")> 0) then
          num_elems = get_final_integer(ctemp1) !return the final integer
          exit read_number_of_elements
       endif
    end do read_number_of_elements
    
!!! allocate memory for element arrays
    if(allocated(elems)) deallocate(elems)
    allocate(elems(num_elems))
    if(allocated(elem_cnct)) deallocate(elem_cnct)
    allocate(elem_cnct(-1:1,0:2,0:num_elems))
    if(allocated(elem_nodes)) deallocate(elem_nodes)
    allocate(elem_nodes(2,num_elems))
    if(allocated(elem_ordrs)) deallocate(elem_ordrs)
    allocate(elem_ordrs(num_ord,num_elems))
    if(allocated(elem_symmetry)) deallocate(elem_symmetry)
    allocate(elem_symmetry(num_elems))
    if(allocated(elem_units_below)) deallocate(elem_units_below)
    allocate(elem_units_below(num_elems))
    if(allocated(elems_at_node)) deallocate(elems_at_node)
    allocate(elems_at_node(num_nodes,0:3))
    if(allocated(elem_field)) deallocate(elem_field)
    allocate(elem_field(num_ne,num_elems))
    if(allocated(elem_direction)) deallocate(elem_direction)
    allocate(elem_direction(3,num_elems))
    if(model_type.eq.'gas_mix')then
       if(allocated(expansile)) deallocate(expansile)
       allocate(expansile(num_elems))
    endif
    
!!! initialise element arrays
    elems = 0
    elem_nodes = 0
    elem_ordrs = 0  ! where the default is that 0==respiratory and 1==conducting
    elem_symmetry = 1
    elem_field = 0.0_dp
    if(model_type.eq.'gas_mix')expansile = .false.
    
    ne=0
    
    read_an_element : do
       !.......read element number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Element")> 0) then
          ne_global = get_final_integer(ctemp1) !return the final integer
          ne=ne+1
          elems(ne)=ne_global
          
          read_element_nodes : do
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "global")> 0) then !found the correct line
                iend=len(ctemp1)
                ibeg=index(ctemp1,":")+1 !get location of first integer in string
                sub_string = adjustl(ctemp1(ibeg:iend)) ! get the characters beyond : remove leading blanks
                i_ss_end=len(sub_string) !get the end location of the sub-string
                ibeg=1
                do nn=1,2
                   iend=index(sub_string," ") !get location of first blank in sub-string
                   read (sub_string(ibeg:iend-1), '(i7)' ) np_global
                   call get_local_node(np_global,np) ! get local node np for global node
                   elem_nodes(nn,ne)=np ! the local node number, not global
                   sub_string = adjustl(sub_string(iend:i_ss_end)) ! get chars beyond blank, remove leading blanks
                end do
                exit read_element_nodes
             endif !index
          end do read_element_nodes
          if(ne.ge.num_elems) exit read_an_element
       endif
       
    end do read_an_element
    
    close(10)
    
    ! calculate the element lengths and directions
    do ne=1,num_elems
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       elem_field(ne_length,ne) = sqrt((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)
       do j=1,3
          elem_direction(j,ne) = (node_xyz(j,np2) - &
               node_xyz(j,np1))/elem_field(ne_length,ne)
       enddo !j
    enddo
    
    call element_connectivity_1d
    call evaluate_ordering

    elem_ordrs(no_type,:) = 1 ! 0 for respiratory, 1 for conducting
    
    call enter_exit(sub_name,2)

  end subroutine define_1d_elements

!!!#############################################################################

  subroutine define_elem_geometry_2d(ELEMFILE,sf_option)
    ! Reads in 2D ipelem file.
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_ELEM_GEOMETRY_2D" :: DEFINE_ELEM_GEOMETRY_2D

    character(len=*) :: ELEMFILE
    character(len=4) :: sf_option
    !     Local Variables
    integer :: ierror,ne,ne_global,nn,np,number_of_elements
    character(len=132) :: ctemp1,readfile
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'define_elem_geometry_2d'
    call enter_exit(sub_name,1)
    
    if(index(ELEMFILE, ".ipelem")> 0) then !full filename is given
       readfile = ELEMFILE
    else ! need to append the correct filename extension
       readfile = trim(ELEMFILE)//'.ipelem'
    endif
    
    open(10, file=readfile, status='old')
    
    read_number_of_elements : do
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "elements")> 0) then
          number_of_elements = get_final_integer(ctemp1) !return the final integer
          exit read_number_of_elements
       endif
    end do read_number_of_elements
    
    num_elems_2d=number_of_elements
    if(.not.allocated(elems_2d)) allocate(elems_2d(num_elems_2d))
    if(.not.allocated(elem_nodes_2d)) allocate(elem_nodes_2d(4,num_elems_2d))
    if(.not.allocated(elem_versn_2d)) allocate(elem_versn_2d(4,num_elems_2d))
    
    ne = 0
    
    read_an_element : do 
       !.......read element number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Element")> 0) then
          ne_global = get_final_integer(ctemp1) !return the final integer
          ne = ne + 1
          elems_2d(ne) = ne_global
          
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
                      elem_versn_2d(nn,ne) = get_final_integer(ctemp1) !return the final integer
                   else
                      elem_versn_2d(nn,ne)= 1
                   endif !nversions
                enddo !nn
                exit read_element_nodes
             endif !index
          end do read_element_nodes
          
          if(ne.ge.number_of_elements) exit read_an_element
       endif
       
    end do read_an_element
    
    close(10)
    
    call element_connectivity_2d
    call line_segments_for_2d_mesh(sf_option)
    
    call enter_exit(sub_name,2)
    
  end subroutine define_elem_geometry_2d
  
!!!#############################################################################

  subroutine define_mesh_geometry_test()
    !*define_mesh_geometry_test:*
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_MESH_GEOMETRY_TEST" :: DEFINE_MESH_GEOMETRY_TEST

    !     Local Variables
    integer :: j,ne,np,np1,np2
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'define_mesh_geometry_test'
    call enter_exit(sub_name,1)
    
    num_elems = 400
    num_nodes = num_elems + 1
    
!!! allocate memory
    if(.not.allocated(nodes)) allocate (nodes(num_nodes))
    if(.not.allocated(node_xyz)) allocate (node_xyz(3,num_nodes))
    if(.not.allocated(node_field)) allocate (node_field(num_nj,num_nodes))
    if(.not.allocated(elems)) allocate(elems(num_elems))
    if(.not.allocated(elem_cnct)) allocate(elem_cnct(-1:1,0:2,0:num_elems))
    if(.not.allocated(elem_nodes)) allocate(elem_nodes(2,num_elems))
    if(.not.allocated(elem_ordrs)) allocate(elem_ordrs(num_ord,num_elems))
    if(.not.allocated(elem_symmetry)) allocate(elem_symmetry(num_elems))
    if(.not.allocated(elem_units_below)) allocate(elem_units_below(num_elems))
    if(.not.allocated(elems_at_node)) allocate(elems_at_node(num_nodes,0:3))
    if(.not.allocated(elem_field)) allocate(elem_field(num_ne,num_elems))
    if(.not.allocated(elem_direction)) allocate(elem_direction(3,num_elems))
    if(.not.allocated(expansile)) allocate(expansile(num_elems))
    
!!! initialise array values
    nodes = 0 !initialise node index values
    node_xyz = 0.0_dp !initialise
    node_field = 0.0_dp !initialise
    elems=0
    elem_nodes=0
    elem_symmetry = 1
    elem_field = 0.0_dp
    expansile = .false.
    
!!! set up node arrays
    nodes(1) = 1
    do np=2,101
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 1.0_dp
    enddo
    
    np=102
    node_xyz(1,np) = node_xyz(1,np-1) - 0.5_dp
    do np=102,151
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) - 0.5_dp
    enddo
    
    np=152
    nodes(np) = np
    node_xyz(1,np) = node_xyz(1,101) + 0.5_dp
    node_xyz(3,np) = node_xyz(3,101) - 0.5_dp
    do np=153,201
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) + 0.5_dp
    enddo
    
    np=202
    nodes(np) = np
    node_xyz(1,np) = node_xyz(1,151) - 0.5_dp
    node_xyz(3,np) = node_xyz(3,151) - 0.5_dp
    do np=203,251
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) - 0.5_dp
    enddo
    
    np=252
    nodes(np) = np
    node_xyz(1,np) = node_xyz(1,151) + 0.5_dp
    node_xyz(3,np) = node_xyz(3,151) - 0.5_dp
    do np=253,301
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) + 0.5_dp
    enddo
    
    np=302
    nodes(np) = np
    node_xyz(1,np) = node_xyz(1,201) - 0.5_dp
    node_xyz(3,np) = node_xyz(3,201) - 0.5_dp
    do np=303,351
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) - 0.5_dp
    enddo
    
    np=352
    nodes(np) = np
    node_xyz(1,np) = node_xyz(1,201) + 0.5_dp
    node_xyz(3,np) = node_xyz(3,201) - 0.5_dp
    do np=353,401
       nodes(np) = np
       node_xyz(3,np) = node_xyz(3,np-1) - 0.5_dp
       node_xyz(1,np) = node_xyz(1,np-1) + 0.5_dp
    enddo
    
!!! set up element arrays
    do ne=1,num_elems
       elems(ne) = ne
       elem_nodes(1,ne) = ne
       elem_nodes(2,ne) = ne+1
    enddo
    
    elem_nodes(1,151) = 101
    elem_nodes(1,201) = 151
    elem_nodes(1,251) = 151
    elem_nodes(1,301) = 201
    elem_nodes(1,351) = 201
    
    elem_field(ne_radius,1:100) = 10.0_dp
    elem_field(ne_radius,101:200) = 5.0_dp
    elem_field(ne_radius,201:400) = sqrt(elem_field(ne_radius,101)**2/2.0_dp)
    
    ! calculate the element lengths and directions
    do ne=1,num_elems
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       elem_field(ne_length,ne) = sqrt((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)
       do j=1,3
          elem_direction(j,ne) = (node_xyz(j,np2) - &
               node_xyz(j,np1))/elem_field(ne_length,ne)
       enddo !j
       ! element volume
       elem_field(ne_vol,ne) = PI * elem_field(ne_radius,ne)**2 * &
            elem_field(ne_length,ne)
       elem_field(ne_a_A,ne) = 1.0_dp ! set default for ratio a/A
    enddo
    
    call element_connectivity_1d
    call evaluate_ordering
    
    call enter_exit(sub_name,2)
    
  end subroutine define_mesh_geometry_test

!!!#############################################################################
  
  subroutine define_node_geometry(NODEFILE)
    !*define_node_geometry:* Reads in an ipnode file to define a tree geometry
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_NODE_GEOMETRY" :: DEFINE_NODE_GEOMETRY
    
    character(len=MAX_FILENAME_LEN), intent(in) :: NODEFILE !Input nodefile
    !     Local Variables
    integer :: i,ierror,np,np_global,num_nodes_temp,num_versions,nv,NJT=0
    real(dp) :: point
    logical :: overwrite = .false. ! initialised
    character(len=300) :: ctemp1,readfile
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'define_node_geometry'
    call enter_exit(sub_name,1)
    
    if(index(NODEFILE, ".ipnode")> 0) then !full filename is given
       readfile = NODEFILE
    else ! need to append the correct filename extension
       readfile = trim(NODEFILE)//'.ipnode'
    endif
    
    open(10, file=readfile, status='old')
    
    if(num_nodes.gt.0) overwrite = .true.
    
    !.....read in the total number of nodes. read each line until one is found
    !.....that has the correct keyword (nodes). then return the integer that is
    !.....at the end of the line
    read_number_of_nodes : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "nodes")> 0) then !keyword "nodes" is found in ctemp1
          num_nodes_temp = get_final_integer(ctemp1) !return the final integer
          exit read_number_of_nodes !exit the named do loop
       endif
    end do read_number_of_nodes
    
    if(.not.overwrite) call allocate_node_arrays(num_nodes_temp) ! don't allocate if just overwriting
    
    !.....read in the number of coordinates
    read_number_of_coords : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "coordinates")> 0) then !keyword "coordinates" is found
          NJT = get_final_integer(ctemp1) !return the final integer
          exit read_number_of_coords !exit the named do loop
       endif
    end do read_number_of_coords
    
    ! note that only the first version of coordinate is currently read in   
    
    !.....read the coordinate, derivative, and version information for each node. 
    np=0
    read_a_node : do !define a do loop name
       !.......read node number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Node")> 0) then
          np_global = get_final_integer(ctemp1) !get node number
          
          np = np+1
          nodes(np) = np_global
          !.......read coordinates
          do i=1,3 ! for the x,y,z coordinates
             num_versions=1
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "versions")> 0) then
                num_versions = get_final_integer(ctemp1)
                if(num_versions > 1)then
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !temporary line
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   point = get_final_real(ctemp1)
                   do nv=2,num_versions
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !temporary line
                      read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   enddo
                else
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   point = get_final_real(ctemp1)
                endif
             else ! no prompting for versions
                point = get_final_real(ctemp1)
             endif
             node_xyz(i,np)=point
          end do !i
       endif !index
       if(np.ge.num_nodes_temp) exit read_a_node
    end do read_a_node
    
    if(.not.overwrite) num_nodes = num_nodes_temp
    
    close(10)
    
    call enter_exit(sub_name,2)
    
  end subroutine define_node_geometry

!!!#############################################################################

  subroutine define_node_geometry_2d(NODEFILE)
    !*define_node_geometry_2d:* Reads in an ipnode file to define surface nodes
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_NODE_GEOMETRY_2D" :: DEFINE_NODE_GEOMETRY_2D
    
    character(len=*),intent(in) :: NODEFILE
    !     Local Variables
    integer :: i,ierror,np,np_global,&
         num_versions,nv
    integer,parameter :: num_derivs = 3
    character(len=132) :: ctemp1,readfile
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'define_node_geometry_2d'
    call enter_exit(sub_name,1)
    
    if(index(NODEFILE, ".ipnode")> 0) then !full filename is given
       readfile = NODEFILE
    else ! need to append the correct filename extension
       readfile = trim(NODEFILE)//'.ipnode'
    endif
    
    open(10, file=readfile, status='old')
    
    !.....read in the total number of nodes. read each line until one is found
    !.....that has the correct keyword (nodes). then return the integer that is
    !.....at the end of the line
    read_number_of_nodes : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "nodes")> 0) then !keyword "nodes" is found in ctemp1
          num_nodes_2d = get_final_integer(ctemp1) !return the final integer
          exit read_number_of_nodes !exit the named do loop
       endif
    end do read_number_of_nodes
    
!!!allocate memory to arrays that require node number
    if(.not.allocated(nodes_2d)) allocate(nodes_2d(num_nodes_2d))
    if(.not.allocated(node_xyz_2d)) allocate(node_xyz_2d(4,10,16,num_nodes_2d))
    if(.not.allocated(node_versn_2d)) allocate(node_versn_2d(num_nodes_2d))
    
    !.....read the coordinate, derivative, and version information for each node. 
    np=0
    read_a_node : do !define a do loop name
       !.......read node number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Node")> 0) then
          np_global = get_final_integer(ctemp1) !get node number
          
          np=np+1
          nodes_2d(np) = np_global
          
          !.......read coordinates
          do i=1,3 ! for the x,y,z coordinates
             num_versions = 0
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "versions")> 0) num_versions = &
                  get_final_integer(ctemp1)
             node_versn_2d(np) = max(1,num_versions) !number of versions for node np
             do nv=1,node_versn_2d(np)
                if(num_versions > 1)then
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1 ! "For version number..."
                endif
                !...........coordinate          
                if(num_versions > 0) &
                     read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                node_xyz_2d(1,nv,i,np) = get_final_real(ctemp1)
                if(num_derivs.ge.1)then
                   !..........derivative 1
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   node_xyz_2d(2,nv,i,np) = get_final_real(ctemp1)
                endif
                if(num_derivs.ge.2)then
                   !..........derivative 2
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   node_xyz_2d(3,nv,i,np) = get_final_real(ctemp1)
                endif
                if(num_derivs.ge.3)then
                   !...........derivative 1&2
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                   node_xyz_2d(4,nv,i,np) = get_final_real(ctemp1)
                endif
                if(num_derivs.ge.4)then
                   write(*,'(''This code is only valid for a surface geometry'')')
                   read(*,*)
                endif
             end do !nv
          end do !i
       endif !index
       if(np.ge.num_nodes_2d) exit read_a_node
    end do read_a_node
    
    close(10)
    
    call enter_exit(sub_name,2)
    
  end subroutine define_node_geometry_2d

!!!#############################################################################

  subroutine define_data_geometry(datafile,is_field,number_of_fields)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_DATA_GEOMETRY" :: DEFINE_DATA_GEOMETRY

!! CAUTION
!! number_of_fields verified for 1

!!! read data points from a file
    
    use arrays,only: dp,data_xyz,data_weight,num_data,field_xyz
    use diagnostics, only: enter_exit
    use indices
    use other_consts, only: MAX_FILENAME_LEN   
!!! dummy arguments
    character(len=*) :: datafile
!!! local variables
    integer :: iend,ierror,length_string,ncount,nj,itemp,number_of_fields
    character(len=132) :: buffer,readfile
    logical,intent(in) :: is_field
    character(len=60) :: sub_name = 'define_data_geometry'

    ! --------------------------------------------------------------------------

    sub_name = 'define_data_geometry'
    call enter_exit(sub_name,1)
    
    !set the counted number of data points to zero
    ncount = 0
    
    !readfile = trim(datafile)//'.ipdata'
    open(10, file=datafile, status='old')
    read(unit=10, fmt="(a)", iostat=ierror) buffer
    
!!! first run through to count the number of data points
    read_line_to_count : do
       read(unit=10, fmt="(a)", iostat=ierror) buffer
       if(ierror<0) exit !ierror<0 means end of file
       ncount = ncount + 1
    end do read_line_to_count
    num_data = ncount
    close (10)
    write(*,'('' Read'',I7,'' data points from file'')') num_data
    
!!! allocate arrays now that we know the size required
    if(allocated(data_xyz)) deallocate(data_xyz)
    if(allocated(data_weight)) deallocate(data_weight)
    allocate(field_xyz(num_data))

    allocate(data_xyz(3,num_data))
    allocate(data_weight(3,num_data))
    
!!! read the data point information
    !readfile = trim(datafile)//'.ipdata'
    open(10, file=datafile, status='old')
    read(unit=10, fmt="(a)", iostat=ierror) buffer
    
    !set the counted number of data points to zero
    ncount = 0
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
       if(is_field) then
          ! read density field
          buffer = adjustl(buffer(iend+1:length_string)) !remove data number from string
          buffer = adjustl(buffer) !remove leading blanks
          length_string = len(buffer) !new length of buffer
          iend=index(buffer," ",.false.)-1 !index returns location of first blank
          read (buffer(1:iend), '(D25.17)') field_xyz(ncount)
       endif
       
       do nj=1,3
          !        ! read weightings
          !        buffer = adjustl(buffer(iend+1:length_string)) !remove data number from string
          !        buffer = adjustl(buffer) !remove leading blanks
          !        length_string = len(buffer) !new length of buffer
          !        iend=index(buffer," ",.false.)-1 !index returns location of first blank
          !        read (buffer(1:iend), '(D25.17)') data_weight(nj,ncount)
          data_weight(nj,ncount)=1.0_dp
       enddo !nj
       
    enddo read_line_of_data
    
    close(10)
    
    call enter_exit(sub_name,2)

  end subroutine define_data_geometry

!!!#############################################################################

  subroutine triangles_from_surface(num_triangles,num_vertices,surface_elems,triangle,vertex_xyz)
    !*triangles_from_surface:* generates a linear surface mesh of triangles
    ! from an existing high order surface mesh. 
    
    integer :: num_triangles,num_vertices
    integer,intent(in) :: surface_elems(:)
    integer,allocatable :: triangle(:,:)
    real(dp),allocatable :: vertex_xyz(:,:)
    ! Local variables
    integer,parameter :: ndiv = 3
    integer :: i,index1,index2,j,ne,nmax_1,nmax_2,num_surfaces, &
         num_tri_vert,nvertex_row,step_1,step_2
    real(dp) :: X(3),xi(3)
    logical :: four_nodes
    character(len=3) :: repeat
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'triangles_from_surface'
    call enter_exit(sub_name,1)
    
    if(.not.allocated(triangle)) allocate(triangle(3,2*num_elems_2d*ndiv**2))
    if(.not.allocated(vertex_xyz)) allocate(vertex_xyz(3,num_elems_2d*(ndiv+1)**2))
    triangle = 0
    vertex_xyz = 0.0_dp
    num_surfaces = count(surface_elems.ne.0)
    num_triangles = 0
    num_vertices = 0
    num_tri_vert = 0 

    do ne=1,num_elems_2d
       four_nodes = .false.
       repeat = '0_0'
       if(elem_nodes_2d(1,ne).eq.elem_nodes_2d(2,ne)) repeat = '1_0'
       if(elem_nodes_2d(1,ne).eq.elem_nodes_2d(3,ne)) repeat = '2_0'
       if(elem_nodes_2d(2,ne).eq.elem_nodes_2d(4,ne)) repeat = '2_1'
       if(elem_nodes_2d(3,ne).eq.elem_nodes_2d(4,ne)) repeat = '1_1'

       select case(repeat)
       case ('0_0')
        
          nmax_1 = ndiv+1 ! ndiv+1 vertices in xi1 direction
          step_1 = 0      ! # of vertices in xi1 is constant
          nmax_2 = ndiv+1 ! ndiv+1 vertices in xi2 direction
          step_2 = 0      ! # of vertices in xi2 is constant
          index1 = 1
          index2 = 2
          four_nodes = .true.
          
       case ('1_0')
          
          nmax_1 = 1      ! start with 1 vertex in xi1 direction
          step_1 = 1      ! increase # of vertices in xi1 with each step in xi2
          nmax_2 = ndiv+1 ! ndiv+1 vertices in xi2 direction
          step_2 = 0      ! # of vertices in xi2 is constant
          index1 = 1
          index2 = 2
          
       case ('1_1')
          
          nmax_1 = ndiv+1 ! start with ndiv+1 vertices in xi1 direction
          step_1 = -1     ! decrease # of vertices in xi1 with each step in xi2
          nmax_2 = ndiv+1 ! ndiv+1 vertices in xi2 direction
          step_2 = 0      ! # of vertices in xi2 is constant
          index1 = 1
          index2 = 2
          
       case ('2_0')
          
          nmax_2 = ndiv+1 ! ndiv+1 vertices in xi1 direction
          step_2 = 0      ! # of vertices in xi1 is constant
          nmax_1 = 1      ! start with 1 vertex in xi2 direction
          step_1 = 1      ! increase # of vertices in xi2 with each step in xi1          
          index2 = 1
          index1 = 2
          
       case ('2_1')
          
          nmax_2 = ndiv+1 ! ndiv+1 vertices in xi1 direction
          step_2 = 0      ! # of vertices in xi1 is constant
          nmax_1 = ndiv+1 ! start with ndiv+1 vertices in xi2 direction
          step_1 = -1     ! decrease # of vertices in xi2 with each step in xi1
          index2 = 1
          index1 = 2
       end select
       
       xi(index2) = 0.0_dp
       do i = 1,nmax_2
          xi(index1) = 0.0_dp
          do j = 1,nmax_1
             num_vertices = num_vertices + 1
             X = coord_at_xi(ne,xi,'hermite')
             vertex_xyz(1:3,num_vertices) = X(1:3)
             if(nmax_1.gt.1) xi(index1) = xi(index1) + 1.0_dp/(nmax_1-1)
             if(i.gt.1.and.j.gt.1)then
                num_triangles = num_triangles + 1
                triangle(1,num_triangles) = num_vertices
                triangle(2,num_triangles) = num_vertices-1
                triangle(3,num_triangles) = nvertex_row+j-1
                if(four_nodes.or.(.not.four_nodes.and.j.lt.nmax_1).or.step_1.eq.-1)then
                   num_triangles = num_triangles + 1
                   triangle(1,num_triangles) = num_vertices
                   triangle(2,num_triangles) = nvertex_row+j
                   triangle(3,num_triangles) = nvertex_row+j-1
                endif
                if(step_1.eq.-1.and.j.eq.nmax_1)then
                   num_triangles = num_triangles + 1
                   triangle(1,num_triangles) = num_vertices
                   triangle(2,num_triangles) = nvertex_row+j+1
                   triangle(3,num_triangles) = nvertex_row+j
                endif
             else if(step_1.eq.-1.and.i.eq.nmax_2.and.j.eq.1)then
                num_triangles = num_triangles + 1
                triangle(1,num_triangles) = num_vertices
                triangle(2,num_triangles) = num_vertices-1
                triangle(3,num_triangles) = num_vertices-2
             endif
          enddo !j
          nvertex_row = num_vertices-nmax_1 !record first vertex # in previous row
          if(nmax_2.gt.1) xi(index2) = xi(index2) + 1.0_dp/(nmax_2-1)
          nmax_1 = nmax_1 + step_1
       enddo !i
    enddo
    
    write(*,'('' Made'',I8,'' triangles to cover'',I6,'' surface elements'')') &
         num_triangles,num_elems_2d
    
    call enter_exit(sub_name,2)
    
  end subroutine triangles_from_surface


!!!#############################################################################
  
  subroutine element_connectivity_2d
    !*element_connectivity_2d:*  Calculates element connectivity in 2D and
    ! stores in array elem_cnct_2d

    ! Local variables
    integer :: ne,ne2,nn,np,noelem2,np_list(4),np_list_2(4)
    integer,parameter :: num_elem_nodes = 4
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'element_connectivity_2d'
    call enter_exit(sub_name,1)
    
    if(.not.allocated(elem_cnct_2d)) allocate(elem_cnct_2d(-2:2,0:10,num_elems_2d))
    if(.not.allocated(elems_at_node_2d)) allocate(elems_at_node_2d(num_nodes_2d,0:10))
    
!!! calculate elems_at_node_2d array: stores the elements that nodes are in
    
    elems_at_node_2d = 0 !initialise all
    
    do ne = 1,num_elems_2d
       do nn = 1,num_elem_nodes
          np = elem_nodes_2d(nn,ne)
          elems_at_node_2d(np,0) = elems_at_node_2d(np,0)+1
          elems_at_node_2d(np,elems_at_node_2d(np,0)) = ne !element that np is in
       enddo !nn
    enddo !noelem
    
!!! calculate elem_cnct_2d array: stores the connectivity of all elements
    
    elem_cnct_2d = 0 !initialise all elem_cnct_2d
    
    do ne = 1,num_elems_2d ! for each of the 2d elements
       np_list(1:4) = elem_nodes_2d(1:4,ne) ! the list of nodes in the element (including repeated)
!!! check the elements attached to the 1st node
       do noelem2 = 1,elems_at_node_2d(np_list(1),0) ! for each element attached to the 1st node
          ne2 = elems_at_node_2d(np_list(1),noelem2) ! attached element number
          if(ne2.ne.ne)then
             np_list_2(1:4) = elem_nodes_2d(1:4,ne2) !list of nodes in attached element
             if(np_list(2).ne.np_list(1))then ! only if first two nodes are not repeated
                if(inlist(np_list(2),np_list_2))then
                   elem_cnct_2d(-2,0,ne) = elem_cnct_2d(-2,0,ne)+1
                   elem_cnct_2d(-2,elem_cnct_2d(-2,0,ne),ne) = ne2 
                endif
             endif
             if(np_list(3).ne.np_list(1))then ! only if the two nodes are not repeated
                if(inlist(np_list(3),np_list_2))then
                   elem_cnct_2d(-1,0,ne) = elem_cnct_2d(-1,0,ne)+1
                   elem_cnct_2d(-1,elem_cnct_2d(-1,0,ne),ne) = ne2 
                endif
             endif
          endif
       enddo
!!! check the elements attached to the 4th node
       do noelem2 = 1,elems_at_node_2d(np_list(4),0) ! for each element attached to the 4th node
          ne2 = elems_at_node_2d(np_list(4),noelem2) ! attached element number
          if(ne2.ne.ne)then
             np_list_2(1:4) = elem_nodes_2d(1:4,ne2) !list of nodes in attached element
             if(np_list(2).ne.np_list(4))then ! only if two nodes are not repeated
                if(inlist(np_list(2),np_list_2))then
                   elem_cnct_2d(1,0,ne) = elem_cnct_2d(1,0,ne)+1
                   elem_cnct_2d(1,elem_cnct_2d(1,0,ne),ne) = ne2 
                endif
             endif
             if(np_list(3).ne.np_list(4))then ! only if the two nodes are not repeated
                if(inlist(np_list(3),np_list_2))then
                   elem_cnct_2d(2,0,ne) = elem_cnct_2d(2,0,ne)+1
                   elem_cnct_2d(2,elem_cnct_2d(2,0,ne),ne) = ne2 
                endif
             endif
          endif
       enddo !noelem2
    enddo
    
    call enter_exit(sub_name,2)
    
  end subroutine element_connectivity_2d

!!!#############################################################################

  subroutine line_segments_for_2d_mesh(sf_option)
    !*line_segments_for_2d_mesh:* set up the line segment arrays for a 2d mesh
    
    character(len=4),intent(in) :: sf_option
    ! Local variables
    integer :: ne,ne_adjacent,ni1,nj,nl,nl_adj,npn(2)
    logical :: MAKE
    logical :: based_on_elems = .true.
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'line_segments_for_2d_mesh'
    call enter_exit(sub_name,1)
    
    if(.not.allocated(elem_lines_2d)) allocate(elem_lines_2d(4,num_elems_2d))
    if(.not.allocated(scale_factors_2d)) allocate(scale_factors_2d(16,num_elems_2d))
    
    elem_lines_2d=0
    num_lines_2d = 0
    
    if(based_on_elems)then
!!! estimate number of lines, for allocating memory to arrays
!!! before setting up arrays, count the number of lines required
       do ne=1,num_elems_2d
          MAKE=.FALSE.
          if(elem_cnct_2d(-1,0,ne) == 0) MAKE=.TRUE. !exterior, make line
          ne_adjacent=elem_cnct_2d(-1,1,ne)
          if(ne_adjacent > 0)then
             if(elem_lines_2d(4,ne_adjacent) == 0) MAKE=.TRUE.
          endif
          if(MAKE) num_lines_2d = num_lines_2d+1
          MAKE=.FALSE.
          if(elem_cnct_2d(-2,0,ne) == 0) MAKE=.TRUE. !exterior, make line
          ne_adjacent=elem_cnct_2d(-2,1,ne)
          if(ne_adjacent > 0)then
             if(elem_lines_2d(2,ne_adjacent) == 0) MAKE=.TRUE.
          endif
          if(MAKE) num_lines_2d=num_lines_2d+1
          num_lines_2d = num_lines_2d+2
          elem_lines_2d(2,ne) = 1 ! at this stage just to tag it for conditional above
          elem_lines_2d(4,ne) = 1 ! at this stage just to tag it for conditional above
       enddo !ne
       
       elem_lines_2d = 0
       
       if(.not.allocated(lines_2d)) allocate(lines_2d(0:num_lines_2d))
       if(.not.allocated(line_versn_2d)) allocate(line_versn_2d(2,3,num_lines_2d))
       if(.not.allocated(lines_in_elem)) allocate(lines_in_elem(0:4,num_lines_2d))
       if(.not.allocated(nodes_in_line)) allocate(nodes_in_line(3,0:3,num_lines_2d))
       if(.not.allocated(arclength)) allocate(arclength(3,num_lines_2d)) 
       lines_in_elem=0
       lines_2d=0
       nodes_in_line=0
       line_versn_2d=0
       num_lines_2d = 0 ! reset to zero for loop below
       
!!! Now run through the same as above, and set up the arrays
       do ne=1,num_elems_2d
          !check whether to make a line
          MAKE=.FALSE.
          if(elem_cnct_2d(-1,0,ne) == 0) MAKE=.TRUE. !exterior, make line
          ne_adjacent=elem_cnct_2d(-1,1,ne)
          if(ne_adjacent.gt.0)then
             if(elem_lines_2d(4,ne_adjacent) == 0) MAKE=.TRUE.
          endif
          if(MAKE)then
             num_lines_2d = num_lines_2d+1
             lines_2d(num_lines_2d) = num_lines_2d !record a new line number
             lines_in_elem(0,num_lines_2d) = lines_in_elem(0,num_lines_2d)+1
             lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d) = ne !line num_lines_2d is in element ne
             elem_lines_2d(3,ne) = num_lines_2d !num_lines_2d is global line # corresponding to local line 3 of ne
             npn(1) = 1
             npn(2) = 3
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
             nodes_in_line(1,0,num_lines_2d)=1 !Xi-direction of line segment num_lines_2d
             do nj=1,3
                nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
                do ni1=1,2
                   line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
                enddo !n
             enddo !nj
          else !get adjacent element line number
             !WARNING:: this only works if all Xi directions are consistent!!!!
             ne_adjacent = elem_cnct_2d(-2,1,ne)
             do nl_adj = 1,4
                nl = elem_lines_2d(nl_adj,ne_adjacent)
                if(nl /= 0)then
                   if(nodes_in_line(2,1,nl) == elem_nodes_2d(1,ne) .and. &
                        nodes_in_line(3,1,nl) == elem_nodes_2d(2,ne))then
                      elem_lines_2d(1,ne) = nl
                   elseif(nodes_in_line(2,1,nl) == elem_nodes_2d(2,ne) .and. &
                        nodes_in_line(3,1,nl) == elem_nodes_2d(1,ne))then
                      elem_lines_2d(1,ne) = nl
                   endif
                endif
             enddo
             !             elem_lines_2d(1,ne)=elem_lines_2d(2,ne_adjacent)
          endif
          
          !*! new:       
          MAKE=.TRUE.
          ne_adjacent=elem_cnct_2d(1,1,ne)
          if(ne_adjacent.gt.0)then
             if(elem_lines_2d(3,ne_adjacent) /= 0) MAKE=.FALSE.
          endif
          
          if(MAKE)then
             num_lines_2d=num_lines_2d+1
             lines_2d(num_lines_2d)=num_lines_2d !record a new line number
             lines_in_elem(0,num_lines_2d)=lines_in_elem(0,num_lines_2d)+1
             lines_in_elem(lines_in_elem(0,num_lines_2d),num_lines_2d)=ne !line num_lines_2d is in element ne
             elem_lines_2d(4,ne)=num_lines_2d !num_lines_2d is global line # corresponding to local line 4 of ne
             npn(1)=2
             npn(2)=4
             nodes_in_line(2,1,num_lines_2d)=elem_nodes_2d(2,ne) !records 1st node in line
             nodes_in_line(3,1,num_lines_2d)=elem_nodes_2d(4,ne) !records 2nd node in line
             nodes_in_line(1,0,num_lines_2d)=2 !Xi-direction of line segment num_lines_2d
             do nj=1,3
                nodes_in_line(1,nj,num_lines_2d)=4 !type of basis function (1 for linear,4 for cubicHermite)
                do ni1=1,2
                   line_versn_2d(ni1,nj,num_lines_2d)=elem_versn_2d(npn(ni1),ne)
                enddo !n
             enddo !nj
          else !get adjacent element line number
             !WARNING:: this only works if all Xi directions are consistent!!!!
             ne_adjacent=elem_cnct_2d(1,1,ne)
             elem_lines_2d(4,ne)=elem_lines_2d(3,ne_adjacent)
          endif
          
          MAKE=.TRUE.
          ne_adjacent=elem_cnct_2d(2,1,ne)
          if(ne_adjacent.gt.0)then
             if(elem_lines_2d(1,ne_adjacent) /= 0) MAKE=.FALSE.
          endif
          
          if(MAKE)then
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
          else !get adjacent element line number
             !WARNING:: this only works if all Xi directions are consistent!!!!
             ne_adjacent=elem_cnct_2d(2,1,ne)
             elem_lines_2d(2,ne)=elem_lines_2d(1,ne_adjacent)
          endif
       enddo !ne
    endif
    
    call calc_scale_factors_2d(sf_option)
    
    call enter_exit(sub_name,2)
    
  end subroutine line_segments_for_2d_mesh

!!!#############################################################################

  subroutine evaluate_ordering()
    !*evaluate_ordering:* calculates generations, Horsfield orders,
    ! Strahler orders for a given tree
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_ORDERING" :: EVALUATE_ORDERING

    ! Local Variables
    integer :: INLETS,ne,ne0,ne2,noelem2,np,np2,nn,num_attach,n_children, &
         n_generation,n_horsfield,OUTLETS,STRAHLER,STRAHLER_ADD,temp1
    LOGICAL :: DISCONNECT,DUPLICATE
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'evaluate_ordering'
    call enter_exit(sub_name,1)
    
!!! Calculate branch generations
    elem_ordrs = 0
    maxgen=1
    do ne=1,num_elems
       ne0=elem_cnct(-1,1,ne) !parent
       if(ne0.NE.0)THEN
          n_generation=elem_ordrs(1,ne0) !parent generation
          if(elem_cnct(1,0,ne0).EQ.1)THEN !single daughter
             elem_ordrs(1,ne)=n_generation + (elem_symmetry(ne)-1)
          else if(elem_cnct(1,0,ne0).GE.2)THEN
             elem_ordrs(1,ne)=n_generation+1
          endif
       else
          elem_ordrs(1,ne)=1 !generation 1
       endif
       maxgen=max(maxgen,elem_ordrs(1,ne))
    enddo !noelem

!!! Calculate the branch orders
    do ne=num_elems,1,-1
       n_horsfield=MAX(elem_ordrs(2,ne),1)
       n_children=elem_cnct(1,0,ne) !number of child branches
       if(n_children.EQ.1)THEN
          if(elem_ordrs(1,elem_cnct(1,1,ne)).EQ.0)  n_children=0
       endif
       STRAHLER=0
       STRAHLER_ADD=1
       if(n_children.GE.2)THEN !branch has two or more daughters
          STRAHLER=elem_ordrs(3,elem_cnct(1,1,ne)) !first daughter
          do noelem2=1,n_children !for all daughters
             ne2=elem_cnct(1,noelem2,ne) !global element # of daughter
             temp1=elem_ordrs(2,ne2) !Horsfield order of daughter
             if(temp1.GT.n_horsfield)then
                n_horsfield=temp1
             endif
             if(elem_ordrs(3,ne2).LT.STRAHLER)THEN
                STRAHLER_ADD=0
             else if(elem_ordrs(3,ne2).GT.STRAHLER)THEN
                STRAHLER_ADD=0
                STRAHLER=elem_ordrs(3,ne2) !highest daughter
             endif
          enddo !noelem2 (ne2)
          n_horsfield=n_horsfield+1 !Horsfield ordering
       else if(n_children.EQ.1)THEN
          ne2=elem_cnct(1,1,ne) !local element # of daughter
          n_horsfield=elem_ordrs(2,ne2)+(elem_symmetry(ne)-1)
          STRAHLER_ADD=elem_ordrs(3,ne2)+(elem_symmetry(ne)-1)
       endif !elem_cnct
       elem_ordrs(2,ne)=n_horsfield !store the Horsfield order
       elem_ordrs(3,ne)=STRAHLER+STRAHLER_ADD !Strahler order
    enddo !noelem
    
!!! Check for disconnected nodes and number of inlets and outlets
    DUPLICATE=.FALSE.
    do ne=1,num_elems
       np=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       if(np.EQ.np2)THEN
          DUPLICATE=.TRUE.
       endif
    enddo
    
    DISCONNECT=.FALSE.
    INLETS=0
    OUTLETS=0
    do np=1,num_nodes
       num_attach=elems_at_node(np,0)
       if(num_attach.EQ.0)THEN
          DISCONNECT=.TRUE.
       elseif(num_attach.EQ.1)THEN
          ne=elems_at_node(np,1)
          if(elem_cnct(1,0,ne).EQ.0) OUTLETS=OUTLETS+1
          if(elem_cnct(-1,0,ne).EQ.0) INLETS=INLETS+1
       elseif(num_attach.GT.3)THEN
          WRITE(*,*) ' Node ',np,' attached to',num_attach,' elements'
       endif
    enddo
    
    call enter_exit(sub_name,2)
    
  end subroutine evaluate_ordering

!!!#############################################################################

  subroutine set_initial_volume(Gdirn,COV,total_volume,Rmax,Rmin)
    !*set_initial_volume:* assigns a volume to terminal units appended on a
    ! tree structure based on an assumption of a linear gradient in the
    ! gravitational direction with max, min, and COV values defined.
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SET_INITIAL_VOLUME" :: SET_INITIAL_VOLUME
    
    integer,intent(in) :: Gdirn
    real(dp),intent(in) :: COV,total_volume,Rmax,Rmin
    !     Local parameters
    integer :: ne,np2,nunit
    real(dp) ::  factor_adjust,max_z,min_z,random_number,range_z,&
         volume_estimate,volume_of_tree,Vmax,Vmin,Xi
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'set_initial_volume'
    call enter_exit(sub_name,1)
    
    volume_estimate = 1.0_dp
    volume_of_tree = 0.0_dp
    
    call volume_of_mesh(volume_estimate,volume_of_tree)
    
    random_number=-1.1_dp
    
    Vmax = Rmax * (total_volume-volume_estimate)/elem_units_below(1)
    Vmin = Rmin * (total_volume-volume_estimate)/elem_units_below(1)
    
!!! for each elastic unit find the maximum and minimum coordinates in the Gdirn direction
    max_z=-1.0e+6_dp
    min_z=1.0e+6_dp
    do nunit=1,num_units
       ne=units(nunit)
       np2=elem_nodes(2,ne)
       max_z=MAX(max_z,node_xyz(Gdirn,np2))
       min_z=MIN(min_z,node_xyz(Gdirn,np2))
    enddo !nunit
    
    range_z=abs(max_z-min_z)
    if(abs(range_z).le.1.0e-5_dp) range_z=1.0_dp
    
!!! for each elastic unit allocate a size based on a gradient in the Gdirn direction, and
!!! perturb by a user-defined COV. This should be calling a random number generator.
    do nunit=1,num_units
       ne=units(nunit)
       np2=elem_nodes(2,ne) !end node
       Xi=(node_xyz(Gdirn,np2)-min_z)/range_z
       random_number=random_number+0.1_dp
       if(random_number.GT.1.0_dp) random_number=-1.1_dp
       unit_field(nu_vol,nunit)=(Vmax*Xi+Vmin*(1.0_dp-Xi))*(1.0_dp+COV*random_number)
       unit_field(nu_vt,nunit)=0.0_dp !initialise the tidal volume to a unit
    enddo !nunit
    
    ! correct unit volumes such that total volume is exactly as specified
    call volume_of_mesh(volume_estimate,volume_of_tree)
    factor_adjust = (total_volume-volume_of_tree)/(volume_estimate-volume_of_tree)
    do nunit=1,num_units
       unit_field(nu_vol,nunit) = unit_field(nu_vol,nunit)*factor_adjust
    enddo
    
    write(*,'('' Number of elements is '',I5)') num_elems
    write(*,'('' Initial volume is '',F6.2,'' L'')') total_volume/1.0e+6_dp
    write(*,'('' Deadspace volume is '',F6.1,'' mL'')') volume_of_tree/1.0e+3_dp
    
    call enter_exit(sub_name,2)
    
  end subroutine set_initial_volume

!!!#############################################################################

  subroutine volume_of_mesh(volume_model,volume_tree)
    !*volume_of_mesh:* calculates the volume of an airway mesh including
    ! conducting and respiratory airways
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_VOLUME_OF_MESH" :: VOLUME_OF_MESH
    
    real(dp) :: volume_model,volume_tree
    !     Local Variables
    integer :: ne,ne0,nunit
    real(dp),allocatable :: vol_anat(:),vol_below(:)
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'volume_of_mesh'
    call enter_exit(sub_name,1)
    
    if(.not.allocated(vol_anat)) allocate(vol_anat(num_elems))
    if(.not.allocated(vol_below)) allocate(vol_below(num_elems))
    
    vol_anat = elem_field(ne_vol,1:num_elems) !initialise to branch volume
    vol_below = elem_field(ne_vol,1:num_elems) !initialise to branch volume

    do nunit = 1,num_units
       ne = units(nunit)
       if(ne.ne.0) vol_below(ne) = vol_below(ne) + unit_field(nu_vol,nunit) !add elastic unit volume
    enddo !nunit

    do ne = num_elems,2,-1
       ne0 = elem_cnct(-1,1,ne)
!!! don't include the respiratory airways (already included via lumped units). multiply
!!! by the element type (0 for respiratory) to account for this
       vol_anat(ne0) = vol_anat(ne0) + dble(elem_symmetry(ne))*dble(elem_ordrs(no_type,ne))*vol_anat(ne)
       vol_below(ne0) = vol_below(ne0) + dble(elem_symmetry(ne))*dble(elem_ordrs(no_type,ne))*vol_below(ne)
    enddo !noelem

    elem_field(ne_vd_bel,:) = vol_anat(:)
    elem_field(ne_vol_bel,:) = vol_below(:)
    volume_model = elem_field(ne_vol_bel,1)
    volume_tree = elem_field(ne_vd_bel,1)

    deallocate(vol_anat)
    deallocate(vol_below)
    
    call enter_exit(sub_name,2)
    
  end subroutine volume_of_mesh

!!!#############################################################################

  subroutine write_geo_file(type, filename)
    !*write_geo_file:* converts a surface mesh (created using make_2d_vessel_from_1d)
    ! into a gmsh formatted mesh and writes to file. 
    ! options on 'type': 1== single layered surface mesh of the vessel wall
    !                    2== double-layered thick-walled volume mesh of vessel wall
    !                    3== volume mesh of vessel lumen
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_WRITE_GEO_FILE" :: WRITE_GEO_FILE

    integer,intent(in) :: type
    character(len=*),intent(in) :: filename
    !     Local parameters
    integer :: j, ncount_loop = 0, ncount_point = 0, ncount_spline = 0, &
         nl_offset,np,np_offset
    integer,parameter :: ifile = 10
    integer,allocatable :: element_spline(:,:),elem_surfaces(:,:)
    real(dp),parameter :: lc0 = 1.0_dp, lc1 = 1.0_dp
    real(dp),allocatable :: node_xyz_offset(:,:)
    character(len=200) :: opfile
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'write_geo_file'
    call enter_exit(sub_name,1)

    opfile = trim(filename)//'.geo'
    open(10, file=opfile, status='replace')
      
    write(ifile,'(''/***********************'')')
    write(ifile,'(''*'')')
    write(ifile,'(''* Conversion of LungSim to GMSH'')')
    write(ifile,'(''*'')')
    write(ifile,'(''***********************/'')')
    
    write(ifile,'(/''lc ='',f8.4,'';'')') lc0
    write(ifile,'(/''sc ='',f8.4,'';'')') lc1
    write(ifile,'(/)')

    allocate(element_spline(4,num_elems_2d*2))
    allocate(elem_surfaces(5,num_elems_2d))
    element_spline = 0
    elem_surfaces = 0
    ncount_spline = 0 
    np_offset = 0

    if(type.eq.1)then
!!! write out a surface mesh that describes a structured vessel surface
       call write_surface_geo(element_spline,elem_surfaces,ifile,ncount_point, &
            ncount_loop,ncount_spline,np_offset)
       
    else if(type.eq.2)then
!!! write out a volume that encloses a thick-walled vessel tree. Make a gmsh .geo file
!!! for the surface of the tree, then copy, scale, and translate to get an 'outer shell'.
!!! Join the inner and outer shells at the entry and exits.

       allocate(node_xyz_offset(3,num_nodes_2d))
       node_xyz_offset = 0.0_dp
       call write_surface_geo(element_spline,elem_surfaces,ifile,ncount_point, &
            ncount_loop,ncount_spline,np_offset)
       call geo_node_offset(node_xyz_offset)

       do np = 1,num_nodes_2d
          forall (j = 1:3) node_xyz_2d(1,1,j,np) = node_xyz_2d(1,1,j,np) &
               + node_xyz_offset(j,np)
       enddo
       np_offset = ncount_point
       nl_offset = ncount_spline
       call write_surface_geo(element_spline,elem_surfaces,ifile,ncount_point, &
            ncount_loop,ncount_spline,np_offset)
       do np = 1,num_nodes_2d
          forall (j = 1:3) node_xyz_2d(1,1,j,np) = node_xyz_2d(1,1,j,np) &
               - node_xyz_offset(j,np)
       enddo
       ! cap the entry and exits
       call geo_entry_exit_cap(element_spline,ifile,ncount_loop, &
            ncount_spline,np_offset,nl_offset)
       deallocate(node_xyz_offset)

    else if(type.eq.3)then
!!! write out a volume mesh for the vessel lumen, where the vessel surface mesh is the
!!! exterior. Make a .gmsh file that includes the vessel surfaces, surfaces that join to a vessel
!!! centreline, and surfaces that 'cap' each vessel segment entry and exit.

       call write_surface_geo(element_spline,elem_surfaces,ifile,ncount_point, &
            ncount_loop,ncount_spline,np_offset)
       call write_3d_geo(element_spline,elem_surfaces,ifile,ncount_point, &
            ncount_loop,ncount_spline)
    endif

    deallocate(element_spline)
    deallocate(elem_surfaces)
    close(ifile)

    call enter_exit(sub_name,2)
    
  end subroutine write_geo_file
  
!!!#############################################################################
  
  function get_final_real(string)
    !*get_final_real:* gets the last real number on a string

    character,intent(in) :: string*(*)
    ! Local parameters
    integer :: ibeg,iend
    real(dp) :: rsign,rtemp,get_final_real
    character :: sub_string*(40)
    
    ! --------------------------------------------------------------------------
    
    iend=len(string) !get the length of the string
    ibeg=index(string,":")+1 !get location of real in string
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond :
    iend=len(sub_string) !get the length of the sub-string
    if(sub_string(1:1).eq.'-')then !check whether negative
       rsign=-1.0_dp
       ibeg=2
    else
       rsign=1.0_dp
       ibeg=1
    endif
   
    read (sub_string(ibeg:iend), * ) rtemp !get real value
    rtemp=rtemp*rsign !apply sign to number
    
    get_final_real=rtemp !return the real value
    
  end function get_final_real

!!!#############################################################################

  subroutine get_final_string(string,rtemp)
    !*get_final_string:* gets the last set of characters surrounded
    !  by whitespace in a string

    character, intent(in) :: string*(*)
    ! Local parameters
    integer :: ibeg,iend
    real(dp), intent(out) :: rtemp
    real(dp) :: rsign
    character :: sub_string*(40)
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'get_final_string'
    call enter_exit(sub_name,1)
    
    iend=len(string) !get the length of the string
    ibeg=index(string,":")+1 !get location of real in string
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond :
    iend=len(sub_string) !get the length of the sub-string
    if(sub_string(1:1).eq.'-')then !check whether negative
       rsign=-1.0_dp
       ibeg=2
    else
       rsign=1.0_dp
       ibeg=1
    endif
    read (sub_string(ibeg:iend), '(D25.17)' ) rtemp !get real value
    rtemp=rtemp*rsign !apply sign to number
    
    call enter_exit(sub_name,2)

  end subroutine get_final_string

!!!#############################################################################

  subroutine get_local_node(np_global,np_local)
    !*get_local_node:* gets the local node number for a given global node

    integer,intent(in) :: np_global
    integer,intent(out) :: np_local
    ! Local parameters
    integer :: np
    logical :: found
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'get_local_node'
    call enter_exit(sub_name,1)
    
    np=1
    found=.false.
    do while (.not.found)
       if(nodes(np).eq.np_global)then
          found=.true.
       elseif(np.gt.num_nodes)then
          found = .true.
          write(*,'('' Global node '',I6,'' not in node list'')') np_global
          read(*,*)
       else
          np=np+1
       endif
    enddo
    
    np_local = np

    call enter_exit(sub_name,2)

  end subroutine get_local_node
  
!!!#############################################################################

  subroutine geo_entry_exit_cap(element_spline,ifile,ncount_loop, &
       ncount_spline,np_offset,nl_offset)

    integer,intent(in) :: element_spline(:,:),ifile,np_offset,nl_offset
    integer :: ncount_loop,ncount_spline
    ! Local variables
    integer :: k,line1,line2,line3,line4,ne,np1,np2
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'geo_entry_exit_cap'
    call enter_exit(sub_name,1)
        
    ne = 1
    do while (ne.le.num_elems_2d)
       
       if(elem_cnct_2d(-2,0,ne).eq.0)then
          do k = 0,3
             np1 = elem_nodes_2d(1,ne+k)
             np2 = np1 + np_offset
             ncount_spline = ncount_spline + 1
             write(10,'(''Line('',I8,'') = {'',I8,'','',I8,''};'')') &
                  ncount_spline,np1,np2
          enddo
          do k = 0,3
             line1 = element_spline(1,ne+k) - nl_offset
             line3 = -element_spline(1,ne+k)
             if(k.lt.3)then
                line2 = ncount_spline + k - 2
                line4 = -(line2 - 1)
             else
                line2 = ncount_spline - 3 ! first new line
                line4 = -ncount_spline
             endif
             ncount_loop = ncount_loop + 1
             write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8, &
                  &'','',i8,'','',i8,''};'')') &
                  ncount_loop, line1, line2, line3, line4
             ncount_loop = ncount_loop + 1
             write(ifile,'(''Surface('',I8,'') = {'',I8,''};'')') &
                  ncount_loop, ncount_loop - 1
          enddo
       endif
       
       if(elem_cnct_2d(2,0,ne).eq.0)then
          do k = 0,3
             np1 = elem_nodes_2d(3,ne+k)
             np2 = np1 + np_offset
             ncount_spline = ncount_spline + 1
             write(10,'(''Line('',I8,'') = {'',I8,'','',I8,''};'')') &
                  ncount_spline,np1,np2
          enddo
          do k = 0,3
             line1 = -element_spline(3,ne+k) - nl_offset
             line3 = -element_spline(3,ne+k)
             if(k.lt.3)then
                line2 = ncount_spline + k - 2
                line4 = -(line2 - 1)
             else
                line2 = ncount_spline - 3 ! first new line
                line4 = -ncount_spline
             endif
             ncount_loop = ncount_loop + 1
             write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','' &
                  &,i8,'','',i8,''};'')') &
                  ncount_loop, line1, line2, line3, line4
             ncount_loop = ncount_loop + 1
             write(ifile,'(''Surface('',I8,'') = {'',I8,''};'')') &
                  ncount_loop, ncount_loop - 1
          enddo
       endif
       
       ne = ne + 4
    enddo

    call enter_exit(sub_name,2)
    
  end subroutine geo_entry_exit_cap

!!!#############################################################################

  subroutine geo_node_offset(node_xyz_offset)
    
    real(dp) :: node_xyz_offset(:,:)
    ! Local variables
    integer:: j,k,ne,np1
    real(dp) :: point_temp(3),point_xyz_centre(3),wall_thickness
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'geo_node_offset'
    call enter_exit(sub_name,1)
        
    ne = 1
    do while (ne.le.num_elems_2d)
       
!!!....nodes at model entry
       if(elem_cnct_2d(-2,0,ne).eq.0)then
          point_xyz_centre(:) = 0.0_dp
          do k = 0,3
             np1 = elem_nodes_2d(1,ne+k)
             forall (j = 1:3) point_xyz_centre(j) = point_xyz_centre(j) + &
                  0.25_dp * node_xyz_2d(1,1,j,np1)
          enddo
          point_temp(1:3) = node_xyz_2d(1,1,1:3,elem_nodes_2d(1,ne)) ! location of first ring node
          wall_thickness = &
               0.2_dp * distance_between_points(point_xyz_centre,point_temp)
          do k = 0,3
             np1 = elem_nodes_2d(1,ne+k)
             point_temp(1:3) = node_xyz_2d(1,1,1:3,np1)
             node_xyz_offset(:,np1) = wall_thickness * &
                  direction_point_to_point(point_xyz_centre,point_temp) 
          enddo  ! k
       endif  ! elem_cnct
       
!!!....nodes at Xi2=1 ends of 'rings'
       point_xyz_centre(:) = 0.0_dp
       do k = 0,3
          np1 = elem_nodes_2d(3,ne+k)
          forall (j = 1:3) point_xyz_centre(j) = point_xyz_centre(j) &
               + 0.25_dp * node_xyz_2d(1,1,j,np1)
       enddo
       point_temp(1:3) = node_xyz_2d(1,1,1:3,elem_nodes_2d(3,ne)) ! location of first ring node
       wall_thickness = 0.2_dp * &
            distance_between_points(point_xyz_centre,point_temp)
       do k = 0,3
          np1 = elem_nodes_2d(3,ne+k)
          point_temp(1:3) = node_xyz_2d(1,1,1:3,np1)
          node_xyz_offset(:,np1) = wall_thickness * &
               direction_point_to_point(point_xyz_centre,point_temp)
       enddo  ! k

!!!....check for crux node
       if(node_versn_2d(elem_nodes_2d(3,ne)).eq.6.or. &
            node_versn_2d(elem_nodes_2d(4,ne)).eq.6)then
          np1 = elem_nodes_2d(3,ne+3) + 1   ! number of crux node
          point_temp(1:3) = node_xyz_2d(1,1,1:3,np1)
          node_xyz_offset(:,np1) =  wall_thickness * &
               direction_point_to_point(point_xyz_centre,point_temp)
       endif

       ne = ne + 4
       
    enddo ! while (ne.le.num_elems_2d)
    
    call enter_exit(sub_name,2)
    
  end subroutine geo_node_offset
  
!!!#############################################################################

  subroutine group_elem_parent_term(ne_parent)
    !*group_elem_parent_term:* group the terminal elements that sit distal to
    !  a given parent element (ne_parent)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GROUP_ELEM_PARENT_TERM" :: GROUP_ELEM_PARENT_TERM
    
    integer,intent(in) :: ne_parent  ! the parent element number
    ! Local Variables
    integer :: ne,ne_count,noelem,num_parents
    integer,allocatable :: templist(:)
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'group_elem_parent_term'
    call enter_exit(sub_name,1)
    
    allocate(templist(num_elems))
    if(.not.allocated(parentlist)) allocate(parentlist(num_elems))
    
    !reset the list of parent elements to zero
    call group_elem_by_parent(ne_parent,templist)
    
    ne_count=count(templist.ne.0)
    num_parents=0
    parentlist=0
    
    do noelem=1,ne_count
       ne = templist(noelem)
       if(elem_cnct(1,0,ne).eq.0)then
          num_parents=num_parents+1
          parentlist(num_parents)=ne
       endif !elem_cnct
    enddo !noelem
    
    deallocate(templist)
    
    call enter_exit(sub_name,2)
    
  end subroutine group_elem_parent_term
  
!!!#############################################################################

  subroutine group_elem_by_parent(ne_parent,elemlist)
    !*group_elem_by_parent:* group elements that sit distal to a given
    ! parent element (ne_parent)

    integer,intent(in) :: ne_parent  ! the parent element number
    integer :: elemlist(:)
    ! Local Variables
    integer :: nt_bns,ne_count,num_nodes,m,n,ne0
    integer,allocatable :: ne_old(:),ne_temp(:)
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name='group_elem_by_parent'
    call enter_exit(sub_name,1)
    
    elemlist = 0
    allocate(ne_old(size(elemlist)))
    allocate(ne_temp(size(elemlist)))
    
    nt_bns=1
    ne_old(1) = ne_parent
    ne_count = 1
    elemlist(ne_count)=ne_parent
    
    do while(nt_bns.ne.0)
       num_nodes=nt_bns
       nt_bns=0
       do m=1,num_nodes
          ne0=ne_old(m) !parent global element number
          do n=1,elem_cnct(1,0,ne0) !for each daughter branch
             nt_bns=nt_bns+1
             ne_temp(nt_bns)=elem_cnct(1,n,ne0)
          enddo !n
       enddo !m
       do n=1,nt_bns
          ne_old(n)=ne_temp(n) !updates list of previous generation element numbers
          ne_count=ne_count+1
          elemlist(ne_count)=ne_temp(n)
       enddo !n
    enddo !while
    
    deallocate(ne_old)
    deallocate(ne_temp)
    
    call enter_exit(sub_name,2)
    
  end subroutine group_elem_by_parent

!!!#############################################################################
  
  subroutine reallocate_node_elem_arrays(num_elems_new,num_nodes_new)
    !*reallocate_node_elem_arrays:* Reallocates the size of geometric
    ! arrays when modifying geometries

    integer,intent(in) :: num_elems_new,num_nodes_new
    ! Local variables
    integer,allocatable :: nodelem_temp(:),enodes_temp(:,:),enodes_temp2(:,:,:)
    real(dp),allocatable :: xyz_temp(:,:),rnodes_temp(:,:)
    logical,allocatable :: exp_temp(:)
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'reallocate_node_elem_arrays'
    call enter_exit(sub_name,1)
    
    allocate(nodelem_temp(num_nodes))
    nodelem_temp = nodes ! copy to temporary array
    deallocate(nodes) !deallocate initially allocated memory
    allocate(nodes(num_nodes_new))
    nodes(1:num_nodes)=nodelem_temp(1:num_nodes)
    deallocate(nodelem_temp) !deallocate the temporary array
    
    allocate(xyz_temp(3,num_nodes))
    xyz_temp=node_xyz
    deallocate(node_xyz)
    allocate(node_xyz(3,num_nodes_new))
    node_xyz(1:3,1:num_nodes)=xyz_temp(1:3,1:num_nodes)
    
    allocate(nodelem_temp(num_elems))
    nodelem_temp = elems ! copy to temporary array
    deallocate(elems) !deallocate initially allocated memory
    allocate(elems(num_elems_new))
    elems(1:num_elems)=nodelem_temp(1:num_elems)
    deallocate(nodelem_temp) !deallocate the temporary array
    
    allocate(enodes_temp(2,num_elems))
    enodes_temp=elem_nodes
    deallocate(elem_nodes)
    allocate(elem_nodes(2,num_elems_new))
    elem_nodes(1:2,1:num_elems)=enodes_temp(1:2,1:num_elems)
    deallocate(enodes_temp)
    
    if(allocated(elem_field).and.num_ne.gt.0)then
       allocate(rnodes_temp(num_ne,num_elems))
       rnodes_temp=elem_field
       deallocate(elem_field)
       allocate(elem_field(num_ne,num_elems_new))
       elem_field(1:num_ne,1:num_elems)=rnodes_temp(1:num_ne,1:num_elems)
       deallocate(rnodes_temp)
       elem_field(1:num_ne,num_elems+1:num_elems_new) = 0.0_dp
    endif
    
    allocate(rnodes_temp(3,num_elems))
    rnodes_temp=elem_direction
    deallocate(elem_direction)
    allocate(elem_direction(3,num_elems_new))
    elem_direction(1:3,1:num_elems)=rnodes_temp(1:3,1:num_elems)
    deallocate(rnodes_temp)
    elem_direction(1:3,num_elems+1:num_elems_new) = 0.0_dp
    
    if(allocated(node_field).and.num_nj.gt.0)then
       allocate(rnodes_temp(num_nj,num_nodes))
       rnodes_temp=node_field
       deallocate(node_field)
       allocate(node_field(num_nj,num_nodes_new))
       node_field(1:num_nj,1:num_nodes)=rnodes_temp(1:num_nj,1:num_nodes)
       deallocate(rnodes_temp)
       node_field(1:num_nj,num_nodes+1:num_nodes_new)=0.0_dp
    endif
    
    allocate(nodelem_temp(num_elems))
    nodelem_temp = elem_symmetry ! copy to temporary array
    deallocate(elem_symmetry) !deallocate initially allocated memory
    allocate(elem_symmetry(num_elems_new))
    elem_symmetry(1:num_elems)=nodelem_temp(1:num_elems)
    deallocate(nodelem_temp) !deallocate the temporary array
    elem_symmetry(num_elems+1:num_elems_new)=1
    
    allocate(enodes_temp2(-1:1,0:2,0:num_elems))
    enodes_temp2=elem_cnct
    deallocate(elem_cnct)
    allocate(elem_cnct(-1:1,0:2,0:num_elems_new))
    elem_cnct(-1:1,0:2,0:num_elems)=enodes_temp2(-1:1,0:2,0:num_elems)
    deallocate(enodes_temp2)
    elem_cnct(-1:1,0:2,num_elems+1:num_elems_new) = 0
    
    allocate(enodes_temp(num_ord,num_elems))
    enodes_temp=elem_ordrs
    deallocate(elem_ordrs)
    allocate(elem_ordrs(num_ord,num_elems_new))
    elem_ordrs(1:num_ord,1:num_elems)=enodes_temp(1:num_ord,1:num_elems)
    deallocate(enodes_temp)
    elem_ordrs(1:num_ord,num_elems+1:num_elems_new) = 0
    
    if(allocated(elem_units_below).and.num_nu.gt.0)then
       allocate(nodelem_temp(num_elems))
       nodelem_temp=elem_units_below
       deallocate(elem_units_below)
       allocate(elem_units_below(num_elems_new))
       elem_units_below(1:num_elems)=nodelem_temp(1:num_elems)
       deallocate(nodelem_temp)
       elem_units_below(num_elems+1:num_elems_new)=0
    endif
    
    allocate(enodes_temp(num_nodes,0:3))
    enodes_temp=elems_at_node
    deallocate(elems_at_node)
    allocate(elems_at_node(num_nodes_new,0:3))
    elems_at_node(1:num_nodes,0:3)=enodes_temp(1:num_nodes,0:3)
    deallocate(enodes_temp)
    elems_at_node(num_nodes+1:num_nodes_new,0:3)=0
    
    if(model_type.eq.'gas_mix')then
       allocate(exp_temp(num_elems))
       exp_temp = expansile
       deallocate(expansile)
       allocate(expansile(num_elems_new))
       expansile(1:num_elems)=exp_temp(1:num_elems)
       deallocate(exp_temp)
       expansile(num_elems+1:num_elems_new)=.false.
    endif
    
    call enter_exit(sub_name,2)
    
  end subroutine reallocate_node_elem_arrays
  
!!!#############################################################################
  
  function get_local_node_f(ndimension,np_global) result(get_local_node)
    
    integer,intent(in) :: ndimension,np_global
    ! Local variables
    integer :: np
    integer :: get_local_node
    logical :: found

    ! --------------------------------------------------------------------------

    found = .false.
    np = 0
    
    select case (ndimension)
       
    case(1)
       do while (.not.found)
          np=np+1
          if(np.gt.num_nodes) then
             found = .true.
             write(*,'('' Warning: local node not found for global node'',I6)') np_global
          endif
          if(np_global.eq.nodes(np))then
             get_local_node = np
             found = .true.
          endif
       enddo
       
    case(2)
       do while (.not.found)
          np=np+1
          if(np.gt.num_nodes_2d) then
             found = .true.
             write(*,'('' Warning: local node not found for global node'',I6)') np_global
             read(*,*)
          endif
          if(np_global.eq.nodes_2d(np))then
             get_local_node = np
             found = .true.
          endif
       enddo
       
    end select
    
  end function get_local_node_f
  
!!!#############################################################################
  
  function get_final_integer(string)
    !*get_final_integer*
    
    character,intent(in) :: string*(*)
    ! Local parameters
    integer :: ibeg,iend,ierror,nsign,ntemp
    character :: sub_string*(40)
    integer :: get_final_integer
    
    ! --------------------------------------------------------------------------
    
    iend=len(string) !get the length of the string
    ibeg=index(string,":")+1 !get location of integer in string, follows ":"
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond ":"
    iend=len(sub_string) !length of the sub-string
    if(sub_string(1:1).eq.'-')then !check for negative sign
       nsign=-1
       ibeg=2
    else
       nsign=1
       ibeg=1
    endif
    
    read (sub_string(ibeg:iend), '(i10)', iostat=ierror ) ntemp !get integer values
    if(ierror.gt.0)then
       !... something wrong with data
       write(*,'(''Data read error'')')
       write(*,'(a)') sub_string(ibeg:iend)
    endif
    ntemp=ntemp*nsign !apply sign to number
    
    get_final_integer=ntemp !return the integer value
    
  end function get_final_integer
  
!!!#############################################################################
  
  subroutine get_four_nodes(ne,string)

    integer, intent(in) :: ne
    character(len=132), intent(in) :: string
    ! Local variables
    integer :: ibeg,iend,i_ss_end,nn,np_global
    character(len=40) :: sub_string
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'get_four_nodes'
    call enter_exit(sub_name,1)

    iend=len(string)
    ibeg=index(string,":")+1 !get location of first integer in string
    sub_string = adjustl(string(ibeg:iend)) ! get the characters beyond : and remove the leading blanks
    i_ss_end=len(sub_string) !get the end location of the sub-string
    
    ibeg=1
    do nn=1,4
       iend=index(sub_string," ") !get location of first blank in sub-string
       read (sub_string(ibeg:iend-1), '(i7)' ) np_global
       sub_string = adjustl(sub_string(iend:i_ss_end)) ! get chars beyond blank, remove leading blanks
       elem_nodes_2d(nn,ne)=get_local_node_f(2,np_global)
    enddo ! nn

    call enter_exit(sub_name,2)
    
  end subroutine get_four_nodes
  
!!!#############################################################################

  subroutine redistribute_mesh_nodes_2d_from_1d

    integer :: i,j,k,ne,nelist(20),ne_adjacent,np,nplist(20),np_adjacent,np_last,num_list, &
         ring1_nodes(4)
    real(dp) :: displace_length,distance_to_crux,distance_to_crux_last,line_length, &
         nedirection(3,20),point1(3),point2(3),point3(3),point4(3),vector(3)
    logical :: continue
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------
    
    sub_name = 'redistribute_mesh_nodes_2d_from_1d'
    call enter_exit(sub_name,1)

    nplist = 0
    nelist = 0

    np = 1
    do while (np.le.num_nodes_2d) ! run through all of the 2d mesh nodes
       ! use the 'front' bifurcation node to identify the crux node (based on the
       ! template structure that we used to set up the bifurcations). Crux node
       ! must be at front_node + 3.
       if(node_versn_2d(np).eq.6)then   ! this is the 'front' node at a bifurcation
          np = np + 3   ! this is the node number for the 'crux' node
          do i = 1,2 ! for each of the two branches that the bifurcation leads to
             ne = elems_at_node_2d(np,2*i-1)   ! get the first (and then third) element that node np (crux) is in
             num_list = 1                      ! count the number of 'rings' of nodes between bifurcations
             nelist(num_list) = ne             ! record the first (and third) element number ne
             np_adjacent = elem_nodes_2d(4,ne) ! node number in the +Xi2 direction (along the branch)
             nplist(num_list) = np_adjacent    ! the node number along the line from one bifurcation to the next
             ! get coordinates for three of the points on the first 'ring' that is in the direction of the
             ! branch. Not straightforward when at a bifurcation.
             point1(1:3) = node_xyz_2d(1,1,1:3,np_adjacent)
             point2(1:3) = node_xyz_2d(1,1,1:3,elem_nodes_2d(3,ne))
             point3(1:3) = node_xyz_2d(1,1,1:3,elem_nodes_2d(3,elem_cnct_2d(-1,1,ne)))
             point4(1:3) = node_xyz_2d(1,1,1:3,np)   ! the location of the crux node
             ! calculate the distance from the crux (point4) to the plane of first adjacent ring
             line_length = distance_from_plane_to_point(point1,point2,point3,point4)
             
             point1(1:3) = node_xyz_2d(1,1,1:3,np-1)   ! the location of the 'back' node of bifurcation
!!             ! calculate the line length from back node to a point on the first ring
!!             distance_to_crux_last = distance_between_points(point1,point2)
             
             continue = .true.
             do while(continue)   ! keep going: note that bifurcation will have > 1 version at nodes
                if(elem_cnct_2d(2,0,ne).eq.0)then ! no adjacent 2d elements in Xi+2 direction
                   continue = .false.
                else
                   ne = elem_cnct_2d(2,1,ne)        ! get the next adjacent element in Xi+2 direction
                   num_list = num_list + 1 ! the number of 'rings' between bifurcations
                   nelist(num_list) = ne ! the element number of the adjacent element
                   np_last = np_adjacent   ! store the previous node number
                   np_adjacent = elem_nodes_2d(4,ne) ! the next node on the line
                   nplist(num_list) = np_adjacent    ! the node number on the line from one bifurcation to the next
                   ! calculate the distance between adjacent rings
                   point1(1:3) = node_xyz_2d(1,1,1:3,np_last)      ! coordinates of node on previous ring
                   point2(1:3) = node_xyz_2d(1,1,1:3,np_adjacent)  ! coordinate of node on current ring
                   line_length = line_length + &
                        distance_between_points(point1,point2)  ! sum distance between rings
                   ! calculate the direction between connected nodes on rings
                   vector(1:3) = node_xyz_2d(1,1,1:3,np_adjacent) - &
                        node_xyz_2d(1,1,1:3,np_last)
                   vector = unit_vector(vector)
                   nedirection(1:3,num_list-1) = vector(1:3)  ! store the direction
                   ! continue until the next bifurcation is detected (nodes with > 1 version)                   
                   if(node_versn_2d(np_adjacent).ne.1) continue = .false.
                endif
             enddo

             line_length = line_length/real(num_list) ! this is the length to redistribute rings to
             
!!!          adjust the location of the nodes in each 'ring'
             do j = 1,num_list - 1   ! only adjust the rings that are between bifns: last 'ring' is actually the next bifn
                
                ! first get the list of nodes in the ring
                ring1_nodes(1) = nplist(j)
                ne_adjacent = elem_cnct_2d(1,1,nelist(j)) ! get the next element in the +Xi1 direction
                do k = 2,4
                   ring1_nodes(k) = elem_nodes_2d(4,ne_adjacent)
                   ne_adjacent = elem_cnct_2d(1,1,ne_adjacent) ! get the next element in the +Xi1 direction
                enddo ! k
                
                ! assume that the direction for adjustment is defined by location of adjacent rings
                vector(1:3) = nedirection(1:3,j)
                
                ! calculate the ring displacement = j*line_length - (distance from crux)
                point1(1:3) = node_xyz_2d(1,1,1:3,ring1_nodes(1))
                point2(1:3) = node_xyz_2d(1,1,1:3,ring1_nodes(2))
                point3(1:3) = node_xyz_2d(1,1,1:3,ring1_nodes(3))
                point4(1:3) = node_xyz_2d(1,1,1:3,np)
                displace_length = real(j) * line_length - &
                     distance_from_plane_to_point(point1,point2,point3,point4)
                
                ! update the location of the four nodes in the current ring
                do k = 1,4
                   node_xyz_2d(1,1,1:3,ring1_nodes(k)) = &
                        node_xyz_2d(1,1,1:3,ring1_nodes(k)) + &
                        vector(1:3) * displace_length
                enddo
             enddo ! j
          enddo ! i
       endif
       np = np + 1 ! increment to check the next node
    enddo
    
    call enter_exit(sub_name,2)
    
  end subroutine redistribute_mesh_nodes_2d_from_1d

!!!#############################################################################
  
  function inlist(item,ilist)
    
    integer :: item,ilist(:)
    ! Local variables
    integer :: n
    logical :: inlist

    ! --------------------------------------------------------------------------

    inlist = .false.
    do n=1,size(ilist)
       if(item == ilist(n)) inlist = .true.
    enddo
    
  end function inlist
  
!!!#############################################################################
  
  function coord_at_xi(ne,xi,basis)
    
    integer,intent(in) :: ne
    real(dp),intent(in) :: xi(:)
    character(len=*),intent(in) :: basis
    ! Local Variables
    integer :: nn,nv
    real(dp) :: phi(4),phi_10(2),phi_11(2),phi_20(2),phi_21(2),x(4,3,4)
    real(dp) :: coord_at_xi(3)

    ! --------------------------------------------------------------------------
    
    select case (basis)
    case('linear')
       forall (nn=1:4) x(1,1:3,nn) = node_xyz_2d(1,1,1:3,elem_nodes_2d(nn,ne))
       phi(1) = (1.0_dp - xi(1))*(1.0_dp - xi(2))
       phi(2) = xi(1)*(1.0_dp - xi(2))
       phi(3) = (1.0_dp - xi(1))*xi(2)
       phi(4) = xi(1)*xi(2)
       
       coord_at_xi(1:3) = phi(1)*x(1,1:3,1)+phi(2)*x(1,1:3,2)+phi(3)* &
            x(1,1:3,3)+phi(4)*x(1,1:3,4)
       
    case('hermite')
       do nn=1,4
          nv = elem_versn_2d(nn,ne)
          x(1:4,1:3,nn) = node_xyz_2d(1:4,nv,1:3,elem_nodes_2d(nn,ne))
       enddo
       phi_10(1) = (2.0_dp*xi(1)-3.0_dp)*xi(1)*xi(1)+1.0_dp  ! 2xi^3-3xi^2+1
       phi_10(2) = (2.0_dp*xi(2)-3.0_dp)*xi(2)*xi(2)+1.0_dp  ! 2xi^3-3xi^2+1
       phi_20(1) = xi(1)*xi(1)*(3.0_dp-2.0_dp*xi(1))         ! -2xi^3+3xi^2
       phi_20(2) = xi(2)*xi(2)*(3.0_dp-2.0_dp*xi(2))         ! -2xi^3+3xi^2
       phi_11(1) = ((xi(1)-2.0_dp)*xi(1)+1.0_dp)*xi(1)       ! xi^3-2xi^2+xi
       phi_11(2) = ((xi(2)-2.0_dp)*xi(2)+1.0_dp)*xi(2)       ! xi^3-2xi^2+xi
       phi_21(1) = xi(1)*xi(1)*(xi(1)-1.0_dp)                ! xi^3-xi^2
       phi_21(2) = xi(2)*xi(2)*(xi(2)-1.0_dp)                ! xi^3-xi^2
       coord_at_xi(1:3) = phi_10(1)*phi_10(2)*x(1,1:3,1) &
            + phi_20(1)*phi_10(2)*x(1,1:3,2) &
            + phi_10(1)*phi_20(2)*x(1,1:3,3) &
            + phi_20(1)*phi_20(2)*x(1,1:3,4) &
            + phi_11(1)*phi_10(2)*x(2,1:3,1) &
            + phi_21(1)*phi_10(2)*x(2,1:3,2) &
            + phi_11(1)*phi_20(2)*x(2,1:3,3) &
            + phi_21(1)*phi_20(2)*x(2,1:3,4) &
            + phi_10(1)*phi_11(2)*x(3,1:3,1) &
            + phi_20(1)*phi_11(2)*x(3,1:3,2) &
            + phi_10(1)*phi_21(2)*x(3,1:3,3) &
            + phi_20(1)*phi_21(2)*x(3,1:3,4) &
            + phi_11(1)*phi_11(2)*x(4,1:3,1) &
            + phi_21(1)*phi_11(2)*x(4,1:3,2) &
            + phi_11(1)*phi_21(2)*x(4,1:3,3) &
            + phi_21(1)*phi_21(2)*x(4,1:3,4)
    end select
    
  end function coord_at_xi
  
!!!#############################################################################

  function get_local_elem(ne_global)

!!! dummy arguments
    integer,intent(in) :: ne_global
!!! local variables
    integer :: ne
    integer :: get_local_elem

    ! --------------------------------------------------------------------------

    do ne=1,num_elems_2d
       if(ne_global.eq.elems_2d(ne)) get_local_elem = ne
    enddo

  end function get_local_elem

!!!#############################################################################

  function get_local_elem_1d(ne_global)

    integer,intent(in) :: ne_global
    ! Local variables
    integer :: ne
    integer :: get_local_elem_1d

    ! --------------------------------------------------------------------------

    get_local_elem_1d = 0
    do ne=1,num_elems
       if(ne_global.eq.elems(ne)) get_local_elem_1d = ne
    enddo

  end function get_local_elem_1d

!!!###########################################################################

  subroutine write_elem_geometry_2d(elemfile)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_WRITE_ELEM_GEOMETRY_2D" :: WRITE_ELEM_GEOMETRY_2D

    character(len=*),intent(in) :: elemfile
    !     Local Variables
    integer :: ne,ne_count,nglobal_list(4),np,nv
    character(len=132) :: writefile
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------

    sub_name = 'write_elem_geometry_2d'
    call enter_exit(sub_name,1)

    writefile = trim(elemfile)//'.ipelem'
    open(10, file=writefile, status='replace')
    
    !.....write the total number of elems
    write(10,'('' CMISS Version 2.1  ipelem File Version 2'')')
    write(10,'('' Heading: 2D surface from 1D centreline'')')
    write(10,'()')
    write(10,'('' The number of elements is [1]:  '',i6)') num_elems_2d
    write(10,'()')

    !    do ne = 1,num_elems_2d
    ne_count = 1
    ne = 0
    do while (ne_count.le.num_elems_2d)
       ne = ne + 1
       if(elems_2d(ne).gt.0)then
          ne_count = ne_count + 1
          write(10,'('' Element number [    1]:  '',i6)')   elems_2d(ne)
          write(10,'('' The number of geometric Xj-coordinates is [3]: 3'')')
          write(10,'('' The basis function type for geometric variable 1 is [1]:  1'')')
          write(10,'('' The basis function type for geometric variable 2 is [1]:  1'')')
          write(10,'('' The basis function type for geometric variable 3 is [1]:  1'')')
          do np = 1,4
             nglobal_list(np) = nodes_2d(elem_nodes_2d(np,ne))
          enddo
          write(10,'('' Enter the 4 global numbers for basis 1: '',4(i6))') &
               nglobal_list(:)
          do np = 1,4
             if(node_versn_2d(elem_nodes_2d(np,ne)).gt.1)then ! has versions
                nv = elem_versn_2d(np,ne)
                write(10,'('' The version number for occurrence  1 of node'' &
                     &,i7,'', njj=1 is [ 1]: '',i3)') nglobal_list(np), nv
                write(10,'('' The version number for occurrence  1 of node'' &
                     &,i7,'', njj=2 is [ 1]: '',i3)') nglobal_list(np), nv
                write(10,'('' The version number for occurrence  1 of node'' &
                     &,i7,'', njj=3 is [ 1]: '',i3)') nglobal_list(np), nv
             endif
          enddo !np
          write(10,'()')
       endif
    enddo ! ne

  end subroutine write_elem_geometry_2d

!!!#############################################################################

  subroutine write_node_geometry_2d(NODEFILE)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_WRITE_NODE_GEOMETRY_2D" :: WRITE_NODE_GEOMETRY_2D

    character(len=*),intent(in) :: NODEFILE
    !     Local Variables
    integer :: i,np,np_count,nv
    character(len=132) :: writefile
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'write_node_geometry_2d'
    call enter_exit(sub_name,1)

    writefile = trim(nodefile)//'.ipnode'
    open(10, file=writefile, status='replace')
    
    !.....write the total number of nodes
    write(10,'('' CMISS Version 1.21 ipnode File Version 2'')')
    write(10,'('' Heading: '')')
    write(10,'()')
    write(10,'('' The number of nodes is [    1]:  '',i6)') num_nodes_2d
    write(10,'('' Number of coordinates [ 3]:  3'')')
    do i=1,3
       write(10,'('' Do you want prompting for different versions of nj='',i1,'' [N]? Y'')') i
    enddo
    do i=1,3
       write(10,'('' The number of derivatives for coordinate '',i1,'' is [0]: 3'')') i
    enddo

    !    do np = 1,num_nodes_2d
    np_count = 1
    np = 0
    do while (np_count.le.num_nodes_2d)
       np = np + 1
       if(nodes_2d(np).gt.0)then
          np_count = np_count + 1
          write(10,'()')
          write(10,'('' Node number [    1]: '',i6)')  nodes_2d(np)
          do i=1,3
             write(10,'('' The number of versions for nj='',i1,'' is [1]:'',i2)')  i,node_versn_2d(np)
             do nv=1,node_versn_2d(np)
                if(node_versn_2d(np).gt.1) write(10,'('' For version number '',i1,'':'')') nv 
                !...........coordinate          
                write(10,'('' The Xj('',i1,'') coordinate is [ 0.00000E+00]: '',f12.5)') &
                     i,node_xyz_2d(1,nv,i,np)
                write(10,'('' The derivative wrt direction 1 is [ 0.00000E+00]: '',f12.5)') &
                     node_xyz_2d(2,nv,i,np)
                write(10,'('' The derivative wrt direction 2 is [ 0.00000E+00]: '',f12.5)') &
                     node_xyz_2d(3,nv,i,np)
                write(10,'('' The derivative wrt directions 1 & 2 is [ 0.00000E+00]: '',f12.5)') &
                     node_xyz_2d(4,nv,i,np)
             enddo !nv
          end do !i
       endif
    enddo
    close(10)

     call enter_exit(sub_name,2)
 
  end subroutine write_node_geometry_2d

!!!#############################################################################

  subroutine write_surface_geo(element_spline,elem_surfaces,ifile,ncount_point, &
       ncount_loop,ncount_spline,np_offset)
    
    integer :: element_spline(:,:),elem_surfaces(:,:),ifile,ncount_point, &
         ncount_loop,ncount_spline,np_offset
    ! Local variables
    integer:: i,j,k,line1,line2,line3,line4,ne,nk,np,np_highest,np1,np2, &
         num_crux_lines,nv1,nv2
    integer,allocatable :: crux_lines(:,:)
    real(dp) :: phi_1_0,phi_1_1,phi_2_0,phi_2_1,point_xyz(3),xidivn(3)
    logical :: repeat
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'write_surface_geo'
    call enter_exit(sub_name,1)

    allocate(crux_lines(num_elems,3))
    
    forall (i=1:3) xidivn(i) = 0.25_dp * i
    
!!!    Make a gmsh 'point' at each of the surface mesh nodes
    write(ifile,'(''/* Points */'')')
    do np = 1,num_nodes_2d
       ncount_point = ncount_point + 1
       write(ifile,'(''Point('',i8,'') = {'',f12.7,'','',f12.7,'','',f12.7,'',lc};'')') &
            ncount_point,node_xyz_2d(1,1,1:3,np)
    enddo
    
!!!    Now work through each 'ring' of four adjoining surface elements. Points are created
!!!    between adjacent nodes in the Xi1 and Xi2 directions.

    element_spline = 0
    num_crux_lines = 0
    ne = 1
    
    do while (ne.le.num_elems_2d)
       
!!!....make intermediate points, lines, surfaces, at Xi2=0 for model entry .........
       if((elem_cnct_2d(-2,0,ne).eq.0) .or. &
            (node_versn_2d(elem_nodes_2d(1,ne)).eq.6.or.&
            node_versn_2d(elem_nodes_2d(2,ne)).eq.6))then
          ! location of points in the Xi+1 direction on the Xi2=0 boundary
          do k = 0,3
             if(elem_cnct_2d(-2,0,ne).eq.0 .or. element_spline(1,ne+k).eq.0)then  ! no line (or points) here already
                np1 = elem_nodes_2d(1,ne+k)
                np2 = elem_nodes_2d(2,ne+k)
                ! check that the line hasn't already been made. If it has, use it.
                repeat = .false.
                do i = 1,num_crux_lines
                   if((np1.eq.crux_lines(i,2).and.np2.eq.crux_lines(i,3)).or. &
                        (np1.eq.crux_lines(i,3).and.np2.eq.crux_lines(i,2)))then
                      repeat = .true.
                      element_spline(1,ne+k) = -crux_lines(i,1)
                   endif
                enddo

                if(.not.repeat)then
                   nv1 = elem_versn_2d(1,ne+k)
                   nv2 = elem_versn_2d(2,ne+k)
                   nk = 2 ! calculate the location in the Xi+1 direction
                   do i = 1,3
                      phi_1_0 = 1.0_dp - 3.0_dp * xidivn(i)**2 + 2.0_dp * xidivn(i)**3
                      phi_1_1 = xidivn(i) * (xidivn(i) - 1.0_dp)**2
                      phi_2_0 = xidivn(i)**2 * (3.0_dp - 2.0_dp * xidivn(i))
                      phi_2_1 = xidivn(i)**2 * (xidivn(i) - 1.0_dp)
                      do j = 1,3
                         point_xyz(j) = phi_1_0 * node_xyz_2d(1,1,j,np1) + &
                              phi_2_0 * node_xyz_2d(1,1,j,np2) &
                              + phi_1_1 * node_xyz_2d(nk,nv1,j,np1) * scale_factors_2d(2,ne+k) &
                              + phi_2_1 * node_xyz_2d(nk,nv2,j,np2) * scale_factors_2d(6,ne+k)
                      enddo
                      
                      ncount_point = ncount_point + 1
                      write(ifile,'(''Point('',i8,'') = {'',f12.7,'','',f12.7,'','',f12.7,'',lc};'')') &
                           ncount_point, point_xyz(1:3)
                   enddo ! i
                   ncount_spline = ncount_spline + 1
                   write(ifile,'(''Spline('',I8,'') = {'',I8,'','',I8,'','' &
                        &,I8,'','',I8,'','',I8,''};'')') ncount_spline, np1+np_offset, &
                        ncount_point-2, ncount_point-1, ncount_point, np2+np_offset
                   element_spline(1,ne+k) = ncount_spline
                   num_crux_lines = num_crux_lines + 1
                   crux_lines(num_crux_lines,1) = ncount_spline
                   crux_lines(num_crux_lines,2) = np1+np_offset
                   crux_lines(num_crux_lines,3) = np2+np_offset
!                   if(elem_cnct_2d(-2,0,ne).ne.0) element_spline(1,elem_cnct_2d(-2,1,ne)) = -ncount_spline
                endif
             endif
          enddo  ! k
       endif  ! elem_cnct etc
       
!!!.......make intermediate points, lines, surfaces, at Xi2=1 for each 'ring' .........
       
       ! location of points in the Xi+1 direction on the Xi2=1 boundary
       np_highest = elem_nodes_2d(3,ne)
       do k = 0,3
          np1 = elem_nodes_2d(3,ne+k)
          np2 = elem_nodes_2d(4,ne+k)
          nv1 = elem_versn_2d(3,ne+k)
          nv2 = elem_versn_2d(4,ne+k)
          nk = 2 
          do i = 1,3
             phi_1_0 = 1.0_dp - 3.0_dp * xidivn(i)**2 + 2.0_dp * xidivn(i)**3
             phi_1_1 = xidivn(i) * (xidivn(i) - 1.0_dp)**2
             phi_2_0 = xidivn(i)**2 * (3.0_dp - 2.0_dp * xidivn(i))
             phi_2_1 = xidivn(i)**2 * (xidivn(i) - 1.0_dp)
             forall (j=1:3) point_xyz(j) = phi_1_0 * node_xyz_2d(1,1,j,np1) &
                  + phi_2_0 * node_xyz_2d(1,1,j,np2) &
                  + phi_1_1 * node_xyz_2d(nk,nv1,j,np1) * scale_factors_2d(10,ne+k) &
                  + phi_2_1 * node_xyz_2d(nk,nv2,j,np2) * scale_factors_2d(14,ne+k) 
             
             ncount_point = ncount_point + 1
             write(ifile,'(''Point('',i8,'') = {'',f12.7,'','',f12.7,'','',f12.7,'',lc};'')') &
                  ncount_point,point_xyz(1:3)
          enddo ! i
          ncount_spline = ncount_spline + 1
          write(ifile,'(''Spline('',I8,'') = {'',I8,'','',I8,'','',I8,'','' &
               & ,I8,'','',I8,''};'')') ncount_spline, np1+np_offset, ncount_point-2, &
               ncount_point-1, ncount_point, np2+np_offset
          element_spline(3,ne+k) = ncount_spline
          if(elem_cnct_2d(2,0,ne+k).gt.0) element_spline(1,elem_cnct_2d(2,1,ne+k)) = ncount_spline
       enddo ! k
       
       ! location of points in the Xi+2 direction on the Xi1=0 boundary
       do k = 0,3
          np1 = elem_nodes_2d(1,ne+k)
          np2 = elem_nodes_2d(3,ne+k)
          nv1 = elem_versn_2d(1,ne+k)
          nv2 = elem_versn_2d(3,ne+k)
          nk = 3
          do i = 1,3
             phi_1_0 = 1.0_dp - 3.0_dp * xidivn(i)**2 + 2.0_dp * xidivn(i)**3
             phi_1_1 = xidivn(i) * (xidivn(i) - 1.0_dp)**2
             phi_2_0 = xidivn(i)**2 * (3.0_dp - 2.0_dp * xidivn(i))
             phi_2_1 = xidivn(i)**2 * (xidivn(i) - 1.0_dp)
             forall (j=1:3) point_xyz(j) = phi_1_0 * node_xyz_2d(1,1,j,np1) &
                  + phi_2_0 * node_xyz_2d(1,1,j,np2) &
                  + phi_1_1 * node_xyz_2d(nk,nv1,j,np1) * scale_factors_2d(3,ne+k) &
                  + phi_2_1 * node_xyz_2d(nk,nv2,j,np2) * scale_factors_2d(11,ne+k)
             
             ncount_point = ncount_point + 1
             write(ifile,'(''Point('',i8,'') = {'',f12.7,'','',f12.7,'','',f12.7,'',lc};'')') &
                  ncount_point,point_xyz(1:3)
          enddo ! i
          ncount_spline = ncount_spline + 1
          write(ifile,'(''Spline('',I8,'') = {'',I8,'','',I8,'','',I8,'','' &
               & ,I8,'','',I8,''};'')') ncount_spline, np1+np_offset, ncount_point-2, &
               ncount_point-1, ncount_point, np2+np_offset
          element_spline(4,ne+k) = ncount_spline
          element_spline(2,elem_cnct_2d(-1,1,ne+k)) = ncount_spline
       enddo ! k
       
       do k = 0,3
          line1 = element_spline(1,ne+k)
          line2 = element_spline(2,ne+k)
          line3 = -element_spline(3,ne+k)
          line4 = -element_spline(4,ne+k)
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','',i8,'','',i8,''};'')') &
               ncount_loop, line1, line2, line3, line4
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Surface('',I8,'') = {'',I8,''};'')') ncount_loop, ncount_loop - 1
          elem_surfaces(3,ne+k) = ncount_loop
       enddo ! k
       
       ne = ne + 4
       
    enddo ! while (ne.le.num_elems_2d)

    deallocate(crux_lines)
    
    call enter_exit(sub_name,2)

  end subroutine write_surface_geo
  
!!!#############################################################################

  subroutine write_3d_geo(element_spline,elem_surfaces,ifile,ncount_point, &
    ncount_loop,ncount_spline)

    integer,intent(in) :: ifile
    integer :: element_spline(:,:),elem_surfaces(:,:),ncount_point,ncount_loop, &
         ncount_spline
    ! Local variables
    integer :: i,j,k,line1,line2,line3,line4,ncount_cap_entry=0,ncount_cap_exit=0, &
         ncount_inner=0,ncount_centre=0,ncount_phys_vol=0,ncount_spline_0, &
         ncount_surface=0,ncount_volume=0,ncount_wall=0,ne,ne_next,np_highest,np1,np2
    integer,allocatable :: centre_points(:),ncap_entry(:),ncap_exit(:), &
         ncentre(:),ninner(:),nphys_vol(:),node_spoke(:,:),nwall(:)
    real(dp) :: point_xyz_centre(3), xidivn(3)
    character(len=60) :: sub_name
    
    ! --------------------------------------------------------------------------
    
    sub_name = 'write_3d_geo'
    call enter_exit(sub_name,1)
        
    allocate(centre_points(num_nodes_2d))
    allocate(ncap_entry(20))  ! assuming could have multiple inlets
    allocate(ncap_exit(num_elems_2d))  
    allocate(ninner(num_elems_2d*2))  
    allocate(ncentre(num_elems_2d))
    allocate(nphys_vol(num_elems_2d))
    allocate(nwall(num_elems_2d))  
    allocate(node_spoke(2,num_nodes_2d))
    node_spoke = 0

    forall (i=1:3) xidivn(i) = 0.25_dp * i

    ne = 1
    do while (ne.le.num_elems_2d)
!!!....the following works on four adjacent surface elements

!!!......... make intermediate points, lines, surfaces, at Xi2=0 for model entry .........
       
       if((elem_cnct_2d(-2,0,ne).eq.0) .or. &
            (node_versn_2d(elem_nodes_2d(1,ne)).eq.6.or.&
            node_versn_2d(elem_nodes_2d(2,ne)).eq.6))then
          ! location of points in the Xi+1 direction on the Xi2=0 boundary
          point_xyz_centre = 0.0_dp
          do k = 0,3
             np1 = elem_nodes_2d(1,ne+k)
             forall (j = 1:3) point_xyz_centre(j) = &
                  point_xyz_centre(j) + 0.25_dp * node_xyz_2d(1,1,j,np1)
          enddo  ! k

          if(elem_cnct_2d(-2,0,ne).eq.0)then
             ! make a point at the centre of the ring
             ncount_point = ncount_point + 1
             write(ifile,'(''Point('',i8,'') = {'',f12.7,'','',f12.7,'','' &
                  &,f12.7,'',lc};'')') ncount_point,point_xyz_centre(1:3)
             
             ! make a 'spoke' from centre of the ring to each surface node
             do k = 0,3
                np1 = elem_nodes_2d(1,ne+k)
                centre_points(np1) = ncount_point   ! record the centre point number for this 'ring'
                ncount_spline = ncount_spline + 1
                write(ifile,'(''Line('',i8,'') = {'',i8,'','',i8,''};'')') &
                     ncount_spline,ncount_point,np1
                node_spoke(1,np1) = ncount_spline
             enddo  ! k

             ! make surfaces at the entry == a 'cap' of four surfaces
             do k = 0,3
                line1 = ncount_spline + k - 3 
                line2 = element_spline(1,ne+k)
                if(k.lt.3)then
                   line3 = -(line1 + 1)
                else
                   line3 = -(ncount_spline - 3)
                endif
                ncount_loop = ncount_loop + 1
                write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','' &
                     &,i8,''};'')') ncount_loop, line1, line2, line3
                ncount_loop = ncount_loop + 1
                write(ifile,'(''Surface('',i8,'') = {'',i8,''};'')') &
                     ncount_loop, ncount_loop - 1
                elem_surfaces(1,ne+k) = ncount_loop
                ncount_cap_entry = ncount_cap_entry + 1
                ncap_entry(ncount_cap_entry) = ncount_loop
             enddo  ! k
          endif
       endif  ! elem_cnct etc
       
!!!......... make intermediate points, lines, surfaces, at Xi2=1 for each 'ring' .........

       ! location of points in the Xi+1 direction on the Xi2=1 boundary
       point_xyz_centre = 0.0_dp
       np_highest = elem_nodes_2d(3,ne)
       ncount_spline_0 = ncount_spline + 1
       do k = 0,3
          np1 = elem_nodes_2d(3,ne+k)
          np2 = elem_nodes_2d(4,ne+k)
          if(node_versn_2d(np1).ge.2.or.node_versn_2d(np2).ge.2)then  ! only for when there is a crux node
             if(np1.gt.np_highest) np_highest = np1
             forall (j = 1:3) point_xyz_centre(j) = point_xyz_centre(j) &
                  + 0.2_dp * node_xyz_2d(1,1,j,np1)
          else
             forall (j = 1:3) point_xyz_centre(j) = point_xyz_centre(j) &
                  + 0.25_dp * node_xyz_2d(1,1,j,np1)
          endif
       enddo ! k

       ncount_point = ncount_point + 1
       if(node_versn_2d(np1).ge.2.or.node_versn_2d(np2).ge.2)then  ! only for when there is a crux node
          np_highest = elem_nodes_2d(2,elem_cnct_2d(1,1,elem_cnct_2d(2,1,ne)))
          forall (j = 1:3) point_xyz_centre(j) = point_xyz_centre(j) &
               + 0.2_dp * node_xyz_2d(1,1,j,np_highest)
          centre_points(np_highest) = ncount_point   ! record the centre point number for this 'ring'
       endif

       write(ifile,'(''Point('',i8,'') = {'',f12.7,'','',f12.7,'','' &
            &,f12.7,'',lc};'')') ncount_point,point_xyz_centre(1:3)

       do k = 0,3
          np1 = elem_nodes_2d(3,ne+k)
          centre_points(np1) = ncount_point   ! record the centre point number for this 'ring'
          ncount_spline = ncount_spline + 1
          write(ifile,'(''Line('',i8,'') = {'',i8,'','',i8,''};'')') &
               ncount_spline,ncount_point,np1
          node_spoke(1,np1) = ncount_spline
       enddo  ! k

       if(node_versn_2d(np1).ge.2.or.node_versn_2d(np2).ge.2)then  ! only for when there is a crux node
          np_highest = elem_nodes_2d(2,elem_cnct_2d(1,1,elem_cnct_2d(2,1,ne)))
          
          if(elems_at_node_2d(np_highest,0).eq.6)then
             do i = 1,elems_at_node_2d(np_highest,0)
                ne_next = elems_at_node_2d(np_highest,i)
                np1 = elem_nodes_2d(1,ne_next)
                if(np1.eq.np_highest) np1 = &
                     elem_nodes_2d(2,elems_at_node_2d(np_highest,i))
                if(node_spoke(1,np1).eq.0)then
                   ! make a spoke from the centre to this node
                   ncount_spline = ncount_spline + 1
                   write(ifile,'(''Line('',i8,'') = {'',i8,'','',i8,''};'')') &
                        ncount_spline,ncount_point,np1
                   node_spoke(1,np1) = ncount_spline
                   elem_surfaces(1,elem_cnct_2d(-2,1,ne_next)) = ncount_loop
                endif
             enddo
          endif
       endif

!!!....make surface elements at the Xi2=1 end
       do k = 0,3
          line1 = node_spoke(1,elem_nodes_2d(3,ne+k))
          line2 = element_spline(3,ne+k)
          line3 = -node_spoke(1,elem_nodes_2d(4,ne+k))
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','',i8,''};'')') &
               ncount_loop, line1, line2, line3
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Surface('',i8,'') = {'',i8,''};'')') &
               ncount_loop, ncount_loop - 1
          elem_surfaces(2,ne+k) = ncount_loop
          if(elem_cnct_2d(2,0,ne+k).ne.0)then
             elem_surfaces(1,elem_cnct_2d(2,1,ne+k)) = ncount_loop
          else
             ! this is an exit 'cap'. store to output as a physical surface
             ncount_cap_exit = ncount_cap_exit + 1
             ncap_exit(ncount_cap_exit) = ncount_loop
          endif
          ncount_inner = ncount_inner + 1
          ninner(ncount_inner) = ncount_loop
       enddo  ! k
             
       if(node_versn_2d(np1).ge.2.or.node_versn_2d(np2).ge.2)then  ! only for crux node
          ncount_spline = ncount_spline + 1
          write(ifile,'(''Line('',i8,'') = {'',i8,'','',i8,''};'')') &
               ncount_spline,ncount_point,np_highest
          node_spoke(1,np_highest) = ncount_spline

          ne_next = elem_cnct_2d(1,1,elem_cnct_2d(2,1,ne))
          line1 = node_spoke(1,elem_nodes_2d(1,ne_next))
          line2 = element_spline(1,ne_next)
          line3 = -node_spoke(1,elem_nodes_2d(2,ne_next))

          ncount_loop = ncount_loop + 1
          write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','',i8,''};'')') &
               ncount_loop, line1, line2, line3
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Surface('',i8,'') = {'',i8,''};'')') &
               ncount_loop, ncount_loop - 1
          ncount_inner = ncount_inner + 1
          ninner(ncount_inner) = ncount_loop

          elem_surfaces(1,ne_next) = ncount_loop
          ne_next = elem_cnct_2d(-2,1,ne_next)
          elem_surfaces(1,ne_next) = ncount_loop

          ne_next = elem_cnct_2d(-1,1,elem_cnct_2d(2,1,elem_cnct_2d(-1,1,ne)))
          line1 = node_spoke(1,elem_nodes_2d(1,ne_next))
          line2 = element_spline(1,ne_next)
          line3 = -node_spoke(1,elem_nodes_2d(2,ne_next))
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','',i8,''};'')') &
               ncount_loop, line1, line2, line3
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Surface('',i8,'') = {'',i8,''};'')') &
               ncount_loop, ncount_loop - 1
          ncount_inner = ncount_inner + 1
          ninner(ncount_inner) = ncount_loop

          elem_surfaces(1,ne_next) = ncount_loop
          ne_next = elem_cnct_2d(-2,1,ne_next)
          elem_surfaces(1,ne_next) = ncount_loop
          
          if(elems_at_node_2d(np_highest,0).eq.6)then
             do i = 1,elems_at_node_2d(np_highest,0)
                ne_next = elems_at_node_2d(np_highest,i)
                if(elem_surfaces(1,ne_next).eq.0)then
                   ! make an extra surface here
                   np1 = elem_nodes_2d(1,ne_next)
                   np2 = elem_nodes_2d(2,ne_next)
                   line1 = node_spoke(1,np1)
                   line2 = element_spline(1,ne_next)
                   line3 = -node_spoke(1,np2)
                   ncount_loop = ncount_loop + 1
                   write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','' &
                        &,i8,'','',i8,''};'')') &
                        ncount_loop, line1, line2, line3
                   ncount_loop = ncount_loop + 1
                   write(ifile,'(''Surface('',i8,'') = {'',i8,''};'')') &
                        ncount_loop, ncount_loop - 1
                   ncount_inner = ncount_inner + 1
                   ninner(ncount_inner) = ncount_loop
                   elem_surfaces(1,ne_next) = ncount_loop
                   ne_next = elem_cnct_2d(-2,1,ne_next)
                   elem_surfaces(1,ne_next) = ncount_loop
                   
                   ne_next = elem_cnct_2d(1,1,ne_next)
                   np1 = elem_nodes_2d(1,ne_next)
                   np2 = elem_nodes_2d(2,ne_next)
                   line1 = node_spoke(1,np1)
                   line2 = element_spline(1,ne_next)
                   line3 = -node_spoke(1,np2)
                   ncount_loop = ncount_loop + 1
                   write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','' &
                        &,i8,'','',i8,''};'')') &
                        ncount_loop, line1, line2, line3
                   ncount_loop = ncount_loop + 1
                   write(ifile,'(''Surface('',i8,'') = {'',i8,''};'')') &
                        ncount_loop, ncount_loop - 1
                   ncount_inner = ncount_inner + 1
                   ninner(ncount_inner) = ncount_loop
                   elem_surfaces(1,ne_next) = ncount_loop
                   ne_next = elem_cnct_2d(-2,1,ne_next)
                   elem_surfaces(1,ne_next) = ncount_loop
                   
                endif
             enddo
          endif
       endif

!!!.........Make line along the centre
       ncount_spline = ncount_spline + 1
       write(ifile,'(''Line('',i8,'') = {'',i8,'','',i8,''};'')') ncount_spline, &
            centre_points(elem_nodes_2d(1,ne)),ncount_point

!!! Make surfaces from the centreline
       do k = 0,3
          line1 = node_spoke(1,elem_nodes_2d(1,ne+k))
          line2 = element_spline(4,ne+k)
          line3 = -node_spoke(1,elem_nodes_2d(3,ne+k))
          line4 = -ncount_spline ! the newest line
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Line Loop('',i8,'') = {'',i8,'','',i8,'','' &
               &,i8,'','',i8,''};'')') &
               ncount_loop, line1, line2, line3, line4
          ncount_loop = ncount_loop + 1
          write(ifile,'(''Surface('',i8,'') = {'',i8,''};'')') &
               ncount_loop, ncount_loop - 1
          ncount_centre = ncount_centre + 1
          ncentre(ncount_centre) = ncount_loop

          elem_surfaces(4,ne+k) = ncount_loop
          if(k.gt.0) elem_surfaces(5,ne+k-1) = ncount_loop
       enddo
       elem_surfaces(5,ne) = elem_surfaces(4,ne+1)
       elem_surfaces(5,ne+3) = elem_surfaces(4,ne)

!!! Make surface loops and volumes
       do k = 0,3
          ncount_volume = ncount_volume + 1
          write(ifile,'(''Surface Loop('',i8,'') = {'',i8,'','' &
               &,i8,'','',i8,'','',i8,'','',i8,''};'')') &
               ncount_volume,elem_surfaces(1:5,ne+k)
            
          ncount_volume = ncount_volume + 1
          write(ifile,'(''Volume('',i8,'') = {'',i8,''};'')') &
               ncount_volume,ncount_volume-1

          ncount_wall = ncount_wall + 1
          nwall(ncount_wall) = elem_surfaces(3,ne+k)
          ncount_phys_vol = ncount_phys_vol + 1
          nphys_vol(ncount_phys_vol) = ncount_volume
       enddo

       ne = ne + 4

    enddo ! while (ne.le.num_elems_2d)

    write(ifile,'(/''/* Physical surface for entry caps */'')')
    write(ifile,'(/''Physical Surface(1) = {'')', advance = "no")
    do i = 1,ncount_cap_entry-1
       write(ifile,'(i6,'','')', advance = "no") ncap_entry(i)
    enddo
    write(ifile,'(i6,''};'')') ncap_entry(ncount_cap_entry)
    
    write(ifile,'(/''/* Physical surface for exit caps */'')')
    write(ifile,'(/''Physical Surface(2) = {'')', advance = "no")
    do i = 1,ncount_cap_exit-1
       write(ifile,'(i6,'','')', advance = "no") ncap_exit(i)
    enddo
    write(ifile,'(i6,''};'')') ncap_exit(ncount_cap_exit)
    
    write(ifile,'(/''/* Physical surface for walls */'')')
    write(ifile,'(/''Physical Surface(3) = {'')', advance = "no")
    do i = 1,ncount_wall-1
       write(ifile,'(i6,'','')', advance = "no") nwall(i)
    enddo
    write(ifile,'(i6,''};'')') nwall(ncount_wall)
    
    write(ifile,'(/''/* Physical surface for centres */'')')
    write(ifile,'(/''Physical Surface(4) = {'')', advance = "no")
    do i = 1,ncount_centre-1
       write(ifile,'(i6,'','')', advance = "no") ncentre(i)
    enddo
    write(ifile,'(i6,''};'')') ncentre(ncount_centre)
    
    write(ifile,'(/''Physical Volume(1) = {'')', advance = "no")
    do i = 1,ncount_phys_vol-1
       write(ifile,'(i6,'','')', advance = "no") nphys_vol(i)
    enddo
    write(ifile,'(i6,''};'')') nphys_vol(ncount_phys_vol)

    write(ifile,'(/)')
    write(ifile,'(''Mesh.Algorithm = 3;'')') ! Anisotropic
    write(ifile,'(''Mesh.Smoothing = 4;'')')
    write(ifile,'(''Mesh.Algorithm3D = 2;'')') ! Netgen

    close(ifile)

    deallocate(centre_points)
    deallocate(ncap_entry)
    deallocate(ncap_exit)
    deallocate(ninner)
    deallocate(ncentre)
    deallocate(nphys_vol)
    deallocate(nwall)
    deallocate(node_spoke)

    call enter_exit(sub_name,2)
    
  end subroutine write_3d_geo

!
!###################################################################################
!
  subroutine splitwords(line, words, nw )
    implicit none
    character(*), intent(in)  :: line
    character(*), intent(out) :: words(:)
    integer,      intent(out) :: nw
    character(len(words)) :: buf( size(words) )
    integer :: k, ios

    nw = 0 ; words(:) = ""

    do k = 1, size(words)
        read( line, *, iostat=ios ) buf( 1 : k )
        if ( ios /= 0 ) exit
        nw = k
        words( 1 : nw ) = buf( 1 : nw )
    enddo

  end subroutine splitwords


!###################################################################################

  subroutine words_to_real( words,ints, ni )
    use arrays,only: dp
    implicit none
    character(*), intent(in)  :: words(:)
    real(dp),     intent(out) :: ints(:)
    integer,      intent(out) :: ni
    integer :: k,  ios
    real(dp) :: val

    ni = 0 ; ints(:) = 0

    do k = 1, size(words)
        read( words( k ), *, iostat=ios ) val
        if ( ios /= 0 ) cycle
        ni = ni + 1
        if ( ni > size(ints) ) stop "size(ints) too small"
        ints( ni ) = val
    enddo

  end subroutine words_to_real

! 
!###########################################################################################
! 
end module geometry
