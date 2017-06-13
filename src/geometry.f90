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
  use other_consts
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
  public define_mesh_geometry_test
  public define_node_geometry
  public define_rad_from_file
  public define_rad_from_geom
  public element_connectivity_1d
  public evaluate_ordering
  public set_initial_volume
  public volume_of_mesh

contains
!
!###################################################################################
!
!*add_mesh:* Reads in an ipmesh file and adds this mesh to the terminal branches of an existing tree geometry
  subroutine add_mesh(AIRWAY_MESHFILE)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ADD_MESH" :: ADD_MESH
    use arrays,only: dp,elems,elem_cnct,elem_direction,elem_field,&
         elem_nodes,elem_ordrs,elem_symmetry,&
         nodes,node_xyz,num_elems,&
         num_nodes,num_units,units
    use indices,only: ne_length,ne_radius,ne_a_A, ne_vol
    use other_consts,only: PI
    use diagnostics, only: enter_exit
    implicit none

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

!
!###################################################################################
!
!*add_matching_mesh:*
  subroutine add_matching_mesh()
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ADD_MATCHING_MESH" :: ADD_MATCHING_MESH
    use arrays,only: dp,elems,elem_cnct,elem_direction,elem_field,&
         elem_nodes,elem_ordrs,elem_symmetry,&
         nodes,node_xyz,num_elems,&
         num_nodes,num_units,units
    use indices,only: ne_length,ne_radius,ne_a_A, ne_vol
    use other_consts,only: PI
    use diagnostics, only: enter_exit
    implicit none
    character(len=60) :: sub_name

    sub_name = 'add_matching_mesh'
    call enter_exit(sub_name,1)


    call enter_exit(sub_name,2)

  end subroutine add_matching_mesh

!
!###################################################################################
!
!*append_units:* Appends terminal units at the end of a tree structure
  subroutine append_units()
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_APPEND_UNITS" :: APPEND_UNITS
    use arrays,only: dp, elem_cnct,elem_symmetry,elem_units_below,&
         num_elems,num_units,units,unit_field
    use indices,only: num_nu
    use diagnostics, only: enter_exit
    implicit none

    integer :: ne,ne0,nu
    character(len=60) :: sub_name

    sub_name = 'append_units'
    call enter_exit(sub_name,1)

    num_units = 0
    DO ne=1,num_elems
       IF(elem_cnct(1,0,ne).eq.0)THEN
          num_units=num_units+1
       ENDIF
    ENDDO

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
    DO ne=1,num_elems
       IF(elem_cnct(1,0,ne).eq.0)THEN
          nu=nu+1
          units(nu)=ne     !Set up units array containing terminals
          elem_units_below(ne)=1
       ENDIF
    ENDDO

    ! count the effective number of elements below each branch
    do ne=num_elems,2,-1
       ne0=elem_cnct(-1,1,ne)
       elem_units_below(ne0) = elem_units_below(ne0) &
            + elem_units_below(ne)*elem_symmetry(ne)
    enddo !ne

    call enter_exit(sub_name,2)

  end subroutine append_units

!
!###################################################################################
!
  !*define_1d_elements:* Reads in an element ipelem file to define a geometry
  subroutine define_1d_elements(ELEMFILE)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_1D_ELEMENTS" :: DEFINE_1D_ELEMENTS

    use arrays,only: dp, elem_direction,elem_field,elems,elem_cnct,elem_nodes,&
         elem_ordrs,elem_symmetry,elems_at_node,elem_units_below,&
         expansile,node_xyz,num_elems,num_nodes
    use indices
    use diagnostics, only: enter_exit
    implicit none

    character(len=MAX_FILENAME_LEN), intent(in) :: ELEMFILE
    !     Local Variables
    integer :: ibeg,iend,ierror,i_ss_end,j,ne,ne_global,&
         nn,np,np1,np2,np_global
    character(LEN=132) :: ctemp1
    character(LEN=40) :: sub_string
    character(len=60) :: sub_name

    sub_name = 'define_1d_elements'
    call enter_exit(sub_name,1)

    open(10, file=ELEMFILE, status='old')

    read_number_of_elements : do
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "elements")> 0) then
          call get_final_integer(ctemp1,num_elems)
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
    if(allocated(expansile)) deallocate(expansile)
    allocate(expansile(num_elems))

!!! initialise element arrays
    elems=0
    elem_nodes=0
    elem_symmetry = 1
    elem_field = 0.0_dp
    expansile = .false.

    ne=0

    read_an_element : do
       !.......read element number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Element")> 0) then
          call get_final_integer(ctemp1,ne_global) !get element number
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
       elem_field(ne_length,ne) = DSQRT((node_xyz(1,np2) - &
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

    call enter_exit(sub_name,2)

  END subroutine define_1d_elements
!
!###################################################################################
!
!*define_mesh_geometry_test:*
  subroutine define_mesh_geometry_test()
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_MESH_GEOMETRY_TEST" :: DEFINE_MESH_GEOMETRY_TEST
    use arrays,only: dp,nodes,node_field,node_xyz,num_nodes,&
         elem_direction,elem_field,elems,elem_cnct,elem_nodes,&
         elem_ordrs,elem_symmetry,elems_at_node,elem_units_below,&
         expansile,node_xyz,num_elems,num_nodes
    use indices
    implicit none

    !     Local Variables
    integer :: j,ne,np,np1,np2

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
    elem_field(ne_radius,201:400) = dsqrt(elem_field(ne_radius,101)**2/2.0_dp)

    ! calculate the element lengths and directions
    do ne=1,num_elems
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       elem_field(ne_length,ne) = DSQRT((node_xyz(1,np2) - &
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

  end subroutine define_mesh_geometry_test
!
!###################################################################################
!
  subroutine define_node_geometry(NODEFILE)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_NODE_GEOMETRY" :: DEFINE_NODE_GEOMETRY

  !*define_node_geometry:* Reads in an ipnode file to define a tree geometry
    use arrays,only: dp,nodes,node_field,node_xyz,num_nodes
    use diagnostics, only: enter_exit
    use indices
    use other_consts, only: MAX_FILENAME_LEN
    implicit none

    character(len=MAX_FILENAME_LEN), intent(in) :: NODEFILE !Input nodefile
    !     Local Variables
    integer :: i,ierror,np,np_global,&
         num_versions,nv,NJT
    character(LEN=132) :: ctemp1
    LOGICAL :: versions
    real(dp) :: point
    character(len=60) :: sub_name

    sub_name = 'define_node_geometry'
    call enter_exit(sub_name,1)

    versions = .TRUE.
    NJT = 0
    open(10, file=NODEFILE, status='old')

    !.....read in the total number of nodes. read each line until one is found
    !.....that has the correct keyword (nodes). then return the integer that is
    !.....at the end of the line
    read_number_of_nodes : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "nodes")> 0) then !keyword "nodes" is found in ctemp1
          call get_final_integer(ctemp1,num_nodes) !return the final integer
          exit read_number_of_nodes !exit the named do loop
       endif
    end do read_number_of_nodes

    if(allocated(nodes)) deallocate (nodes)
    allocate (nodes(num_nodes))
    if(allocated(node_xyz)) deallocate (node_xyz)
    allocate (node_xyz(3,num_nodes))
    if(allocated(node_field)) deallocate (node_field)
    allocate (node_field(num_nj,num_nodes))
    nodes = 0 !initialise node index values
    node_xyz = 0.0_dp !initialise
    node_field = 0.0_dp !initialise

    !.....read in the number of coordinates
    read_number_of_coords : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "coordinates")> 0) then !keyword "coordinates" is found
          call get_final_integer(ctemp1,NJT) !return the final integer
          exit read_number_of_coords !exit the named do loop
       endif
    end do read_number_of_coords

    !.....check whether versions are prompted (>1)
    read_versions : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "different")> 0) then !keyword "different" is found
          if(index(ctemp1, " N")> 0) then !keyword " N" is found
             versions=.false.
          endif
          exit read_versions !exit the named do loop
       endif
    end do read_versions

!!! WARNING :: following should be in general code
    ! note that only the first version of coordinate is currently read in

    !.....read the coordinate, derivative, and version information for each node.
    np=0
    read_a_node : do !define a do loop name
       !.......read node number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Node")> 0) then
          call get_final_integer(ctemp1,np_global) !get node number
          np=np+1
          nodes(np)=np_global
          !.......read coordinates and derivatives
          do i=1,NJT ! for the NJT coordinates
             !...........coordinate
             num_versions=1
             if(versions)then
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                call get_final_integer(ctemp1,num_versions)
             endif
             if(num_versions > 1)then
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !temporary line
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                call get_final_real(ctemp1,point)
                do nv=2,num_versions
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !temporary line
                   read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                enddo
             else
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                call get_final_real(ctemp1,point)
             endif
             node_xyz(i,np)=point
          end do !i

       endif !index
       if(np.ge.num_nodes) exit read_a_node
    end do read_a_node

    close(10)

    call enter_exit(sub_name,2)

  END subroutine define_node_geometry
!
!###################################################################################
!
  subroutine define_rad_from_file(FIELDFILE, radius_type_in)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_RAD_FROM_FILE" :: DEFINE_RAD_FROM_FILE

  !*define_rad_from_file:* reads in a radius field associated with an aiway tree
! and assigns radius information to each element, also calculates volume of each
! element
    use arrays,only: dp,elem_field,elem_cnct,elem_nodes,&
         elems_at_node,num_elems,num_nodes,node_field
    use indices,only: ne_a_A,ne_length,ne_radius,ne_vol,&
      ne_radius_in,ne_radius_out

    use diagnostics, only: enter_exit
    implicit none

    character(len=MAX_FILENAME_LEN), intent(in) :: FIELDFILE
    character(len=MAX_STRING_LEN), optional ::  radius_type_in
    !     Local Variables
    integer :: ierror,ne,np,np1,np2,np_global,surround
    character(len=MAX_STRING_LEN) ::  radius_type
    character(LEN=132) :: ctemp1
    LOGICAL :: versions
    real(dp) :: constrict,radius
    character(len=60) :: sub_name

    sub_name = 'define_rad_from_file'
    call enter_exit(sub_name,1)
    
    versions = .TRUE.
    if(present(radius_type_in))then
      radius_type = radius_type_in
    else
      radius_type = 'no_taper'
    endif

!!! note that 'constrict' should not be used here (so is set to 1.0).
!!! this should be specified and used as part of a simulation, not when
!!! reading in airway geometry
    constrict = 1.0_dp
    open(10, file=FIELDFILE, status='old')

    !.....check whether versions are prompted (>1)
    read_versions : do !define a do loop name
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1 !read a line into ctemp1
       if(index(ctemp1, "different")> 0) then !keyword "different" is found
          if(index(ctemp1, " N")> 0) then !keyword " N" is found
             versions=.false.
          endif
          exit read_versions !exit the named do loop
       endif
    end do read_versions

    np = 0
    !.....read the coordinate, derivative, and version information for each node.
    read_a_node : do !define a do loop name
       !.......read node number
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       if(index(ctemp1, "Node")> 0) then
          call get_final_integer(ctemp1,np_global) !get global node number
          ! find the corresponding local node number
          call get_local_node(np_global,np) ! get local node np for global node
          surround=elems_at_node(np,0)         !get number of surrounding elems
          ne=elems_at_node(np,1)  !First element at this node
          if(surround==1)then !only one element at this node so either a terminal or inlet
             if(radius_type.eq.'taper')then !inlet radius needs to be defined
               if(elem_cnct(-1,0,ne).eq.0)then!Inlet as it has no parent need to set up radius into this vessel
                 read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                 read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                 if(index(ctemp1, "version number")>0) then
                    read(unit=10, fmt="(a)", iostat=ierror) ctemp1
                 endif
                 if(index(ctemp1, "value")> 0) then
                    call get_final_real(ctemp1,radius)
                    elem_field(ne_radius_in,ne)=constrict*radius
                 endif
               endif
             endif
             if(elem_cnct(-1,0,ne).eq.0)cycle      !No parent therefore inlet. Skip and go to the next
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "version number")>0) then
                read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             endif
             if(index(ctemp1, "value")> 0) then
                call get_final_real(ctemp1,radius)
                 if(radius_type.eq.'taper')then
                   elem_field(ne_radius_out,ne)=constrict*radius
                 else
                  elem_field(ne_radius,ne)=constrict*radius
                 endif
             endif
          elseif(surround.gt.1)then !Terminal airway - use first radius
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             read(unit=10, fmt="(a)", iostat=ierror) ctemp1
             if(index(ctemp1, "value")> 0) then
                call get_final_real(ctemp1,radius)
                  if(radius_type.eq.'taper')then
                    elem_field(ne_radius_out,ne)=constrict*radius
                  else
                    elem_field(ne_radius,ne)=constrict*radius
                  endif
             endif
          endif
       endif !index
       if(np.ge.num_nodes) exit read_a_node
    end do read_a_node

!If airway type calculate element volume
    ! calculate the element volumes
    do ne=1,num_elems
       if(radius_type.eq.'taper')then
         if(elem_cnct(-1,0,ne).ne.0)then !radius in is radius of upstream vessel
            elem_field(ne_radius_in,ne)=elem_field(ne_radius_out,elem_cnct(-1,1,ne))
         endif
         elem_field(ne_radius,ne)=(elem_field(ne_radius_in,ne)+elem_field(ne_radius_out,ne))/2
       endif
       !       ne_global=elems(noelem)
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       ! element volume
       elem_field(ne_vol,ne) = PI * elem_field(ne_radius,ne)**2 * &
            elem_field(ne_length,ne)
       elem_field(ne_a_A,ne) = 1.0_dp ! set default for ratio a/A
    enddo
    call enter_exit(sub_name,2)

  END subroutine define_rad_from_file
!
!##################################################################################
!
!*define_rad_from_geom:* Defines vessel or airway radius based on their geometric structure
  subroutine define_rad_from_geom(ORDER_SYSTEM, CONTROL_PARAM, START_FROM, START_RAD, group_type_in, group_option_in)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_DEFINE_RAD_FROM_GEOM" :: DEFINE_RAD_FROM_GEOM
    use arrays,only: dp,num_elems,elem_field,elem_ordrs,maxgen,elem_cnct
    use indices
    use diagnostics, only: enter_exit
    implicit none
   character(LEN=100), intent(in) :: ORDER_SYSTEM,START_FROM
   character(LEN=100), optional :: group_type_in, group_option_in
   real(dp), intent(in) :: CONTROL_PARAM,START_RAD
   !Input options ORDER_SYSTEM=STRAHLER (CONTROL_PARAM=RDS), HORSFIELD (CONTROL_PARAM=RDH)

   !Local variables
   character(LEN=100) :: group_type, group_options
   integer :: ne_min,ne_max,nindex,ne,n_max_ord,n,ne_start,&
      inlet_count
   real(dp) :: radius
   character(len=60) :: sub_name

   sub_name = 'define_rad_from_geom'
    call enter_exit(sub_name,1)
    !define list of elements you are going to operate on
    if(present(group_type_in))then
      group_type = group_type_in
    else!default to all
      group_type='all'
    endif
    if(group_type.eq.'all')then
       ne_min=1
       ne_max=num_elems
    elseif(group_type.eq.'efield')then

    endif

    !Define start element
    if(START_FROM.eq.'inlet')then
      inlet_count=0
      do ne=ne_min,ne_max
         if(elem_cnct(-1,0,ne).eq.0)then
           inlet_count=inlet_count+1
           ne_start=ne
         endif
         if(inlet_count.gt.1)then
            WRITE(*,*) ' More than one inlet in this group, using last found, ne = ',ne
         endif
      enddo
    else!element number defined
       read (START_FROM,'(I10)') ne_start
    endif

    !Strahler and Horsfield ordering system
    if(ORDER_SYSTEM(1:5).EQ.'strah')THEN
      nindex=no_sord !for Strahler ordering
    else if(ORDER_SYSTEM(1:5).eq.'horsf')then
      nindex = no_hord !for Horsfield ordering
    endif

    ne=ne_start
    n_max_ord=elem_ordrs(nindex,ne)
    elem_field(ne_radius,ne)=START_RAD

    do ne=1,num_elems
     radius=10.0_dp**(log10(CONTROL_PARAM)*dble(elem_ordrs(nindex,ne)-n_max_ord)&
        +log10(START_RAD))
     elem_field(ne_radius,ne)=radius
     elem_field(ne_radius_in,ne)=radius
     elem_field(ne_radius_out,ne)=radius
    enddo

    call enter_exit(sub_name,2)

  END subroutine define_rad_from_geom
!
!###########################################################################
!
!*element_connectivity_1d:*  Calculates element connectivity in 1D and stores in elelem_cnct
  subroutine element_connectivity_1d()
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ELEMENT_CONNECTIVITY_1D" :: ELEMENT_CONNECTIVITY_1D
    use arrays,only: elem_cnct,elem_nodes,elems_at_node,num_elems,num_nodes
    use diagnostics, only: enter_exit
    implicit none
    !     Local Variables
    integer :: ne,ne2,nn,noelem,np,np2,np1
    integer,parameter :: NNT=2
    character(len=60) :: sub_name

    sub_name = 'element_connectivity_1d'
    call enter_exit(sub_name,1)

    elem_cnct = 0 !initialise

    ! calculate elems_at_node array: stores the elements that nodes are in
    elems_at_node(1:num_nodes,0) = 0 !initialise number of adjacent elements

    DO ne=1,num_elems
       DO nn=1,2
          np=elem_nodes(nn,ne)
          elems_at_node(np,0)=elems_at_node(np,0)+1
          elems_at_node(np,elems_at_node(np,0))=ne ! local element that np is in
       ENDDO !nn
    ENDDO !noelem

    ! calculate elem_cnct array: stores the connectivity of all elements

    elem_cnct=0 !initialise all elem_cnct

    DO ne=1,num_elems
       !     ne_global=elems(noelem)
       IF(NNT == 2) THEN !1d
          np1=elem_nodes(1,ne) !first local node
          np2=elem_nodes(2,ne) !second local node
          DO noelem=1,elems_at_node(np2,0)
             ne2=elems_at_node(np2,noelem)
             IF(ne2 /= ne)THEN
                elem_cnct(-1,0,ne2)=elem_cnct(-1,0,ne2)+1
                elem_cnct(-1,elem_cnct(-1,0,ne2),ne2)=ne !previous element
                elem_cnct(1,0,ne)=elem_cnct(1,0,ne)+1
                elem_cnct(1,elem_cnct(1,0,ne),ne)=ne2
             ENDIF !ne2
          ENDDO !noelem2

       ENDIF
    ENDDO

    call enter_exit(sub_name,2)

  END subroutine element_connectivity_1d

!
!###################################################################################
!
!*evaluate_ordering:* calculates generations, Horsfield orders, Strahler orders for a given tree
  subroutine evaluate_ordering()
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_ORDERING" :: EVALUATE_ORDERING
    use arrays,only: elem_cnct,elem_nodes,elem_ordrs,elem_symmetry,&
         elems_at_node,num_elems,num_nodes,maxgen
    use diagnostics, only: enter_exit
    implicit none

    integer :: INLETS,ne,ne0,ne2,noelem2,np,np2, &
         num_attach,n_children,n_generation, &
         n_horsfield,OUTLETS,STRAHLER,STRAHLER_ADD,temp1
    LOGICAL :: DISCONNECT,DUPLICATE
    character(len=60) :: sub_name
    sub_name = 'evaluate_ordering'
    call enter_exit(sub_name,1)

    !Calculate generations, Horsfield orders, Strahler orders
    !.....Calculate branch generations

    maxgen=1
    DO ne=1,num_elems
       ne0=elem_cnct(-1,1,ne) !parent
       IF(ne0.NE.0)THEN
          n_generation=elem_ordrs(1,ne0) !parent generation
          IF(elem_cnct(1,0,ne0).EQ.1)THEN !single daughter
             elem_ordrs(1,ne)=n_generation + (elem_symmetry(ne)-1)
          ELSE IF(elem_cnct(1,0,ne0).GE.2)THEN
             elem_ordrs(1,ne)=n_generation+1
          ENDIF
       ELSE
          elem_ordrs(1,ne)=1 !generation 1
       ENDIF
       maxgen=max(maxgen,elem_ordrs(1,ne))
    ENDDO !noelem

    !.....Calculate the branch orders
    DO ne=num_elems,1,-1
       n_horsfield=MAX(elem_ordrs(2,ne),1)
       n_children=elem_cnct(1,0,ne) !number of child branches
       IF(n_children.EQ.1)THEN
          IF(elem_ordrs(1,elem_cnct(1,1,ne)).EQ.0)  n_children=0
       ENDIF
       STRAHLER=0
       STRAHLER_ADD=1
       IF(n_children.GE.2)THEN !branch has two or more daughters
          STRAHLER=elem_ordrs(3,elem_cnct(1,1,ne)) !first daughter
          DO noelem2=1,n_children !for all daughters
             ne2=elem_cnct(1,noelem2,ne) !global element # of daughter
             temp1=elem_ordrs(2,ne2) !Horsfield order of daughter
             IF(temp1.GT.n_horsfield) n_horsfield=temp1
             IF(elem_ordrs(3,ne2).LT.STRAHLER)THEN
                STRAHLER_ADD=0
             ELSE IF(elem_ordrs(3,ne2).GT.STRAHLER)THEN
                STRAHLER_ADD=0
                STRAHLER=elem_ordrs(3,ne2) !highest daughter
             ENDIF
          ENDDO !noelem2 (ne2)
          n_horsfield=n_horsfield+1 !Horsfield ordering
       ELSE IF(n_children.EQ.1)THEN
          ne2=elem_cnct(1,1,ne) !local element # of daughter
          n_horsfield=elem_ordrs(2,ne2)+(elem_symmetry(ne)-1)
          STRAHLER_ADD=elem_ordrs(3,ne2)+(elem_symmetry(ne)-1)
       ENDIF !elem_cnct
       elem_ordrs(2,ne)=n_horsfield !store the Horsfield order
       elem_ordrs(3,ne)=STRAHLER+STRAHLER_ADD !Strahler order
    ENDDO !noelem

    !       Check for disconnected nodes and number of inlets and outlets
    DUPLICATE=.FALSE.
    DO ne=1,num_elems
       np=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       IF(np.EQ.np2)THEN
          DUPLICATE=.TRUE.
       ENDIF
    ENDDO

    DISCONNECT=.FALSE.
    INLETS=0
    OUTLETS=0
    DO np=1,num_nodes
       num_attach=elems_at_node(np,0)
       IF(num_attach.EQ.0)THEN
          DISCONNECT=.TRUE.
       ELSEIF(num_attach.EQ.1)THEN
          ne=elems_at_node(np,1)
          IF(elem_cnct(1,0,ne).EQ.0) OUTLETS=OUTLETS+1
         IF(elem_cnct(-1,0,ne).EQ.0) INLETS=INLETS+1
       ELSEIF(num_attach.GT.3)THEN
          WRITE(*,*) ' Node ',np,' attached to',num_attach,' elements'
       ENDIF
    ENDDO

    call enter_exit(sub_name,2)

  end subroutine evaluate_ordering
!
!###################################################################################
!
!>*set_initial_volume:* assigns a volume to terminal units appended on a tree structure
!>based on an assumption of a linear gradient in the gravitational direction with max
!> min and COV values defined.
  subroutine set_initial_volume(Gdirn,COV,total_volume,Rmax,Rmin)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SET_INITIAL_VOLUME" :: SET_INITIAL_VOLUME

    use arrays,only: dp,elem_nodes,elem_units_below,&
         node_xyz,num_elems,num_units,units,unit_field
    use indices,only: nu_vol,nu_vt
    use diagnostics, only: enter_exit
    implicit none

    !     Parameter List
    integer,intent(in) :: Gdirn
    real(dp),intent(in) :: COV,total_volume,Rmax,Rmin
    !     Local parameters
    integer :: ne,np2,nunit
    real(dp) ::  factor_adjust,max_z,min_z,random_number,range_z,&
         volume_estimate,volume_of_tree,Vmax,Vmin,Xi
    character(len=60) :: sub_name

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

    range_z=DABS(max_z-min_z)
    if(DABS(range_z).le.1.0e-5_dp) range_z=1.0_dp

!!! for each elastic unit allocate a size based on a gradient in the Gdirn direction, and
!!! perturb by a user-defined COV. This should be calling a random number generator.
    do nunit=1,num_units
       ne=units(nunit)
       np2=elem_nodes(2,ne) !end node
       Xi=(node_xyz(Gdirn,np2)-min_z)/range_z
       random_number=random_number+0.1_dp
       IF(random_number.GT.1.0_dp) random_number=-1.1_dp
       unit_field(nu_vol,nunit)=(Vmax*Xi+Vmin*(1.0_dp-Xi))*(1.0_dp+COV*random_number)
       unit_field(nu_vt,nunit)=0.0_dp !initialise the tidal volume to a unit
    enddo !nunit

    ! correct unit volumes such that total volume is exactly as specified
    call volume_of_mesh(volume_estimate,volume_of_tree)
    factor_adjust = (total_volume-volume_of_tree)/(volume_estimate-volume_of_tree)
    do nunit=1,num_units
       unit_field(nu_vol,nunit) = unit_field(nu_vol,nunit)*factor_adjust
    enddo

    write(*,'('' Number of elements is'',I5)') num_elems
    write(*,'('' Initial volume is'',F6.2,'' L'')') total_volume/1.0e+6_dp
    write(*,'('' Deadspace volume is'',F6.1,'' mL'')') volume_of_tree/1.0e+3_dp

    call enter_exit(sub_name,2)

  end subroutine set_initial_volume

!
!###################################################################################
!
!*volume_of_mesh:* calculates the volume of an airway mesh including conducting and respiratory airways
  subroutine volume_of_mesh(volume_model,volume_tree)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_VOLUME_OF_MESH" :: VOLUME_OF_MESH
    use arrays,only: dp,elem_cnct,elem_field,elem_symmetry,&
         num_elems,num_units,units,unit_field
    use indices,only: ne_vol,nu_vol
    use diagnostics, only: enter_exit
    implicit none

    real(dp) :: volume_model,volume_tree
    !     Local Variables
    integer :: ne,ne0,nunit
    real(dp),allocatable :: vol_anat(:),vol_below(:)
    character(len=60) :: sub_name

    sub_name = 'volume_of_mesh'
    call enter_exit(sub_name,1)

    if(.not.allocated(vol_anat)) allocate(vol_anat(num_elems))
    if(.not.allocated(vol_below)) allocate(vol_below(num_elems))

    vol_anat = elem_field(ne_vol,1:num_elems) !initialise to branch volume
    vol_below = elem_field(ne_vol,1:num_elems) !initialise to branch volume

    do nunit=1,num_units
       ne=units(nunit)
       vol_below(ne) = vol_below(ne) + unit_field(nu_vol,nunit) !add elastic unit volume
    enddo !nunit

    do ne=num_elems,2,-1
       ne0=elem_cnct(-1,1,ne)
       vol_anat(ne0) = vol_anat(ne0) + DBLE(elem_symmetry(ne))*vol_anat(ne)
       vol_below(ne0) = vol_below(ne0) + DBLE(elem_symmetry(ne))*vol_below(ne)
    enddo !noelem

    volume_model = vol_below(1)
    volume_tree = vol_anat(1)

    deallocate(vol_anat)
    deallocate(vol_below)

    call enter_exit(sub_name,2)

  end subroutine volume_of_mesh

!
!#####################################################################
!
!*reallocate_node_elem_arrays:* Reallocates the size of arrays when modifying geometries
  subroutine reallocate_node_elem_arrays(num_elems_new,num_nodes_new)
    use arrays,only: dp,elems,elem_cnct,elem_direction,elem_field,&
         elem_ordrs,elem_nodes,&
         elem_symmetry,elem_units_below,elems_at_node,expansile,&
         nodes,node_field,node_xyz,num_elems,num_nodes
    use indices
    use diagnostics, only: enter_exit
    implicit none

!!! Parameters
    integer,intent(in) :: num_elems_new,num_nodes_new

!!! Local variables
    integer,allocatable :: nodelem_temp(:),enodes_temp(:,:),enodes_temp2(:,:,:)
    real(dp),allocatable :: xyz_temp(:,:),rnodes_temp(:,:)
    logical,allocatable :: exp_temp(:)
    character(len=60) :: sub_name

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

    allocate(rnodes_temp(num_ne,num_elems))
    rnodes_temp=elem_field
    deallocate(elem_field)
    allocate(elem_field(num_ne,num_elems_new))
    elem_field(1:num_ne,1:num_elems)=rnodes_temp(1:num_ne,1:num_elems)
    deallocate(rnodes_temp)
    elem_field(1:num_ne,num_elems+1:num_elems_new) = 0.0_dp

    allocate(rnodes_temp(3,num_elems))
    rnodes_temp=elem_direction
    deallocate(elem_direction)
    allocate(elem_direction(3,num_elems_new))
    elem_direction(1:3,1:num_elems)=rnodes_temp(1:3,1:num_elems)
    deallocate(rnodes_temp)
    elem_direction(1:3,num_elems+1:num_elems_new) = 0.0_dp

    allocate(rnodes_temp(num_nj,num_nodes))
    rnodes_temp=node_field
    deallocate(node_field)
    allocate(node_field(num_nj,num_nodes_new))
    node_field(1:num_nj,1:num_nodes)=rnodes_temp(1:num_nj,1:num_nodes)
    deallocate(rnodes_temp)
    node_field(1:num_nj,num_nodes+1:num_nodes_new)=0.0_dp

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

    allocate(nodelem_temp(num_elems))
    nodelem_temp=elem_units_below
    deallocate(elem_units_below)
    allocate(elem_units_below(num_elems_new))
    elem_units_below(1:num_elems)=nodelem_temp(1:num_elems)
    deallocate(nodelem_temp)
    elem_units_below(num_elems+1:num_elems_new)=0

    allocate(enodes_temp(num_nodes,0:3))
    enodes_temp=elems_at_node
    deallocate(elems_at_node)
    allocate(elems_at_node(num_nodes_new,0:3))
    elems_at_node(1:num_nodes,0:3)=enodes_temp(1:num_nodes,0:3)
    deallocate(enodes_temp)
    elems_at_node(num_nodes+1:num_nodes_new,0:3)=0

    allocate(exp_temp(num_elems))
    exp_temp = expansile
    deallocate(expansile)
    allocate(expansile(num_elems_new))
    expansile(1:num_elems)=exp_temp(1:num_elems)
    deallocate(exp_temp)
    expansile(num_elems+1:num_elems_new)=.false.

    call enter_exit(sub_name,2)

  end subroutine reallocate_node_elem_arrays

!
!###################################################################################
!
!>get_final_real
  subroutine get_final_real(string,rtemp)
    use arrays,only: dp
    implicit none
    character, intent(in) :: string*(132)
    integer :: ibeg,iend
    real(dp), intent(out) :: rtemp
    real(dp) :: rsign
    character :: sub_string*(40)

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

  end subroutine get_final_real

!
!###################################################################################
!
!*get_final_string*
  subroutine get_final_string(string,rtemp)
    use arrays,only: dp
    implicit none
    character, intent(in) :: string*(132)
    integer :: ibeg,iend
    real(dp), intent(out) :: rtemp
    real(dp) :: rsign
    character :: sub_string*(40)

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

  end subroutine get_final_string

!
!###################################################################################
!
!*get_local_node*
  subroutine get_local_node(np_global,np_local)
    use arrays,only: nodes,num_nodes
    implicit none

    integer,intent(in) :: np_global
    integer,intent(out) :: np_local

    integer :: np
    logical :: found

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
    return

  end subroutine get_local_node

!
!###################################################################################
!
!*get_final_integer*
  subroutine get_final_integer(string,num)
    implicit none
    character,intent(in) :: string*(132)
    integer,intent(out) :: num
    integer :: ibeg,iend,nsign,ntemp
    character :: sub_string*(40)

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
    read (sub_string(ibeg:iend), '(i10)' ) ntemp !get integer values
    ntemp=ntemp*nsign !apply sign to number

    num=ntemp !return the integer value

  end subroutine get_final_integer

!
! ##########################################################################      
!

  function inlist(item,ilist)
!!! dummy arguments
    integer :: item,ilist(:)
! local variables
    integer :: n
    logical :: inlist

    inlist = .false.
    do n=1,size(ilist)
       if(item == ilist(n)) inlist = .true.
    enddo

  end function inlist

!
!###########################################################################################
!
end module geometry

