module growtree
  !
  !*Brief Description:* This module generates a volume-filling tree within a bounding surface.
  !
  !*LICENSE:*
  !
  !
  !
  !*Full Description:*
  ! A volume-filling tree is generated using a recursive algorithm that creates branches from the
  ! ends of an initial tree of 1D elements, towards the centre of mass of sets of 'seed points'
  ! that fill a bounding surface (referred to as the 'host').  Seeds are grouped with the closest
  ! candidate parent branch, then each set of seeds is split in two by a plane that is defined by
  ! the parent branch and the centre of mass of the seeds. Two new branches are created that 
  ! originate at the end of the parent branch, and are directed towards the centres of mass of the
  ! two subsets of seed points.
  
  ! Based on Tawhai et al. (2000) 'Generation of an anatomically based three-dimensional model 
  ! of the conducting airways', Ann. Biomed. Eng. 28(7): 793-802 and its modifications in 
  ! Tawhai et al. (2004) 'CT-based geometry analysis and finite element models of the human and 
  ! ovine bronchial tree' J. Appl. Physiol. 97(6): 2310-21.
  
  use arrays
  !use diagnostics, only: enter_exit      !! for diagnostics at enter/exit of subroutines
  use geometry
  use indices
  use other_consts   !! pi 
  use mesh_utilities   !! general functions for geometric/mesh calculations
  use math_utilities   !! general utility functions for sorting etc
  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private
  public grow_tree,smooth_1d_tree
  
contains
  
  !###############################################################
  !  
  !*adjust_branch_angle:* Adjusts child branch angle so not larger than a max value, 
  ! keeping child in-plane. Checks the branching angle between a parent and child branch. 
  ! If the branch angle is greater than the branch angle limit, then the branch angle is
  ! reduced to the limit value, such that the daughter branch remains in the original
  ! branching plane.
  ! 
  subroutine adjust_branch_angle(Nth,ne,np1,np2,np,angle_max,angle_min)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ADJUST_BRANCH_ANGLE" :: ADJUST_BRANCH_ANGLE
    use diagnostics, only: enter_exit
    implicit none
    
    integer :: Nth,ne,np1,np2,np
    real(dp) :: angle_max,angle_min

    !Local variables
    real(dp) :: a_min,a_lim,angle,angle_sibling,length,LU,LV,U(3),V(3),W(3)
    real(dp),parameter :: loose_tol = 1.0e-4_dp
    
    character(len=60) :: sub_name

    sub_name = 'adjust_branch_angle'
    call enter_exit(sub_name,1)

    length = elem_field(ne_length,ne)

    a_lim = angle_max*PI/180.0_dp ! maximum branching angle, in radians
    a_min = angle_min*PI/180.0_dp ! minimum branching angle, in radians
    
    U(1:3) = node_xyz(1:3,np1)-node_xyz(1:3,np2) ! direction of parent
    V(1:3) = node_xyz(1:3,np)-node_xyz(1:3,np1)  ! direction of this branch
    W(1:3) = V(1:3) ! will store direction if no angle change
    LU = vector_length(U)
    LV = vector_length(V)
    angle = angle_btwn_vectors(U,V)
    U = unit_vector(U)
    V = unit_vector(V)
    W = unit_vector(W)

    if(abs(angle).gt.a_lim)then !reduce angle
       if(abs(dabs(angle)-pi).lt.loose_tol)then !reduce angle from 180 degrees
          if(Nth.eq.2)then
             V(1:3) = elem_direction(1:3,ne-1) !direction of sibling
          else if(Nth.eq.1)then
             V(1:3) = elem_direction(1:3,ne+1) !direction of sibling
          endif
          angle_sibling = angle_btwn_vectors(U,V)
          W = vector_for_angle_limit(U,V,a_lim,angle_sibling+a_lim)
       else
          W = vector_for_angle_limit(U,V,a_lim,angle-a_lim)
       endif
       node_xyz(1:3,np) = node_xyz(1:3,np1)+W(1:3)*0.5_dp*length
       elem_field(ne_length,ne) = 0.5_dp*length
       elem_direction(1:3,ne) = W(1:3)

    elseif(abs(angle).lt.a_min)then
       if(abs(angle).lt.loose_tol)then !increase angle from 0 degrees
          if(Nth.eq.2)then
             V(1:3) = elem_direction(1:3,ne-1) !direction of sibling
          else if(Nth.eq.1)then
             V(1:3) = elem_direction(1:3,ne+1) !direction of sibling
          endif
          angle_sibling = angle_btwn_vectors(U,V) !angle between branch and sibling
          W = vector_for_angle_limit(U,V,a_min,angle_sibling+a_min)
       else
          W = vector_for_angle_limit(U,V,a_min,angle-a_min)
       endif
       node_xyz(1:3,np) = node_xyz(1:3,np1)+W(1:3)*0.5_dp*length
       elem_field(ne_length,ne) = 0.5_dp*length
       elem_direction(1:3,ne) = W(1:3)
       
    endif !abs(angle).gt.a_lim
    
    call enter_exit(sub_name,2)

  end subroutine adjust_branch_angle
  
  
  !###############################################################
  !
  !*branch_to_cofm*: Creates a new branch towards the cofm of a set of seed points. 
  ! Used in a volume-filling branching method to create a branch that runs from a
  ! defined point to some fraction along a line towards the centre of mass of a 
  ! collection of seed points.
  !
  subroutine branch_to_cofm(map_seed_to_elem,nen,np1,COFM,branch_fraction,length_limit,&
    length_parent,shortest_length,candidate_xyz,make_branch)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_BRANCH_TO_COFM" :: BRANCH_TO_COFM

    use diagnostics, only: enter_exit
    implicit none
    
    integer :: map_seed_to_elem(*),nen,np1
    real(dp) :: COFM(3),branch_fraction,length_limit,length_parent,&
         shortest_length,candidate_xyz(3)
    logical :: make_branch

    !Local variables
    integer :: N,NCLOSEST(100),nd,nsp,NUM_CLOSEST,NUM_ND,number_of_points
    real(dp) :: CLOSEST(100),DIST,L_COFM,LENGTH,MIN_DIST,VECTOR(3)
    
    character(len=60) :: sub_name

    sub_name = 'branch_to_cofm'
    call enter_exit(sub_name,1)

    number_of_points = 0
    do nd = 1,num_data
       if(map_seed_to_elem(nd).eq.nen) number_of_points=number_of_points+1
    enddo
    candidate_xyz(1:3) = node_xyz(1:3,np1) + branch_fraction*(COFM(1:3)-node_xyz(1:3,np1))
    VECTOR(1:3) = COFM(1:3)-node_xyz(1:3,np1)
    LENGTH = distance_between_points(candidate_xyz,node_xyz(1,np1))
    L_COFM = vector_length(VECTOR)
    VECTOR = unit_vector(VECTOR)
    
    if(LENGTH.ge.LENGTH_LIMIT)then !the branch will not be terminal
       make_branch = .TRUE.
    else
       if(elem_ordrs(no_gen,nen).lt.12)then
          make_branch = .true.
          candidate_xyz(1:3) = node_xyz(1:3,np1)+VECTOR(1:3)*0.33_dp*length_parent
       else
          make_branch = .FALSE.
       endif
    endif !LENGTH.ge.L_LIM

    !*** Position branch
    MIN_DIST = 1.0e+6_dp
    NUM_ND = 0
    NUM_CLOSEST = 1
    NCLOSEST(1) = 1
    CLOSEST(1) = 1.0e+6_dp
    do nd = 1,num_data
       nsp = map_seed_to_elem(nd) !space # that random point belongs to
       if(nsp.eq.nen)then !random point belongs to this element space
          dist = distance_between_points(data_xyz(1,nd),candidate_xyz)
          if(DIST.lt.MIN_DIST) MIN_DIST = DIST
          NUM_ND = NUM_ND+1
          if(.NOT.make_branch)then !remove closest data points
             if(DIST.lt.CLOSEST(NUM_CLOSEST))then !store this data point
                if(NUM_CLOSEST.lt.1)then 
                   NUM_CLOSEST = NUM_CLOSEST+1 !increment number of closest
                endif
                CLOSEST(NUM_CLOSEST) = DIST !store distance
                NCLOSEST(NUM_CLOSEST) = nd !store data point number
             endif !DIST
             call sort_real_list(NUM_CLOSEST,CLOSEST,NCLOSEST) !sort into ascending
          endif !NOT.make_branch
       endif !nsp
    enddo !nd
    
    if(LENGTH.lt.0.5_dp*length_parent)then
       candidate_xyz(1:3) = node_xyz(1:3,np1)+VECTOR(1:3)*0.5_dp*length_parent
    endif
    
    if(LENGTH.lt.shortest_length)then
       candidate_xyz(1:3) = node_xyz(1:3,np1)+VECTOR(1:3)*shortest_length
    endif
    
    if(.NOT.make_branch)then !remove the closest data points
       do N = 1,NUM_CLOSEST
          nd = NCLOSEST(N)
          map_seed_to_elem(nd) = 0
       enddo
    endif !NOT.make_branch
    
    call enter_exit(sub_name,2)

  end subroutine branch_to_cofm
  

  !###############################################################
  !
  !*calculate_seed_cofm:* calculates centre of mass of a list of seed
  ! points by averaging their coordinates and returns result in 'cofm'
  !
  subroutine calculate_seed_cofm(map_seed_to_elem,nen,COFM)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CALCULATE_SEED_COFM" :: CALCULATE_SEED_COFM

    use diagnostics, only: enter_exit
    implicit none
    
    integer :: map_seed_to_elem(*),nen
    real(dp) :: COFM(3)
    
    !Local variables
    integer :: DAT,nd,nsp
    
    character(len=60) :: sub_name

    sub_name = 'calculate_seed_cofm'
    call enter_exit(sub_name,1)
    
    DAT = 0
    cofm = 0.0_dp
    
    do nd=1,num_data
       nsp=map_seed_to_elem(nd) !the space # for the nd-th data point
       if(nsp.eq.nen)then
          DAT=DAT+1
          COFM(1:3)=COFM(1:3)+data_xyz(1:3,nd)
       endif
    enddo !nd
    if(DAT.ne.0) COFM(1:3) = COFM(1:3)/DAT !centre of mass

    call enter_exit(sub_name,2)
    
  end subroutine calculate_seed_cofm
 
 
  !###############################################################
  !  
  !*check_branch_rotation_plane:* limits the angle between branching planes
  ! to a maximum (user-defined) value, and makes sure branch remains internal
  ! to the host volume
  !
  subroutine check_branch_rotation_plane(map_seed_to_elem,ne,&
       ne_grnd_parent,ne_parent,local_parent_temp,num_next_parents,&
       np,np1,np2,np3,num_terminal,rotation_limit)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CHECK_BRANCH_ROTATION_PLANE" :: CHECK_BRANCH_ROTATION_PLANE

    use diagnostics, only: enter_exit
    implicit none
    
    integer :: map_seed_to_elem(*),ne,ne_grnd_parent,ne_parent, &
         num_next_parents,np,np1,np2,np3,num_terminal
    ! np == end node; np1 == np_start; np2 == np_prnt_start; np3 == np_grnd_start
    integer :: local_parent_temp(*)
    real(dp),intent(in) :: rotation_limit
    
    !Local variables
    integer :: COUNT,nd_min,ne_other,nes,np4,offset
    double precision :: ANGLE,length,ROT_ANGLE,candidate_xyz(3)
    logical :: INTERNAL
    
    character(len=60) :: sub_name

    sub_name = 'check_branch_rotation_plane'
    call enter_exit(sub_name,1)

    ! find the appropriate other point for calculating the branching plane
    if(elem_cnct(1,1,ne_grnd_parent).eq.ne_parent)then
       if(elem_cnct(1,2,ne_grnd_parent).eq.0)then
          ne_other = elem_cnct(-1,1,ne_grnd_parent)
          np4 = elem_nodes(1,ne_other)
       else
          ne_other = elem_cnct(1,2,ne_grnd_parent)
          np4 = elem_nodes(2,ne_other)
       endif
    else
       ne_other = elem_cnct(1,1,ne_grnd_parent)
       np4 = elem_nodes(2,ne_other)
    endif
    
    ROT_ANGLE=ROTATION_LIMIT*PI/180.0_dp
    call CHECK_ROTATION_ANGLE(ne,np4,np3,np2,np1,np,np-1,rotation_limit)
    angle = rotation_angle(np2,np1,np4,np,np-1)
    INTERNAL=.FALSE.
    COUNT=0
    do while(.NOT.INTERNAL.and.COUNT.lt.2)
       candidate_xyz(1:3) = node_xyz(1:3,np-1)
       !                  call CHECK_POINT_INTERNAL(IBT,Ido,INP,NBJ,NDLIST,NEELEM, &
       !                       NEP,NHOST,NKJE,np-1,np1,NPF,NPNE,NVJE,elem_cnct,SE, &
       !                       XA,XE,XIP,XP,candidate_xyz,data_xyz,INTERNAL,.FALSE.)
       internal=.true.
       if(.not.internal)then ! halve the length, and halve the angle from parent
          length = distance_between_points(node_xyz(1,np1),node_xyz(1,np-1))
!          node_xyz(1:3,np-1) = node_xyz(1:3,np1) + 0.5_dp*length*elem_direction(1:3,np-1)
          node_xyz(1:3,np-1) = node_xyz(1:3,np1) + 0.5_dp*length*elem_direction(1:3,ne-1)
          candidate_xyz(1:3) = node_xyz(1:3,np-1)
          call reduce_branch_angle(np1,np2,np-1,candidate_xyz,0.5_dp) ! reduces the branch angle by a half
       endif !internal
       COUNT=COUNT+1
    enddo
    if(.NOT.INTERNAL)then
       !...............Remove the branch from the list of next generation parents
       offset=0
       do nes=1,num_next_parents
          if(local_parent_temp(nes).eq.ne-1) offset=1
          local_parent_temp(nes)=local_parent_temp(nes+offset)
       enddo
       num_next_parents=num_next_parents-1
       num_terminal=num_terminal+1
       
       !...............Remove the closest data point to the end of the branch
       nd_min = closest_seed_to_node(map_seed_to_elem,np-1)
       
       map_seed_to_elem(nd_min)=0
       
    endif
    
    INTERNAL=.FALSE.
    COUNT=0
    do while(.NOT.INTERNAL.and.COUNT.lt.2)
       candidate_xyz(1:3) = node_xyz(1:3,np)
       !                  call CHECK_POINT_INTERNAL(IBT,Ido,INP,NBJ,NDLIST,NEELEM, &
       !                       NEP,NHOST,NKJE,np,np1,NPF,NPNE,NVJE,elem_cnct,SE,XA, &
       !                       XE,XIP,XP,candidate_xyz,data_xyz,INTERNAL,.FALSE.)
       internal=.true.
       if(.not.internal)then
          length = distance_between_points(node_xyz(1,np1),node_xyz(1,np))
!          node_xyz(1:3,np) = node_xyz(1:3,np1)+0.5_dp*length*elem_direction(1:3,np)
          node_xyz(1:3,np) = node_xyz(1:3,np1)+0.5_dp*length*elem_direction(1:3,ne)
          call reduce_branch_angle(np1,np2,np,candidate_xyz,0.5_dp) ! reduces the branch angle by a half
       endif
       COUNT=COUNT+1
    enddo !while
    if(.NOT.INTERNAL)then
       !...............Remove the branch from the list of next generation parents
       offset=0
       do nes=1,num_next_parents
          if(local_parent_temp(nes).eq.ne) offset=1
          local_parent_temp(nes)=local_parent_temp(nes+offset)
       enddo
       num_next_parents=num_next_parents-1
       num_terminal=num_terminal+1
       !...............Remove the closest data point to the end of the branch
       nd_min = closest_seed_to_node(map_seed_to_elem,np-1)
       
       map_seed_to_elem(nd_min)=0
       
    endif
    
    call enter_exit(sub_name,2)

  end subroutine check_branch_rotation_plane

  
  !##################################################
  !
  !*check_rotation_angle:* adjusts the branch locations such that the angle between
  ! branching planes is less than a user-defined maximum. 
  ! Calculates using quaternions. For angle ROTATION_ANGLE and unit
  ! vector a,b,c , calculate rotation matrix for arbitrary point.
  !
  subroutine check_rotation_angle(ne,np00,np0,np1,np2,np3,np4,rotation_limit)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CHECK_ROTATION_ANGLE" :: CHECK_ROTATION_ANGLE
    
    use diagnostics, only: enter_exit
    implicit none

    integer,intent(in) :: ne,np00,np0,np1,np2,np3,np4
    real(dp),intent(in) :: rotation_limit

    ! Local variables
    integer :: IT,ITMAX
    real(dp) :: ANGLE0,ANGLE,AXIS(3),DIRECTION(3),NRML(3), &
         NRML_PARENT(3),ROT_ANGLE,U(3),V(3),Q0,Q1,Q2,Q3,Q(3,3),X(3), &
         ANGLE_BETWEEN,length
    logical :: COMPLETE
    
    character(len=60) :: sub_name

    sub_name = 'check_rotation_angle'
    call enter_exit(sub_name,1)

    ITmAX=10
    
    !...parent branching plane, cross-product of parent and grand-parent
    U(1:3)=node_xyz(1:3,np2)-node_xyz(1:3,np1) !parent
    V(1:3)=node_xyz(1:3,np1)-node_xyz(1:3,np0) !grandparent
    U = unit_vector(U)
    V = unit_vector(U)
    NRML_PARENT = cross_product(U,V) !calculate branching plane
    NRML_PARENT = unit_vector(NRML_PARENT)
    
    !...current branching plane
    U(1:3)=node_xyz(1:3,np3)-node_xyz(1:3,np2) !branch
    V(1:3)=node_xyz(1:3,np4)-node_xyz(1:3,np2) !sibling
    U = unit_vector(U)
    V = unit_vector(U)
    NRML = cross_product(U,V) !calculate branching plane
    NRML = unit_vector(NRML)
    
    !...angle between branching planes
    ANGLE=angle_btwn_vectors(NRML,NRML_PARENT)
    ANGLE_BETWEEN=ANGLE
    if(ANGLE_BETWEEN.gt.PI/2.0_dp)then
       ANGLE_BETWEEN=ANGLE_BETWEEN-PI
       ANGLE=ANGLE_BETWEEN
    endif
    !      ANGLE=PI/2.0_dp-ANGLE
    ANGLE=PI/2.0_dp+ANGLE
    
    if(abs(ANGLE_BETWEEN).gt.ROTATION_LIMIT.and.abs(ANGLE_BETWEEN) &
         .lt.PI/2.0_dp-ROTATION_LIMIT)then
       if(ANGLE.lt.0.0_dp)then
          ROT_ANGLE=-(ANGLE+ROTATION_LIMIT)
       else
          ROT_ANGLE=-(ANGLE-ROTATION_LIMIT)
       endif
       
       ANGLE0=ANGLE
       
       !...if the difference in angles is not correct, rotate branches
       AXIS(1:3)=node_xyz(1:3,np2)-node_xyz(1:3,np1)
       AXIS = unit_vector(AXIS)
       
       Q0=DCOS(ROT_ANGLE/2.0_dp)
       Q1=DSIN(ROT_ANGLE/2.0_dp)*AXIS(1)
       Q2=DSIN(ROT_ANGLE/2.0_dp)*AXIS(2)
       Q3=DSIN(ROT_ANGLE/2.0_dp)*AXIS(3)
       
       Q(1,1)=Q0**2+Q1**2-Q2**2-Q3**2
       Q(1,2)=2*(Q1*Q2-Q0*Q3)
       Q(1,3)=2*(Q1*Q3+Q0*Q2)
       Q(2,1)=2*(Q2*Q1+Q0*Q3)
       Q(2,2)=Q0**2-Q1**2+Q2**2-Q3**2
       Q(2,3)=2*(Q2*Q3-Q0*Q1)
       Q(3,1)=2*(Q3*Q1-Q0*Q2)
       Q(3,2)=2*(Q3*Q2+Q0*Q1)
       Q(3,3)=Q0**2-Q1**2-Q2**2+Q3**2
       
!       X(1:3) = elem_direction(1:3,np3) ! unit vector
       X(1:3) = elem_direction(1:3,ne) ! unit vector
       length = distance_between_points(node_xyz(1,np3),node_xyz(1,np2))
       
       node_xyz(1:3,np3)=node_xyz(1:3,np2)+length*(Q(1:3,1)*X(1)+Q(1:3,2) &
            *X(2)+Q(1:3,3)*X(3))
       DIRECTION(1:3)=(node_xyz(1:3,np3)-node_xyz(1:3,np2))
       DIRECTION = unit_vector(DIRECTION)
!       elem_direction(1:3,np3)=DIRECTION(1:3)
       elem_direction(1:3,ne)=DIRECTION(1:3)

!       X(1:3)=elem_direction(1:3,np4) !unit vector
       X(1:3)=elem_direction(1:3,ne-1) !unit vector
       length = distance_between_points(node_xyz(1,np3),node_xyz(1,np2))

       node_xyz(1:3,np4)=node_xyz(1:3,np2)+length*(Q(1:3,1)*X(1)+Q(1:3,2) &
            *X(2)+Q(1:3,3)*X(3))
       DIRECTION(1:3)=(node_xyz(1:3,np4)-node_xyz(1:3,np2))
       DIRECTION = unit_vector(DIRECTION)
       elem_direction(1:3,ne-1)=DIRECTION(1:3)
       
       U(1:3)=elem_direction(1:3,ne) !direction of a branch
       V(1:3)=elem_direction(1:3,ne-1) !direction of its sibling
       U = unit_vector(U)
       V = unit_vector(V)
       NRML = cross_product(U,V) !calculate branching plane
       NRML = unit_vector(NRML)
       
       !...angle between branching planes
       ANGLE=angle_btwn_vectors(NRML,NRML_PARENT)
       
       ! should find that 90degrees minus new angle is within the limit range
       if(abs(ANGLE_BETWEEN).gt.ROTATION_LIMIT.and.abs(ANGLE_BETWEEN) &
            .lt.PI/2.0_dp-ROTATION_LIMIT)then
          COMPLETE=.TRUE.
       else
          COMPLETE=.FALSE.
       endif
       do while(.NOT.COMPLETE)
          IT=IT+1
          ANGLE_BETWEEN=ANGLE
          ANGLE=PI/2.0_dp-ANGLE
          if(ANGLE.lt.0.0_dp)then
             ROT_ANGLE=-(ANGLE+ROTATION_LIMIT)
          else
             ROT_ANGLE=-(ANGLE-ROTATION_LIMIT)
          endif
          
          Q0=DCOS(ROT_ANGLE/2.0_dp)
          Q1=DSIN(ROT_ANGLE/2.0_dp)*AXIS(1)
          Q2=DSIN(ROT_ANGLE/2.0_dp)*AXIS(2)
          Q3=DSIN(ROT_ANGLE/2.0_dp)*AXIS(3)
          
          Q(1,1)=Q0**2+Q1**2-Q2**2-Q3**2
          Q(1,2)=2*(Q1*Q2-Q0*Q3)
          Q(1,3)=2*(Q1*Q3+Q0*Q2)
          Q(2,1)=2*(Q2*Q1+Q0*Q3)
          Q(2,2)=Q0**2-Q1**2+Q2**2-Q3**2
          Q(2,3)=2*(Q2*Q3-Q0*Q1)
          Q(3,1)=2*(Q3*Q1-Q0*Q2)
          Q(3,2)=2*(Q3*Q2+Q0*Q1)
          Q(3,3)=Q0**2-Q1**2-Q2**2+Q3**2
          
!          X(1:3)=elem_direction(1:3,np3) !unit vector
          X(1:3)=elem_direction(1:3,ne) !unit vector
          length = distance_between_points(node_xyz(1,np3),node_xyz(1,np2))
          
          node_xyz(1:3,np3)=node_xyz(1:3,np2)+length*(Q(1:3,1)*X(1)+Q(1:3,2) &
               *X(2)+Q(1:3,3)*X(3))
          DIRECTION(1:3)=(node_xyz(1:3,np3)-node_xyz(1:3,np2))
          DIRECTION = unit_vector(DIRECTION)
          elem_direction(1:3,ne)=DIRECTION(1:3)

          X(1:3)=elem_direction(1:3,ne-1) !unit vector
          length = distance_between_points(node_xyz(1,np3),node_xyz(1,np2))
          node_xyz(1:3,np4)=node_xyz(1:3,np2)+length*(Q(1:3,1)*X(1)+Q(1:3,2) &
               *X(2)+Q(1:3,3)*X(3))
          DIRECTION(1:3)=(node_xyz(1:3,np4)-node_xyz(1:3,np2))
          DIRECTION = unit_vector(DIRECTION)
          elem_direction(1:3,ne-1)=DIRECTION(1:3)
          
          U(1:3)=elem_direction(1:3,ne) !direction of a branch
          V(1:3)=elem_direction(1:3,ne-1) !direction of its sibling
          U = unit_vector(U)
          V = unit_vector(V)
          NRML = cross_product(U,V) !calculate branching plane
          NRML = unit_vector(NRML)
          
          !...angle between branching planes
          ANGLE=angle_btwn_vectors(NRML,NRML_PARENT)
          
          if(abs(ROTATION_LIMIT-abs(PI/2.0_dp-ANGLE)).LE.0.001_dp)then
             COMPLETE=.TRUE.
          endif
          
          if(IT.gt.ITMAX)then
             WRITE(*,*) 'WARNING!!!! rotation angle = ',ANGLE*180.0_dp/PI
          endif
          
       enddo !do while not found
       
       !.......Alternate calculation for the rotation angle
       angle = rotation_angle(np1,np2,np00,np3,np4)
    else
       !        write(*,*) 'Not',ANGLE_BETWEEN*180.0_dp/PI
    endif
    
    call enter_exit(sub_name,2)

  end subroutine check_rotation_angle

  
  !###############################################################
  !
  !*create_new_node:* sets up arrays for a new mesh node and element.
  !
  subroutine create_new_node(ne,ne_start,np,np_start,MAKE)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CREATE_NEW_NODE" :: CREATE_NEW_NODE

    use diagnostics, only: enter_exit
    implicit none
    
    integer :: ne,ne_start,np,np_start
    logical :: MAKE

    !Local variables    
    character(len=60) :: sub_name

    sub_name = 'create_new_node'
    call enter_exit(sub_name,1)
    
    if(MAKE)then
       ne=ne+1
       elems(ne) = ne ! store global element number
       elem_nodes(1,ne) = np_start
       elems_at_node(np_start,0)=elems_at_node(np_start,0)+1
       elems_at_node(np_start,elems_at_node(np_start,0))=ne
       
       np = np+1
       nodes(np) = np
       elems_at_node(np,0) = 0 !initialise
       elem_nodes(2,ne) = np !end node of new element
       elems_at_node(np,0) = elems_at_node(np,0)+1
       elems_at_node(np,elems_at_node(np,0)) = ne
       
       elem_cnct(1,0,ne)=0 !initialise number of proximal branches
       if(ne_start.ne.0)then
          elem_cnct(-1,0,ne)=1
          elem_cnct(-1,elem_cnct(-1,0,ne),ne)=ne_start
          elem_cnct(1,0,ne_start)=elem_cnct(1,0,ne_start)+1
          elem_cnct(1,elem_cnct(1,0,ne_start),ne_start)=ne
       endif
    endif
    
    call enter_exit(sub_name,2)

  end subroutine create_new_node
  

  !###############################################################
  !
  !*group_seeds_with_branch:* groups a set of seed points with the 
  ! closest candidate parent branches. reassigns data (seed) points
  ! to the closest ending of branches in the current generation.
  !
  subroutine group_seeds_with_branch(map_array,num_next_parents,num_seeds_from_elem, &
       num_terminal,local_parent,DISTANCE_LIMIT,FIRST)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GROUP_SEEDS_WITH_BRANCH" :: GROUP_SEEDS_WITH_BRANCH

    use diagnostics, only: enter_exit
    implicit none
    
    integer :: num_next_parents,local_parent(:),map_array(:),num_seeds_from_elem(*),&
         num_terminal
    real(dp),intent(in) :: DISTANCE_LIMIT
    logical :: FIRST
    
    !Local variables
    integer :: i,n,m,nd,nd_min,ne,n_elm_temp,ne_min,noelem,np,np_temp
    integer :: size_map
    integer,allocatable :: map_array_copy(:),my_closest(:)
    real(dp) :: dist,min_dist

    character(len=60) :: sub_name

    sub_name = 'group_seeds_with_branch'
    call enter_exit(sub_name,1)

    size_map = size(map_array)
    allocate(my_closest(size_map))
    
    if(first)then
!!! for the first seed point allocation, need to make sure that every terminal branch has 
!!! at least two seed points. find the closest two data points to each terminal node, and
!!! assign to terminals. the 'my_closest' array is used to indicate that the points have 
!!! already been allocated (when calculating the distribution of points below). Also need
!!! to ensure that the two closest seed points are unique (assigned to only one terminal).
!       n_closest = 0
!       do i=1,2
!          do N=1,num_next_parents
!             ne_min=local_parent(N)
!             np_temp=elem_nodes(2,ne_min)
!             min_dist = 1.0e+10_dp
!             nd_min = 0
!             do nd = 1,num_data
!                if(map_array(nd).ne.0)then
!                   if(.not.inlist(nd,my_closest))then
!                      dist = distance_between_points(data_xyz(1,nd),node_xyz(1,np_temp))
!                      if(dist.lt.min_dist)then
!                         nd_min = nd
!                         min_dist = dist
!                      endif !dist
!                   endif
!                endif
!             enddo
!             n_closest = n_closest + 1
!             my_closest(n_closest) = nd_min
!             map_array(nd_min) = ne_min
!          enddo
!       enddo !i

    else !use the data groupings from previous
       allocate(map_array_copy(size_map))
       map_array_copy(1:size_map) = map_array(1:size_map)
       do n=1,num_next_parents
          ne_min = local_parent(n)
          np_temp = elem_nodes(2,ne_min)
          MIN_DIST=1.0e+10_dp
          do nd=1,num_data
             if(map_array(nd).eq.ne_min)then ! was associated with this element
                dist = distance_between_points(data_xyz(1,nd),node_xyz(1,np_temp))
                if(dist.lt.min_dist)then
                   nd_min = nd
                   min_dist = dist
                endif !DIST
             endif
          enddo
          my_closest(N) = nd_min
          map_array(nd_min) = ne_min
       enddo
    endif

    do nd = 1,num_data            ! for all seed/data points
       if(map_array(nd).ne.0)then ! the data point is still in use
          if(.not.inlist(nd,my_closest))then
             MIN_DIST=1.0e+10_dp     ! initialise the minimum (closest) distance
             do noelem = 1,num_next_parents ! for each parent in the next branch generation
                ne = local_parent(noelem)
                np = elem_nodes(2,ne)
                dist = distance_between_points(data_xyz(1,nd),node_xyz(1,np))
                if(DIST.lt.MIN_DIST)then
                   ne_min = ne
                   MIN_DIST=DIST
                endif
             enddo
             if(first)then
                map_array(nd)=ne_min
             else
                if(MIN_DIST.lt.DISTANCE_LIMIT)then !keep seed points
                   map_array(nd)=ne_min
                else
                   map_array(nd)=0 !too far from branch ends, so discard
                endif
             endif
          endif
       endif
    enddo
    
    num_seeds_from_elem(1:num_elems) = 0 !initialise the count of nd
    do nd=1,num_data
       if(map_array(nd).ne.0)then
          ne_min = map_array(nd)
          num_seeds_from_elem(ne_min) = num_seeds_from_elem(ne_min)+1
       endif !map_array
    enddo !nd
    
!!! If there is only 0 or 1 seed point grouped with an element then set it as a
!!! terminal and remove a single seed point. Also involves modifying the local list of parents.
    
    if(.not.first)then
       N_ELM_TEMP=num_next_parents
       do N=1,num_next_parents
          ne_min=local_parent(N)
          if(num_seeds_from_elem(ne_min).eq.0)then !find closest point to end node
             nd_min = my_closest(N)
             map_array(nd_min)=0
             N_ELM_TEMP=N_ELM_TEMP-1
             local_parent(N)=0
             num_terminal=num_terminal+1

          else if(num_seeds_from_elem(ne_min).eq.1)then
             do nd=1,num_data
                if(map_array(nd).eq.ne_min)then
                   map_array(nd)=0
                   local_parent(N)=0
                   N_ELM_TEMP=N_ELM_TEMP-1
                   num_terminal=num_terminal+1
                endif
             enddo !nd
             
          endif !num_seeds_from_elem
       enddo !N
       
       do N=1,num_next_parents
          if(local_parent(N).eq.0)then
             I=0
             do while((N+I.lt.num_next_parents).and.(local_parent(N+I).eq.0))
                I=I+1
             enddo
             do M=N,num_next_parents-I
                local_parent(M)=local_parent(M+I)
             enddo !M
          endif !local_parent
       enddo !N
       num_next_parents = N_ELM_TEMP
       
       call sort_integer_list(num_next_parents,local_parent)
       
    endif
    
    deallocate(my_closest)
    if(.not.first) deallocate(map_array_copy)

    call enter_exit(sub_name,2)

  end subroutine group_seeds_with_branch
  

  !###############################################################
  !
  !*grow_tree:* the main growing subroutine (public). Genertes a volume-filling
  ! tree into a closed surface.
  !
  subroutine grow_tree(parent_ne,surface_elems,angle_max,angle_min,&
       branch_fraction,length_limit,shortest_length,rotation_limit)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GROW_TREE" :: GROW_TREE

    use diagnostics, only: enter_exit
    implicit none

    integer,intent(in)  :: parent_ne                ! list of end branch elements to grow from
    integer,intent(in)  :: surface_elems(:)         ! list of surface elements defining the host region
    real(dp),intent(in) :: angle_max                ! maximum branch angle with parent; in degrees
    real(dp),intent(in) :: angle_min                ! minimum branch angle with parent; in degrees
    real(dp),intent(in) :: branch_fraction          ! fraction of distance (to COFM) to branch
    real(dp),intent(in) :: length_limit             ! minimum length of a generated branch (shorter == terminal)
    real(dp),intent(in) :: shortest_length          ! length that short branches are reset to (shortest in model)
    real(dp),intent(in) :: rotation_limit           ! maximum angle of rotation of branching plane

    !Local variables
    integer,allocatable :: local_parent(:)          ! stores current generation of local parent elements 
    integer,allocatable :: local_parent_temp(:)     ! temporary storage of next generation of local parent elems
    integer,allocatable :: map_seed_to_elem(:)      ! records current elem associated w. data points
    integer,allocatable :: map_seed_to_space(:)     ! records initial elem associated w. data points (the 'space')
    integer,allocatable :: num_seeds_from_elem(:)   ! records # of seeds currently grouped with an elem
    integer,allocatable :: triangle(:,:)

    integer :: i,j,kount,M,N,nd,nd_min,ne,ne_grnd_parent,ne_parent,ne_stem,&
         noelem_parent,np,np_start,np_prnt_start,np_grnd_start,num_seeds_in_space,num_next_parents, &
         num_parents,num_triangles,num_vertices,num_elems_new,num_nodes_new,num_terminal
    
    real(dp),allocatable :: vertex_xyz(:,:)
    real(dp),dimension(3) :: COFM,candidate_xyz
    real(dp) :: distance_limit = 300.0_dp,length_parent

    logical :: make_branch,enough_points,first_group,internal, &
         limit_branching_angle = .true., &  ! option to restrict branch angle
         limit_branching_plane = .false.    ! option to restrict angle between branching planes

    character(len=60) :: sub_name

    sub_name = 'grow_tree'
    call enter_exit(sub_name,1)

    call triangles_from_surface(num_triangles,num_vertices,surface_elems,triangle,vertex_xyz)

!!! We can estimate the number of elements in the generated model based on the 
!!! number of data (seed) points. i.e. N = 2*N_data - 1. So the total number of 
!!! elements following tree generation will be ~ num_elems + 2*num_data. Use this estimate
!!! to increase the node and element arrays.
    num_elems_new = num_elems + 2*num_data + 100
    num_nodes_new = num_nodes + 2*num_data + 100
    call reallocate_node_elem_arrays(num_elems_new,num_nodes_new)
    
!!! Allocate memory for temporary arrays (need a more intelligent way of estimating size!)
    allocate(local_parent_temp(num_elems_new))
    allocate(local_parent(num_elems_new))
    allocate(num_seeds_from_elem(num_elems_new))
    allocate(map_seed_to_elem(num_data))
    allocate(map_seed_to_space(num_data))
    
!!! Initialise local_parent to the list of parent elements, and num_parents (current
!!! number of parent branches) to the number of parent branches.
    local_parent(1:size(parentlist)) = parentlist(1:size(parentlist))
    num_parents = count(parentlist.ne.0) !initial number of 'terminal' parent branches

    NUM_SEEDS_FROM_ELEM = 0
    num_next_parents = num_parents
    
!!! Calculate the initial grouping of data points with terminal elements
!!! this defines the 'space' with which each seed is associated
!!! For a single parent, all seed points will initially be mapped to
!!! it; for multiple parents 'group_seeds_with_branch' used to be called to calculate
!!! the closest parent end-point to each seed point. This has been replaced by splitting
!!! seed points using the orthogonal to branching planes of the upper tree.
    map_seed_to_space(1:num_data) = parentlist(1) !#! this is done for the new-style growing (full grow per terminal)
    if(num_parents.gt.1)then
       first_group = .true.
!       call group_seeds_with_branch(map_seed_to_space,num_next_parents,num_seeds_from_elem,&
       !            num_terminal,local_parent,500.0_dp,first_group)
       call split_seed_points_initial(map_seed_to_space,parent_ne)
    endif !parentlist.gt.1
    first_group = .false.

    WRITE(*,'(''  parent  #seeds  #terminal'')')
    
    ! Set initial values for local and global nodes and elements
    ne = num_elems !initialise mesh global element #
    np = num_nodes !initialise mesh global node #
    
!!! loop over the initial parent list (the initial conditions/terminal elements for growing)
!!! growing is done into 'spaces', where each 'space' is the initial grouping of seed
!!! points with the closest terminal branch. This grouping can be manipulated to take into account
!!! the size of the initial terminal branches. e.g. smaller diameter --> smaller set of seeds
    do noelem_parent = 1,num_parents
       ne_stem = parentlist(noelem_parent) ! the 'stem' parent element for the 'space'
       map_seed_to_elem = 0 ! initialise the seed mapping array
       num_seeds_in_space = 0 !initialise the number of seed points in the 'space'
       do nd = 1,num_data ! for all of the seed points (stored in data_xyz array)
          if(map_seed_to_space(nd).eq.ne_stem)then ! for the points in this space
             map_seed_to_elem(nd) = ne_stem ! record the current element associated with seed point nd
             num_seeds_in_space = num_seeds_in_space+1 ! count number of seed points in the space
          endif
       enddo

       num_next_parents = 1 ! initialise the number of current local parent branches
       local_parent(1) = ne_stem ! first local parent branch is the 'stem' branch 
       num_terminal = 0 ! initialise the number of definite terminal branches
       
!!! bifurcating distributive algorithm
       do while(num_next_parents.ne.0) !while still some parent branches with seed points
          num_parents = num_next_parents ! update the number of current local parent branches
          num_next_parents = 0 ! reset the number of local parent branches in next generation
          
          do M = 1,num_parents ! for each of the current local parent branches
             ne_parent = local_parent(M) !parent element #
             ! Calculate centre of mass of current seed point set
             call calculate_seed_cofm(map_seed_to_elem,ne_parent,COFM)
             
             ne_grnd_parent = elem_cnct(-1,1,ne_parent) !grandparent global element #
             np_start = elem_nodes(2,ne_parent) !parent global end node #
             np_prnt_start = elem_nodes(1,ne_parent) !parent global start node #
             np_grnd_start = elem_nodes(1,ne_grnd_parent) !grandparent global start node #
             
             length_parent = elem_field(ne_length,ne_parent)
             
!!! Split each set of seed points using the plane defined by the
!!! parent branch and the centre of mass. Seed points get associated with NEW elements (ne+1,ne+2)
             call split_seed_points(map_seed_to_elem,ne_parent,ne,np_start,&
                  np_prnt_start,np_grnd_start,COFM,enough_points)
             
!!! check whether enough seed points remaining in BOTH seed groups for branching to be done
!!! (note: this could be improved to continue branching in one set of seeds)
             if(enough_points)then
                do N = 1,2 !for each of the two new branches
                   ! Set up arrays for new element and node 
                   ! after create_new_node the current element == ne and current node == np
                   call create_new_node(ne,ne_parent,np,np_start,.TRUE.)
                   ! find the centre of mass of seed points
                   call calculate_seed_cofm(map_seed_to_elem,ne,COFM)
                   ! Generate a branch directed towards the centre of mass. Returns location
                   ! of end node in candidate_xyz (adjusted below based on length and shape criteria)         
                   call branch_to_cofm(map_seed_to_elem,ne,np_start,&
                        COFM,branch_fraction,length_limit,length_parent,shortest_length,&
                        candidate_xyz,make_branch)
                   node_xyz(1:3,np) = candidate_xyz(1:3) ! the new node location is as returned by 'branch_to_cofm'
                   call calc_branch_direction(ne) ! calculate direction of the new branch
                   elem_field(ne_length,ne) = distance_between_points(node_xyz(1,np_start),node_xyz(1,np))
                   ! Check whether this is a new parent branch or a terminal branch
                   if(make_branch.and.enough_points)then ! meets all criteria for continuing branching
                      num_next_parents = num_next_parents+1 ! increment the number of next parents
                      local_parent_temp(num_next_parents) = ne !records the elements that are parents 
                   else ! this is a terminal branch
                      num_terminal = num_terminal+1 ! increment the number of terminal branches
                   endif
                enddo !N (for both new branches)

                if(limit_branching_angle)then
                   ! Check that the branch angles are not too large or too small
                   ! Correct such that branches stay in the original branching plane
                   call limit_branch_angles(ne,ne_parent,np,&
                        np_prnt_start,np_start,angle_max,angle_min)

                   internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np),&
                        vertex_xyz)
                   if(.not.internal)then ! halve the length, and make terminal
                      elem_field(ne_length,ne) = 0.5_dp*distance_between_points&
                           (node_xyz(1,np_start),node_xyz(1,np))
                      node_xyz(1:3,np) = node_xyz(1:3,np_start) + 0.5_dp*&
                           elem_field(ne_length,ne)*elem_direction(1:3,ne)
                      internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np),&
                           vertex_xyz)
                      kount = 0
                      do while(.not.internal) 
                         kount = kount+1
                         call shorten_branch_and_children(ne_parent)
                         internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np),&
                              vertex_xyz)
                         if(kount.ge.3)then
                            call shorten_branch_and_children(elem_cnct(-1,1,ne_parent))
                            internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np),&
                                 vertex_xyz)
                         endif
                         if(kount.gt.5.and.(.not.internal))then
                            write(*,'('' WARNING: element'',i6,'' not internal'')') ne
                            internal = .true.
                         endif
                      enddo
                      if(inlist(ne,local_parent_temp))then ! set to be terminal
                         local_parent_temp(num_next_parents) = 0
                         num_next_parents = num_next_parents-1 ! decrement the number of next parents
                         num_terminal = num_terminal+1
                         nd_min = closest_seed_to_node_in_group(map_seed_to_elem,ne,np) ! closest seed point
                         map_seed_to_elem(nd_min) = 0 ! remove seed point from list
                      endif
                   endif !.not.internal

                   internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np-1),&
                        vertex_xyz)
                   if(.not.internal)then ! halve the length, and halve the angle from parent
                      elem_field(ne_length,ne-1) = 0.5_dp*distance_between_points&
                           (node_xyz(1,np_start),node_xyz(1,np-1))
                      node_xyz(1:3,np-1) = node_xyz(1:3,np_start) + 0.5_dp*&
                           elem_field(ne_length,ne-1)*elem_direction(1:3,ne-1)
                      internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np-1),&
                           vertex_xyz)
                      do while(.not.internal)
                         kount = kount+1
                         call shorten_branch_and_children(ne_parent)
                         internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np),&
                              vertex_xyz)
                         if(kount.ge.3)then
                            call shorten_branch_and_children(elem_cnct(-1,1,ne_parent))
                            internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np),&
                                 vertex_xyz)
                         endif
                         if(kount.gt.4.and.(.not.internal))then
                            write(*,'('' WARNING: element'',i6,'' not internal'')') ne-1
                            internal = .true.
                         endif
                      enddo
                      if(inlist(ne-1,local_parent_temp))then ! set to be terminal
                         do i = 1,num_next_parents
                            if(local_parent_temp(i).eq.ne-1)then
                               do j = i,num_next_parents-1
                                  local_parent_temp(j) = local_parent_temp(i+1)
                               enddo
                            endif
                         enddo
                         num_next_parents = num_next_parents-1 ! decrement the number of next parents
                         num_terminal = num_terminal+1
                         nd_min = closest_seed_to_node_in_group(map_seed_to_elem,ne-1,np-1) ! closest seed point
                         map_seed_to_elem(nd_min) = 0 ! remove seed point from list
                      endif
                   endif ! .not.internal

                endif

                if(limit_branching_plane)then
                   ! Check the angle of rotation between child and parent branching planes.
                   ! If absolute angle is larger than a user-specified limit, adjust to be
                   ! the limit value. This is used for making sure that the CFD geometry turns out ok.
                   call check_branch_rotation_plane(map_seed_to_elem,ne,ne_grnd_parent,ne_parent, &
                        local_parent_temp,num_next_parents, &
                        np,np_start,np_prnt_start,np_grnd_start,num_terminal,rotation_limit)
                   internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np),vertex_xyz)
                endif

             else
                write(*,*) 'terminal, not enough points',ne !!! never happens!!!
                read(*,*)
                ! Not enough seed points in the set during the split.
                ! Find the closest seed point to node np_start, and remove from seeds
                num_terminal = num_terminal+1 ! increment number of terminal branches
                nd_min = closest_seed_to_node(map_seed_to_elem,np_start) ! closest seed point
                map_seed_to_elem(nd_min) = 0 ! remove seed point from list
                map_seed_to_space(nd_min) = ne ! record the element to data point mapping
             endif
          enddo ! for all current parent branches
          ! Copy the temporary list of branches to local_parent. These become the
          ! parent elements for the next branching   
          local_parent(1:num_next_parents) = local_parent_temp(1:num_next_parents)
          ! Regroup the seed points with the closest current parent 
          call group_seeds_with_branch(map_seed_to_elem,num_next_parents,num_seeds_from_elem,&
               num_terminal,local_parent,DISTANCE_LIMIT,.FALSE.)
       enddo ! while still parent branches

       write(*,'(I7,I8,I9)') ne_parent,num_seeds_in_space,num_terminal

    enddo ! for each initial parent
    
!!! set new total numbers of nodes and elements
    num_nodes=np !highest node # in nr
    num_elems=ne !highest element # in nr
    
!!! update the tree connectivity
    call element_connectivity_1d
!!! calculate branch generations and orders
    call evaluate_ordering
!!! deallocate temporary arrays
    deallocate(vertex_xyz)
    deallocate(triangle)
    deallocate(local_parent_temp)
    deallocate(local_parent)
    deallocate(map_seed_to_elem)
    deallocate(map_seed_to_space)
    deallocate(num_seeds_from_elem)
    
    call enter_exit(sub_name,2)

  end subroutine grow_tree

  
  !###############################################################
  !
  !*limit_branch_angles:* checks both new branches at once to make sure the angles
  ! to parent branch not too large or too small. Both branches are checked to make sure
  ! they are not co-linear with parent. One branch is checked for its angle, and then 
  ! the other is made to be in-plane with the parent and sibling branch.
  !
  subroutine limit_branch_angles(ne,ne_parent,np,&
       np_prnt_start,np_start,angle_max,angle_min)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_LIMIT_BRANCH_ANGLES" :: LIMIT_BRANCH_ANGLES

    use diagnostics, only: enter_exit
    implicit none
    
    integer,intent(in) :: ne,ne_parent,np,np_prnt_start,np_start
    real(dp) :: angle_max,angle_min
    
    !Local variables
    integer :: elem_check,elem_in_plane,node_check,node_in_plane
    real(dp),dimension(3) :: B,NRML,U,V,W
    real(dp),dimension(3,3) :: A
    real(dp) :: angle_uw,angle_uv,angle_vw,length,length_w
    
    character(len=60) :: sub_name

    sub_name = 'limit_branch_angles'
    call enter_exit(sub_name,1)

    elem_check = ne-1  ! the default element to check angle for is the '1st' branch
    node_check = np-1  ! node associated with ne-1
    elem_in_plane = ne ! the default element to make in-plane is the '2nd' branch
    node_in_plane = np ! node associated with ne
    ! check whether default element is co-linear with parent element
    u(1:3) = elem_direction(1:3,ne_parent)  ! parent branch direction
    v(1:3) = elem_direction(1:3,ne-1)       ! new branch direction
    length = elem_field(ne_length,ne-1)
    nrml = cross_product(u,v) ! calculate normal
    if(abs(nrml(1))+abs(nrml(2))+abs(nrml(3)).lt.1.0e-4)then !co-linear branches
       ! check whether the other branch is co-linear with parent
       v(1:3) = elem_direction(1:3,ne)
       nrml = cross_product(u,v) ! calculate normal
       if(abs(nrml(1))+abs(nrml(2))+abs(nrml(3)).lt.1.0e-4)then !co-linear points
          write(*,*) 'warning: both branches co-linear with parent'
          read(*,*)
       else
          elem_check = ne
          node_check = np
          elem_in_plane = ne-1
          node_in_plane = np-1
       endif
    endif
    call adjust_branch_angle(1,elem_check,np_start,np_prnt_start,node_check, &
         angle_max,angle_min)
       
    
!!! for the other branch, make sure it is in-plane with parent and sibling branch
    u(1:3) = elem_direction(1:3,ne_parent)  ! parent branch direction
    v(1:3) = elem_direction(1:3,elem_check)       ! sibling branch direction
    w(1:3) = elem_direction(1:3,elem_in_plane)         ! current vector
    length_w = distance_between_points(node_xyz(1,np),node_xyz(1,np_start))
    nrml = cross_product(u,v) ! calculate normal to the parent and sibling branching plane
    nrml = unit_vector(nrml)
    !.....Adjust the direction of the second branch s.t. in plane
    angle_uw = angle_btwn_vectors(u,w)
    angle_uv = angle_btwn_vectors(u,v)
    angle_vw = angle_uw+angle_uv
!!! set up system of linear equations
    A(1,1:3) = u(1:3)     ! dotprod parent and new element
    A(2,1:3) = nrml(1:3)  ! dotprod normal and new element
    A(3,1:3) = v(1:3)     ! dotprod sibling and new element
    
    B(1) = DCOS(angle_uw) ! same angle with parent
    B(2) = 0.0_dp         ! in-plane with parent and sibling
    B(3) = DCOS(angle_vw) ! same angle with sibling
    
    w = mesh_a_x_eq_b(A,B)
    w = unit_vector(w)
    
    node_xyz(1:3,node_in_plane) = node_xyz(1:3,np_start)+length_w*w(1:3)
    elem_direction(1:3,elem_in_plane) = W(1:3)
    call adjust_branch_angle(2,elem_in_plane,np_start,np_prnt_start,node_in_plane,angle_max,angle_min)
    
    call enter_exit(sub_name,2)

  end subroutine limit_branch_angles


  !##################################################
  !
  !*reduce_branch_angle:* calculates the direction of a branch for a given branch angle

  subroutine reduce_branch_angle(np1,np2,np,candidate_xyz,factor)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_REDUCE_BRANCH_ANGLE" :: REDUCE_BRANCH_ANGLE

    use diagnostics, only: enter_exit
    implicit none
    
    integer,intent(in) :: np1,np2,np
    real(dp),intent(in) :: factor
    real(dp) :: candidate_xyz(3)

    !Local variables
    real(dp) :: angle,LV,U(3),V(3),W(3)
    
    character(len=60) :: sub_name

    sub_name = 'reduce_branch_angle'
    call enter_exit(sub_name,1)
    
    U(1:3)=node_xyz(1:3,np1)-node_xyz(1:3,np2) !direction of parent
    V(1:3)=candidate_xyz(1:3)-node_xyz(1:3,np1) !direction of this branch
    LV = vector_length(V) 
    U = unit_vector(U)
    V = unit_vector(V)
    angle = angle_btwn_vectors(U,V)

    W = vector_for_angle_limit(U,V,angle,factor*angle) 
    
    node_xyz(1:3,np) = node_xyz(1:3,np1) + W(1:3)*LV ! use original length
    candidate_xyz(1:3) = node_xyz(1:3,np) ! adjust the candidate node location
    
    call enter_exit(sub_name,2)

  end subroutine reduce_branch_angle
  

  !###############################################################
  !
  !*shorten_branch_and_children:* shorten the specified branch, its children,
  ! and their children
  !
  subroutine shorten_branch_and_children(ne)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SHORTEN_BRANCH_AND_CHILDREN" :: SHORTEN_BRANCH_AND_CHILDREN

    use diagnostics, only: enter_exit
    implicit none
    
    integer,intent(in) :: ne

    !Local variables
    integer :: i,j,ne1,ne2,np0,np1,np2
    character(len=60) :: sub_name = 'shorten_branch_and_children'

    call enter_exit(sub_name,1)

    np0 = elem_nodes(1,ne)
    np1 = elem_nodes(2,ne)

    elem_field(ne_length,ne) = 0.8_dp*elem_field(ne_length,ne)
    node_xyz(1:3,np1) = node_xyz(1:3,np0)+elem_field(ne_length,ne)*elem_direction(1:3,ne)

    np0 = np1
    do i=1,elem_cnct(1,0,ne) ! for each child
       ne1 = elem_cnct(1,i,ne)
       np1 = elem_nodes(2,ne1)
       elem_field(ne_length,ne1) = 0.8_dp*elem_field(ne_length,ne1)
       node_xyz(1:3,np1) = node_xyz(1:3,np0)+elem_field(ne_length,ne1)*elem_direction(1:3,ne1)
       do j = 1,elem_cnct(1,0,ne1)
          ne2 = elem_cnct(1,j,ne1)
          np2 = elem_nodes(2,ne2)
          elem_field(ne_length,ne2) = 0.8_dp*elem_field(ne_length,ne2)
          node_xyz(1:3,np2) = node_xyz(1:3,np1)+elem_field(ne_length,ne2)*elem_direction(1:3,ne2)
       enddo
    enddo

    call enter_exit(sub_name,2)

  end subroutine shorten_branch_and_children


  !###############################################################
  !
  !*smooth_1d_tree:* smooth a tree geometry by placing the end of a parent branch
  ! at the average of the parent end and child end coordinates. This is used to
  ! improve the topology of generated trees, minimising the impact of 'odd' branching
  !
  subroutine smooth_1d_tree(num_elem_start,length_limit)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SMOOTH_1D_TREE" :: SMOOTH_1D_TREE

    use diagnostics, only: enter_exit
    implicit none

    integer,intent(in) :: num_elem_start
    real(dp),intent(in) :: length_limit

    integer :: n,ne,ne1,ne2,np,np0,np1,np2,n_smoothing_steps = 2
    real(dp) :: new_xyz(3)
    character(len=60) :: sub_name

    sub_name = 'smooth_1d_tree'
    call enter_exit(sub_name,1)

    do n = 1,n_smoothing_steps
       do ne = num_elems,num_elem_start,-1
          if(elem_cnct(1,0,ne).eq.0)then ! terminal, check branch length
             if(elem_field(ne_length,ne).lt.0.75_dp*length_limit)then
                elem_field(ne_length,ne) = 0.75_dp*length_limit
                np1 = elem_nodes(1,ne) ! the start node
                np2 = elem_nodes(2,ne) ! the end node
                node_xyz(:,np2) = node_xyz(:,np1) + elem_direction(:,ne)*0.75_dp*length_limit
             else if(elem_field(ne_length,ne).gt.1.5_dp*length_limit)then
                elem_field(ne_length,ne) = 1.5_dp*length_limit
                np1 = elem_nodes(1,ne) ! the start node
                np2 = elem_nodes(2,ne) ! the end node
                node_xyz(:,np2) = node_xyz(:,np1) + elem_direction(:,ne)*1.5_dp*length_limit
             endif
          else
             if(elem_cnct(1,0,ne).eq.2)then
                ne1 = elem_cnct(1,1,ne)
                ne2 = elem_cnct(1,2,ne)
                np0 = elem_nodes(1,ne)
                np  = elem_nodes(2,ne)
                np1 = elem_nodes(2,ne1)
                np2 = elem_nodes(2,ne2)
                new_xyz(:) = node_xyz(:,np0)*0.5_dp + node_xyz(:,np1)*0.25_dp + node_xyz(:,np2)*0.25_dp
                node_xyz(:,np) = new_xyz(:)
             endif
          endif
       enddo
    enddo
    
    call enter_exit(sub_name,2)

  end subroutine smooth_1d_tree


  !###############################################################
  !
  !*split_seed_points:* divides a set of seed points into two subsets 
  ! using the plane that contains the parent branch and the seed point
  ! centre of mass. Decides which side of a plane a seed point is on by calculating
  ! the distance between two parallel planes: one which is defined by
  ! the parent and grandparent branch, and the other which contains a
  ! seed point.
  !
  subroutine split_seed_points(map_seed_to_elem,ne1,ne,np1,np2,np3,COFM,enough_points)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SPLIT_SEED_POINTS" :: SPLIT_SEED_POINTS
    
    use diagnostics, only: enter_exit
    implicit none
    
    integer :: map_seed_to_elem(*),ne,ne1,np1,np2,np3
    real(dp) :: COFM(3)
    logical :: enough_points

    !Local variables
    integer :: DAT1,DAT2,nd,ND1_1ST,ND2_1ST,nsp,NPOINTS
    real(dp) :: DIST,NORML(4),P(3),Q(3),R(3)
    logical :: COLINEAR
    
    character(len=60) :: sub_name

    sub_name = 'split_seed_points'
    call enter_exit(sub_name,1)
    
    R = COFM ! split based on cofm and branch
    P(1:3) = node_xyz(1:3,np2) ! point at start of parent branch
    Q(1:3) = node_xyz(1:3,np1) ! point at end of parent branch

!!! check whether the centre of mass and the parent start & end branches 
!!! are co-linear. if so, will need to use 'aunt' branch for split    
    colinear = check_colinear_points(P,Q,R)
    if(colinear) R(1:3) = node_xyz(1:3,np3) !split based on parent and aunt
    call make_plane_from_3points(NORML,1,P,Q,R) !calculate plane
    
    NPOINTS=0
    DAT1=0
    DAT2=0
    ND1_1ST=0
    ND2_1ST=0
    do nd=1,num_data
       nsp=map_seed_to_elem(nd) !space # that random point belongs to
       if(nsp.eq.ne1)then !random point belongs to this element space
          NPOINTS=NPOINTS+1
          dist = -scalar_product_3(norml,data_xyz(1,nd)) - norml(4) ! distance between two planes
          if(dist.ge.0.0_dp)then
             if(dat1.eq.0) nd1_1st = nd
             DAT1=DAT1+1
             map_seed_to_elem(nd)=ne+1
          else if(DIST.lt.0.0_dp)then
             if(DAT2.eq.0) ND2_1ST=nd
             DAT2=DAT2+1
             map_seed_to_elem(nd)=ne+2
          endif
       endif
    enddo !nd
    
    if(dat1.eq.0.and.dat2.eq.0)then
       enough_points = .false.
       write(*,'('' Zero seed points associated with parent'',I6)') ne1
       read(*,*)
    else
       if(dat1.eq.0)then
          map_seed_to_elem(nd2_1st) = ne+1
          dat1 = dat1+1
          dat2 = dat2-1
       elseif(dat2.eq.0)then
          map_seed_to_elem(nd1_1st) = ne+2
          dat2 = dat2+1
          dat1 = dat1-1
       endif
       
       enough_points = .true.    
    endif

    call enter_exit(sub_name,2)

  end subroutine split_seed_points
  
  
  !###############################################################
  !
  !*split_seed_points_initial:* divides a set of seed points into N subsets 
  ! to match N terminal branches, using the plane that is orthogonal to the branching plane
  ! of child branches, and that passes mid-way between child branches.
  !
  subroutine split_seed_points_initial(map_array,ne_stem)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SPLIT_SEED_POINTS_INITIAL" :: SPLIT_SEED_POINTS_INITIAL

    use diagnostics, only: enter_exit
    implicit none
    
    integer :: map_array(:),ne_stem

    !Local variables
    integer :: local_parent(20),local_parent_temp(20),M,nd,ne_parent, &
         ne1,ne2,np0,np1,np2,nsp,num_points,num_next_parents,num_parents
    real(dp) :: DIST,dist_p1,dist_p2,NORML(4),P(3),Q(3),R(3)
    
    character(len=60) :: sub_name

    sub_name = 'split_seed_points_initial'
    call enter_exit(sub_name,1)

    ne_parent = ne_stem
    do while(elem_cnct(1,0,ne_parent).eq.1)
       ne_parent = elem_cnct(1,1,ne_parent) ! get the next element in a refined branch
    enddo
    map_array(:) = ne_parent ! initialise that all seed points map to the stem branch

    num_next_parents = 1
    local_parent(1) = ne_parent
    
    do while(num_next_parents.ne.0) !while still some parent branches with seed points
       num_parents = num_next_parents ! update the number of current local parent branches
       num_next_parents = 0 ! reset the number of local parent branches in next generation
       do M = 1,num_parents ! for each of the current local parent branches
          ne_parent = local_parent(M) !parent element #
          do while(elem_cnct(1,0,ne_parent).eq.1)
             ne_parent = elem_cnct(1,1,ne_parent) ! get the next element in a refined branch
          enddo
          np0 = elem_nodes(2,ne_parent)
          ne1 = elem_cnct(1,1,ne_parent)
          do while(elem_cnct(1,0,ne1).eq.1)
             ne1 = elem_cnct(1,1,ne1) ! get the next element in a refined branch
          enddo
          ne2 = elem_cnct(1,2,ne_parent)
          do while(elem_cnct(1,0,ne2).eq.1)
             ne2 = elem_cnct(1,1,ne2) ! get the next element in a refined branch
          enddo
          np1 = elem_nodes(2,ne1)
          np2 = elem_nodes(2,ne2)
          P(:) = node_xyz(:,np0) ! point at end of parent branch
          Q(:) = node_xyz(:,np1) ! point at end of child1 branch
          R(:) = node_xyz(:,np2) ! point at end of child2 branch
          call make_plane_from_3points(NORML,1,P,Q,R) !calculate plane
          
          P(:) = node_xyz(:,np0) ! point at end of parent branch
          Q(:) = 0.5_dp*(node_xyz(:,np1)+node_xyz(:,np2))
          R(1:3) = Q(1:3) + NORML(1:3)
          call make_plane_from_3points(NORML,1,P,Q,R) !calculate plane
!!! NORML is now the plane between the child branches    
          
          dist_p1 = -scalar_product_3(norml,node_xyz(:,np1)) - norml(4) ! distance between two planes
          dist_p2 = -scalar_product_3(norml,node_xyz(:,np2)) - norml(4) ! distance between two planes

          do nd = 1,num_data
             nsp = map_array(nd) !space # that random point belongs to
             if(nsp.eq.ne_parent)then !random point belongs to this element space
                dist = -scalar_product_3(norml,data_xyz(1,nd)) - norml(4) ! distance between two planes
                if(dist.ge.0.0_dp.and.dist_p1.ge.0.0_dp)then
                   map_array(nd) = ne1
                else if(dist.ge.0.0_dp.and.dist_p1.lt.0.0_dp)then
                   map_array(nd) = ne2
                else if(dist.le.0.0_dp.and.dist_p2.le.0.0_dp)then
                   map_array(nd) = ne2
                else if(dist.le.0.0_dp.and.dist_p2.gt.0.0_dp)then
                   map_array(nd) = ne1
                endif
             endif
          enddo !nd

          num_points = count(map_array.eq.ne1)
          if(num_points.eq.0)then
             write(*,'('' Warning: number of points for element'',i6,'' is zero'')') ne1
             write(*,'('' Press enter to continue; however the code is likely to fail'')')
          endif
          num_points = count(map_array.eq.ne2)
          if(num_points.eq.0)then
             write(*,'('' Warning: number of points for element'',i6,'' is zero'')') ne2
             write(*,'('' Press enter to continue; however the code is likely to fail'')')
          endif
       
          if(elem_cnct(1,0,ne1).ne.0)then
             num_next_parents = num_next_parents+1
             local_parent_temp(num_next_parents) = ne1
          endif
          if(elem_cnct(1,0,ne2).ne.0)then
             num_next_parents = num_next_parents+1
             local_parent_temp(num_next_parents) = ne2
          endif
       enddo ! num_parents
       local_parent(1:num_next_parents) = local_parent_temp(1:num_next_parents)
    enddo 

    call enter_exit(sub_name,2)

  end subroutine split_seed_points_initial
  
  
  !###############################################################
  !
  !*closest_seed_to_node:* returns the closest seed point to a given branch node
  !
  function closest_seed_to_node(map_seed_to_elem,np)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CLOSEST_SEED_TO_NODE" :: CLOSEST_SEED_TO_NODE

    use diagnostics, only: enter_exit
    implicit none
    
    integer,intent(in) :: map_seed_to_elem(*),np
    
    !Local variables
    integer :: nd
    integer :: closest_seed_to_node
    real(dp) :: distance,min_distance
    
    min_distance = 1.0e+10_dp
    do nd=1,num_data
       if(map_seed_to_elem(nd).ne.0)then
          distance = distance_between_points(data_xyz(1,nd),node_xyz(1,np))
          if(distance.lt.min_distance)then
             closest_seed_to_node = nd
             min_distance = distance
          endif !DIST
       endif !map_seed_to_elem
    enddo !nd
    
  end function closest_seed_to_node
  

  !###############################################################
  !
  !*closest_seed_to_node_in_group:* finds the closest seed point to a node
  ! that is in a group of seed points currently associated with a specific element
  !
  function closest_seed_to_node_in_group(map_seed_to_elem,ne,np)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CLOSEST_SEED_TO_NODE_IN_GROUP" :: CLOSEST_SEED_TO_NODE_IN_GROUP
    
    integer,intent(in) :: map_seed_to_elem(*),ne,np
    
    !Local variables
    integer :: nd
    integer :: closest_seed_to_node_in_group
    real(dp) :: distance,min_distance
    
    min_distance = 1.0e+10_dp
    do nd=1,num_data
       if(map_seed_to_elem(nd).eq.ne)then
          distance = distance_between_points(data_xyz(1,nd),node_xyz(1,np))
          if(distance.lt.min_distance)then
             closest_seed_to_node_in_group = nd
             min_distance = distance
          endif !DIST
       endif !map_seed_to_elem
    enddo !nd
    
  end function closest_seed_to_node_in_group
  

  !###############################################################
  !
  !*rotation_angle:* calculates angle between two branching planes
  !
  function rotation_angle(np1,np2,np3,np4,np5)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ROTATION_ANGLE" :: ROTATION_ANGLE
    
    integer,intent(in) :: np1,np2,np3,np4,np5

    !Local variables
    real(dp) :: norm_1(4),norm_2(4)
    real(dp) :: rotation_angle
    
    call make_plane_from_3points(norm_1,2,node_xyz(1,np1),node_xyz(1,np2),node_xyz(1,np3))
    call make_plane_from_3points(norm_2,2,node_xyz(1,np2),node_xyz(1,np4),node_xyz(1,np5))
    rotation_angle = angle_btwn_vectors(norm_1,norm_2)
    
  end function rotation_angle
  

  !###############################################################
  !
  !*vector_for_angle_limit:* Calculates the new direction
  ! of an element, given current and target angles. Specifically, solves small
  ! system of equations to get new direction of branch (w), such that the branch
  ! remains in-plane (n.w = 0), the angle of w with u is defined 
  ! (as u.w = cos(angle_with_u)), and the angle of w with original direction (v) 
  ! is defined (as v.w = cos(angle_with_v)).
  !
  function vector_for_angle_limit(U,V,angle_with_u,angle_with_v)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_VECTOR_FOR_ANGLE_LIMIT" :: VECTOR_FOR_ANGLE_LIMIT

    real(dp),intent(in) :: U(*),V(*),angle_with_u,angle_with_v

    !Local variables
    real(dp) :: A(3,3),N_UV(3),VECTOR(3),W(3)
    real(dp) :: vector_for_angle_limit(3)

    N_UV = cross_product(U,V)     ! calculate the normal to the vectors U and V
    N_UV = unit_vector(N_UV) ! unit vector for normal
    A(1,1:3) = N_UV(1:3)
    A(2,1:3) = U(1:3)
    A(3,1:3) = V(1:3)
       
    VECTOR(1) = 0.0_dp
    VECTOR(2) = DCOS(angle_with_u)
    VECTOR(3) = DCOS(angle_with_v)
       
    w = mesh_a_x_eq_b(A,VECTOR)
    vector_for_angle_limit = unit_vector(W)

  end function vector_for_angle_limit
  

  !###############################################################
  
end module growtree

!###############################################################
  
