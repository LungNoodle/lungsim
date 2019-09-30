MODULE growtree
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

  USE arrays
  !use diagnostics, only: enter_exit      !! for diagnostics at enter/exit of subroutines
  USE geometry
  USE indices
  USE other_consts   !! pi
  USE mesh_utilities   !! general functions for geometric/mesh calculations
  USE math_utilities   !! general utility functions for sorting etc
  IMPLICIT NONE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  PRIVATE
  PUBLIC grow_tree,smooth_1d_tree

CONTAINS

  !###############################################################
  !
  !*adjust_branch_angle:* Adjusts child branch angle so not larger than a max value,
  ! keeping child in-plane. Checks the branching angle between a parent and child branch.
  ! If the branch angle is greater than the branch angle limit, then the branch angle is
  ! reduced to the limit value, such that the daughter branch remains in the original
  ! branching plane.
  !
  SUBROUTINE adjust_branch_angle(Nth,ne,np1,np2,np,angle_max,angle_min)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ADJUST_BRANCH_ANGLE" :: ADJUST_BRANCH_ANGLE
    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER :: Nth,ne,np1,np2,np
    REAL(dp) :: angle_max,angle_min

    !Local variables
    REAL(dp) :: a_min,a_lim,angle,angle_sibling,length,LU,LV,U(3),V(3),W(3)
    REAL(dp),PARAMETER :: loose_tol = 1.0e-4_dp

    CHARACTER(len=60) :: sub_name

    sub_name = 'adjust_branch_angle'
    CALL enter_exit(sub_name,1)

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

    IF(ABS(angle).GT.a_lim)THEN !reduce angle
       IF(ABS(dabs(angle)-pi).LT.loose_tol)THEN !reduce angle from 180 degrees
          IF(Nth.EQ.2)THEN
             V(1:3) = elem_direction(1:3,ne-1) !direction of sibling
          ELSE IF(Nth.EQ.1)THEN
             V(1:3) = elem_direction(1:3,ne+1) !direction of sibling
          ENDIF
          angle_sibling = angle_btwn_vectors(U,V)
          W = vector_for_angle_limit(U,V,a_lim,angle_sibling+a_lim)
       ELSE
          W = vector_for_angle_limit(U,V,a_lim,angle-a_lim)
       ENDIF
       node_xyz(1:3,np) = node_xyz(1:3,np1)+W(1:3)*0.5_dp*length
       elem_field(ne_length,ne) = 0.5_dp*length
       elem_direction(1:3,ne) = W(1:3)

    ELSEIF(ABS(angle).LT.a_min)THEN
       IF(ABS(angle).LT.loose_tol)THEN !increase angle from 0 degrees
          IF(Nth.EQ.2)THEN
             V(1:3) = elem_direction(1:3,ne-1) !direction of sibling
          ELSE IF(Nth.EQ.1)THEN
             V(1:3) = elem_direction(1:3,ne+1) !direction of sibling
          ENDIF
          angle_sibling = angle_btwn_vectors(U,V) !angle between branch and sibling
          W = vector_for_angle_limit(U,V,a_min,angle_sibling+a_min)
       ELSE
          W = vector_for_angle_limit(U,V,a_min,angle-a_min)
       ENDIF
       node_xyz(1:3,np) = node_xyz(1:3,np1)+W(1:3)*0.5_dp*length
       elem_field(ne_length,ne) = 0.5_dp*length
       elem_direction(1:3,ne) = W(1:3)

    ENDIF !abs(angle).gt.a_lim

    CALL enter_exit(sub_name,2)

  END SUBROUTINE adjust_branch_angle


  !###############################################################
  !
  !*branch_to_cofm*: Creates a new branch towards the cofm of a set of seed points.
  ! Used in a volume-filling branching method to create a branch that runs from a
  ! defined point to some fraction along a line towards the centre of mass of a
  ! collection of seed points.
  !
  SUBROUTINE branch_to_cofm(map_seed_to_elem,nen,np1,COFM,branch_fraction,length_limit,&
       length_parent,shortest_length,candidate_xyz,make_branch)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_BRANCH_TO_COFM" :: BRANCH_TO_COFM

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER :: map_seed_to_elem(*),nen,np1
    REAL(dp) :: COFM(3),branch_fraction,length_limit,length_parent,&
         shortest_length,candidate_xyz(3)
    LOGICAL :: make_branch

    !Local variables
    INTEGER :: N,NCLOSEST(100),nd,nsp,NUM_CLOSEST,NUM_ND,number_of_points
    REAL(dp) :: CLOSEST(100),DIST,L_COFM,LENGTH,MIN_DIST,VECTOR(3)

    CHARACTER(len=60) :: sub_name

    sub_name = 'branch_to_cofm'
    CALL enter_exit(sub_name,1)

    number_of_points = 0
    DO nd = 1,num_data
       IF(map_seed_to_elem(nd).EQ.nen) number_of_points=number_of_points+1
    ENDDO
    candidate_xyz(1:3) = node_xyz(1:3,np1) + branch_fraction*(COFM(1:3)-node_xyz(1:3,np1))
    VECTOR(1:3) = COFM(1:3)-node_xyz(1:3,np1)
    LENGTH = distance_between_points(candidate_xyz,node_xyz(1,np1))
    L_COFM = vector_length(VECTOR)
    VECTOR = unit_vector(VECTOR)

    IF(LENGTH.GE.LENGTH_LIMIT)THEN !the branch will not be terminal
       make_branch = .TRUE.
    ELSE
       IF(elem_ordrs(no_gen,nen).LT.12)THEN
          make_branch = .TRUE.
          candidate_xyz(1:3) = node_xyz(1:3,np1)+VECTOR(1:3)*0.33_dp*length_parent
       ELSE
          make_branch = .FALSE.
       ENDIF
    ENDIF !LENGTH.ge.L_LIM

    !*** Position branch
    MIN_DIST = 1.0e+6_dp
    NUM_ND = 0
    NUM_CLOSEST = 1
    NCLOSEST(1) = 1
    CLOSEST(1) = 1.0e+6_dp
    DO nd = 1,num_data
       nsp = map_seed_to_elem(nd) !space # that random point belongs to
       IF(nsp.EQ.nen)THEN !random point belongs to this element space
          dist = distance_between_points(data_xyz(1,nd),candidate_xyz)
          IF(DIST.LT.MIN_DIST) MIN_DIST = DIST
          NUM_ND = NUM_ND+1
          IF(.NOT.make_branch)THEN !remove closest data points
             IF(DIST.LT.CLOSEST(NUM_CLOSEST))THEN !store this data point
                IF(NUM_CLOSEST.LT.1)THEN
                   NUM_CLOSEST = NUM_CLOSEST+1 !increment number of closest
                ENDIF
                CLOSEST(NUM_CLOSEST) = DIST !store distance
                NCLOSEST(NUM_CLOSEST) = nd !store data point number
             ENDIF !DIST
             CALL sort_real_list(NUM_CLOSEST,CLOSEST,NCLOSEST) !sort into ascending
          ENDIF !NOT.make_branch
       ENDIF !nsp
    ENDDO !nd

    IF(LENGTH.LT.0.5_dp*length_parent)THEN
       candidate_xyz(1:3) = node_xyz(1:3,np1)+VECTOR(1:3)*0.5_dp*length_parent
    ENDIF

    IF(LENGTH.LT.shortest_length)THEN
       candidate_xyz(1:3) = node_xyz(1:3,np1)+VECTOR(1:3)*shortest_length
    ENDIF

    IF(.NOT.make_branch)THEN !remove the closest data points
       DO N = 1,NUM_CLOSEST
          nd = NCLOSEST(N)
          map_seed_to_elem(nd) = 0
       ENDDO
    ENDIF !NOT.make_branch

    CALL enter_exit(sub_name,2)

  END SUBROUTINE branch_to_cofm


  !###############################################################
  !
  !*calculate_seed_cofm:* calculates centre of mass of a list of seed
  ! points by averaging their coordinates and returns result in 'cofm'
  !
  SUBROUTINE calculate_seed_cofm(map_seed_to_elem,nen,COFM)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CALCULATE_SEED_COFM" :: CALCULATE_SEED_COFM

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER :: map_seed_to_elem(*),nen
    REAL(dp) :: COFM(3)

    !Local variables
    INTEGER :: DAT,nd,nsp

    CHARACTER(len=60) :: sub_name

    sub_name = 'calculate_seed_cofm'
    CALL enter_exit(sub_name,1)

    DAT = 0
    cofm = 0.0_dp

    DO nd=1,num_data
       nsp=map_seed_to_elem(nd) !the space # for the nd-th data point
       IF(nsp.EQ.nen)THEN
          DAT=DAT+1
          COFM(1:3)=COFM(1:3)+data_xyz(1:3,nd)
       ENDIF
    ENDDO !nd
    IF(DAT.NE.0) COFM(1:3) = COFM(1:3)/DAT !centre of mass

    CALL enter_exit(sub_name,2)

  END SUBROUTINE calculate_seed_cofm


  !###############################################################
  !
  !*check_branch_rotation_plane:* limits the angle between branching planes
  ! to a maximum (user-defined) value, and makes sure branch remains internal
  ! to the host volume
  !
  SUBROUTINE check_branch_rotation_plane(map_seed_to_elem,ne,&
       ne_grnd_parent,ne_parent,local_parent_temp,num_next_parents,&
       np,np1,np2,np3,num_terminal,rotation_limit)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CHECK_BRANCH_ROTATION_PLANE" :: CHECK_BRANCH_ROTATION_PLANE

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER :: map_seed_to_elem(*),ne,ne_grnd_parent,ne_parent, &
         num_next_parents,np,np1,np2,np3,num_terminal
    ! np == end node; np1 == np_start; np2 == np_prnt_start; np3 == np_grnd_start
    INTEGER :: local_parent_temp(*)
    REAL(dp),INTENT(in) :: rotation_limit

    !Local variables
    INTEGER :: COUNT,nd_min,ne_other,nes,np4,offset
    DOUBLE PRECISION :: ANGLE,length,ROT_ANGLE,candidate_xyz(3)
    LOGICAL :: INTERNAL

    CHARACTER(len=60) :: sub_name

    sub_name = 'check_branch_rotation_plane'
    CALL enter_exit(sub_name,1)

    ! find the appropriate other point for calculating the branching plane
    IF(elem_cnct(1,1,ne_grnd_parent).EQ.ne_parent)THEN
       IF(elem_cnct(1,2,ne_grnd_parent).EQ.0)THEN
          ne_other = elem_cnct(-1,1,ne_grnd_parent)
          np4 = elem_nodes(1,ne_other)
       ELSE
          ne_other = elem_cnct(1,2,ne_grnd_parent)
          np4 = elem_nodes(2,ne_other)
       ENDIF
    ELSE
       ne_other = elem_cnct(1,1,ne_grnd_parent)
       np4 = elem_nodes(2,ne_other)
    ENDIF

    ROT_ANGLE=ROTATION_LIMIT*PI/180.0_dp
    CALL CHECK_ROTATION_ANGLE(ne,np4,np3,np2,np1,np,np-1,rotation_limit)
    angle = rotation_angle(np2,np1,np4,np,np-1)
    INTERNAL=.FALSE.
    COUNT=0
    DO WHILE(.NOT.INTERNAL.AND.COUNT.LT.2)
       candidate_xyz(1:3) = node_xyz(1:3,np-1)
       !                  call CHECK_POINT_INTERNAL(IBT,Ido,INP,NBJ,NDLIST,NEELEM, &
       !                       NEP,NHOST,NKJE,np-1,np1,NPF,NPNE,NVJE,elem_cnct,SE, &
       !                       XA,XE,XIP,XP,candidate_xyz,data_xyz,INTERNAL,.FALSE.)
       internal=.TRUE.
       IF(.NOT.internal)THEN ! halve the length, and halve the angle from parent
          length = distance_between_points(node_xyz(1,np1),node_xyz(1,np-1))
          !          node_xyz(1:3,np-1) = node_xyz(1:3,np1) + 0.5_dp*length*elem_direction(1:3,np-1)
          node_xyz(1:3,np-1) = node_xyz(1:3,np1) + 0.5_dp*length*elem_direction(1:3,ne-1)
          candidate_xyz(1:3) = node_xyz(1:3,np-1)
          CALL reduce_branch_angle(np1,np2,np-1,candidate_xyz,0.5_dp) ! reduces the branch angle by a half
       ENDIF !internal
       COUNT=COUNT+1
    ENDDO
    IF(.NOT.INTERNAL)THEN
       !...............Remove the branch from the list of next generation parents
       offset=0
       DO nes=1,num_next_parents
          IF(local_parent_temp(nes).EQ.ne-1) offset=1
          local_parent_temp(nes)=local_parent_temp(nes+offset)
       ENDDO
       num_next_parents=num_next_parents-1
       num_terminal=num_terminal+1

       !...............Remove the closest data point to the end of the branch
       nd_min = closest_seed_to_node(map_seed_to_elem,np-1)

       map_seed_to_elem(nd_min)=0

    ENDIF

    INTERNAL=.FALSE.
    COUNT=0
    DO WHILE(.NOT.INTERNAL.AND.COUNT.LT.2)
       candidate_xyz(1:3) = node_xyz(1:3,np)
       !                  call CHECK_POINT_INTERNAL(IBT,Ido,INP,NBJ,NDLIST,NEELEM, &
       !                       NEP,NHOST,NKJE,np,np1,NPF,NPNE,NVJE,elem_cnct,SE,XA, &
       !                       XE,XIP,XP,candidate_xyz,data_xyz,INTERNAL,.FALSE.)
       internal=.TRUE.
       IF(.NOT.internal)THEN
          length = distance_between_points(node_xyz(1,np1),node_xyz(1,np))
          !          node_xyz(1:3,np) = node_xyz(1:3,np1)+0.5_dp*length*elem_direction(1:3,np)
          node_xyz(1:3,np) = node_xyz(1:3,np1)+0.5_dp*length*elem_direction(1:3,ne)
          CALL reduce_branch_angle(np1,np2,np,candidate_xyz,0.5_dp) ! reduces the branch angle by a half
       ENDIF
       COUNT=COUNT+1
    ENDDO !while
    IF(.NOT.INTERNAL)THEN
       !...............Remove the branch from the list of next generation parents
       offset=0
       DO nes=1,num_next_parents
          IF(local_parent_temp(nes).EQ.ne) offset=1
          local_parent_temp(nes)=local_parent_temp(nes+offset)
       ENDDO
       num_next_parents=num_next_parents-1
       num_terminal=num_terminal+1
       !...............Remove the closest data point to the end of the branch
       nd_min = closest_seed_to_node(map_seed_to_elem,np-1)

       map_seed_to_elem(nd_min)=0

    ENDIF

    CALL enter_exit(sub_name,2)

  END SUBROUTINE check_branch_rotation_plane


  !##################################################
  !
  !*check_rotation_angle:* adjusts the branch locations such that the angle between
  ! branching planes is less than a user-defined maximum.
  ! Calculates using quaternions. For angle ROTATION_ANGLE and unit
  ! vector a,b,c , calculate rotation matrix for arbitrary point.
  !
  SUBROUTINE check_rotation_angle(ne,np00,np0,np1,np2,np3,np4,rotation_limit)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CHECK_ROTATION_ANGLE" :: CHECK_ROTATION_ANGLE

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER,INTENT(in) :: ne,np00,np0,np1,np2,np3,np4
    REAL(dp),INTENT(in) :: rotation_limit

    ! Local variables
    INTEGER :: IT,ITMAX
    REAL(dp) :: ANGLE0,ANGLE,AXIS(3),DIRECTION(3),NRML(3), &
         NRML_PARENT(3),ROT_ANGLE,U(3),V(3),Q0,Q1,Q2,Q3,Q(3,3),X(3), &
         ANGLE_BETWEEN,length
    LOGICAL :: COMPLETE

    CHARACTER(len=60) :: sub_name

    sub_name = 'check_rotation_angle'
    CALL enter_exit(sub_name,1)

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
    IF(ANGLE_BETWEEN.GT.PI/2.0_dp)THEN
       ANGLE_BETWEEN=ANGLE_BETWEEN-PI
       ANGLE=ANGLE_BETWEEN
    ENDIF
    !      ANGLE=PI/2.0_dp-ANGLE
    ANGLE=PI/2.0_dp+ANGLE

    IF(ABS(ANGLE_BETWEEN).GT.ROTATION_LIMIT.AND.ABS(ANGLE_BETWEEN) &
         .LT.PI/2.0_dp-ROTATION_LIMIT)THEN
       IF(ANGLE.LT.0.0_dp)THEN
          ROT_ANGLE=-(ANGLE+ROTATION_LIMIT)
       ELSE
          ROT_ANGLE=-(ANGLE-ROTATION_LIMIT)
       ENDIF

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
       IF(ABS(ANGLE_BETWEEN).GT.ROTATION_LIMIT.AND.ABS(ANGLE_BETWEEN) &
            .LT.PI/2.0_dp-ROTATION_LIMIT)THEN
          COMPLETE=.TRUE.
       ELSE
          COMPLETE=.FALSE.
       ENDIF
       DO WHILE(.NOT.COMPLETE)
          IT=IT+1
          ANGLE_BETWEEN=ANGLE
          ANGLE=PI/2.0_dp-ANGLE
          IF(ANGLE.LT.0.0_dp)THEN
             ROT_ANGLE=-(ANGLE+ROTATION_LIMIT)
          ELSE
             ROT_ANGLE=-(ANGLE-ROTATION_LIMIT)
          ENDIF

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

          IF(ABS(ROTATION_LIMIT-ABS(PI/2.0_dp-ANGLE)).LE.0.001_dp)THEN
             COMPLETE=.TRUE.
          ENDIF

          IF(IT.GT.ITMAX)THEN
             WRITE(*,*) 'WARNING!!!! rotation angle = ',ANGLE*180.0_dp/PI
          ENDIF

       ENDDO !do while not found

       !.......Alternate calculation for the rotation angle
       angle = rotation_angle(np1,np2,np00,np3,np4)
    ELSE
       !        write(*,*) 'Not',ANGLE_BETWEEN*180.0_dp/PI
    ENDIF

    CALL enter_exit(sub_name,2)

  END SUBROUTINE check_rotation_angle


  !###############################################################
  !
  !*create_new_node:* sets up arrays for a new mesh node and element.
  !
  SUBROUTINE create_new_node(ne,ne_start,np,np_start,MAKE)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CREATE_NEW_NODE" :: CREATE_NEW_NODE

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER :: ne,ne_start,np,np_start
    LOGICAL :: MAKE

    !Local variables
    CHARACTER(len=60) :: sub_name

    sub_name = 'create_new_node'
    CALL enter_exit(sub_name,1)

    IF(MAKE)THEN
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
       IF(ne_start.NE.0)THEN
          elem_cnct(-1,0,ne)=1
          elem_cnct(-1,elem_cnct(-1,0,ne),ne)=ne_start
          elem_cnct(1,0,ne_start)=elem_cnct(1,0,ne_start)+1
          elem_cnct(1,elem_cnct(1,0,ne_start),ne_start)=ne
       ENDIF
    ENDIF

    CALL enter_exit(sub_name,2)

  END SUBROUTINE create_new_node


  !###############################################################
  !
  !*group_seeds_with_branch:* groups a set of seed points with the
  ! closest candidate parent branches. reassigns data (seed) points
  ! to the closest ending of branches in the current generation.
  !
  SUBROUTINE group_seeds_with_branch(map_array,num_next_parents,num_seeds_from_elem, &
       num_terminal,local_parent,DISTANCE_LIMIT,FIRST)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GROUP_SEEDS_WITH_BRANCH" :: GROUP_SEEDS_WITH_BRANCH

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER :: num_next_parents,local_parent(:),map_array(:),num_seeds_from_elem(*),&
         num_terminal
    REAL(dp),INTENT(in) :: DISTANCE_LIMIT
    LOGICAL :: FIRST

    !Local variables
    INTEGER :: i,n,m,nd,nd_min,ne,n_elm_temp,ne_min,noelem,np,np_temp
    INTEGER :: size_map
    INTEGER,ALLOCATABLE :: map_array_copy(:),my_closest(:)
    REAL(dp) :: dist,min_dist

    CHARACTER(len=60) :: sub_name

    sub_name = 'group_seeds_with_branch'
    CALL enter_exit(sub_name,1)

    size_map = SIZE(map_array)
    ALLOCATE(my_closest(size_map))

    IF(first)THEN
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

    ELSE !use the data groupings from previous
       ALLOCATE(map_array_copy(size_map))
       map_array_copy(1:size_map) = map_array(1:size_map)
       DO n=1,num_next_parents
          ne_min = local_parent(n)
          np_temp = elem_nodes(2,ne_min)
          MIN_DIST=1.0e+10_dp
          DO nd=1,num_data
             IF(map_array(nd).EQ.ne_min)THEN ! was associated with this element
                dist = distance_between_points(data_xyz(1,nd),node_xyz(1,np_temp))
                IF(dist.LT.min_dist)THEN
                   nd_min = nd
                   min_dist = dist
                ENDIF !DIST
             ENDIF
          ENDDO
          my_closest(N) = nd_min
          map_array(nd_min) = ne_min
       ENDDO
    ENDIF

    DO nd = 1,num_data            ! for all seed/data points
       IF(map_array(nd).NE.0)THEN ! the data point is still in use
          IF(.NOT.inlist(nd,my_closest))THEN
             MIN_DIST=1.0e+10_dp     ! initialise the minimum (closest) distance
             DO noelem = 1,num_next_parents ! for each parent in the next branch generation
                ne = local_parent(noelem)
                np = elem_nodes(2,ne)
                dist = distance_between_points(data_xyz(1,nd),node_xyz(1,np))
                IF(DIST.LT.MIN_DIST)THEN
                   ne_min = ne
                   MIN_DIST=DIST
                ENDIF
             ENDDO
             IF(first)THEN
                map_array(nd)=ne_min
             ELSE
                IF(MIN_DIST.LT.DISTANCE_LIMIT)THEN !keep seed points
                   map_array(nd)=ne_min
                ELSE
                   map_array(nd)=0 !too far from branch ends, so discard
                ENDIF
             ENDIF
          ENDIF
       ENDIF
    ENDDO

    num_seeds_from_elem(1:num_elems) = 0 !initialise the count of nd
    DO nd=1,num_data
       IF(map_array(nd).NE.0)THEN
          ne_min = map_array(nd)
          num_seeds_from_elem(ne_min) = num_seeds_from_elem(ne_min)+1
       ENDIF !map_array
    ENDDO !nd

!!! If there is only 0 or 1 seed point grouped with an element then set it as a
!!! terminal and remove a single seed point. Also involves modifying the local list of parents.

    IF(.NOT.first)THEN
       N_ELM_TEMP=num_next_parents
       DO N=1,num_next_parents
          ne_min=local_parent(N)
          IF(num_seeds_from_elem(ne_min).EQ.0)THEN !find closest point to end node
             nd_min = my_closest(N)
             map_array(nd_min)=0
             N_ELM_TEMP=N_ELM_TEMP-1
             local_parent(N)=0
             num_terminal=num_terminal+1

          ELSE IF(num_seeds_from_elem(ne_min).EQ.1)THEN
             DO nd=1,num_data
                IF(map_array(nd).EQ.ne_min)THEN
                   map_array(nd)=0
                   local_parent(N)=0
                   N_ELM_TEMP=N_ELM_TEMP-1
                   num_terminal=num_terminal+1
                ENDIF
             ENDDO !nd

          ENDIF !num_seeds_from_elem
       ENDDO !N

       DO N=1,num_next_parents
          IF(local_parent(N).EQ.0)THEN
             I=0
             DO WHILE((N+I.LT.num_next_parents).AND.(local_parent(N+I).EQ.0))
                I=I+1
             ENDDO
             DO M=N,num_next_parents-I
                local_parent(M)=local_parent(M+I)
             ENDDO !M
          ENDIF !local_parent
       ENDDO !N
       num_next_parents = N_ELM_TEMP

       CALL sort_integer_list(num_next_parents,local_parent)

    ENDIF

    DEALLOCATE(my_closest)
    IF(.NOT.first) DEALLOCATE(map_array_copy)

    CALL enter_exit(sub_name,2)

  END SUBROUTINE group_seeds_with_branch


  !###############################################################
  !
  !*grow_tree:* the main growing subroutine (public). Genertes a volume-filling
  ! tree into a closed surface.
  !
  SUBROUTINE grow_tree(parent_ne,surface_elems,angle_max,angle_min,&
       branch_fraction,length_limit,shortest_length,rotation_limit,to_export,filename)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GROW_TREE" :: GROW_TREE

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER,INTENT(in)  :: parent_ne                ! list of end branch elements to grow from
    INTEGER,INTENT(in)  :: surface_elems(:)         ! list of surface elements defining the host region
    REAL(dp),INTENT(in) :: angle_max                ! maximum branch angle with parent; in degrees
    REAL(dp),INTENT(in) :: angle_min                ! minimum branch angle with parent; in degrees
    REAL(dp),INTENT(in) :: branch_fraction          ! fraction of distance (to COFM) to branch
    REAL(dp),INTENT(in) :: length_limit             ! minimum length of a generated branch (shorter == terminal)
    REAL(dp),INTENT(in) :: shortest_length          ! length that short branches are reset to (shortest in model)
    REAL(dp),INTENT(in) :: rotation_limit           ! maximum angle of rotation of branching plane
    LOGICAL,INTENT(in) :: to_export                 ! option to export terminal element mapping to datapoints
    CHARACTER(len=*),INTENT(in) :: filename

    !Local variables
    INTEGER,ALLOCATABLE :: local_parent(:)          ! stores current generation of local parent elements
    INTEGER,ALLOCATABLE :: local_parent_temp(:)     ! temporary storage of next generation of local parent elems
    INTEGER,ALLOCATABLE :: map_seed_to_elem(:)      ! records current elem associated w. data points
    INTEGER,ALLOCATABLE :: map_seed_to_space(:)     ! records initial elem associated w. data points (the 'space')
    INTEGER,ALLOCATABLE :: num_seeds_from_elem(:)   ! records # of seeds currently grouped with an elem
    INTEGER,ALLOCATABLE :: triangle(:,:)
    CHARACTER(len=100) :: writefile

    INTEGER :: i,j,kount,M,N,nd,nd_min,ne,ne_grnd_parent,ne_parent,ne_stem,&
         noelem_parent,np,np_start,np_prnt_start,np_grnd_start,num_seeds_in_space,num_next_parents, &
         num_parents,num_triangles,num_vertices,num_elems_new,num_nodes_new,num_terminal

    REAL(dp),ALLOCATABLE :: vertex_xyz(:,:)
    REAL(dp),DIMENSION(3) :: COFM,candidate_xyz
    REAL(dp) :: distance_limit = 300.0_dp,length_parent

    LOGICAL :: make_branch,enough_points,first_group,internal, &
         limit_branching_angle = .TRUE., &  ! option to restrict branch angle
         limit_branching_plane = .FALSE.    ! option to restrict angle between branching planes

    CHARACTER(len=60) :: sub_name

    sub_name = 'grow_tree'
    CALL enter_exit(sub_name,1)


    IF(to_export)THEN
!!! export vertices as nodes
       writefile = TRIM(filename)//'.txt'
       OPEN(40, file = writefile, status='replace')
       WRITE(40,'('' Data point number          Terminal element number'')')
    ENDIF



    CALL triangles_from_surface(num_triangles,num_vertices,surface_elems,triangle,vertex_xyz)

!!! We can estimate the number of elements in the generated model based on the
!!! number of data (seed) points. i.e. N = 2*N_data - 1. So the total number of
!!! elements following tree generation will be ~ num_elems + 2*num_data. Use this estimate
!!! to increase the node and element arrays.
    num_elems_new = num_elems + 2*num_data + 100
    num_nodes_new = num_nodes + 2*num_data + 100
    CALL reallocate_node_elem_arrays(num_elems_new,num_nodes_new)

!!! Allocate memory for temporary arrays (need a more intelligent way of estimating size!)
    ALLOCATE(local_parent_temp(num_elems_new))
    ALLOCATE(local_parent(num_elems_new))
    ALLOCATE(num_seeds_from_elem(num_elems_new))
    ALLOCATE(map_seed_to_elem(num_data))
    ALLOCATE(map_seed_to_space(num_data))

!!! Initialise local_parent to the list of parent elements, and num_parents (current
!!! number of parent branches) to the number of parent branches.
    local_parent(1:SIZE(parentlist)) = parentlist(1:SIZE(parentlist))
    num_parents = COUNT(parentlist.NE.0) !initial number of 'terminal' parent branches

    NUM_SEEDS_FROM_ELEM = 0
    num_next_parents = num_parents


!!! Calculate the initial grouping of data points with terminal elements
!!! this defines the 'space' with which each seed is associated
!!! For a single parent, all seed points will initially be mapped to
!!! it; for multiple parents 'group_seeds_with_branch' used to be called to calculate
!!! the closest parent end-point to each seed point. This has been replaced by splitting
!!! seed points using the orthogonal to branching planes of the upper tree.
    map_seed_to_space(1:num_data) = parentlist(1) !#! this is done for the new-style growing (full grow per terminal)
    IF(num_parents.GT.1)THEN
       first_group = .TRUE.
       !       call group_seeds_with_branch(map_seed_to_space,num_next_parents,num_seeds_from_elem,&
       !            num_terminal,local_parent,500.0_dp,first_group)
       CALL split_seed_points_initial(map_seed_to_space,parent_ne)
    ENDIF !parentlist.gt.1
    first_group = .FALSE.

    WRITE(*,'(''  parent  #seeds  #terminal'')')

    ! Set initial values for local and global nodes and elements
    ne = num_elems !initialise mesh global element #
    np = num_nodes !initialise mesh global node #

!!! loop over the initial parent list (the initial conditions/terminal elements for growing)
!!! growing is done into 'spaces', where each 'space' is the initial grouping of seed
!!! points with the closest terminal branch. This grouping can be manipulated to take into account
!!! the size of the initial terminal branches. e.g. smaller diameter --> smaller set of seeds


    DO noelem_parent = 1,num_parents
       ne_stem = parentlist(noelem_parent) ! the 'stem' parent element for the 'space'
       map_seed_to_elem = 0 ! initialise the seed mapping array
       num_seeds_in_space = 0 !initialise the number of seed points in the 'space'
       DO nd = 1,num_data ! for all of the seed points (stored in data_xyz array)
          IF(map_seed_to_space(nd).EQ.ne_stem)THEN ! for the points in this space
             map_seed_to_elem(nd) = ne_stem ! record the current element associated with seed point nd
             num_seeds_in_space = num_seeds_in_space+1 ! count number of seed points in the space
          ENDIF
       ENDDO

       num_next_parents = 1 ! initialise the number of current local parent branches
       local_parent(1) = ne_stem ! first local parent branch is the 'stem' branch
       num_terminal = 0 ! initialise the number of definite terminal branches

!!! bifurcating distributive algorithm
       DO WHILE(num_next_parents.NE.0) !while still some parent branches with seed points
          num_parents = num_next_parents ! update the number of current local parent branches
          num_next_parents = 0 ! reset the number of local parent branches in next generation

          DO M = 1,num_parents ! for each of the current local parent branches
             ne_parent = local_parent(M) !parent element #
             ! Calculate centre of mass of current seed point set
             CALL calculate_seed_cofm(map_seed_to_elem,ne_parent,COFM)

             ne_grnd_parent = elem_cnct(-1,1,ne_parent) !grandparent global element #
             np_start = elem_nodes(2,ne_parent) !parent global end node #
             np_prnt_start = elem_nodes(1,ne_parent) !parent global start node #
             np_grnd_start = elem_nodes(1,ne_grnd_parent) !grandparent global start node #

             length_parent = elem_field(ne_length,ne_parent)

!!! Split each set of seed points using the plane defined by the
!!! parent branch and the centre of mass. Seed points get associated with NEW elements (ne+1,ne+2)
             CALL split_seed_points(map_seed_to_elem,ne_parent,ne,np_start,&
                  np_prnt_start,np_grnd_start,COFM,enough_points)

!!! check whether enough seed points remaining in BOTH seed groups for branching to be done
!!! (note: this could be improved to continue branching in one set of seeds)
             IF(enough_points)THEN
                DO N = 1,2 !for each of the two new branches
                   ! Set up arrays for new element and node
                   ! after create_new_node the current element == ne and current node == np
                   CALL create_new_node(ne,ne_parent,np,np_start,.TRUE.)
                   ! find the centre of mass of seed points
                   CALL calculate_seed_cofm(map_seed_to_elem,ne,COFM)
                   ! Generate a branch directed towards the centre of mass. Returns location
                   ! of end node in candidate_xyz (adjusted below based on length and shape criteria)
                   CALL branch_to_cofm(map_seed_to_elem,ne,np_start,&
                        COFM,branch_fraction,length_limit,length_parent,shortest_length,&
                        candidate_xyz,make_branch)
                   node_xyz(1:3,np) = candidate_xyz(1:3) ! the new node location is as returned by 'branch_to_cofm'
                   CALL calc_branch_direction(ne) ! calculate direction of the new branch
                   elem_field(ne_length,ne) = distance_between_points(node_xyz(1,np_start),node_xyz(1,np))
                   ! Check whether this is a new parent branch or a terminal branch
                   IF(make_branch.AND.enough_points)THEN ! meets all criteria for continuing branching
                      num_next_parents = num_next_parents+1 ! increment the number of next parents
                      local_parent_temp(num_next_parents) = ne !records the elements that are parents
                   ELSE ! this is a terminal branch
                      num_terminal = num_terminal+1 ! increment the number of terminal branches
                   ENDIF
                ENDDO !N (for both new branches)

                IF(limit_branching_angle)THEN
                   ! Check that the branch angles are not too large or too small
                   ! Correct such that branches stay in the original branching plane
                   CALL limit_branch_angles(ne,ne_parent,np,&
                        np_prnt_start,np_start,angle_max,angle_min)

                   internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np),&
                        vertex_xyz)
                   IF(.NOT.internal)THEN ! halve the length, and make terminal
                      elem_field(ne_length,ne) = 0.5_dp*distance_between_points&
                           (node_xyz(1,np_start),node_xyz(1,np))
                      node_xyz(1:3,np) = node_xyz(1:3,np_start) + 0.5_dp*&
                           elem_field(ne_length,ne)*elem_direction(1:3,ne)
                      internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np),&
                           vertex_xyz)
                      kount = 0
                      DO WHILE(.NOT.internal)
                         kount = kount+1
                         CALL shorten_branch_and_children(ne_parent)
                         internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np),&
                              vertex_xyz)
                         IF(kount.GE.3)THEN
                            CALL shorten_branch_and_children(elem_cnct(-1,1,ne_parent))
                            internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np),&
                                 vertex_xyz)
                         ENDIF
                         IF(kount.GT.5.AND.(.NOT.internal))THEN
                            WRITE(*,'('' WARNING: element'',i6,'' not internal'')') ne
                            internal = .TRUE.
                         ENDIF
                      ENDDO
                      IF(inlist(ne,local_parent_temp))THEN ! set to be terminal
                         local_parent_temp(num_next_parents) = 0
                         num_next_parents = num_next_parents-1 ! decrement the number of next parents
                         num_terminal = num_terminal+1
                         nd_min = closest_seed_to_node_in_group(map_seed_to_elem,ne,np) ! closest seed point
                         map_seed_to_elem(nd_min) = 0 ! remove seed point from list
                         map_seed_to_space(nd_min) = ne ! recording element number

                         IF(to_export) THEN
                            WRITE(40,*) nd_min,ne
                         ENDIF

                      ENDIF
                   ENDIF !.not.internal

                   internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np-1),&
                        vertex_xyz)
                   IF(.NOT.internal)THEN ! halve the length, and halve the angle from parent
                      elem_field(ne_length,ne-1) = 0.5_dp*distance_between_points&
                           (node_xyz(1,np_start),node_xyz(1,np-1))
                      node_xyz(1:3,np-1) = node_xyz(1:3,np_start) + 0.5_dp*&
                           elem_field(ne_length,ne-1)*elem_direction(1:3,ne-1)
                      internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np-1),&
                           vertex_xyz)
                      DO WHILE(.NOT.internal)
                         kount = kount+1
                         CALL shorten_branch_and_children(ne_parent)
                         internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np),&
                              vertex_xyz)
                         IF(kount.GE.3)THEN
                            CALL shorten_branch_and_children(elem_cnct(-1,1,ne_parent))
                            internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np),&
                                 vertex_xyz)
                         ENDIF
                         IF(kount.GT.4.AND.(.NOT.internal))THEN
                            WRITE(*,'('' WARNING: element'',i6,'' not internal'')') ne-1
                            internal = .TRUE.
                         ENDIF
                      ENDDO
                      IF(inlist(ne-1,local_parent_temp))THEN ! set to be terminal
                         DO i = 1,num_next_parents
                            IF(local_parent_temp(i).EQ.ne-1)THEN
                               DO j = i,num_next_parents-1
                                  local_parent_temp(j) = local_parent_temp(i+1)
                               ENDDO
                            ENDIF
                         ENDDO
                         num_next_parents = num_next_parents-1 ! decrement the number of next parents
                         num_terminal = num_terminal+1
                         nd_min = closest_seed_to_node_in_group(map_seed_to_elem,ne-1,np-1) ! closest seed point
                         map_seed_to_elem(nd_min) = 0 ! remove seed point from list
                         map_seed_to_space(nd_min) = ne ! recording element number

                         IF(to_export) THEN
                            WRITE(40,*) nd_min,ne
                         ENDIF
                      ENDIF
                   ENDIF ! .not.internal

                ENDIF

                IF(limit_branching_plane)THEN
                   ! Check the angle of rotation between child and parent branching planes.
                   ! If absolute angle is larger than a user-specified limit, adjust to be
                   ! the limit value. This is used for making sure that the CFD geometry turns out ok.
                   CALL check_branch_rotation_plane(map_seed_to_elem,ne,ne_grnd_parent,ne_parent, &
                        local_parent_temp,num_next_parents, &
                        np,np_start,np_prnt_start,np_grnd_start,num_terminal,rotation_limit)
                   internal = point_internal_to_surface(num_vertices,triangle,node_xyz(1:3,np),vertex_xyz)
                ENDIF

             ELSE
                WRITE(*,*) 'terminal, not enough points',ne !!! never happens!!!
                READ(*,*)
                ! Not enough seed points in the set during the split.
                ! Find the closest seed point to node np_start, and remove from seeds
                num_terminal = num_terminal+1 ! increment number of terminal branches
                nd_min = closest_seed_to_node(map_seed_to_elem,np_start) ! closest seed point
                map_seed_to_elem(nd_min) = 0 ! remove seed point from list
                map_seed_to_space(nd_min) = ne ! record the element to data point mapping
             ENDIF
          ENDDO ! for all current parent branches
          ! Copy the temporary list of branches to local_parent. These become the
          ! parent elements for the next branching
          local_parent(1:num_next_parents) = local_parent_temp(1:num_next_parents)
          ! Regroup the seed points with the closest current parent
          CALL group_seeds_with_branch(map_seed_to_elem,num_next_parents,num_seeds_from_elem,&
               num_terminal,local_parent,DISTANCE_LIMIT,.FALSE.)

       ENDDO ! while still parent branches

       WRITE(*,'(I7,I8,I9)') ne_parent,num_seeds_in_space,num_terminal

    ENDDO ! for each initial parent
    CLOSE(40)


!!! set new total numbers of nodes and elements
    num_nodes=np !highest node # in nr
    num_elems=ne !highest element # in nr

!!! update the tree connectivity
    CALL element_connectivity_1d
!!! calculate branch generations and orders
    CALL evaluate_ordering
!!! deallocate temporary arrays
    DEALLOCATE(vertex_xyz)
    DEALLOCATE(triangle)
    DEALLOCATE(local_parent_temp)
    DEALLOCATE(local_parent)
    DEALLOCATE(map_seed_to_elem)
    DEALLOCATE(map_seed_to_space)
    DEALLOCATE(num_seeds_from_elem)

    CALL enter_exit(sub_name,2)
  END SUBROUTINE grow_tree


  !###############################################################
  !
  !*limit_branch_angles:* checks both new branches at once to make sure the angles
  ! to parent branch not too large or too small. Both branches are checked to make sure
  ! they are not co-linear with parent. One branch is checked for its angle, and then
  ! the other is made to be in-plane with the parent and sibling branch.
  !
  SUBROUTINE limit_branch_angles(ne,ne_parent,np,&
       np_prnt_start,np_start,angle_max,angle_min)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_LIMIT_BRANCH_ANGLES" :: LIMIT_BRANCH_ANGLES

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER,INTENT(in) :: ne,ne_parent,np,np_prnt_start,np_start
    REAL(dp) :: angle_max,angle_min

    !Local variables
    INTEGER :: elem_check,elem_in_plane,node_check,node_in_plane
    REAL(dp),DIMENSION(3) :: B,NRML,U,V,W
    REAL(dp),DIMENSION(3,3) :: A
    REAL(dp) :: angle_uw,angle_uv,angle_vw,length,length_w

    CHARACTER(len=60) :: sub_name

    sub_name = 'limit_branch_angles'
    CALL enter_exit(sub_name,1)

    elem_check = ne-1  ! the default element to check angle for is the '1st' branch
    node_check = np-1  ! node associated with ne-1
    elem_in_plane = ne ! the default element to make in-plane is the '2nd' branch
    node_in_plane = np ! node associated with ne
    ! check whether default element is co-linear with parent element
    u(1:3) = elem_direction(1:3,ne_parent)  ! parent branch direction
    v(1:3) = elem_direction(1:3,ne-1)       ! new branch direction
    length = elem_field(ne_length,ne-1)
    nrml = cross_product(u,v) ! calculate normal
    IF(ABS(nrml(1))+ABS(nrml(2))+ABS(nrml(3)).LT.1.0e-4)THEN !co-linear branches
       ! check whether the other branch is co-linear with parent
       v(1:3) = elem_direction(1:3,ne)
       nrml = cross_product(u,v) ! calculate normal
       IF(ABS(nrml(1))+ABS(nrml(2))+ABS(nrml(3)).LT.1.0e-4)THEN !co-linear points
          WRITE(*,*) 'warning: both branches co-linear with parent'
          READ(*,*)
       ELSE
          elem_check = ne
          node_check = np
          elem_in_plane = ne-1
          node_in_plane = np-1
       ENDIF
    ENDIF
    CALL adjust_branch_angle(1,elem_check,np_start,np_prnt_start,node_check, &
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
    CALL adjust_branch_angle(2,elem_in_plane,np_start,np_prnt_start,node_in_plane,angle_max,angle_min)

    CALL enter_exit(sub_name,2)

  END SUBROUTINE limit_branch_angles


  !##################################################
  !
  !*reduce_branch_angle:* calculates the direction of a branch for a given branch angle

  SUBROUTINE reduce_branch_angle(np1,np2,np,candidate_xyz,factor)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_REDUCE_BRANCH_ANGLE" :: REDUCE_BRANCH_ANGLE

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER,INTENT(in) :: np1,np2,np
    REAL(dp),INTENT(in) :: factor
    REAL(dp) :: candidate_xyz(3)

    !Local variables
    REAL(dp) :: angle,LV,U(3),V(3),W(3)

    CHARACTER(len=60) :: sub_name

    sub_name = 'reduce_branch_angle'
    CALL enter_exit(sub_name,1)

    U(1:3)=node_xyz(1:3,np1)-node_xyz(1:3,np2) !direction of parent
    V(1:3)=candidate_xyz(1:3)-node_xyz(1:3,np1) !direction of this branch
    LV = vector_length(V)
    U = unit_vector(U)
    V = unit_vector(V)
    angle = angle_btwn_vectors(U,V)

    W = vector_for_angle_limit(U,V,angle,factor*angle)

    node_xyz(1:3,np) = node_xyz(1:3,np1) + W(1:3)*LV ! use original length
    candidate_xyz(1:3) = node_xyz(1:3,np) ! adjust the candidate node location

    CALL enter_exit(sub_name,2)

  END SUBROUTINE reduce_branch_angle


  !###############################################################
  !
  !*shorten_branch_and_children:* shorten the specified branch, its children,
  ! and their children
  !
  SUBROUTINE shorten_branch_and_children(ne)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SHORTEN_BRANCH_AND_CHILDREN" :: SHORTEN_BRANCH_AND_CHILDREN

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER,INTENT(in) :: ne

    !Local variables
    INTEGER :: i,j,ne1,ne2,np0,np1,np2
    CHARACTER(len=60) :: sub_name = 'shorten_branch_and_children'

    CALL enter_exit(sub_name,1)

    np0 = elem_nodes(1,ne)
    np1 = elem_nodes(2,ne)

    elem_field(ne_length,ne) = 0.8_dp*elem_field(ne_length,ne)
    node_xyz(1:3,np1) = node_xyz(1:3,np0)+elem_field(ne_length,ne)*elem_direction(1:3,ne)

    np0 = np1
    DO i=1,elem_cnct(1,0,ne) ! for each child
       ne1 = elem_cnct(1,i,ne)
       np1 = elem_nodes(2,ne1)
       elem_field(ne_length,ne1) = 0.8_dp*elem_field(ne_length,ne1)
       node_xyz(1:3,np1) = node_xyz(1:3,np0)+elem_field(ne_length,ne1)*elem_direction(1:3,ne1)
       DO j = 1,elem_cnct(1,0,ne1)
          ne2 = elem_cnct(1,j,ne1)
          np2 = elem_nodes(2,ne2)
          elem_field(ne_length,ne2) = 0.8_dp*elem_field(ne_length,ne2)
          node_xyz(1:3,np2) = node_xyz(1:3,np1)+elem_field(ne_length,ne2)*elem_direction(1:3,ne2)
       ENDDO
    ENDDO

    CALL enter_exit(sub_name,2)

  END SUBROUTINE shorten_branch_and_children


  !###############################################################
  !
  !*smooth_1d_tree:* smooth a tree geometry by placing the end of a parent branch
  ! at the average of the parent end and child end coordinates. This is used to
  ! improve the topology of generated trees, minimising the impact of 'odd' branching
  !
  SUBROUTINE smooth_1d_tree(num_elem_start,length_limit)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SMOOTH_1D_TREE" :: SMOOTH_1D_TREE

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER,INTENT(in) :: num_elem_start
    REAL(dp),INTENT(in) :: length_limit

    INTEGER :: n,ne,ne1,ne2,np,np0,np1,np2,n_smoothing_steps = 2
    REAL(dp) :: new_xyz(3)
    CHARACTER(len=60) :: sub_name

    sub_name = 'smooth_1d_tree'
    CALL enter_exit(sub_name,1)

    DO n = 1,n_smoothing_steps
       DO ne = num_elems,num_elem_start,-1
          IF(elem_cnct(1,0,ne).EQ.0)THEN ! terminal, check branch length
             IF(elem_field(ne_length,ne).LT.0.75_dp*length_limit)THEN
                elem_field(ne_length,ne) = 0.75_dp*length_limit
                np1 = elem_nodes(1,ne) ! the start node
                np2 = elem_nodes(2,ne) ! the end node
                node_xyz(:,np2) = node_xyz(:,np1) + elem_direction(:,ne)*0.75_dp*length_limit
             ELSE IF(elem_field(ne_length,ne).GT.1.5_dp*length_limit)THEN
                elem_field(ne_length,ne) = 1.5_dp*length_limit
                np1 = elem_nodes(1,ne) ! the start node
                np2 = elem_nodes(2,ne) ! the end node
                node_xyz(:,np2) = node_xyz(:,np1) + elem_direction(:,ne)*1.5_dp*length_limit
             ENDIF
          ELSE
             IF(elem_cnct(1,0,ne).EQ.2)THEN
                ne1 = elem_cnct(1,1,ne)
                ne2 = elem_cnct(1,2,ne)
                np0 = elem_nodes(1,ne)
                np  = elem_nodes(2,ne)
                np1 = elem_nodes(2,ne1)
                np2 = elem_nodes(2,ne2)
                new_xyz(:) = node_xyz(:,np0)*0.5_dp + node_xyz(:,np1)*0.25_dp + node_xyz(:,np2)*0.25_dp
                node_xyz(:,np) = new_xyz(:)
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    CALL enter_exit(sub_name,2)

  END SUBROUTINE smooth_1d_tree


  !###############################################################
  !
  !*split_seed_points:* divides a set of seed points into two subsets
  ! using the plane that contains the parent branch and the seed point
  ! centre of mass. Decides which side of a plane a seed point is on by calculating
  ! the distance between two parallel planes: one which is defined by
  ! the parent and grandparent branch, and the other which contains a
  ! seed point.
  !
  SUBROUTINE split_seed_points(map_seed_to_elem,ne1,ne,np1,np2,np3,COFM,enough_points)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SPLIT_SEED_POINTS" :: SPLIT_SEED_POINTS

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER :: map_seed_to_elem(*),ne,ne1,np1,np2,np3
    REAL(dp) :: COFM(3)
    LOGICAL :: enough_points

    !Local variables
    INTEGER :: DAT1,DAT2,nd,ND1_1ST,ND2_1ST,nsp,NPOINTS
    REAL(dp) :: DIST,NORML(4),P(3),Q(3),R(3)
    LOGICAL :: COLINEAR

    CHARACTER(len=60) :: sub_name

    sub_name = 'split_seed_points'
    CALL enter_exit(sub_name,1)

    R = COFM ! split based on cofm and branch
    P(1:3) = node_xyz(1:3,np2) ! point at start of parent branch
    Q(1:3) = node_xyz(1:3,np1) ! point at end of parent branch

!!! check whether the centre of mass and the parent start & end branches
!!! are co-linear. if so, will need to use 'aunt' branch for split
    colinear = check_colinear_points(P,Q,R)
    IF(colinear) R(1:3) = node_xyz(1:3,np3) !split based on parent and aunt
    CALL make_plane_from_3points(NORML,1,P,Q,R) !calculate plane

    NPOINTS=0
    DAT1=0
    DAT2=0
    ND1_1ST=0
    ND2_1ST=0
    DO nd=1,num_data
       nsp=map_seed_to_elem(nd) !space # that random point belongs to
       IF(nsp.EQ.ne1)THEN !random point belongs to this element space
          NPOINTS=NPOINTS+1
          dist = -scalar_product_3(norml,data_xyz(1,nd)) - norml(4) ! distance between two planes
          IF(dist.GE.0.0_dp)THEN
             IF(dat1.EQ.0) nd1_1st = nd
             DAT1=DAT1+1
             map_seed_to_elem(nd)=ne+1
          ELSE IF(DIST.LT.0.0_dp)THEN
             IF(DAT2.EQ.0) ND2_1ST=nd
             DAT2=DAT2+1
             map_seed_to_elem(nd)=ne+2
          ENDIF
       ENDIF
    ENDDO !nd

    IF(dat1.EQ.0.AND.dat2.EQ.0)THEN
       enough_points = .FALSE.
       WRITE(*,'('' Zero seed points associated with parent'',I6)') ne1
       READ(*,*)
    ELSE
       IF(dat1.EQ.0)THEN
          map_seed_to_elem(nd2_1st) = ne+1
          dat1 = dat1+1
          dat2 = dat2-1
       ELSEIF(dat2.EQ.0)THEN
          map_seed_to_elem(nd1_1st) = ne+2
          dat2 = dat2+1
          dat1 = dat1-1
       ENDIF

       enough_points = .TRUE.
    ENDIF

    CALL enter_exit(sub_name,2)

  END SUBROUTINE split_seed_points


  !###############################################################
  !
  !*split_seed_points_initial:* divides a set of seed points into N subsets
  ! to match N terminal branches, using the plane that is orthogonal to the branching plane
  ! of child branches, and that passes mid-way between child branches.
  !
  SUBROUTINE split_seed_points_initial(map_array,ne_stem)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SPLIT_SEED_POINTS_INITIAL" :: SPLIT_SEED_POINTS_INITIAL

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER :: map_array(:),ne_stem

    !Local variables
    INTEGER :: local_parent(20),local_parent_temp(20),M,nd,ne_parent, &
         ne1,ne2,np0,np1,np2,nsp,num_points,num_next_parents,num_parents
    REAL(dp) :: DIST,dist_p1,dist_p2,NORML(4),P(3),Q(3),R(3)

    CHARACTER(len=60) :: sub_name

    sub_name = 'split_seed_points_initial'
    CALL enter_exit(sub_name,1)

    ne_parent = ne_stem
    DO WHILE(elem_cnct(1,0,ne_parent).EQ.1)
       ne_parent = elem_cnct(1,1,ne_parent) ! get the next element in a refined branch
    ENDDO
    map_array(:) = ne_parent ! initialise that all seed points map to the stem branch

    num_next_parents = 1
    local_parent(1) = ne_parent

    DO WHILE(num_next_parents.NE.0) !while still some parent branches with seed points
       num_parents = num_next_parents ! update the number of current local parent branches
       num_next_parents = 0 ! reset the number of local parent branches in next generation
       DO M = 1,num_parents ! for each of the current local parent branches
          ne_parent = local_parent(M) !parent element #
          DO WHILE(elem_cnct(1,0,ne_parent).EQ.1)
             ne_parent = elem_cnct(1,1,ne_parent) ! get the next element in a refined branch
          ENDDO
          np0 = elem_nodes(2,ne_parent)
          ne1 = elem_cnct(1,1,ne_parent)
          DO WHILE(elem_cnct(1,0,ne1).EQ.1)
             ne1 = elem_cnct(1,1,ne1) ! get the next element in a refined branch
          ENDDO
          ne2 = elem_cnct(1,2,ne_parent)
          DO WHILE(elem_cnct(1,0,ne2).EQ.1)
             ne2 = elem_cnct(1,1,ne2) ! get the next element in a refined branch
          ENDDO
          np1 = elem_nodes(2,ne1)
          np2 = elem_nodes(2,ne2)
          P(:) = node_xyz(:,np0) ! point at end of parent branch
          Q(:) = node_xyz(:,np1) ! point at end of child1 branch
          R(:) = node_xyz(:,np2) ! point at end of child2 branch
          CALL make_plane_from_3points(NORML,1,P,Q,R) !calculate plane

          P(:) = node_xyz(:,np0) ! point at end of parent branch
          Q(:) = 0.5_dp*(node_xyz(:,np1)+node_xyz(:,np2))
          R(1:3) = Q(1:3) + NORML(1:3)
          CALL make_plane_from_3points(NORML,1,P,Q,R) !calculate plane
!!! NORML is now the plane between the child branches

          dist_p1 = -scalar_product_3(norml,node_xyz(:,np1)) - norml(4) ! distance between two planes
          dist_p2 = -scalar_product_3(norml,node_xyz(:,np2)) - norml(4) ! distance between two planes

          DO nd = 1,num_data
             nsp = map_array(nd) !space # that random point belongs to
             IF(nsp.EQ.ne_parent)THEN !random point belongs to this element space
                dist = -scalar_product_3(norml,data_xyz(1,nd)) - norml(4) ! distance between two planes
                IF(dist.GE.0.0_dp.AND.dist_p1.GE.0.0_dp)THEN
                   map_array(nd) = ne1
                ELSE IF(dist.GE.0.0_dp.AND.dist_p1.LT.0.0_dp)THEN
                   map_array(nd) = ne2
                ELSE IF(dist.LE.0.0_dp.AND.dist_p2.LE.0.0_dp)THEN
                   map_array(nd) = ne2
                ELSE IF(dist.LE.0.0_dp.AND.dist_p2.GT.0.0_dp)THEN
                   map_array(nd) = ne1
                ENDIF
             ENDIF
          ENDDO !nd

          num_points = COUNT(map_array.EQ.ne1)
          IF(num_points.EQ.0)THEN
             WRITE(*,'('' Warning: number of points for element'',i6,'' is zero'')') ne1
             WRITE(*,'('' Press enter to continue; however the code is likely to fail'')')
          ENDIF
          num_points = COUNT(map_array.EQ.ne2)
          IF(num_points.EQ.0)THEN
             WRITE(*,'('' Warning: number of points for element'',i6,'' is zero'')') ne2
             WRITE(*,'('' Press enter to continue; however the code is likely to fail'')')
          ENDIF

          IF(elem_cnct(1,0,ne1).NE.0)THEN
             num_next_parents = num_next_parents+1
             local_parent_temp(num_next_parents) = ne1
          ENDIF
          IF(elem_cnct(1,0,ne2).NE.0)THEN
             num_next_parents = num_next_parents+1
             local_parent_temp(num_next_parents) = ne2
          ENDIF
       ENDDO ! num_parents
       local_parent(1:num_next_parents) = local_parent_temp(1:num_next_parents)
    ENDDO

    CALL enter_exit(sub_name,2)

  END SUBROUTINE split_seed_points_initial


  !###############################################################
  !
  !*closest_seed_to_node:* returns the closest seed point to a given branch node
  !
  FUNCTION closest_seed_to_node(map_seed_to_elem,np)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CLOSEST_SEED_TO_NODE" :: CLOSEST_SEED_TO_NODE

    USE diagnostics, ONLY: enter_exit
    IMPLICIT NONE

    INTEGER,INTENT(in) :: map_seed_to_elem(*),np

    !Local variables
    INTEGER :: nd
    INTEGER :: closest_seed_to_node
    REAL(dp) :: distance,min_distance

    min_distance = 1.0e+10_dp
    DO nd=1,num_data
       IF(map_seed_to_elem(nd).NE.0)THEN
          distance = distance_between_points(data_xyz(1,nd),node_xyz(1,np))
          IF(distance.LT.min_distance)THEN
             closest_seed_to_node = nd
             min_distance = distance
          ENDIF !DIST
       ENDIF !map_seed_to_elem
    ENDDO !nd

  END FUNCTION closest_seed_to_node


  !###############################################################
  !
  !*closest_seed_to_node_in_group:* finds the closest seed point to a node
  ! that is in a group of seed points currently associated with a specific element
  !
  FUNCTION closest_seed_to_node_in_group(map_seed_to_elem,ne,np)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CLOSEST_SEED_TO_NODE_IN_GROUP" :: CLOSEST_SEED_TO_NODE_IN_GROUP

    INTEGER,INTENT(in) :: map_seed_to_elem(*),ne,np

    !Local variables
    INTEGER :: nd
    INTEGER :: closest_seed_to_node_in_group
    REAL(dp) :: distance,min_distance

    min_distance = 1.0e+10_dp
    DO nd=1,num_data
       IF(map_seed_to_elem(nd).EQ.ne)THEN
          distance = distance_between_points(data_xyz(1,nd),node_xyz(1,np))
          IF(distance.LT.min_distance)THEN
             closest_seed_to_node_in_group = nd
             min_distance = distance
          ENDIF !DIST
       ENDIF !map_seed_to_elem
    ENDDO !nd

  END FUNCTION closest_seed_to_node_in_group


  !###############################################################
  !
  !*rotation_angle:* calculates angle between two branching planes
  !
  FUNCTION rotation_angle(np1,np2,np3,np4,np5)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_ROTATION_ANGLE" :: ROTATION_ANGLE

    INTEGER,INTENT(in) :: np1,np2,np3,np4,np5

    !Local variables
    REAL(dp) :: norm_1(4),norm_2(4)
    REAL(dp) :: rotation_angle

    CALL make_plane_from_3points(norm_1,2,node_xyz(1,np1),node_xyz(1,np2),node_xyz(1,np3))
    CALL make_plane_from_3points(norm_2,2,node_xyz(1,np2),node_xyz(1,np4),node_xyz(1,np5))
    rotation_angle = angle_btwn_vectors(norm_1,norm_2)

  END FUNCTION rotation_angle


  !###############################################################
  !
  !*vector_for_angle_limit:* Calculates the new direction
  ! of an element, given current and target angles. Specifically, solves small
  ! system of equations to get new direction of branch (w), such that the branch
  ! remains in-plane (n.w = 0), the angle of w with u is defined
  ! (as u.w = cos(angle_with_u)), and the angle of w with original direction (v)
  ! is defined (as v.w = cos(angle_with_v)).
  !
  FUNCTION vector_for_angle_limit(U,V,angle_with_u,angle_with_v)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_VECTOR_FOR_ANGLE_LIMIT" :: VECTOR_FOR_ANGLE_LIMIT

    REAL(dp),INTENT(in) :: U(*),V(*),angle_with_u,angle_with_v

    !Local variables
    REAL(dp) :: A(3,3),N_UV(3),VECTOR(3),W(3)
    REAL(dp) :: vector_for_angle_limit(3)

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

  END FUNCTION vector_for_angle_limit


  !###############################################################

END MODULE growtree

!###############################################################

