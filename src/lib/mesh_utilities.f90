module mesh_utilities

!!! Subroutines and functions for general calculations. Not specific to any
!!! particular application, although it is expected these will generally be used 
!!! for calculations to do with a mesh (not fields or solutions). 
!!! Any function that is used by more than one module should appear in here. 
!!! ALL subroutines and functions in this module are public.


  use other_consts
  use arrays, only: dp,zero_tol

  implicit none

  private

  public  area_between_three_points,area_between_two_vectors,calc_branch_direction,&
       angle_btwn_points,angle_btwn_vectors,calc_scale_factors_2d,check_colinear_points,cross_product,&
       distance_between_points,make_plane_from_3points,mesh_a_x_eq_b,ph3,pl1,&
       point_internal_to_surface,scalar_product_3,scalar_triple_product,scale_mesh,&
       unit_norm_to_plane_two_vectors,unit_norm_to_three_points,unit_vector,&
       vector_length,volume_internal_to_surface

contains

!!! list of subroutines

  ! calc_branch_direction
  ! ..... calculates direction of branch and stores in elem_direction

  ! calc_scale_factors_2d
  ! ..... calculate the scale factors for a 2d mesh

  ! make_plane_from_3points
  ! ..... finds the equation of a plane in 3D and a vector normal to the plane from three
  ! ..... non-colinear points.

  ! scale_mesh
  ! ..... multiply mesh (coordinates and derivatives) by a constant

!!! list of functions

  ! angle_btwn_points
  ! .... returns the angle between three points

  ! angle_btwn_vectors
  ! .... returns the angle between two vectors
 
  ! check_colinear_points 
  ! .... returns true or false for whether 3 points are colinear

  ! cross_product ***
  ! .... returns the vector cross product of A*B

  ! distance_between_points ***
  ! .... calculates the distance between two arbitrary points

  ! mesh_a_x_eq_b ***
  ! .... solves a small matrix system

  ! scalar_product_3 ***
  ! .... dot product of two 3x1 vectors

  ! unit_vector ***
  ! .... Calculates the unit vector for an arbitrary 3x1 vector

  ! vector_length ***
  ! .... Calculates the length of a 3x1 vector

!!!#####################################################################
  
  subroutine calc_branch_direction(ne)
    
!!! calculates the direction of element ne and stores in elem_direction    
    
    use arrays,only: elem_direction,elem_nodes,node_xyz
    
    integer,intent(in) :: ne
    
    integer :: np_end,np_start
    real(dp) :: length
    
    np_start = elem_nodes(1,ne)
    np_end = elem_nodes(2,ne)
    
    length = distance_between_points(node_xyz(1,np_end),node_xyz(1,np_start))
    elem_direction(1:3,ne)=(node_xyz(1:3,np_end)-node_xyz(1:3,np_start))/length

  end subroutine calc_branch_direction

!!! ##########################################################################      

  subroutine calc_scale_factors_2d(sf_option)

!!! calculates the arclengths and scale factors for 2d surface elements,
!!! stores in scale_factors_2d
  
    use arrays,only: arclength,elem_lines_2d,elem_nodes_2d,lines_2d,lines_in_elem,&
         line_versn_2d,nodes_in_line,node_xyz_2d,num_elems_2d,num_lines_2d,scale_factors_2d

    character(len=4),intent(in) :: sf_option
!!! local variables
    integer,parameter :: num_deriv = 4
    integer :: ido(num_deriv,2),IG(4),it,IT_count,ITMAX=20,k,N,NAE,ne,&
         ng,NGA=4,NI1(3),ni,ni2,nj,nk,nk2,nl,nn,nn2,NNK,no_nl,&
         np,ns,nv,NNL(2,4)
    real(dp) :: DA,SUM1,SUM2,SUM3,SUM4,W,WG_LOCAL(10),XA_LOCAL(4,3),XI,&
         XIGG(10),XN_LOCAL(2,3,4)
    logical :: FOUND
    
    XIGG = [0.6_dp,0.2113248654051_dp,0.7886751345948_dp,0.1127016653792_dp,&
         0.6_dp,0.8872983346207_dp,0.0694318442029_dp,0.3300094782075_dp,&
         0.6699905217924_dp,0.9305681557970_dp]
    WG_LOCAL = [1.0_dp,0.6_dp,0.6_dp,0.2777777777778_dp,0.4444444444444_dp,&
         0.2777777777778_dp,0.1739274225687_dp,0.3260725774313_dp,&
         0.3260725774313_dp,0.1739274225687_dp]
    IG = [0,1,3,6]
    ido = reshape ([1,2,1,2,1,1,2,2],shape(ido))
    NI1 = [1,2,1]
    NNL = reshape([1,2,3,4,1,3,2,4],shape(NNL))
    
    if(.not.allocated(scale_factors_2d)) allocate(scale_factors_2d(16,num_elems_2d))

    select case (sf_option)
    case ('unit')
       scale_factors_2d = 1.0_dp
    case('arcl')
       do no_nl=1,num_lines_2d !loop over global lines
          nl=lines_2d(no_nl)
          ni = nodes_in_line(1,0,nl)
          do n = 1,2                  !for each node on the line
             np = nodes_in_line(n+1,1,nl)       !np1,np2
             do nj = 1,3
                nv = line_versn_2d(N,nj,nl)
                XN_LOCAL(1,nj,n) = node_xyz_2d(1,nv,nj,np)
                ne = lines_in_elem(1,nl)
                ni2 = 1+MOD(ni,2)
                nn = 1
                FOUND =.FALSE.
                do WHILE((nn.LE.4).AND.(.NOT.FOUND))
                   if(np.EQ.elem_nodes_2d(nn,ne))then
                      FOUND=.TRUE.
                   else
                      nn=nn+1
                   endif
                enddo
                do nk=2,4           !dxi1, dxi2, d2xi1xi2
                   if(IDO(nk,ni).EQ.2.AND.IDO(nk,ni2).EQ.1) then
                      XN_LOCAL(2,nj,n) = node_xyz_2d(nk,nv,nj,np)
                   endif
                enddo !nk
             enddo !nj
          enddo !n
          
          SUM2=0.0_dp
          do ng=1,NGA
             XI=XIGG(IG(NGA)+ng)
             W=WG_LOCAL(IG(NGA)+ng)
             do nj=1,3
                do k=1,2
                   XA_LOCAL(k,nj)=PL1(1,k,XI)*XN_LOCAL(1,nj,1) &
                        +PL1(2,k,XI)*XN_LOCAL(1,nj,2)
                enddo
             enddo
             
             SUM1=XA_LOCAL(2,1)**2+XA_LOCAL(2,2)**2+XA_LOCAL(2,3)**2
             SUM2=SUM2+W*DSQRT(SUM1)
          enddo !ng
          
          arclength(1:3,nl)=SUM2
          
          it=0
          iterative_loop : do
             it=it+1
             IT_count=it
             SUM3=0.0_dp
             SUM4=0.0_dp
             do ng=1,NGA
                XI=XIGG(IG(NGA)+ng)
                W=WG_LOCAL(IG(NGA)+ng)
                do nj=1,3
                   do k=1,2
                      XA_LOCAL(k,nj)=0.0_dp
                      do n=1,2
                         XA_LOCAL(k,nj)=XA_LOCAL(k,nj)+ &
                              PH3(n,1,k,XI)*XN_LOCAL(1,nj,n) &
                              +PH3(n,2,k,XI)*XN_LOCAL(2,nj,n)*arclength(n,nl)
                      enddo
                   enddo
                   XA_LOCAL(3,nj)=0.0_dp
                   do n=1,2
                      XA_LOCAL(3,nj)=XA_LOCAL(3,nj)+ &
                           PH3(n,2,2,XI)*XN_LOCAL(2,nj,n)
                   enddo
                   XA_LOCAL(4,nj)=0.0_dp
                   do n=1,2
                      XA_LOCAL(4,nj)=XA_LOCAL(4,nj)+ &
                           PH3(n,2,1,XI)*XN_LOCAL(2,nj,n)
                   enddo
                enddo
                SUM1=XA_LOCAL(2,1)**2+XA_LOCAL(2,2)**2+XA_LOCAL(2,3)**2
                SUM2=0.0_dp
                do nj=1,3
                   SUM2=SUM2+XA_LOCAL(2,nj)*XA_LOCAL(3,nj)
                enddo                  !nj
                SUM3=SUM3+W*DSQRT(SUM1)
                if(SUM1.GT.1.0e-6_dp) SUM4=SUM4+W*SUM2/DSQRT(SUM1)
             enddo                     !ng
             DA=-(arclength(3,nl)-SUM3)/(1.0_dp-SUM4)
             if(DABS(DA).GT.1.0e+6_dp) then
                arclength(3,nl)=1.0_dp
                exit iterative_loop
             endif
             
             arclength(3,nl) = arclength(3,nl)+DA      !is new arclength
             arclength(1:2,nl) = arclength(3,nl)
             
             if(it.eq.ITMAX) exit iterative_loop
             
          enddo iterative_loop      !iteration
       enddo !loop over lines
       
       scale_factors_2d = 1.0_dp !initialise
       
       do ne=1,num_elems_2d
          do NAE=1,4
             nl = elem_lines_2d(NAE,ne)
             if(nl /= 0)then
                ni = nodes_in_line(1,0,nl)
                ni2 = NI1(ni+1)
                do N=1,2
                   nn=NNL(N,NAE)
                   ns=0
                   do nn2=1,nn-1
                      do nk2=1,num_deriv
                         ns=ns+1
                      enddo
                   enddo
                   do nk=2,num_deriv
                      if(IDO(nk,ni2).EQ.1) then
                         scale_factors_2d(nk+ns,ne) = arclength(N,nl)
                         if(DABS(scale_factors_2d(nk+ns,ne)).LT.1.0e-6_dp) scale_factors_2d(nk+ns,ne) = 1.0_dp
                      endif
                   enddo !nk
                enddo !N=1,2
             endif
          enddo !NAE (nl)
          
          NNK=0
          do ns=1,4
             scale_factors_2d(NNK+4,ne)=scale_factors_2d(NNK+2,ne)*scale_factors_2d(NNK+3,ne)
             NNK=NNK+4
          enddo !nn
          
       enddo !noelem (ne)
    end select

  end subroutine calc_scale_factors_2d
  
!!!###############################################################
  
  subroutine make_plane_from_3points(NORML,NORMALTYPE,POINT1,POINT2,POINT3)
    
    !###    make_plane_from_3points finds the equation of a plane in three
    !###    dimensions and a vector normal to the plane from three
    !###    non-colinear points.
    !###    NORMALTYPE=1 for raw normal and plane equation
    !###    NORMALTYPE=2 for unit normal and plane equation
    !###    The coefficients represent aX + bY + cZ + d = 0
    !###    NORML(1)=a,NORML(2)=b,NORML(3)=c,NORML(4)=d
    
    
    !     Parameter list
    integer :: NORMALTYPE
    real(dp) :: POINT1(3),POINT2(3),POINT3(3),NORML(4)
    !     Local variables
    real(dp) :: DifF1(3),DifF2(3),NORMSIZE
    logical :: COLINEAR
    
    
    ! Check for colinearity
    COLINEAR=.FALSE.
    colinear = check_colinear_points(POINT1,POINT2,POINT3)
    if(.NOT.COLINEAR) then
       DifF1(1:3)=POINT2(1:3)-POINT1(1:3)
       DifF2(1:3)=POINT2(1:3)-POINT3(1:3)
       
       NORML(1)=(DifF1(2)*DifF2(3))-(DifF1(3)*DifF2(2))
       NORML(2)=(DifF1(3)*DifF2(1))-(DifF1(1)*DifF2(3))
       NORML(3)=(DifF1(1)*DifF2(2))-(DifF1(2)*DifF2(1))
       
       if(NORMALTYPE.EQ.2) then
          NORMSIZE = vector_length(NORML)
          NORML(1:3)=NORML(1:3)/NORMSIZE
       endif

       NORML(4) = -scalar_product_3(NORML,POINT1)
       
    else !Colinear
       
       WRITE(*,*) ' COLINEAR points in make_plane_from_3points '
       NORML = 0.0_dp
    endif
  end subroutine make_plane_from_3points
  
!!!##################################################

  subroutine scale_mesh(scaling,type)

    use arrays,only: node_xyz,node_xyz_2d,scale_factors_2d

    real(dp),intent(in) :: scaling
    character(len=2),intent(in) :: type

    select case(type)
    case('1d')
       node_xyz = node_xyz * scaling
    case('2d')
       node_xyz_2d = node_xyz_2d * scaling
       node_xyz_2d(4,:,:,:) = 0.0_dp
    end select
    
    scale_factors_2d = 1.0_dp

  end subroutine scale_mesh

!!!##################################################
  
  function area_between_two_vectors(vect_a,vect_b)
    
    !### 
    
    real(dp),intent(in) :: vect_a(3),vect_b(3)
    real(dp) :: cross(3)
    real(dp) :: area_between_two_vectors
    
    ! area = 1/2 x magnitude of the cross-product of vectors a and b
    
    cross = cross_product(vect_a,vect_b)
    area_between_two_vectors = 0.5_dp * sqrt(dot_product(cross,cross))
    
  end function area_between_two_vectors
  

!!!##################################################
  
  function area_between_three_points(point_a,point_b,point_c)
    
    !### 
    
    real(dp),intent(in) :: point_a(3),point_b(3),point_c(3)
    real(dp) :: norm(3),vect_a(3),vect_b(3)
    real(dp) :: area_between_three_points
    
    ! area = 1/2 x magnitude of the cross-product of vectors a and b
    
    vect_a(1:3) = point_a(1:3) - point_b(1:3)
    vect_b(1:3) = point_a(1:3) - point_c(1:3)
    norm = cross_product(vect_a,vect_b)
    area_between_three_points = 0.5_dp * sqrt(dot_product(norm,norm))
    
  end function area_between_three_points
  

!!! ##########################################################################      

  function ph3(I,J,K,XI)
    
!!! dummy arguments
    integer :: I,I_J_K,J,K
    real(dp) :: XI
!!! local variables
    real(dp) :: ph3
    
    ! K is 1,2, or 3; J is 1 or 2; I is 1 or 2
    
    I_J_K = 100*I + 10*J + K
    
    select case(I_J_K)
    case(111) !i=1,j=1,k=1
       PH3=(2.0_dp*XI-3.0_dp)*XI*XI+1.0_dp  ! 2xi^3-3xi^2+1
    case(121) !i=1,j=2,k=1
       PH3=((XI-2.0_dp)*XI+1.0_dp)*XI      ! xi^3-2xi^2+xi
    case(211) !i=2,j=1,k=1
       PH3=XI*XI*(3.0_dp-2.0_dp*XI)        ! -2xi^3+3xi^2
    case(221) !i=2,j=2,k=1
       PH3=XI*XI*(XI-1.0_dp)              ! xi^3-xi^2
    case(112) !i=1,j=1,k=2
       PH3=6.0_dp*XI*(XI-1.0_dp)           ! 6xi^2-6xi
    case(122) !i=1,j=2,k=2
       PH3=(3.0_dp*XI-4.0_dp)*XI+1.0_dp     ! 3xi^2-4xi+1
    case(212) !i=2,j=1,k=2
       PH3=6.0_dp*XI*(1.0_dp-XI)           ! -6xi^2+6xi
    case(222) !i=2,j=2,k=2
       PH3=XI*(3.0_dp*XI-2.0_dp)           ! 3xi^2-2xi
    case(113) !i=1,j=1,k=3
       PH3=12.0_dp*XI-6.0_dp               ! 12xi-6
    case(123) !i=1,j=2,k=3
       PH3=6.0_dp*XI-4.0_dp                ! 6xi-4
    case(213) !i=2,j=1,k=3
       PH3=6.0_dp-12.0_dp*XI               ! -12xi+6
    case(223) !i=2,j=2,k=3
       PH3=6.0_dp*XI-2.0_dp                ! 6xi-2
    end select
    
  end function ph3

!!! ##########################################################################      

  function pl1(I,K,XI)
    
!!! dummy arguments
    integer :: I,I_K,K
    real(dp) :: XI
!!! local variables
    real(dp) :: pl1
    
    I_K = 10*I + K
    
    select case(I_K)
    case(11) !i=1,k=1
       PL1=1.0_dp-XI
    case(21) !i=2,k=1
       PL1=XI
    case(12) !i=1,k=2
       PL1=-1.0_dp
    case(22) !i=2,k=2
       PL1=1.0_dp
    case(30 :) !k=3
       PL1=0.0_dp
    end select
    
    return
  end function pl1

!!!##################################################
  
  function unit_norm_to_plane_two_vectors(vect_a,vect_b)
    
    real(dp),intent(in) :: vect_a(3),vect_b(3)
    real(dp) :: magnitude,norm(3)
    real(dp) :: unit_norm_to_plane_two_vectors(3)
    
    norm = cross_product(vect_a,vect_b)
    magnitude = sqrt(dot_product(norm,norm))
    unit_norm_to_plane_two_vectors = norm/magnitude
    
  end function unit_norm_to_plane_two_vectors
  

!!!##################################################
  
  function unit_norm_to_three_points(point_a,point_b,point_c)

    real(dp),intent(in) :: point_a(3),point_b(3),point_c(3)
    real(dp) :: magnitude,norm(3),vect_a(3),vect_b(3)
    real(dp) :: unit_norm_to_three_points(3)

    vect_a(1:3) = point_a(1:3) - point_b(1:3)
    vect_b(1:3) = point_a(1:3) - point_c(1:3)
    norm = cross_product(vect_a,vect_b)
    magnitude = sqrt(dot_product(norm,norm))
    unit_norm_to_three_points = norm/magnitude

  end function unit_norm_to_three_points

!!!##################################################
  

  function angle_btwn_points(A,B,C)
    
    !###    calculates the angle between three points
    
    real(dp),intent(in) :: A(3),B(3),C(3)

    real(dp) :: U(3),V(3)
    real(dp) :: angle_btwn_points

    U = A - B
    V = C - B
    angle_btwn_points = angle_btwn_vectors(U,V)
        
  end function angle_btwn_points
  
!!!##################################################

  function angle_btwn_vectors(U,V)
    
    !###    ANGLE calculates the angle between two vectors
    
    real(dp),intent(in) :: U(3),V(3)

    real(dp) :: ANGLE,angle_btwn_vectors,N_U(3),N_V(3)
    
    N_U = unit_vector(U)
    N_V = unit_vector(V)

    ANGLE = scalar_product_3(N_U,N_V)
    ANGLE = max(-1.0_dp,ANGLE)
    ANGLE = min(1.0_dp,ANGLE)
    ANGLE = acos(ANGLE)

    angle_btwn_vectors=ANGLE
    
  end function angle_btwn_vectors
  
!!!###############################################################
  
  function check_colinear_points(POINT1,POINT2,POINT3)
    
    !###    check_colinear_points checks whether two vectors are colinear.
     
    !     Parameter list
    real(dp) :: POINT1(3),POINT2(3),POINT3(3)
    !     Local variables
    real(dp) :: ERR1(3),ERR2(3),LU,LV,U(3),V(3)
    logical :: check_colinear_points
    
    
    check_colinear_points =.FALSE.
    U(1:3)=POINT2(1:3)-POINT1(1:3)
    V(1:3)=POINT3(1:3)-POINT1(1:3)
    LU = vector_length(U)
    LV = vector_length(V)
    ! If 2 of the points are the same then LU and LV
    ! can be zero causing div by zero below and resulting in
    ! the wrong answer (on Linux) 
    if((dabs(LU)>zero_tol).AND.(dabs(LV)>zero_tol)) then
       ERR1(1:3)=DABS(U(1:3)/LU-V(1:3)/LV)
       ERR2(1:3)=DABS(U(1:3)/LU+V(1:3)/LV)
       if((ERR1(1).LE.ZERO_TOL.AND.ERR1(2).LE.ZERO_TOL.AND.ERR1(3).LE. &
            ZERO_TOL).OR.(ERR2(1).LE.ZERO_TOL.AND.ERR2(2).LE.ZERO_TOL.AND. &
            ERR2(3).LE.ZERO_TOL)) check_colinear_points=.TRUE.
    else
       check_colinear_points=.TRUE.
    endif
  end function check_colinear_points
  

!!!###############################################################
  
  function cross_product(A,B)
    
    !###  cross_product returns the vector cross product of A*B in C.
    
    !     Parameter List
    real(dp),intent(in) :: A(3),B(3)

    real(dp) :: cross_product(3)
    
    cross_product(1) = A(2)*B(3)-A(3)*B(2)
    cross_product(2) = A(3)*B(1)-A(1)*B(3)
    cross_product(3) = A(1)*B(2)-A(2)*B(1)
    
  end function cross_product
  
!!!###############################################################
  
  function scalar_triple_product(A,B,C)
    
    !###  scalar_triple_product returns A.(BxC)
    
    !     Parameter List
    real(dp),intent(in) :: A(3),B(3),C(3)

    real(dp) :: scalar_triple_product
    
    scalar_triple_product = A(1)*(B(2)*C(3)-B(3)*C(2)) + &
         A(2)*(B(3)*C(1)-B(1)*C(3)) + A(3)*(B(1)*C(2)-B(2)*C(1))
    
  end function scalar_triple_product
  
!!!###############################################################
  
  function distance_between_points(point1, point2)
    
    !###    calculates the distance between two arbitrary points
    
    real(dp),intent(in) :: point1(3),point2(3)
    integer :: i
    real(dp) :: distance_between_points
    
    distance_between_points = 0.0_dp
    do i=1,3
       distance_between_points = distance_between_points + (point1(i)-point2(i))**2
    enddo
    distance_between_points = dsqrt(distance_between_points)
    
  end function distance_between_points
  
!!!###############################################################
  
  function mesh_a_x_eq_b(MATRIX,VECTOR)
    
    real(dp) :: MATRIX(3,3),VECTOR(3)
    !Local variables
    integer :: i,j,k,pivot_row
    real(dp) :: A(3,4),max,pivot_value,SOLUTION(3),TEMP(4)
    real(dp) :: mesh_a_x_eq_b(3)
    
    
    A(1:3,1:3) = MATRIX(1:3,1:3)
    A(1:3,4) = VECTOR(1:3)
    do k=1,2
       max=0.0_dp
       do i=k,3
          if(DABS(A(i,k)).GT.max)then
             max=DABS(A(i,k))
             pivot_row=i
          endif
       enddo !i
       if(pivot_row.ne.k)then
          do j=1,4
             TEMP(j)=A(k,j)
             A(k,j)=A(pivot_row,j)
             A(pivot_row,j)=TEMP(j)
          enddo !j
       endif
       pivot_value = A(k,k)
       A(k,1:4) = A(k,1:4)/pivot_value
       do i=k+1,3
          do j=k+1,4
             A(i,j) = A(i,j)-A(i,k)*A(k,j)
          enddo
          A(i,k) = 0.0_dp
       enddo
    enddo !N
    A(3,4) = A(3,4)/A(3,3)
    A(2,4) = A(2,4)-A(3,4)*A(2,3)
    A(1,4) = A(1,4)-A(3,4)*A(1,3)-A(2,4)*A(1,2)

    SOLUTION(1:3) = A(1:3,4)
    
    mesh_a_x_eq_b = solution

  end function mesh_a_x_eq_b
  
!!!##################################################
  
  function scalar_product_3(A,B)
    
    !### calculates scalar product of two vectors A,B of length 3.
    
    real(dp),intent(in) :: A(*),B(*)

    integer :: i
    real(dp) :: scalar_product_3
    
    scalar_product_3 = 0.0_dp
    do i=1,3
       scalar_product_3 = scalar_product_3 + A(i)*B(i)
    enddo
    
  end function scalar_product_3
  
!!!###############################################################
  
  function unit_vector(A)
    
    !###  Calculates the unit vector for an arbitrary 3x1 vector 
    
    real(dp),intent(in) :: A(*)
    real(dp) :: length_a,unit_vector(3)

    length_a = vector_length(A)
    if(length_a.gt.1.0e-6_dp)then
       unit_vector(1:3) = A(1:3)/length_a
    else
       WRITE(*,*) ' >>WARNING: Cannot normalise a zero length vector'
       WRITE(*,*) ' We recommend debugging, but hit enter to continue'
       read(*,*)
    endif

  end function unit_vector

!!!##################################################
  
  function vector_length(A)
    
    !###  Calculates the length of a 3x1 vector 
    
    real(dp),intent(in) :: A(*)
    real(dp) :: vector_length
    integer :: i
    
    vector_length = 0.0_dp
    do i=1,3
       vector_length = vector_length + A(i)*A(i)
    enddo
    vector_length = dsqrt(vector_length)
    
  end function vector_length
  
!!!###############################################################

  function volume_internal_to_surface(triangles,vertex_xyz)

    ! calculates the volume enclosed by a list of surface elements

    integer,intent(in) :: triangles(:,:)
    real(dp),intent(in) :: vertex_xyz(:,:)
    real(dp) :: volume_internal_to_surface

!!! Local Variables
    integer :: ntri,num_triangles
    real(dp) :: volume,V1(3),V2(3),V3(3),P4(3)

    num_triangles = count(triangles(:,:).ne.0)/3

    P4 = sum(vertex_xyz,dim=2)/size(vertex_xyz,dim=2)

    volume = 0.0_dp

    do ntri = 1,num_triangles
       V1(1:3) = P4(1:3) - vertex_xyz(1:3,triangles(1,ntri))
       V2(1:3) = P4(1:3) - vertex_xyz(1:3,triangles(2,ntri))
       V3(1:3) = P4(1:3) - vertex_xyz(1:3,triangles(3,ntri))
       volume = volume + abs(scalar_triple_product(V1,V2,V3))
    enddo

    volume_internal_to_surface = volume/6.0_dp

  end function volume_internal_to_surface

!!!###############################################################

  function point_internal_to_surface(num_vertices,triangles,point_xyz,vertex_xyz)
!!! Cast a line in positive x-direction from each data point and 
!!! then work out how many triangular elements it crosses. If even it is in the 
!!! shape and if odd it is outside the shape

    integer,intent(in) :: num_vertices,triangles(:,:)
    real(dp),intent(in) :: point_xyz(3),vertex_xyz(:,:)
    logical :: point_internal_to_surface

!!! Local Variables
    integer :: i,ncrossed,ntri,num_triangles
    real(dp) :: area,area_triangle,cofm_surfaces(3),denominator,&
         norm_v(3),point(3),P1(3),P2(3),P3(3),u
    real(dp),parameter :: dist_tol = 1.0e-4_dp, user_tol = 1.0e-14_dp
    logical :: cross_any

    num_triangles = count(triangles(:,:).ne.0)/3

    forall (i=1:3) cofm_surfaces(i) = sum(vertex_xyz(i,1:num_vertices))/num_vertices
!    write(*,*) 'cofm',cofm_surfaces

! check whether the line that joins the centre of mass of the surface mesh and the point
! in question crosses ANY face. If it does, then point not inside.

    cross_any = .false.
    ncrossed = 0

    do ntri = 1,num_triangles
       P1(1:3) = vertex_xyz(1:3,triangles(1,ntri))
       P2(1:3) = vertex_xyz(1:3,triangles(2,ntri))
       P3(1:3) = vertex_xyz(1:3,triangles(3,ntri))
       norm_v = unit_norm_to_three_points(P1,P2,P3) ! unit normal to triangle plane
       ! u = (a*x1+b*y1+c*z1+d)/(a*(x1-x2)+b*(y1-y2)+c*(z1-z2))
       denominator = norm_v(1)*(point_xyz(1)-cofm_surfaces(1)) + &
            norm_v(2)*(point_xyz(2)-cofm_surfaces(2)) + &
            norm_v(3)*(point_xyz(3)-cofm_surfaces(3))
       ! denominator is zero for line parallel to plane
       if(abs(denominator).gt.user_tol)then
          ! calculate the distance of the surface point from point_xyz
          u = (dot_product(norm_v,point_xyz)-dot_product(norm_v,P1))/denominator
          if(u.ge.0.0_dp.and.u.le.1.0_dp)then ! POTENTIALLY crosses. Test further (angle)
             point = point_xyz + u*(cofm_surfaces-point_xyz) ! projection to surface
             area = area_between_two_vectors(P1-point,P2-point)+ &
                  area_between_two_vectors(P1-point,P3-point)+area_between_two_vectors(P2-point,P3-point)
             area_triangle = area_between_two_vectors(P1-P2,P1-P3)
             if(abs(area_triangle-area).lt.dist_tol)then
                cross_any = .true.
                ncrossed = ncrossed + 1
             endif
          endif
       endif
    enddo
    
    if(.not.cross_any)then
       point_internal_to_surface = .true.
    else
       if(ncrossed.eq.2)then
          point_internal_to_surface = .true.
       else
          point_internal_to_surface = .false.
       endif
    endif

  end function point_internal_to_surface



end module mesh_utilities
