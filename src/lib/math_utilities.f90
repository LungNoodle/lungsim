module math_utilities
!*Brief Description:* This module contains mathematical and numerical utilities
!
!*LICENSE:*
!
!
!
!*Full Description:*
!
!
  use arrays
  use diagnostics
  use other_consts
  use precision

  implicit none
  private

  public ax_cr,diagonal_pointer_cr,ilu_cr,lus_cr,mult_givens,rearrange_cr,bessel_complex
  public sort_integer_list
  public sort_real_list
  public angle_btwn_vectors,check_vectors_same,cross_product,inlist
  public mesh_a_x_eq_b,scalar_product_3,scalar_triple_product
  public unit_vector,vector_length


contains

  function angle_btwn_vectors(U,V)
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

!
!###########################################################################
!
  function check_vectors_same(vector1, vector2)
    real(dp) :: vector1(3),vector2(3)

    real(dp) :: norm_v1(3),norm_v2(3),u(3),v(3)
    logical :: check_vectors_same

    check_vectors_same = .false.
    norm_v1 = unit_vector(vector1)
    norm_v2 = unit_vector(vector2)
    u(1:3) = norm_v1(1:3) - norm_v2(1:3)
    v(1:3) = norm_v1(1:3) + norm_v2(1:3)

    if((abs(u(1))+abs(u(2))+abs(u(3)).lt.zero_tol).or. &
         (abs(v(1))+abs(v(2))+abs(v(3)).lt.zero_tol)) &
         check_vectors_same = .true.
  end function check_vectors_same

!
!###########################################################################
!
  function cross_product(A,B)
    real(dp),intent(in) :: A(3),B(3)
    real(dp) :: cross_product(3)

    cross_product(1) = A(2)*B(3)-A(3)*B(2)
    cross_product(2) = A(3)*B(1)-A(1)*B(3)
    cross_product(3) = A(1)*B(2)-A(2)*B(1)
  end function cross_product

!
!###########################################################################
!
  function inlist(item,ilist)
    integer :: item,ilist(:)
    integer :: n
    logical :: inlist

    inlist = .false.
    do n=1,size(ilist)
       if(item == ilist(n)) inlist = .true.
    enddo
  end function inlist

!
!###########################################################################
!
  function mesh_a_x_eq_b(MATRIX,VECTOR)
    real(dp) :: MATRIX(3,3),VECTOR(3)

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

!
!###########################################################################
!
  function scalar_product_3(A,B)
    real(dp),intent(in) :: A(*),B(*)

    integer :: i
    real(dp) :: scalar_product_3

    scalar_product_3 = 0.0_dp
    do i=1,3
       scalar_product_3 = scalar_product_3 + A(i)*B(i)
    enddo
  end function scalar_product_3

!
!###########################################################################
!
  function scalar_triple_product(A,B,C)
    real(dp),intent(in) :: A(3),B(3),C(3)
    real(dp) :: scalar_triple_product

    scalar_triple_product = A(1)*(B(2)*C(3)-B(3)*C(2)) + &
         A(2)*(B(3)*C(1)-B(1)*C(3)) + A(3)*(B(1)*C(2)-B(2)*C(1))
  end function scalar_triple_product

!
!###########################################################################
!
  function unit_vector(A)
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

!
!###########################################################################
!
  function vector_length(A)
    real(dp),intent(in) :: A(*)
    real(dp) :: vector_length
    integer :: i

    vector_length = 0.0_dp
    do i=1,3
       vector_length = vector_length + A(i)*A(i)
    enddo
    vector_length = dsqrt(vector_length)
  end function vector_length

  subroutine bessel_complex(z, bessel0, bessel1)
    use, intrinsic :: ieee_arithmetic
    implicit none

    complex(dp), intent(in)  :: z
    complex(dp), intent(out) :: bessel0, bessel1

    real(dp)    :: a(12), a1(10), b(12)
    real(dp)    :: a0, absz, scale
    complex(dp) :: cr, z1, ca, zr, zwork
    integer     :: k, k0

    ! -----------------------
    ! User-tunable safety caps
    ! -----------------------
    real(dp), parameter :: Z_SMALL   = 1.0e-12_dp     ! small-|z| threshold
    real(dp), parameter :: Z_ABS_MAX = 200.0_dp       ! cap |z| (keeps series/asymptotic stable)
    real(dp), parameter :: Z_RE_MAX  = 50.0_dp        ! cap Re(z) to avoid exp overflow / solver blow-up
    real(dp), parameter :: REL_TOL   = 1.0e-15_dp

    ! Coefficients (unchanged from your original)
    a = (/ &
      0.125e00_dp,            7.03125e-02_dp, &
      7.32421875e-02_dp,      1.1215209960938e-01_dp, &
      2.2710800170898e-01_dp, 5.7250142097473e-01_dp, &
      1.7277275025845e00_dp,  6.0740420012735e00_dp, &
      2.4380529699556e01_dp,  1.1001714026925e02_dp, &
      5.5133589612202e02_dp,  3.0380905109224e03_dp /)

    a1 = (/ &
      0.125e00_dp,            0.2109375e00_dp, &
      1.0986328125e00_dp,     1.1775970458984e01_dp, &
      2.1461706161499e002_dp, 5.9511522710323e03_dp, &
      2.3347645606175e05_dp,  1.2312234987631e07_dp, &
      8.401390346421e08_dp,   7.2031420482627e10_dp /)

    b = (/ &
     -0.375e00_dp,           -1.171875e-01_dp, &
     -1.025390625e-01_dp,     -1.4419555664063e-01_dp, &
     -2.7757644653320e-01_dp, -6.7659258842468e-01_dp, &
     -1.9935317337513e00_dp, -6.8839142681099e00_dp, &
     -2.7248827311269e01_dp, -1.2159789187654e02_dp, &
     -6.0384407670507e02_dp, -3.3022722944809e03_dp /)

    ! -----------------------
    ! 0) Input sanitisation
    ! -----------------------
    if (.not. ieee_is_finite(real(z)) .or. .not. ieee_is_finite(aimag(z))) then
      bessel0 = cmplx(0.0_dp, 0.0_dp, kind=dp)
      bessel1 = cmplx(0.0_dp, 0.0_dp, kind=dp)
      return
    end if

    ! -----------------------
    ! 1) Optional clipping
    !    (prevents solver-destroying growth from exp(Re(z)))
    ! -----------------------
    zwork = z

    ! Clip real part (prevents exp overflow and insane growth)
    if (real(zwork) > Z_RE_MAX) zwork = cmplx(Z_RE_MAX, aimag(zwork), kind=dp)
    if (real(zwork) < -Z_RE_MAX) zwork = cmplx(-Z_RE_MAX, aimag(zwork), kind=dp)

    ! Clip magnitude (keeps powers zr**k and series stable)
    absz = abs(zwork)
    if (absz > Z_ABS_MAX .and. absz > 0.0_dp) then
      scale = Z_ABS_MAX / absz
      zwork = zwork * scale
      absz  = Z_ABS_MAX
    end if

    a0 = absz

    ! -----------------------
    ! 2) Robust small-|z| handling
    ! -----------------------
    if (a0 <= max(zero_tol, Z_SMALL)) then
      bessel0 = cmplx(1.0_dp, 0.0_dp, kind=dp)  ! I0(0)=1
      bessel1 = 0.5_dp * zwork                  ! I1(z) ~ z/2 near 0
      return
    end if

    ! Preserve original symmetry rule based on ORIGINAL z (not clipped zwork)
    if (real(z) < 0.0_dp) then
      z1 = -zwork
    else
      z1 =  zwork
    end if

    ! -----------------------
    ! 3) Main evaluation
    ! -----------------------
    if (a0 <= 18.0_dp) then
      ! ---- Power series (safe region) ----

      ! I0(z) series
      bessel0 = cmplx(1.0_dp, 0.0_dp, kind=dp)
      cr      = cmplx(1.0_dp, 0.0_dp, kind=dp)
      do k = 1, 50
        cr = 0.25_dp * cr * (z1*z1) / real(k*k, dp)
        bessel0 = bessel0 + cr
        if (abs(cr) < REL_TOL * max(1.0_dp, abs(bessel0))) exit
      end do

      ! I1(z) series (match z1 consistently)
      bessel1 = cmplx(1.0_dp, 0.0_dp, kind=dp)
      cr      = cmplx(1.0_dp, 0.0_dp, kind=dp)
      do k = 1, 50
        cr = 0.25_dp * cr * (z1*z1) / real(k*(k+1), dp)
        bessel1 = bessel1 + cr
        if (abs(cr) < REL_TOL * max(1.0_dp, abs(bessel1))) exit
      end do
      bessel1 = 0.5_dp * z1 * bessel1

    else
      ! ---- Asymptotic branch (large |z|) ----

      if (a0 < 35.0_dp) then
        k0 = 12
      else if (a0 < 50.0_dp) then
        k0 = 9
      else
        k0 = 7
      end if

      ! Guard against tiny z1 (paranoia)
      if (abs(z1) <= max(zero_tol, Z_SMALL)) then
        bessel0 = cmplx(1.0_dp, 0.0_dp, kind=dp)
        bessel1 = 0.5_dp * z1
        return
      end if

      ! If exp(z1) would overflow, DO NOT inject huge().
      ! Return bounded values so your network doesn't explode.
      if (real(z1) > log(huge(1.0_dp)) - 2.0_dp) then
        bessel0 = cmplx(0.0_dp, 0.0_dp, kind=dp)
        bessel1 = cmplx(0.0_dp, 0.0_dp, kind=dp)
        return
      end if

      zr = 1.0_dp / z1
      ca = exp(z1) / sqrt(2.0_dp*pi*z1)

      bessel0 = cmplx(1.0_dp, 0.0_dp, kind=dp)
      do k = 1, k0
        bessel0 = bessel0 + a(k) * (zr ** k)
      end do
      bessel0 = ca * bessel0

      bessel1 = cmplx(1.0_dp, 0.0_dp, kind=dp)
      do k = 1, k0
        bessel1 = bessel1 + b(k) * (zr ** k)
      end do
      bessel1 = ca * bessel1
    end if

    ! Preserve original sign correction rule
    if (real(z) < 0.0_dp) then
      bessel1 = -bessel1
    end if

    ! -----------------------
    ! 4) Final sanitisation (NEVER return NaN/Inf)
    ! -----------------------
    if (.not. ieee_is_finite(real(bessel0)) .or. .not. ieee_is_finite(aimag(bessel0))) then
      bessel0 = cmplx(0.0_dp, 0.0_dp, kind=dp)
    end if
    if (.not. ieee_is_finite(real(bessel1)) .or. .not. ieee_is_finite(aimag(bessel1))) then
      bessel1 = cmplx(0.0_dp, 0.0_dp, kind=dp)
    end if

  end subroutine bessel_complex


!
!###########################################################################
!
!*ax_cr:* Computes A*x for a matrix stored in sparse compressed row form
  subroutine ax_cr ( n, ia, ja, a, x, w )
    implicit none

    integer ( kind = 4 ) n !the order of the system
    integer ( kind = 4 ) ia(*) !ia(n+1) row indices
    integer ( kind = 4 ) ja(*) !ja(nz_num) column indices
    real ( kind = 8 ) a(*) !a(nz_num) Matrix values
    real ( kind = 8 ) x(*) !x(n) Vector to be multiplied by A
    real ( kind = 8 ) w(*) !w(n) Value of A*x

    integer ( kind = 4 ) i
    integer ( kind = 4 ) k1
    integer ( kind = 4 ) k2

    w(1:n) = 0.0_dp

    do i = 1, n
       k1 = ia(i)
       k2 = ia(i+1) - 1
       w(i) = w(i) + dot_product ( a(k1:k2), x(ja(k1:k2)) )
    end do

    return
  end subroutine ax_cr
!
!##############################################################################
!
! *ILU_CR:* computes the incomplete LU factorization of a matrix. For a matrix
! stored in compressed row format.
    !    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
    !    of each row.
    !    Output, real ( kind = 8 ) L(NZ_NUM), the ILU factorization of A.
  subroutine ilu_cr ( n, nz_num, ia, ja, a, ua, l )
    integer ( kind = 4 ) n
    integer ( kind = 4 ) nz_num
    integer ( kind = 4 ) ia(*) !ia(n+1)
    integer ( kind = 4 ) ja(*) !ja(nz_num)
    real ( kind = 8 ) a(*) !a(nz_num)
    integer ( kind = 4 ) ua(*) !ua(n)
    real ( kind = 8 ) l(*) !l(nz_num)

    integer ( kind = 4 ) i
    integer ( kind = 4 ) iw(n)
    integer ( kind = 4 ) j
    integer ( kind = 4 ) jj
    integer ( kind = 4 ) jrow
    integer ( kind = 4 ) jw
    integer ( kind = 4 ) k
    real ( kind = 8 ) tl


    !  Copy A.
    l(1:nz_num) = a(1:nz_num)

    do i = 1, n ! for each row, up to max number of rows
       !  IW points to the nonzero entries in row I.
       iw(1:n) = -1
       do k = ia(i), ia(i+1) - 1 !for each
          iw(ja(k)) = k
       end do
       do j = ia(i), ia(i+1) - 1
          jrow = ja(j)
          if ( i <= jrow ) then
             exit
          end if
          tl = l(j) * l(ua(jrow))
          l(j) = tl
          do jj = ua(jrow) + 1, ia(jrow+1) - 1
             jw = iw(ja(jj))
             if ( jw /= -1 ) then
                l(jw) = l(jw) - tl * l(jj)
             end if
          end do
       end do
       ua(i) = j
       if ( jrow /= i ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'ILU_CR - Fatal error!'
          write ( *, '(a)' ) '  JROW ~= I'
          write ( *, '(a,i8)' ) '  JROW = ', jrow
          write ( *, '(a,i8)' ) '  I    = ', i
          stop
       end if
       if ( abs(l(j)) .le. zero_tol ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'ILU_CR - Fatal error!'
          write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
          write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
          stop
       end if
       l(j) = 1.0_dp / l(j)
    end do

    l(ua(1:n)) = 1.0_dp / l(ua(1:n))

    return
  end subroutine ilu_cr
!
!##############################################################################
!
!*DIAGONAL_POINTER_CR:* finds diagonal entries in a sparse compressed row matrix.
    !    The array UA can be used to locate the diagonal elements of the matrix.
    !    It is assumed that every row of the matrix includes a diagonal element,
    !    and that the elements of each row have been ascending sorted.
subroutine diagonal_pointer_cr ( n, ia, ja, ua )
    integer ( kind = 4 ) n
    integer ( kind = 4 ) ia(*) !ia(n+1)
    integer ( kind = 4 ) ja(*) !ja(nz_num)
    integer ( kind = 4 ) ua(*) !ua(n)

    integer ( kind = 4 ) i
    integer ( kind = 4 ) k

    ua(1:n) = -1

    do i = 1, n
       do k = ia(i), ia(i+1) - 1
          if ( ja(k) == i ) then
             ua(i) = k
          end if
       end do
    end do
    return
  end subroutine diagonal_pointer_cr

  !*****************************************************************************80

  subroutine lus_cr ( n, ia, ja, l, ua, r, z )
!!! LUS_CR applies the incomplete LU preconditioner.
    !    The linear system M * Z = R is solved for Z.  M is the incomplete
    !    LU preconditioner matrix, and R is a vector supplied by the user.
    !    So essentially, we're solving L * U * Z = R.
    !    Input, integer ( kind = 4 ) UA(N), the index of the diagonal element
    !    of each row.
    !    Input, real ( kind = 8 ) R(N), the right hand side.
    !    Output, real ( kind = 8 ) Z(N), the solution of the system M * Z = R.
    implicit none

    integer ( kind = 4 ) n
    integer ( kind = 4 ) ia(*) !ia(n+1)
    integer ( kind = 4 ) ja(*) !ja(nz_num)
    real ( kind = 8 ) l(*) !l(nz_num)
    integer ( kind = 4 ) ua(*) !ua(n)
    real ( kind = 8 ) r(*) !r(n)

    integer ( kind = 4 ) i
    integer ( kind = 4 ) j
    real ( kind = 8 ) w(n)
    real ( kind = 8 ) z(n)

    !  Copy R in.
    w(1:n) = r(1:n)

    !  Solve L * w = w where L is unit lower triangular.
    do i = 2, n
       do j = ia(i), ua(i) - 1
          w(i) = w(i) - l(j) * w(ja(j))
       end do
    end do

    !  Solve U * w = w, where U is upper triangular.
    do i = n, 1, -1
       do j = ua(i) + 1, ia(i+1) - 1
          w(i) = w(i) - l(j) * w(ja(j))
       end do
       w(i) = w(i) / l(ua(i))
    end do

    !  Copy Z out.
    z(1:n) = w(1:n)

    return
  end subroutine lus_cr

  !*****************************************************************************80
  subroutine mult_givens ( c, s, k, g )
!!! MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
    !    In order to make it easier to compare this code with the Original C,
    !    the vector indexing is 0-based.
    !    Input, real ( kind = 8 ) C, S, the cosine and sine of a Givens
    !    rotation.
    !
    !    Input, integer ( kind = 4 ) K, indicates the location of the first
    !    vector entry.
    !
    !    Input/output, real ( kind = 8 ) G(1:K+1), the vector to be modified.
    !    On output, the Givens rotation has been applied to entries G(K) and G(K+1).

    implicit none

    real ( kind = 8 ) c
    real ( kind = 8 ) s
    integer ( kind = 4 ) k
    real ( kind = 8 ) g(*) !g(1:k+1)

    real ( kind = 8 ) g1
    real ( kind = 8 ) g2

    g1 = c * g(k) - s * g(k+1)
    g2 = s * g(k) + c * g(k+1)

    g(k)   = g1
    g(k+1) = g2

    return
  end subroutine mult_givens


  !*****************************************************************************80
  subroutine rearrange_cr ( n, ia, ja, a )
!!! REARRANGE_CR sorts a sparse compressed row matrix.
    !    This routine guarantees that the entries in the CR matrix
    !    are properly sorted.
    !
    !    After the sorting, the entries of the matrix are rearranged in such
    !    a way that the entries of each column are listed in ascending order
    !    of their column values.
    !    Input, integer ( kind = 4 ) N, the order of the system.
    !
    !    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
    !
    !    Input, integer ( kind = 4 ) IA(N+1), the compressed row indices.
    !
    !    Input/output, integer ( kind = 4 ) JA(NZ_NUM), the column indices.
    !    On output, these may have been rearranged by the sorting.
    !
    !    Input/output, real ( kind = 8 ) A(NZ_NUM), the matrix values.  On output,
    !    the matrix values may have been moved somewhat because of the sorting.
    !
    implicit none

    integer ( kind = 4 ) n
    integer ( kind = 4 ) ia(*) !ia(n+1)
    integer ( kind = 4 ) ja(*) !ja(nz_num)
    real ( kind = 8 ) a(*) !a(nz_num)

    integer ( kind = 4 ) i
    integer ( kind = 4 ) i4temp
    integer ( kind = 4 ) k
    integer ( kind = 4 ) l
    real ( kind = 8 ) r8temp

    do i = 1, n

       do k = ia(i), ia(i+1) - 2
          do l = k + 1, ia(i+1) - 1

             if ( ja(l) < ja(k) ) then
                i4temp = ja(l)
                ja(l)  = ja(k)
                ja(k)  = i4temp

                r8temp = a(l)
                a(l)   = a(k)
                a(k)   = r8temp
             end if

          end do
       end do

    end do

    return
  end subroutine rearrange_cr


!######################################################################

  !
  !*sort_integer_list:* sorts a list of integer values into a non-decreasing order.
  ! sorts N integer IDATA values into a non-decreasing sequence using IHEAPSORT
  ! (N>50) or ISHELLSORT (N>50) and then  removes all duplicates from the list. On
  ! exit N contains the number of unique elements in the list.
  !
  subroutine sort_integer_list(N,IDATA)
    integer :: IDATA(:),N

    !Local Variables
    integer :: count1,index,itemp,nolist
    logical :: continue

    character(len=60) :: sub_name

    sub_name = 'sort_integer_list'
    call enter_exit(sub_name,1)

    !order the array non-decreasing
    do nolist=2,N
       count1=0
       continue=.true.
       do while(continue)
          if(IDATA(nolist-count1).lt.IDATA(nolist-count1-1))then
             itemp=IDATA(nolist-1)
             IDATA(nolist-1)=IDATA(nolist)
             IDATA(nolist)=itemp
             count1=count1+1
             if(nolist-count1-1.eq.0) continue=.false.
          else
             continue=.false.
          endif
       enddo !while
    enddo !N

    !eliminate duplicate entries
    index=0
    do nolist=2,N
       if(IDATA(nolist).eq.IDATA(nolist-1)) then
          index=index+1
       else
          IDATA(nolist-index)=IDATA(nolist)
       endif
    enddo !nolist

    N=N-index

    call enter_exit(sub_name,2)

  end subroutine sort_integer_list

!!!#########################################################################
  !*sort_real_list:* sorts a list of real values into a non-decreasing order
  ! using a bubble sort algorithm.

  subroutine sort_real_list(n,RDATA,INDEX)

    integer :: INDEX(*),n
    real(dp) :: RDATA(*)

    !Local Variables
    integer :: FLAG,i,ITEMP,j,k
    real(dp) :: TEMP

    character(len=60) :: sub_name

    sub_name = 'sort_real_list'
    call enter_exit(sub_name,1)

    if(N.LE.1) then
    else
       FLAG=n
       do i=1,n
          k=FLAG-1
          FLAG=0
          do j=1,k
             if(RDATA(j).gt.RDATA(j+1)) then
                TEMP=RDATA(j)
                RDATA(j)=RDATA(j+1)
                RDATA(j+1)=TEMP
                ITEMP=INDEX(j)
                INDEX(j)=INDEX(j+1)
                INDEX(j+1)=ITEMP
                FLAG=j
             endif
          enddo
          if(FLAG.eq.0) then
             write(*,*) 'warning in rsort'
          endif
       enddo
    endif

    call enter_exit(sub_name,2)

  end subroutine sort_real_list


end module math_utilities
