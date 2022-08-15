module math_utilities
!*Brief Description:* This module contains solvers required for lung problems
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


contains

subroutine bessel_complex(z,bessel0,bessel1)

  complex(dp), intent(in) :: z
  complex(dp), intent(out) :: bessel0,bessel1

  real(dp) :: a(12),a1(10),b(12)
  real(dp) :: a0
  complex(dp) :: ci,cr,z1,ca,zr
  integer :: k,k0
!  complex ( kind = 8 ) ca
!  complex ( kind = 8 ) cb
!  complex ( kind = 8 ) ci
!  complex ( kind = 8 ) cr
!  complex ( kind = 8 ) cs
!  complex ( kind = 8 ) ct
!  complex ( kind = 8 ) cw
!  integer ( kind = 4 ) k
!  integer ( kind = 4 ) k0
!  real ( kind = 8 ) pi
!  real ( kind = 8 ) w0
!  complex ( kind = 8 ) z
!  complex ( kind = 8 ) z1
!  complex ( kind = 8 ) z2
!  complex ( kind = 8 ) zr
!  complex ( kind = 8 ) zr2
  a = (/ &
    0.125e00_dp,           7.03125e-02_dp,&
    7.32421875e-02_dp,      1.1215209960938e-01_dp,&
    2.2710800170898e-01_dp, 5.7250142097473e-01_dp,&
    1.7277275025845e00_dp, 6.0740420012735e00_dp,&
    2.4380529699556e01_dp, 1.1001714026925e02_dp,&
    5.5133589612202e02_dp, 3.0380905109224e03_dp /)
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

!
  ci = cmplx (0.0_dp,1.0_dp,8)
  a0 = abs (z)
  z1 = z
!
  if(abs(a0).le.zero_tol)then
    bessel0 = cmplx(1.0_dp,0.0_dp,8)
    bessel1 = cmplx(0.0_dp,0.0_dp,8)
  endif

  if(real(z).lt.0.0_dp) then
    z1 = -z
  endif
!
  if( a0 <= 18.0_dp ) then

    bessel0 =cmplx(1.0_dp,0.0_dp,8)
    cr = cmplx(1.0_dp,0.0_dp,8)
    do k = 1,50
      cr = 0.25_dp*cr* z1**2/k**2
      bessel0 = bessel0+cr
      if (abs (cr/bessel0).LT.1.0e-15_dp) then
        exit
      endif
    enddo

    bessel1 =cmplx(1.0_dp,0.0_dp,8)
    cr = cmplx(1.0_dp,0.0_dp,8)
    do k = 1,50
      cr = 0.25_dp*cr*z**2/(k*(k+1))
      bessel1 = bessel1+cr
      if (abs (cr/bessel1).LT.1.0e-15_dp) then
        exit
      endif
    enddo

    bessel1 = 0.5_dp*z1*bessel1

  else

    if ( a0 < 35.0_dp ) then
      k0 = 12
    else if ( a0 < 50.0_dp ) then
      k0 = 9
    else
      k0 = 7
    end if

    ca = exp(z1)/sqrt(2.0_dp*pi*z1)
    bessel0 = cmplx(1.0_dp,0.0_dp,8)
    zr = 1.0_dp/z1
    do k = 1,k0
      bessel0 = bessel0 + a(k) * zr ** k
    enddo
    bessel0 = ca * bessel0
    bessel1 = cmplx (1.0_dp,0.0_dp,8)
    do k = 1,k0
      bessel1 =bessel1+b(k)*zr**k
    end do
    bessel1 = ca * bessel1
  endif


  if ( real (z).lt.0.0_dp)then
    bessel1 = - bessel1
  endif

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
