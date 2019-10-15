module utils_c

  implicit none

  public c_to_f_string, strcpy, strncpy

contains

!!!######################################################################
  FUNCTION f_c_string (F_STRING) RESULT (C_STRING)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR
    use other_consts, only: MAX_FILENAME_LEN
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: F_STRING
    CHARACTER(LEN=1,KIND=C_CHAR) :: C_STRING(MAX_FILENAME_LEN+1)
    character(len=1, kind=c_char), allocatable :: p_c_str(:)

    INTEGER                      :: N, I, test

    test = 0
    N = MAX_FILENAME_LEN
    DO I = N, 1, -1
      C_STRING(I) = F_STRING(I:I)
      if ((test == 0) .and. (f_string(I:I) /= ' ')) then
        c_string(I+1) = C_NULL_CHAR
        test = 1
      end if
    END DO
    C_STRING(N + 1) = C_NULL_CHAR

  END FUNCTION f_c_string

  function c_to_f_string(s) result(str)
    use iso_c_binding, only: c_char, c_null_char
    implicit none
    character(kind=c_char,len=1), intent(in) :: s(*)
    character(len=:), allocatable :: str
    integer i, nchars
    i = 1
    do
       if (s(i) == c_null_char) exit
       i = i + 1
    end do
    nchars = i - 1  ! Exclude null character from Fortran string
    allocate(character(LEN=255) :: str)
    str = transfer(s(1:nchars), str)
  end function c_to_f_string

  subroutine strcpy(dest, src)
    use iso_c_binding, only: c_char, c_null_char
    implicit none

    character(kind=c_char,len=1), intent(in) :: src(:)
    character(len=*), intent(out) :: dest
    !Local variables
    integer :: i,length

    IF(LEN(dest)>=SIZE(src,1)) THEN
      length=SIZE(src,1)
    ELSE
      length=LEN(dest)
    ENDIF
    dest=""
    DO i=1,length
      IF(src(i) == c_null_char) THEN
        EXIT
      ELSE
        dest(i:i) = src(i)
      ENDIF
    ENDDO !i

  end subroutine strcpy

  subroutine strncpy(dest, src, src_length)
    use iso_c_binding, only: c_char, c_null_char, c_f_pointer, c_ptr
    implicit none

    integer, intent(in) :: src_length
    type(c_ptr), value, intent(in) :: src
    character(len=*), intent(out) :: dest
    CHARACTER(LEN=1,KIND=C_CHAR), POINTER :: src_c_chars(:)

    !Local variables
    integer :: i,length

    call c_f_pointer(src, src_c_chars, [src_length])
    IF(LEN(dest)>=src_length) THEN
      length=src_length
    ELSE
      length=LEN(dest)
    ENDIF
    dest=""
    DO i=1,length
      IF(src_c_chars(i) == c_null_char) THEN
        EXIT
      ELSE
        dest(i:i) = src_c_chars(i)
      ENDIF
    ENDDO !i

  end subroutine strncpy

end module utils_c
