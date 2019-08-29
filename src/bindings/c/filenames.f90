module filenames_c

contains
  !###################################################################################
  subroutine read_geometry_evaluate_flow_c() bind(C, name="read_geometry_evaluate_flow_c")
    use filenames, only : read_geometry_evaluate_flow
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_read_geometry_evaluate_flow()
#else
    call read_geometry_evaluate_flow()
#endif

  end subroutine read_geometry_evaluate_flow_c

  !###################################################################################
  subroutine read_geometry_main_c() bind(C, name="read_geometry_main_c")
    use filenames, only : read_geometry_main
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_read_geometry_main()
#else
    call read_geometry_main()
#endif

  end subroutine read_geometry_main_c

  !###################################################################################
  subroutine get_filename_c(label, label_len, filename) bind(C, name="get_filename_c")

    use iso_c_binding, only: c_ptr, c_f_pointer, c_char, c_int, c_null_char
    use utils_c, only: strncpy, f_c_string
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    use filenames, only : get_filename
    implicit none

    integer(c_int),intent(in) :: label_len
    type(c_ptr), value, intent(in) :: label
    !type(c_ptr), value, intent(in) :: filename
    character(c_char), dimension(MAX_FILENAME_LEN + 1), intent(inout) :: filename
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: label_f
    character(len=1,kind=c_char), pointer :: filename_c_chars(:)
    integer :: I, N, test

    call strncpy(label_f, label, label_len)

    filename_f = get_filename(label_f)
    test = 0
    N = MAX_FILENAME_LEN
    DO I = N, 1, -1
      filename(I) = filename_f(I:I)
      if ((test == 0) .and. (filename_f(I:I) /= ' ')) then
        filename(I+1) = C_NULL_CHAR
        test = 1
      end if
    END DO
    filename(MAX_FILENAME_LEN + 1) = C_NULL_CHAR

  end subroutine get_filename_c
!
!#####################################################################################################
!
end module filenames_c
