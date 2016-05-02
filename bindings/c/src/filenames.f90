module filenames_c

contains
  !###################################################################################
  subroutine read_geometry_evaluate_flow_c() bind(C, name="read_geometry_evaluate_flow_c")
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_READ_GEOMETRY_EVALUATE_FLOW_C" :: READ_GEOMETRY_EVALUATE_FLOW_C
    use filenames, only : read_geometry_evaluate_flow
    implicit none

    call read_geometry_evaluate_flow()

  end subroutine read_geometry_evaluate_flow_c

  !###################################################################################
  subroutine read_geometry_main_c() bind(C, name="read_geometry_main_c")
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_READ_GEOMETRY_MAIN_C" :: READ_GEOMETRY_MAIN_C
  use filenames, only : read_geometry_main
  implicit none

  call read_geometry_main()

  end subroutine read_geometry_main_c

  !###################################################################################
  subroutine get_filename_c(label, label_len, filename) bind(C, name="get_filename_c")
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"DLL_GET_FILENAME_C" :: READ_GET_FILENAME_C
    use iso_c_binding, only: c_ptr, c_f_pointer, c_char, c_loc
    use utils_c, only: strncpy, f_c_string
    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
    use filenames, only : get_filename
    implicit none

    integer,intent(in) :: label_len
    type(c_ptr), value, intent(in) :: label
    type(c_ptr), value, intent(in) :: filename
    character(len=MAX_FILENAME_LEN) :: filename_f
    character(len=MAX_STRING_LEN) :: label_f
    character(len=1,kind=c_char), pointer :: filename_c_chars(:)

    call strncpy(label_f, label, label_len)

    filename_f = get_filename(label_f)
    call c_f_pointer(filename,filename_c_chars,[MAX_FILENAME_LEN + 1])
    filename_c_chars = f_c_string(filename_f)

  end subroutine get_filename_c
!
!#####################################################################################################
!
end module filenames_c
