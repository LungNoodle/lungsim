module growtree_c
  implicit none
  private

contains
! 
!#########################################################################
! 
!*growtree:* the main growing subroutine. Generates a volume-filling tree into a closed surface.
  subroutine grow_tree_c(parent_ne, surface_elems, angle_max, angle_min, branch_fraction, length_limit,&
shortest_length, rotation_limit, to_export, filename, filename_len) bind(C, name="grow_tree_c")

    use arrays,only: dp
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use growtree,only: grow_tree
    implicit none

    integer,intent(in) :: parent_ne
    integer,intent(in) :: surface_elems(:)
    real(dp),intent(in) :: angle_max
    real(dp),intent(in) :: angle_min
    real(dp),intent(in) :: branch_fraction
    real(dp),intent(in) :: length_limit
    real(dp),intent(in) :: shortest_length
    real(dp),intent(in) :: rotation_limit
    logical,intent(in) :: to_export
    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: filename
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, filename, filename_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_grow_tree(parent_ne, surface_elems, angle_max, angle_min, branch_fraction, length_limit,&
shortest_length, rotation_limit, to_export, filename_f)
#else
    call grow_tree(parent_ne, surface_elems, angle_max, angle_min, branch_fraction, length_limit,&
shortest_length, rotation_limit, to_export, filename_f)
#endif

  end subroutine grow_tree_c

! 
!#########################################################################
! 

end module growtree_c
