module growtree_c
  implicit none
  private

contains

  !
  !###################################################################################
  !
  ! the main growing subroutine. Generates a volume-filling tree into a closed surface.
  subroutine grow_tree_c(surface_elems_len, surface_elems, parent_ne, supernumerary_ne, &
       angle_max, angle_min, branch_fraction, length_limit, shortest_length, rotation_limit, &
       to_export, filename, filename_len, grouping, grouping_len) bind(C, name="grow_tree_c")
    
    use arrays,only: dp
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use growtree,only: grow_tree
    implicit none
    
    integer,intent(in) :: surface_elems_len
    integer,intent(in) :: surface_elems(surface_elems_len)
    integer,intent(in) :: parent_ne
    integer,intent(in) :: supernumerary_ne
    real(dp),intent(in) :: angle_max
    real(dp),intent(in) :: angle_min
    real(dp),intent(in) :: branch_fraction
    real(dp),intent(in) :: length_limit
    real(dp),intent(in) :: shortest_length
    real(dp),intent(in) :: rotation_limit
    logical,intent(in) :: to_export
    integer,intent(in) :: filename_len
    integer,intent(in) :: grouping_len
    type(c_ptr), value, intent(in) :: filename, grouping
    character(len=MAX_FILENAME_LEN) :: filename_f, grouping_f
    
    call strncpy(filename_f, filename, filename_len)
    call strncpy(grouping_f, grouping, grouping_len)

    call grow_tree(surface_elems, parent_ne, supernumerary_ne, angle_max, angle_min, branch_fraction, length_limit,&
         shortest_length, rotation_limit, to_export, filename_f, grouping_f)

  end subroutine grow_tree_c

  !
  !###################################################################################
  !
  ! option to smooth branching in a generated tree
  subroutine smooth_1d_tree_c(num_elem_start, length_limit) bind(C, name="smooth_1d_tree_c")
    
    use arrays,only: dp
    use growtree,only: smooth_1d_tree
    implicit none
    
    integer,intent(in) :: num_elem_start
    real(dp),intent(in) :: length_limit

    call smooth_1d_tree(num_elem_start, length_limit)

  end subroutine smooth_1d_tree_c
! 
!#########################################################################
! 

end module growtree_c
