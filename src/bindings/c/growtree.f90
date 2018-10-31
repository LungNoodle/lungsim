module growtree_c
  implicit none
  private

contains
! 
!#########################################################################
! 
!*growtree:* the main growing subroutine. Generates a volume-filling tree into a closed surface.
  subroutine grow_tree_c(parent_ne, surface_elems, angle_max, angle_min, branch_fraction, length_limit,&
shortest_length, rotation_limit) bind(C, name="grow_tree_c")

    use arrays,only: dp
    use iso_c_binding, only: c_ptr
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

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_grow_tree(parent_ne, surface_elems, angle_max, angle_min, branch_fraction, length_limit,&
shortest_length, rotation_limit)
#else
    call grow_tree(parent_ne, surface_elems, angle_max, angle_min, branch_fraction, length_limit,&
shortest_length, rotation_limit)
#endif

  end subroutine grow_tree_c

! 
!#########################################################################
! 

end module growtree_c
