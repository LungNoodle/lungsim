module growtree_c
  implicit none
  private

contains

  !
  !###################################################################################
  !
  ! the main growing subroutine. Generates a volume-filling tree into a closed surface.
  subroutine grow_tree_c(parent_elems_len, parent_elems, &
       angle_max, angle_min, branch_fraction, length_limit, shortest_length, rotation_limit, &
       grouping, grouping_len) bind(C, name="grow_tree_c")
    
    use arrays,only: dp
    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use growtree,only: grow_tree
    implicit none
    
    integer,intent(in) :: parent_elems_len
    integer,intent(in) :: parent_elems(parent_elems_len)
    real(dp),intent(in) :: angle_max
    real(dp),intent(in) :: angle_min
    real(dp),intent(in) :: branch_fraction
    real(dp),intent(in) :: length_limit
    real(dp),intent(in) :: shortest_length
    real(dp),intent(in) :: rotation_limit
    integer,intent(in) :: grouping_len
    type(c_ptr), value, intent(in) :: grouping
    character(len=MAX_FILENAME_LEN) :: grouping_f
    
    call strncpy(grouping_f, grouping, grouping_len)

    call grow_tree(parent_elems, angle_max, angle_min, branch_fraction, length_limit,&
         shortest_length, rotation_limit, grouping_f)

  end subroutine grow_tree_c

  !
  !###################################################################################
  !
  ! export the seed-to-terminal mapping produced by the most recent tree growth.
  subroutine export_terminal_seed_mapping_c(filename, filename_len) &
       bind(C, name="export_terminal_seed_mapping_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use growtree,only: export_terminal_seed_mapping
    implicit none

    integer,intent(in) :: filename_len
    type(c_ptr), value, intent(in) :: filename
    character(len=MAX_FILENAME_LEN) :: filename_f

    call strncpy(filename_f, filename, filename_len)
    call export_terminal_seed_mapping(filename_f)

  end subroutine export_terminal_seed_mapping_c

  !
  !###################################################################################
  !
  ! option to smooth branching in a generated tree
  subroutine smooth_1d_tree_c(num_elem_start, n_smoothing_steps) bind(C, name="smooth_1d_tree_c")
    
    use arrays,only: dp
    use growtree,only: smooth_1d_tree
    implicit none
    
    integer,intent(in) :: num_elem_start, n_smoothing_steps

    call smooth_1d_tree(num_elem_start, n_smoothing_steps)

  end subroutine smooth_1d_tree_c
  ! 
  !#########################################################################
  !
  subroutine align_terminals_with_seeds_c(seed_dist) bind(C, name="align_terminals_with_seeds_c")

    use arrays,only: dp
    use growtree,only: align_terminals_with_seeds
    implicit none

    real(dp),intent(in) :: seed_dist

    call align_terminals_with_seeds(seed_dist)

  end subroutine align_terminals_with_seeds_c
  ! 
  !#########################################################################
  ! 

end module growtree_c
