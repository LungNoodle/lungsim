module test_geometry
  use testdrive, only : new_unittest, unittest_type, error_type, check
  implicit none
  private

  public :: collect_geometry

contains

subroutine collect_geometry(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("matching mesh with noncontiguous element labels", &
      test_matching_mesh_noncontiguous_labels) &
    ]

end subroutine collect_geometry

subroutine test_matching_mesh_noncontiguous_labels(error)
  use arrays, only: elem_cnct, elem_nodes, elems, nodes, num_elems, num_nodes
  use geometry, only: add_matching_mesh, append_units, define_1d_elements, &
    define_node_geometry
  use indices, only: define_problem_type
  use other_consts, only: MAX_FILENAME_LEN
  implicit none

  type(error_type), allocatable, intent(out) :: error
  character(len=MAX_FILENAME_LEN) :: node_file, elem_file, problem_type

  node_file = 'matching_mesh_gap.exnode'
  elem_file = 'matching_mesh_gap.exelem'
  problem_type = 'perfusion'
  call write_test_mesh(node_file, elem_file)
  call define_problem_type(problem_type)
  call define_node_geometry(node_file)
  call define_1d_elements(elem_file)
  call append_units()
  call add_matching_mesh()

  call check(error, num_nodes, 6, 'matching mesh node count')
  if (allocated(error)) return
  call check(error, num_elems, 5, 'matching mesh element count')
  if (allocated(error)) return

  ! The arterial labels contain a gap: 1, 3. The matching veins and terminal
  ! connector must receive fresh global labels without using those labels as
  ! local array indices.
  call check(error, elems, [1, 3, 4, 5, 6], 'global element labels')
  if (allocated(error)) return
  call check(error, nodes, [1, 2, 3, 4, 5, 6], 'global node labels')
  if (allocated(error)) return

  call check(error, elem_nodes(:,3), [5, 4], 'first matching vein nodes')
  if (allocated(error)) return
  call check(error, elem_nodes(:,4), [6, 5], 'second matching vein nodes')
  if (allocated(error)) return
  call check(error, elem_nodes(:,5), [3, 6], 'terminal connector nodes')
  if (allocated(error)) return

  call check(error, elem_cnct(1,1,2), 5, 'artery-to-connector connectivity')
  if (allocated(error)) return
  call check(error, elem_cnct(-1,1,4), 5, 'connector-to-vein connectivity')

end subroutine test_matching_mesh_noncontiguous_labels

subroutine write_test_mesh(node_file, elem_file)
  character(len=*), intent(in) :: node_file, elem_file
  integer :: unit

  open(newunit=unit, file=node_file, status='replace', action='write')
  write(unit,'(a)') 'Group name: gap_test'
  write(unit,'(a)') '#Fields=1'
  write(unit,'(a)') '1) coordinates, coordinate, rectangular cartesian, #Components=3'
  write(unit,'(a)') ' x. Value index=1, #Derivatives=0'
  write(unit,'(a)') ' y. Value index=2, #Derivatives=0'
  write(unit,'(a)') ' z. Value index=3, #Derivatives=0'
  write(unit,'(a)') 'Node: 1'
  write(unit,'(a)') ' 0.0 0.0 0.0'
  write(unit,'(a)') 'Node: 2'
  write(unit,'(a)') ' 1.0 0.0 0.0'
  write(unit,'(a)') 'Node: 3'
  write(unit,'(a)') ' 2.0 0.0 0.0'
  close(unit)

  open(newunit=unit, file=elem_file, status='replace', action='write')
  write(unit,'(a)') 'Group name: gap_test'
  write(unit,'(a)') 'Shape. Dimension=1'
  write(unit,'(a)') 'Element: 1 0 0'
  write(unit,'(a)') ' Nodes:'
  write(unit,'(a)') ' 1 2'
  write(unit,'(a)') 'Element: 3 0 0'
  write(unit,'(a)') ' Nodes:'
  write(unit,'(a)') ' 2 3'
  close(unit)

end subroutine write_test_mesh

end module test_geometry
