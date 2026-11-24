module test_micro_remodel
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use precision, only: dp
  use capillaryflow, only: micro_remodel_radius,select_remodeling_side
  implicit none
  private
  public :: collect_micro_remodel

contains

subroutine collect_micro_remodel(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)
  testsuite = [ &
    new_unittest('grade_zero', test_grade_zero), &
    new_unittest('side_selection', test_side_selection), &
    new_unittest('diameter_bands', test_diameter_bands), &
    new_unittest('side_optional', test_side_optional) &
  ]
end subroutine collect_micro_remodel

subroutine test_grade_zero(error)
  type(error_type), allocatable, intent(out) :: error
  real(dp) :: radius,prune_scale
  radius=80.0e-6_dp
  call micro_remodel_radius(75.0e-6_dp,radius,0.0_dp,'both','arterial',prune_scale)
  call check(error,abs(radius-80.0e-6_dp)<1.0e-12_dp)
  if (allocated(error)) return
  call check(error,abs(prune_scale-1.0_dp)<1.0e-12_dp)
end subroutine test_grade_zero

subroutine test_side_selection(error)
  type(error_type), allocatable, intent(out) :: error
  real(dp) :: radius,prune_scale,arterial_radius
  radius=80.0e-6_dp
  call micro_remodel_radius(75.0e-6_dp,radius,60.0_dp,'venous','arterial',prune_scale)
  call check(error,abs(radius-80.0e-6_dp)<1.0e-12_dp)
  if (allocated(error)) return
  radius=80.0e-6_dp
  call micro_remodel_radius(75.0e-6_dp,radius,60.0_dp,'arterial','arterial',prune_scale)
  arterial_radius=radius
  call check(error,arterial_radius<75.0e-6_dp)
  if (allocated(error)) return
  radius=80.0e-6_dp
  call micro_remodel_radius(75.0e-6_dp,radius,60.0_dp,'both','arterial',prune_scale)
  call check(error,abs(radius-arterial_radius)<1.0e-12_dp)
  if (allocated(error)) return
  call check(error,prune_scale>1.0_dp)
end subroutine test_side_selection

subroutine test_diameter_bands(error)
  type(error_type), allocatable, intent(out) :: error
  real(dp) :: radius,prune_scale
  ! 20 um diameter: eligible for effective pruning, outside narrowing.
  radius=10.0e-6_dp
  call micro_remodel_radius(10.0e-6_dp,radius,60.0_dp,'venous','venous',prune_scale)
  call check(error,abs(radius-10.0e-6_dp)<1.0e-12_dp)
  if (allocated(error)) return
  call check(error,prune_scale>1.0_dp)
  if (allocated(error)) return
  ! 600 um diameter: hypertrophy only; zero pressure change means no narrowing.
  radius=300.0e-6_dp
  call micro_remodel_radius(300.0e-6_dp,radius,60.0_dp,'arterial','arterial',prune_scale)
  call check(error,abs(radius-300.0e-6_dp)<1.0e-12_dp)
  if (allocated(error)) return
  call check(error,abs(prune_scale-1.0_dp)<1.0e-12_dp)
end subroutine test_diameter_bands

subroutine test_side_optional(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=8) :: side
  ! Healthy runs need no side, and any side given is ignored.
  call select_remodeling_side(0.0_dp,'simple_tree',selected_side=side)
  call check(error,side=='none')
  if (allocated(error)) return
  call select_remodeling_side(0.0_dp,'full_plus_ladder','none',side)
  call check(error,side=='none')
  if (allocated(error)) return
  call select_remodeling_side(0.0_dp,'full_plus_ladder','arterial',side)
  call check(error,side=='none')
  if (allocated(error)) return
  call select_remodeling_side(40.0_dp,'full_plus_ladder','both',side)
  call check(error,side=='both')
end subroutine test_side_optional

end module test_micro_remodel
