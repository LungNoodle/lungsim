module test_micro_remodel
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use precision, only: dp
  use capillaryflow, only: micro_remodel_radius,micro_prune_alive,select_remodeling_side
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
    new_unittest('side_optional', test_side_optional), &
    new_unittest('prune_cascade', test_prune_cascade) &
  ]
end subroutine collect_micro_remodel

subroutine test_grade_zero(error)
  type(error_type), allocatable, intent(out) :: error
  real(dp) :: radius
  radius=80.0e-6_dp
  call micro_remodel_radius(75.0e-6_dp,radius,0.0_dp,'both','arterial')
  call check(error,abs(radius-80.0e-6_dp)<1.0e-12_dp)
end subroutine test_grade_zero

subroutine test_side_selection(error)
  type(error_type), allocatable, intent(out) :: error
  real(dp) :: radius,arterial_radius
  radius=80.0e-6_dp
  call micro_remodel_radius(75.0e-6_dp,radius,60.0_dp,'venous','arterial')
  call check(error,abs(radius-80.0e-6_dp)<1.0e-12_dp)
  if (allocated(error)) return
  radius=80.0e-6_dp
  call micro_remodel_radius(75.0e-6_dp,radius,60.0_dp,'arterial','arterial')
  arterial_radius=radius
  call check(error,arterial_radius<75.0e-6_dp)
  if (allocated(error)) return
  radius=80.0e-6_dp
  call micro_remodel_radius(75.0e-6_dp,radius,60.0_dp,'both','arterial')
  call check(error,abs(radius-arterial_radius)<1.0e-12_dp)
end subroutine test_side_selection

subroutine test_diameter_bands(error)
  type(error_type), allocatable, intent(out) :: error
  real(dp) :: radius
  ! 20 um diameter: below the narrowing band, unchanged.
  radius=10.0e-6_dp
  call micro_remodel_radius(10.0e-6_dp,radius,60.0_dp,'venous','venous')
  call check(error,abs(radius-10.0e-6_dp)<1.0e-12_dp)
  if (allocated(error)) return
  ! 600 um diameter: hypertrophy only; zero distension means no change.
  radius=300.0e-6_dp
  call micro_remodel_radius(300.0e-6_dp,radius,60.0_dp,'arterial','arterial')
  call check(error,abs(radius-300.0e-6_dp)<1.0e-12_dp)
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

subroutine test_prune_cascade(error)
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: ngen=9
  real(dp) :: rad_a(ngen),rad_v(ngen),alive(ngen)
  integer :: i,k
  ! Example ladder: terminal artery 400 um diameter, vein 380 um, tapering to 20/18 um.
  do i=1,ngen
    rad_a(i)=200.0e-6_dp+(10.0e-6_dp-200.0e-6_dp)*i/ngen
    rad_v(i)=190.0e-6_dp+(9.0e-6_dp-190.0e-6_dp)*i/ngen
  enddo
  ! Healthy or grade <= 20: nothing pruned.
  call micro_prune_alive(rad_a,rad_v,ngen,20.0_dp,'arterial',alive,k)
  call check(error,all(abs(alive-1.0_dp)<1.0e-12_dp) .and. k==0)
  if (allocated(error)) return
  ! Grade 30 (limit 300 um, p_f=1/16): arterioles are 358, 316, 273 um, so gen 3 is first.
  call micro_prune_alive(rad_a,rad_v,ngen,30.0_dp,'arterial',alive,k)
  call check(error,k==3)
  if (allocated(error)) return
  call check(error,all(abs(alive(1:2)-1.0_dp)<1.0e-12_dp))
  if (allocated(error)) return
  call check(error,all(abs(alive(3:ngen)-(1.0_dp-1.0_dp/16.0_dp))<1.0e-12_dp))
  if (allocated(error)) return
  ! Grade 60 (limit 500 um, p_f=1/4): gen 1 qualifies, so the whole ladder loses 1/4.
  call micro_prune_alive(rad_a,rad_v,ngen,60.0_dp,'both',alive,k)
  call check(error,k==1 .and. all(abs(alive-0.75_dp)<1.0e-12_dp))
  if (allocated(error)) return
  ! Venous side uses the venule diameters: 340, 299.6 um, so gen 2 is first.
  call micro_prune_alive(rad_a,rad_v,ngen,30.0_dp,'venous',alive,k)
  call check(error,k==2)
  if (allocated(error)) return
  call micro_prune_alive(rad_a,rad_v,ngen,60.0_dp,'none',alive,k)
  call check(error,all(abs(alive-1.0_dp)<1.0e-12_dp) .and. k==0)
end subroutine test_prune_cascade

end module test_micro_remodel
