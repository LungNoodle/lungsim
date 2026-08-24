module test_lymphatic_parameters
  use testdrive, only : new_unittest, unittest_type, error_type, check
  implicit none
  private

  public :: collect_lymphatic_parameters

contains

subroutine collect_lymphatic_parameters(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("test_update_lymphatics", test_update_lymphatics) &
    ]

end subroutine collect_lymphatic_parameters

subroutine test_update_lymphatics(error)
  use precision, only: dp
  use parameter_types, only: lymphatic_params, update_lymphatics
  implicit none

  type(error_type), allocatable, intent(out) :: error

  call check(error, lymphatic_params%lung_mass_g, 639.0_dp)
  if (allocated(error)) return

  call update_lymphatics('lung_mass_g', 1.5_dp)
  call check(error, lymphatic_params%lung_mass_g, 1.5_dp)
  if (allocated(error)) return

  call update_lymphatics('lymphatic_surface_area_ratio', 2.0_dp)
  call check(error, lymphatic_params%lymphatic_density, 2.0_dp)
  if (allocated(error)) return

  call update_lymphatics('integration_steps_per_transit', 48.0_dp)
  call check(error, lymphatic_params%integration_steps_per_transit, 48)

  ! Restore defaults so this test does not alter later model tests.
  call update_lymphatics('lung_mass_g', 639.0_dp)
  call update_lymphatics('lymphatic_density', 1.0_dp)
  call update_lymphatics('integration_steps_per_transit', 96.0_dp)

end subroutine test_update_lymphatics

end module test_lymphatic_parameters
