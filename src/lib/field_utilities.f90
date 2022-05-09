module field_utilities
  !*Brief Description:* This module contains all the subroutines that perform general operations
  !on fields
  !*LICENSE:*
  !
  !
  !
  !*Full Description:*
  !More info on what the module does if necessary
  !
  use arrays
  use diagnostics
  use indices
  use other_consts
  use precision
  
  implicit none
  
  !Module parameters
  
  !Module types
  
  !Module variables
  
  !Interfaces
  private
  public scale_flow_to_inlet
  
contains

!!!#############################################################################
  
  subroutine scale_flow_to_inlet(inlet_flow,VorQ)
    !*scale_flow_field:* Scales a flow field to an 'inlet flow' value (real units).

    real(dp),intent(in) :: inlet_flow
    character(len=1), intent(in) :: VorQ
    ! Local variables
    real(dp) :: ratio
    character(len=60) :: sub_name

    ! --------------------------------------------------------------------------

    sub_name = 'scale_flow_to_inlet'
    call enter_exit(sub_name,1)

    if(VorQ.eq.'V')then
      if(abs(elem_field(ne_Vdot,1)).gt.zero_tol)then
         ratio = inlet_flow/elem_field(ne_Vdot,1)
         unit_field(nu_Vdot0,1:num_units) = unit_field(nu_Vdot0,1:num_units)*ratio
         elem_field(ne_Vdot,1:num_elems) = elem_field(ne_Vdot,1:num_elems)*ratio
         elem_field(ne_dvdt,1:num_elems) = elem_field(ne_dvdt,1:num_elems)*ratio
      else
         write(*,'('' Cannot scale to zero flow'')')
      endif
    elseif(VorQ.eq.'Q')then
      if(abs(elem_field(ne_Qdot,1)).gt.zero_tol)then
        ratio = inlet_flow/elem_field(ne_Qdot,1)
        unit_field(nu_perf,1:num_units) = unit_field(nu_perf,1:num_units)*ratio
        elem_field(ne_Qdot,1:num_elems) = elem_field(ne_Qdot,1:num_elems)*ratio
      else
        write(*,'('' Cannot scale to zero flow'')')
      endif
    else
     print*, ' WARNING: Not a ventilation or perfusion scaling, no calculations done'
    endif

    call enter_exit(sub_name,2)
    
  end subroutine scale_flow_to_inlet

!!!#############################################################################

end module field_utilities
