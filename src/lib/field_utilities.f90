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
  use other_consts
  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private
  public scale_flow_to_inlet

contains

!
!*scale_flow_field* Scales a flow field to an 'inlet flow' value (real units).
  subroutine scale_flow_to_inlet(inlet_flow,VorQ)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SCALE_FLOW_TO_INLET" :: SCALE_FLOW_TO_INLET
    use arrays,only: dp,elem_field,num_elems,num_units,unit_field,zero_tol
    use indices,only: ne_dvdt,ne_Vdot,nu_Vdot0,ne_Qdot,nu_perf
    use diagnostics, only: enter_exit

    real(dp),intent(in) :: inlet_flow
    character(len=1), intent(in) :: VorQ
    real(dp) :: ratio
    character(len=60) :: sub_name

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

end module field_utilities
