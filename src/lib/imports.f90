module imports
!*Brief Description:* This module contains all the subroutines required to
!import fields, previous model results, etc.
!*LICENSE:*
!
!
!
!*Full Description:*
!
!
  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private 
  public import_ventilation
  public import_perfusion

contains
!
!##############################################################################
!
!>*import_ventilation:* This subroutine reads in the results of a ventilation model that
! has been saved in an exelem format as a single flow field (elements listed with
! ventilation as field values).
 subroutine import_ventilation(FLOWFILE)
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_IMPORT_VENTILATION" :: IMPORT_VENTILATION
   use arrays,only: dp,elem_field,num_elems,num_units,units,unit_field,zero_tol,elem_cnct
   use geometry,only: get_final_real
   use indices
   use ventilation,only: sum_elem_field_from_periphery
   use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
   use diagnostics, only: enter_exit

   character(len=MAX_FILENAME_LEN),intent(in) :: FLOWFILE
   !local variables
   integer :: ierror,ne,nunit
   character(LEN=132) :: ctemp1,exfile
   real(dp) :: flow,flow_unit,maxflow

   character(len=60) :: sub_name

   sub_name = 'import_ventilation'
   call enter_exit(sub_name,1)

   print *, 'Reading in ventilation results'
   call import_exelemfield(FLOWFILE,ne_Vdot)
   do nunit = 1,num_units
     ne = units(nunit)
     if(elem_field(ne_Vdot,ne).lt.0.0_dp) elem_field(ne_Vdot,ne) = zero_tol
     unit_field(nu_Vdot0,nunit) = elem_field(ne_Vdot,ne)
   enddo

!!! sum the fields up the tree
   call sum_elem_field_from_periphery(ne_Vdot) !sum the air flows recursively UP the tree
   maxflow = elem_field(ne_Vdot,1)


   call enter_exit(sub_name,2)
 end subroutine import_ventilation

!
!###########################################################################################
!
!>*import_perfusion:* This subroutine reads in the results of a ventilation model that
! has been saved in an exelem format as a single flow field (elements listed with
! ventilation as field values).
 subroutine import_perfusion(FLOWFILE)
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_IMPORT_PERFUSION" :: IMPORT_PERFUSION
   use arrays,only: dp,elem_field,num_elems,num_units,units,unit_field,zero_tol,elem_cnct
   use geometry,only: get_final_real
   use indices
   use ventilation,only: sum_elem_field_from_periphery
   use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
   use diagnostics, only: enter_exit

   character(len=MAX_FILENAME_LEN),intent(in) :: FLOWFILE
   !local variables
   integer :: ierror,ne,nunit
   character(LEN=132) :: ctemp1,exfile
   real(dp) :: flow,flow_unit,maxflow

   character(len=60) :: sub_name

   sub_name = 'import_perfusion'
   call enter_exit(sub_name,1)

   print *, 'Reading in perfusion results'
   call import_exelemfield(FLOWFILE,ne_Qdot)
   do nunit = 1,num_units
     ne = units(nunit)
     if(elem_field(ne_Qdot,ne).lt.0.0_dp) elem_field(ne_Qdot,ne) = zero_tol
     unit_field(nu_perf,nunit) = elem_field(ne_Qdot,ne)
   enddo

!!! sum the fields up the tree
   call sum_elem_field_from_periphery(ne_Qdot) !sum the air flows recursively UP the tree
   maxflow = elem_field(ne_Qdot,1)

   call enter_exit(sub_name,2)
 end subroutine import_perfusion

!
!##############################################################################
!
!>*import_exelemfield:* This subroutine reads in the content of an exelem field file (1 field)
 subroutine import_exelemfield(FLOWFILE,field_no)
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_IMPORT_EXELEMFIELD" :: IMPORT_EXELEMFIELD
   use arrays,only: dp,elem_field,num_elems,num_units,units,unit_field,zero_tol,elem_cnct
   use geometry,only: get_final_real
   use indices
   use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
   use diagnostics, only: enter_exit

   character(len=MAX_FILENAME_LEN),intent(in) :: FLOWFILE
   integer, intent(in) :: field_no
   !local variables
   integer :: ierror,ne,nunit
   character(LEN=132) :: ctemp1,exfile
   real(dp) :: flow,flow_unit,maxflow

   character(len=60) :: sub_name

   sub_name = 'import_exelemfield'
   call enter_exit(sub_name,1)

   open(10, file=FLOWFILE, status='old')
   ne = 0
   read_elem_flow : do !define a do loop name
     !.......read element flow
     read(unit=10, fmt="(a)", iostat=ierror) ctemp1
     if(index(ctemp1, "Values:")> 0) then
       ne = ne+1
       read(unit=10, fmt="(a)", iostat=ierror) ctemp1
       call get_final_real(ctemp1,flow)
       if(flow.lt.0.0_dp) flow = zero_tol
         elem_field(field_no,ne) = flow! read it in
       end if
       if(ne.ge.num_elems) exit read_elem_flow
     end do read_elem_flow

   close(10)

    call enter_exit(sub_name,2)
 end subroutine import_exelemfield

end module imports
