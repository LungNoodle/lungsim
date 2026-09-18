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
  use arrays
  use diagnostics
  use geometry
  use indices
  use other_consts
  use ventilation
  
  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private
  public import_capillary
  public import_terminal
  public import_ventilation
  public import_perfusion

contains
!
!###########################################################################################
!
!> Read mean pressure, transit time, and sheet surface area from micro_flow_unit.out.
 subroutine import_capillary(micro_unit_file)
   !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_IMPORT_CAPILLARY" :: IMPORT_CAPILLARY

   character(len=MAX_FILENAME_LEN),intent(in) :: micro_unit_file
   integer :: count_units,io_unit,ios,ne,nunit
   real(dp) :: Pin,Pout,xyz(3),before_sa(4),TOTAL_SHEET_SA,TT_TOTAL,after_tt(2)

   character(len=60) :: sub_name

   sub_name = 'import_capillary'
   call enter_exit(sub_name,1)

   open(newunit=io_unit, file=trim(micro_unit_file), status='old', action='read', iostat=ios)
   if (ios /= 0) error stop 'Unable to open capillary result file'

   count_units = 0
   do
      read(io_unit, *, iostat=ios) ne,xyz,Pin,Pout,before_sa,TOTAL_SHEET_SA,TT_TOTAL,after_tt
      if (ios < 0) exit
      if (ios > 0) error stop 'Unable to parse capillary result file'

      nunit = int(elem_field(ne_unit,elem_cnct(-1,1,ne)))
      if (nunit < 1 .or. nunit > num_units) error stop 'Invalid unit mapping in capillary result file'
      unit_field(nu_capillary_press,nunit) = (Pin+Pout)/2.0_dp
      unit_field(nu_tt,nunit) = TT_TOTAL
      unit_field(nu_sa,nunit) = TOTAL_SHEET_SA
      count_units = count_units + 1
   enddo

   close(io_unit)

   if (count_units /= num_units) then
      write(*,'(''WARNING: capillary result has '',i0,'' units; geometry has '',i0)') &
           count_units,num_units
   endif

   call enter_exit(sub_name,2)

 end subroutine import_capillary
!
!##############################################################################
!
!>*import_ventilation:* This subroutine reads in the results of a ventilation model that
! has been saved in an exelem format as a single flow field (elements listed with
! ventilation as field values).
 subroutine import_ventilation(FLOWFILE)

   character(len=MAX_FILENAME_LEN),intent(in) :: FLOWFILE
   !local variables
   integer :: ierror,ne,nunit
   character(LEN=132) :: ctemp1,exfile
   real(dp) :: flow,flow_unit,maxflow

   character(len=60) :: sub_name

   sub_name = 'import_ventilation'
   call enter_exit(sub_name,1)

   if(.not.allocated(gasex%Vdot)) allocate(gasex%Vdot(num_units))
   gasex%Vdot = 0.0_dp
   
   print *, 'Reading in ventilation results'
   call import_exelemfield(FLOWFILE,ne_Vdot)
   do nunit = 1,num_units
     ne = units(nunit)
     if(elem_field(ne_Vdot,ne).lt.0.0_dp) elem_field(ne_Vdot,ne) = zero_tol
     unit_field(nu_Vdot0,nunit) = elem_field(ne_Vdot,ne)
     gasex%Vdot(nunit) = elem_field(ne_Vdot,ne)
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

   character(len=MAX_FILENAME_LEN),intent(in) :: FLOWFILE
   !local variables
   integer :: ierror,ne,nunit
   character(LEN=132) :: ctemp1,exfile
   real(dp) :: flow,flow_unit,maxflow

   character(len=60) :: sub_name

   sub_name = 'import_perfusion'
   call enter_exit(sub_name,1)

   if(.not.allocated(gasex%Qdot)) allocate(gasex%Qdot(num_units))
   gasex%Qdot = 0.0_dp
   
   print *, 'Reading in perfusion results'
   call import_exelemfield(FLOWFILE,ne_Qdot)
   do nunit = 1,num_units
     ne = units(nunit)
     if(elem_field(ne_Qdot,ne).lt.0.0_dp) elem_field(ne_Qdot,ne) = zero_tol
     unit_field(nu_perf,nunit) = elem_field(ne_Qdot,ne)
     gasex%Qdot(nunit) = elem_field(ne_Qdot,ne)
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
       flow = get_final_real(ctemp1)
       if(flow.lt.0.0_dp) flow = zero_tol
         elem_field(field_no,ne) = flow! read it in
       end if
       if(ne.ge.num_elems) exit read_elem_flow
     end do read_elem_flow

   close(10)

    call enter_exit(sub_name,2)
 end subroutine import_exelemfield

!
!##############################################################################
!
!> Import terminal-element number and minimum/maximum elastic recoil pressure.
 subroutine import_terminal(EXFILE)
   !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_IMPORT_TERMINAL" :: IMPORT_TERMINAL

   character(len=MAX_FILENAME_LEN),intent(in) :: EXFILE
   integer :: count_units,field_label(10),i,ibeg,iend,ierror,io_unit,nunit,n_fields,terminal_ne
   integer,allocatable :: unit_by_element(:)
   real(dp) :: rtemp
   character(len=132) :: ctemp,label
   character(len=300) :: readfile

   character(len=60) :: sub_name

   sub_name = 'import_terminal'
   call enter_exit(sub_name,1)

   if (index(EXFILE, '.exnode') > 0) then
      readfile = EXFILE
   else
      readfile = trim(EXFILE)//'.exnode'
   endif

   open(newunit=io_unit, file=trim(readfile), status='old', action='read', iostat=ierror)
   if (ierror /= 0) error stop 'Unable to open terminal ventilation result file'

   allocate(unit_by_element(num_elems))
   unit_by_element = 0
   do i = 1,num_units
      unit_by_element(units(i)) = i
   enddo

   n_fields = 0
   ierror = 0
   read_field_labels : do
      read(io_unit, fmt='(a)', iostat=ierror) ctemp
      if (ierror /= 0) error stop 'Terminal ventilation result has no node data'
      if (index(ctemp, 'Node:') > 0) exit read_field_labels
      if (index(ctemp, ') ') > 0) then
         ibeg = index(ctemp, ') ')+1
         iend = index(ctemp, ',')-1
         label = adjustl(ctemp(ibeg:iend))
         n_fields = n_fields + 1
         if (n_fields > size(field_label)) error stop 'Too many terminal fields'
         field_label(n_fields) = 0

         select case(trim(label))
         case('max_Pe')
            field_label(n_fields) = nu_Pe_max
         case('min_Pe')
            field_label(n_fields) = nu_Pe_min
         end select
      endif
   enddo read_field_labels

   count_units = 0
   do
      if (index(ctemp, 'Node:') > 0) then
         do i = 1,3
            read(io_unit, *, iostat=ierror) rtemp
            if (ierror /= 0) exit
         enddo

         read(io_unit, *, iostat=ierror) rtemp
         if (ierror /= 0) exit
         terminal_ne = nint(rtemp)
         if (terminal_ne < 1 .or. terminal_ne > num_elems) &
              error stop 'Terminal result contains an invalid element number'
         nunit = unit_by_element(terminal_ne)
         if (nunit == 0) error stop 'Terminal result does not match perfusion geometry'
         count_units = count_units + 1

         do i = 3,n_fields
            read(io_unit, *, iostat=ierror) rtemp
            if (ierror /= 0) exit
            if (field_label(i) > 0) unit_field(field_label(i),nunit) = rtemp
         enddo
      endif
      if (ierror /= 0) exit
      read(io_unit, fmt='(a)', iostat=ierror) ctemp
      if (ierror /= 0) exit
   enddo

   close(io_unit)
   deallocate(unit_by_element)
   if (count_units /= num_units) then
      write(*,'(''WARNING: terminal result has '',i0,'' units; geometry has '',i0)') count_units,num_units
   endif

   call enter_exit(sub_name,2)

 end subroutine import_terminal

end module imports
