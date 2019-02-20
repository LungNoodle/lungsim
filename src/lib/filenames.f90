!> \file
!> \author Merryn Tawhai, Alys Clark
!> \brief This module handles all read and write filenames.
!>
!> \section LICENSE
!>
!>
!> Contributor(s):
!>

!> This module handles all read and write filenames.
module filenames
  use other_consts, only : MAX_FILENAME_LEN
  implicit none
 !Module parameters

  CHARACTER(LEN=MAX_FILENAME_LEN) :: AIRWAY_ELEMFILE,AIRWAY_NODEFILE,AIRWAY_FIELDFILE, &
       AIRWAY_EXNODEFILE,AIRWAY_EXELEMFILE,AIRWAYFIELD_EXELEMFILE,&
       AIRWAY_MESHFILE, &
       ARTERY_ELEMFILE,ARTERY_NODEFILE,ARTERY_FIELDFILE, &
       ARTERY_EXNODEFILE,ARTERY_EXELEMFILE,ARTERYFIELD_EXELEMFILE,&
       ARTERY_MESHFILE, &
       TERMINAL_EXNODEFILE, MAIN_GEOMETRY_FILE, &
       FLOW_GEOMETRY_FILE,MAIN_PARAMETER_FILE,FLOW_PARAMETER_FILE,&
       FLOW_EXELEMFILE, FLOW_EXNODEFILE, FLOW_RADIUS_EXELEM, EMPTY_FILENAME

public AIRWAY_ELEMFILE,AIRWAY_NODEFILE,AIRWAY_FIELDFILE, &
       AIRWAY_EXNODEFILE,AIRWAY_EXELEMFILE,AIRWAYFIELD_EXELEMFILE,&
       AIRWAY_MESHFILE, &
       ARTERY_ELEMFILE,ARTERY_NODEFILE,ARTERY_FIELDFILE, &
       ARTERY_EXNODEFILE,ARTERY_EXELEMFILE,ARTERYFIELD_EXELEMFILE,&
       ARTERY_MESHFILE, &
       TERMINAL_EXNODEFILE, MAIN_GEOMETRY_FILE, &
       FLOW_GEOMETRY_FILE,MAIN_PARAMETER_FILE,FLOW_PARAMETER_FILE,&
       FLOW_EXELEMFILE, FLOW_EXNODEFILE, FLOW_RADIUS_EXELEM
  !Module types

  !Module variables

  !Interfaces

  private
  public read_geometry_main
  public read_geometry_evaluate_flow
  public get_filename

contains
!
!###################################################################################
!
!> reads in output filenames typically used to analyse and visualise ventilation model results
  subroutine read_geometry_evaluate_flow()
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_READ_GEOMETRY_EVALUATE_FLOW" :: READ_GEOMETRY_EVALUATE_FLOW

    use diagnostics, only: enter_exit
    implicit none

    ! Input related variables
    character(len=255) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios
    integer :: line
    character(len=60) :: sub_name

    sub_name = 'read_geometry_evaluate_flow'
    call enter_exit(sub_name,1)
    
    ios = 0
    line = 0

    open(fh, file='Parameters/geometry_evaluate_flow.txt')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, '    ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          select case (label)
          case ('exelem')
             read(buffer, *, iostat=ios) AIRWAYFIELD_EXELEMFILE
          case ('exnode')
             read(buffer, *, iostat=ios) TERMINAL_EXNODEFILE
          case ('flowexelem')
             read(buffer, *, iostat=ios) FLOW_EXELEMFILE
          case ('flowexnode')
             read(buffer, *, iostat=ios) FLOW_EXNODEFILE
          case ('flowradiusexelem')
             read(buffer, *, iostat=ios) FLOW_RADIUS_EXELEM
          case default
             !print *, 'Skipping invalid label at line', line
          end select
       end if
    end do

    close(fh)
    call enter_exit(sub_name,2)

  end subroutine read_geometry_evaluate_flow

  !###################################################################################

  subroutine read_geometry_main()
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_READ_GEOMETRY_MAIN" :: READ_GEOMETRY_MAIN

    use diagnostics, only: enter_exit
    implicit none

    ! Input related variables
    character(len=255) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios
    integer :: line
    character(len=60) :: sub_name
!
! ###########################################################################
!
!> Reads in input and output filenames required to generate and export a geometry
    sub_name = 'read_geometry_main'
    call enter_exit(sub_name,1)

    ios = 0
    line = 0
    open(fh, file='Parameters/geometry_main.txt')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, '    ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          select case (label)
          case ('airway_ipnode')
             read(buffer, *, iostat=ios) AIRWAY_NODEFILE
          case ('airway_ipelem')
             read(buffer, *, iostat=ios) AIRWAY_ELEMFILE
          case ('airway_ipfiel')
             read(buffer, *, iostat=ios) AIRWAY_FIELDFILE
          case ('airway_ipmesh')
             read(buffer, *, iostat=ios) AIRWAY_MESHFILE
          case ('airway_exnode')
             read(buffer, *, iostat=ios) AIRWAY_EXNODEFILE
          case ('airway_exelem')
             read(buffer, *, iostat=ios) AIRWAY_EXELEMFILE
          case ('artery_ipnode')
             read(buffer, *, iostat=ios) ARTERY_NODEFILE
          case ('artery_ipelem')
             read(buffer, *, iostat=ios) ARTERY_ELEMFILE
          case ('artery_ipfiel')
             read(buffer, *, iostat=ios) ARTERY_FIELDFILE
          case ('artery_ipmesh')
             read(buffer, *, iostat=ios) ARTERY_MESHFILE
          case ('artery_exnode')
             read(buffer, *, iostat=ios) ARTERY_EXNODEFILE
          case ('artery_exelem')
             read(buffer, *, iostat=ios) ARTERY_EXELEMFILE
          case default
             !print *, 'Skipping invalid label at line', line
          end select
       end if
    end do

    close(fh)
    call enter_exit(sub_name,2)

  end subroutine read_geometry_main
!
!#####################################################################################################
!
  function get_filename(label) result(str)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_GET_FILENAME" :: GET_FILENAME
    use other_consts, only: MAX_FILENAME_LEN, MAX_STRING_LEN
    implicit none

    character(len=MAX_STRING_LEN), intent(in) :: label
    character(len=MAX_FILENAME_LEN) :: str

    select case (label)
    case ('airway_ipnode')
       str = AIRWAY_NODEFILE
    case ('airway_ipelem')
       str = AIRWAY_ELEMFILE
    case ('airway_ipfiel')
       str = AIRWAY_FIELDFILE
    case ('airway_ipmesh')
       str = AIRWAY_MESHFILE
    case ('airway_exnode')
       str = AIRWAY_EXNODEFILE
    case ('airway_exelem')
       str = AIRWAY_EXELEMFILE
    case ('artery_ipnode')
       str = ARTERY_NODEFILE
    case ('artery_ipelem')
       str = ARTERY_ELEMFILE
    case ('artery_ipfiel')
       str = ARTERY_FIELDFILE
    case ('artery_ipmesh')
       str = ARTERY_MESHFILE
    case ('artery_exnode')
       str = ARTERY_EXNODEFILE
    case ('artery_exelem')
       str = ARTERY_EXELEMFILE
    case ('exelem')
       str =  AIRWAYFIELD_EXELEMFILE
    case ('exnode')
       str =  TERMINAL_EXNODEFILE
    case ('flowexelem')
       str =  FLOW_EXELEMFILE
    case ('flowexnode')
       str =  FLOW_EXNODEFILE
    case ('flowradiusexelem')
       str =  FLOW_RADIUS_EXELEM
    case default
       print *, 'Umm I dont know this label, sorry'
       str = EMPTY_FILENAME
    end select

  end function get_filename
end module filenames
