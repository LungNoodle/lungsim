
========================
Fortran style guidelines
========================

Here are some examples of fortran programming style to follow 

Preamble
========
::

    module module name
    !*Brief Description:* A one line descriptor of what the module does. 
    !
    !*LICENSE:*
    !
    !
    !
    !*Full Description:*
    !More info on what the module does if necessary
    !
      use anyothermodulesused
      implicit none
  
      !Module parameters
  
      !Module types

      !Module variables

      !Interfaces
      private list1
      public list2
      public list3

    contains
    !
    !##############################################################################
    !

Indenting
=========
Please indent using two spaces.

Subroutines
===========
Include subroutines in alphabetical order.

Separate subroutines as follows
::

    !
    !##############################################################################
    !

Begin each subroutine with its name and a brief description of its functionality
::

    !*subname:* What the subroutine does

Subroutine contents:
::

    subroutine subname(things,you,pass)
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SUBNAME: SUBNAME
      use module1,only: stuff,from,module1 
      use module2
      use diagnostics, only: enter_exit

      integer, intent(in):: things
      integer, intent(out) :: you
      real(dp), intent(inout)  :: pass
      !local variables
      integer :: other,stuff

      character(len=60) :: sub_name

      sub_name = 'subname'
      call enter_exit(sub_name,1)

      !BODY OF CODE

      call enter_exit(sub_name,2)
    end subroutine subname

Functions
=========
Include functions in alphabetical order below subroutines.
