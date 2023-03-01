
=================
Wrapping a module
=================

When we want to make modules or functions/subroutines available from Python we cannot use the Fortran code directly, we need to wrap it in some way.  Because it is not a straightforward task to wrap the Aether library, which is written in Fortran, directly into Python we use a tool called Simplified Wrapping and Interface Generator (SWIG).  SWIG performs a lot of the work transforming Python data structures; integers, floats, strings, etc.. into C or C++ for us.  However since Aether is written in Fortran we still have a bit of work to do, namely making the Fortran code available from C.  It is then our task to wrap the Aether library into a C interface and set up SWIG to complete the process.

This documentation describes how we create C bindings to the Fortran code and then how to get SWIG to wrap the code into Python.

C bindings
==========

There are three files associated with a single module in the Aether library that are required for making the C bindings.  These files are:

   - A C header file
   - A C implementation file
   - A Fortran implementation file

The C header file is used by SWIG to generate the wrapped functions.  The C implementation file starts the procedure of converting the C variables into Fortran and the Fortran implementation file completes the conversion of the C variables and calls the functions in the Aether library.

If we are adding a new module follow the instructions in :doc:`adding_modules` and then return here to implement the wrapping.

C header file
-------------

The C header file contains the definitions of all the functions that are public from the current module.  The less functions made public the better, once a function is public and other developers or users have started using it changes to the function will break their code.  Keep this in mind when determining whether a function should be available in Python or not.

When a new C header file is created there is a little bit of housekeeping to perform.  We must add a little C to avoid compilation issues.  The code we require is::

    #ifndef AETHER_MODULE_NAME_H
    #define AETHER_MODULE_NAME_H

    #include "symbol_export.h"

    /* Declaration of public functions goes here. */
    
    #endif /* AETHER_MODULE_NAME_H */

where in the above snippet *MODULE_NAME* is replaced with the actual name of the module written in capital letters.

The #ifndef, #define, and #endif directives make sure we don't try and declare our functions multiple times.  The #include "symbol_export.h" declares a macro that we will use to control the visibility of the public functions.

Once the C header file has got it's house in order we can start to declare the public functions we are going to use from Python.  Below are a few examples, first a selection of subroutines with different arguments from the Aether library and then their equivalent declaration in C.

Fortran subroutines/functions::

    subroutine set_diagnostics_on(state)
    function get_ne_radius() result(res)
    subroutine add_mesh(AIRWAY_MESHFILE)

equivalent C function declarations::

    SHO_PUBLIC void set_diagnostics_on(int state);
    SHO_PUBLIC int get_ne_radius();
    SHO_PUBLIC void add_mesh(const char *AIRWAY_MESHFILE);

where the *SHO_PUBLIC* macro makes the functions visible to other programs.  

We can see here that the arguments and the return types are comparable across languages the naming of the functions in the two different languages are exacty the same.  This is by design, we want to keep the usage of the code as similar as possible across the different languages.  This similarity will flow through to Python  when SWIG creates the bindings.

The table below shows a few comparison of Fortran and C datatypes:

=========  ===================
Fortran    C
=========  ===================
integer    long int or int
---------  -------------------
logical    long int or int
---------  -------------------
real       float
---------  -------------------
real*8     double
---------  -------------------
character  char *
=========  ===================

Keep this table in mind when working out which data type best matches your arguments when creating C bindings for you own functions.

C implementation file
=====================

The C implemention file is where we start to handover to the Fortran language.
The first task in creating a C implementation file is to include the associated header file for our example module *module_name* we would include the header file at the very top of the file like so::

    #include module_name.h

After the corresponding header file is included we can include any system headers that we may need, for instance if we are working with characters or strings we would include the *string.h* header file **after** *module_name.h*.

The second task is to declare the C form of the Fortran binding function.
We give this function the same name as the original function but with '_c' appended to it.  For the examples above this would look like::

    void set_diagnostics_on_c(int *state);
    int get_ne_radius_c();
    void add_mesh_c(const char *AIRWAY_MESHFILE, int *filename_len);

There are two things here that we have to take care of; The first is that all arguments in Fortran are passed by reference and not by value.
Thus C must pass Fortran arguments as a pointer, The second is that strings or character arrays are dealt with quite differently in C and Fortran.
We will explain more as we go further.

The third task we must perform is the implementation of the C function that calls the corresponding Fortran function that we have just declared (but not yet implemented).
Let's look at the implementation of our example functions::

    void set_diagnostics_on(int state)
    {
      set_diagnostics_on_c(&state);
    }

    int get_ne_radius()
    {
      return get_ne_radius_c();
    }

    void add_mesh(const char *AIRWAY_MESHFILE)
    {
      int filename_len = strlen(AIRWAY_MESHFILE);
      add_mesh_c(AIRWAY_MESHFILE, &filename_len);
    }

In *set_diagnostics_on* we simply pass the argument *state* by reference to the corresponding Fortran function.
The returned integer from *get_ne_radius_c* is already ready for us to use in C.
The only real work we have to do here is calculate the length of the string we are passing to the Fortran function as Fortran character arrays have no notion of a termination character to signal the end of a string.
There are other differences but we can make use of some utility functions to hide most of the differences from us.

The standard we are using for adding the length of string argument is to add it directly after the string argument in the function argument list.
It then follows that if we have more than one string argument or mixed string and value arguments then the string argument is always followed by it's length argument. 

Fortran implementation file
---------------------------

The Fortran implmentation file is where the majority of the work is done.
We have to tell compilers what to bind the Fortran function name to so the C compiler can locate the function when linking.
We also have to implment the conversion from C char pointers to Fortran character arrays.

The first task we have to do is setup the module, for our example module *module_name* we would write the following::

    module module_name_c

    implicit none
    private

    contains

    ! module subroutines declared here

    end module module_name_c

The second task is to implement the subroutine that will call into the corresponding subroutine in the Aether library that we are binding.
For our example functions we have::

    !
    !###################################################################################
    !
      subroutine set_diagnostics_on_c(state) bind(C, name="set_diagnostics_on_c")
        use diagnostics, only: set_diagnostics_on
        implicit none

        logical, intent(in) :: state

        call set_diagnostics_on(state)

      end subroutine set_diagnostics_on_c

    !
    !###################################################################################
    !
      function get_ne_radius_c() result(res) bind(C, name="get_ne_radius_c")

        use indices, only: get_ne_radius
        implicit none
        integer :: res

        res = get_ne_radius()

      end function get_ne_radius_c

    !
    !###################################################################################
    !
      subroutine add_mesh_c(AIRWAY_MESHFILE, filename_len) bind(C, name="add_mesh_c")
        use iso_c_binding, only: c_ptr
        use utils_c, only: strncpy
        use other_consts, only: MAX_FILENAME_LEN
        use geometry, only: add_mesh
        implicit none

        integer,intent(in) :: filename_len
        type(c_ptr), value, intent(in) :: AIRWAY_MESHFILE
        character(len=MAX_FILENAME_LEN) :: filename_f

        call strncpy(filename_f, AIRWAY_MESHFILE, filename_len)

        call add_mesh(filename_f)

      end subroutine add_mesh_c

We can see that on the function/subroutine declaration we have added the *bind(C)* attribute.
This attribute tells the compiler that this symbol must be operable with C.
With this attribute we also set the name of the symbol that we want to be able to find from C.
This name matches the name of the function we declared at the top of the C implementation file.

The last thing we need to consider is the way that C string is dealt with in *add_mesh_c*.
We have to be careful when converting from C to Fortran but we can make use of the *strncpy* utility to make life easier.
When we are converting strings from C to Fortran, like in this situation, we will replicate what is shown here and accept that it works.
The details of the conversion are not going to be explained here.

SWIG interface
==============

When creating a new module we need to create an interface file so that SWIG creates a corresponding module in the target language.
The interface file is typically very simple but we can add some directives in this file to help map from C to the target language and vice versa.
In the simplest case we just describe the interface using the C header file.
For our example module *module_name* the interface file looks like the following::

    %module(package="aether") module_name
    %include symbol_export.h
    %include module_name.h

    %{
    #include "module_name.h"
    %}

Here we declare the package that we want this module to belong to (*aether* in this case) and the name of the module.
Then we define the files that SWIG needs to create the bindings from and lastly a C part that defines the header files that are required for compilation.
