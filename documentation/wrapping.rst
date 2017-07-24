
=================
Wrapping a Module
=================

When we want to make modules or functions/subroutines available from Python we need to wrap the Fortran code.  The method we follow uses SWIG to help automate at least part of this process.

 - Need to create a virtual environment for each build configuration.
 - Need to make functions public that we want to use from say the Python interface.  This is because the C interface makes all functions private by default this means we cannot call or use these functions when we are using the library.  So when the Python interface is using the library it does not have access to private functions and any attempt at trying to use them from Python will cause an error.
  - In the bindings we have four files associated with a fortran module
    - a fortran file
    - a c header file
    - a c implementation file
    - a swig interface file
    
    
    - Declare the function in the c header file with the correct arguments, type and order.
    - Define your function in the c source file with the correct arguments, type and order.
    - Add the wrapper function in the corresponding fortran file
    
    


