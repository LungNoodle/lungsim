
%module(package="aether") tree_wave_propagation
%include symbol_export.h
%include tree_wave_propagation.h

%{
#include "tree_wave_propagation.h"
%}

%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"

%init %{
  import_array();
%}

%apply (int DIM1, double* IN_ARRAY1) {(int len1, double* vec1),
                                      (int len2, double* vec2)}


%rename (evaluate_wave_propagation) my_evaluate_wave_propagation;

%exception my_evaluate_wave_propagation {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{
void my_evaluate_wave_propagation(double a0, int len1, double* vec1, int len2, double* vec2) {
    if (len1 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "Arrays of lengths (%d,%d) given",
                     len1, len2);
        return;
    }
    evaluate_wave_propagation(a0, len1, vec1, vec2);
}
%}
