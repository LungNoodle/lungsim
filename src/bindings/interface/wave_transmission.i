
%module(package="aether") wave_transmission
%include symbol_export.h
%include wave_transmission.h

%{
#include "wave_transmission.h"
%}

%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"

%init %{
  import_array();
%}

%apply (int DIM1, double* IN_ARRAY1) {(int len1, double* vec1),
                                      (int len2, double* vec2),
                                      (int len3, double* vec3),
									  (int len4, double* vec4)}


%rename (evaluate_wave_transmission) my_evaluate_wave_transmission;

%exception my_evaluate_wave_transmission {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{
void my_evaluate_wave_transmission(int grav_dirn, double grav_factor,int n_time, double heartrate, double a0, int len1, double* vec1, int len2, double* vec2, int len3, double* vec3, int len4, double* vec4, int cap_model, double remodeling_grade, const char *bc_type, const char *lobe_imped) {
    if (len1 != len2) {
        PyErr_Format(PyExc_ValueError,
                     "Arrays of lengths (%d,%d) given",
                     len1, len2);
        return;
    }
    evaluate_wave_transmission(grav_dirn, grav_factor, n_time, heartrate, a0, len1, vec1, vec2, len3, vec3, len4, vec4, cap_model, remodeling_grade, bc_type, lobe_imped);
}
%}
