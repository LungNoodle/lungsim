%module(package="aether") geometry
%include symbol_export.h

%{
#include "geometry.h"
%}

// define_rad_from_file has an optional argument that C cannot replicate,
// so we use SWIG to override with a C++ version that can.
void define_rad_from_file(const char *FIELDFILE, const char *radius_type="no_taper");
%include geometry.h

