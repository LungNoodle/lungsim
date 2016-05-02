%module(package="aether") geometry
%include symbol_export.h

%{
#include "geometry.h"
%}

void define_rad_from_file(const char *FIELDFILE, const char *radius_type="no_taper");
%include geometry.h

