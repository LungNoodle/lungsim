
%module(package="aether") surface_fitting
%include symbol_export.h
%include surface_fitting.h

%{
#include "surface_fitting.h"
%}

void fit_surface_geometry(int niterations, const char *fitting_file);
void initialise_fit_mesh();


