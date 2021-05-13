%module(package="aether") geometry
  %include symbol_export.h
  
  %typemap(in) (int elemlist_len, int elemlist[]) {
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    SWIG_fail;
  }
  $1 = PyList_Size($input);
  $2 = (int *) malloc(($1)*sizeof(int));
  for (i = 0; i < $1; i++) {
    PyObject *o = PyList_GetItem($input, i);
    if (!PyInt_Check(o)) {
      free($2);
      PyErr_SetString(PyExc_ValueError, "List items must be integers");
      SWIG_fail;
    }
    $2[i] = PyInt_AsLong(o);
  }
 }

%typemap(freearg) (int elemlist_len, int elemlist[]) {
  if ($2) free($2);
 }

%{
#include "geometry.h"
  %}

// define_rad_from_file has an optional argument that C cannot replicate,
// so we use SWIG to override with a C++ version that can.
void define_elem_geometry_2d(const char *ELEMFILE, const char *sf_option="arcl");
void define_node_geometry_2d(const char *NODEFILE);
void import_node_geometry_2d(const char *NODEFILE);
void write_elem_geometry_2d(const char *ELEMFILE);
void write_geo_file(int ntype, const char *GEOFILE);
void write_node_geometry_2d(const char *NODEFILE);
void define_rad_from_file(const char *FIELDFILE, const char *radius_type="no_taper");
void define_rad_from_geom(const char *ORDER_SYSTEM, double CONTROL_PARAM, const char *START_FROM, double START_RAD, const char *group_type_in="all", const char *group_option_in="dummy");


%include geometry.h

