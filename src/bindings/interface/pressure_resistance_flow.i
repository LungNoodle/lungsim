
%module(package="aether") pressure_resistance_flow
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
  if (!PyLong_Check(o)) {
    free($2);
    PyErr_SetString(PyExc_ValueError, "List items must be integers");
    SWIG_fail;
  }
  $2[i] = (int) PyLong_AsLong(o);
}
}

%typemap(in) (int elemlist2_len, int elemlist2[]) {
int i;
if (!PyList_Check($input)) {
  PyErr_SetString(PyExc_ValueError, "Expecting a list");
  SWIG_fail;
}
$1 = PyList_Size($input);
$2 = (int *) malloc(($1)*sizeof(int));
for (i = 0; i < $1; i++) {
  PyObject *o = PyList_GetItem($input, i);
  if (!PyLong_Check(o)) {
    free($2);
    PyErr_SetString(PyExc_ValueError, "List items must be integers");
    SWIG_fail;
  }
  $2[i] = (int) PyLong_AsLong(o);
}
}

%typemap(freearg) (int elemlist_len, int elemlist[]) {
if ($2) free($2);
}

%typemap(freearg) (int elemlist2_len, int elemlist2[]) {
if ($2) free($2);
}

%typemap(arginit) (int element_ids_len, int element_ids[]) {
$2 = NULL;
}

%typemap(in) (int element_ids_len, int element_ids[]) {
int i;
if (!PyList_Check($input)) {
  PyErr_SetString(PyExc_ValueError, "Expecting a list of element IDs");
  SWIG_fail;
}
$1 = PyList_Size($input);
$2 = (int *) malloc(($1)*sizeof(int));
for (i = 0; i < $1; i++) {
  PyObject *o = PyList_GetItem($input, i);
  if (!PyLong_Check(o)) {
    PyErr_SetString(PyExc_ValueError, "Element IDs must be integers");
    SWIG_fail;
  }
  $2[i] = (int) PyLong_AsLong(o);
}
}

%typemap(arginit) (int flow_values_len, double flow_values[]) {
$2 = NULL;
}

%typemap(in) (int flow_values_len, double flow_values[]) {
int i;
if (!PyList_Check($input)) {
  PyErr_SetString(PyExc_ValueError, "Expecting a list of flow values");
  SWIG_fail;
}
$1 = PyList_Size($input);
$2 = (double *) malloc(($1)*sizeof(double));
for (i = 0; i < $1; i++) {
  PyObject *o = PyList_GetItem($input, i);
  if (!PyFloat_Check(o) && !PyLong_Check(o)) {
    PyErr_SetString(PyExc_ValueError, "Flow values must be numbers");
    SWIG_fail;
  }
  $2[i] = PyFloat_AsDouble(o);
}
}

%typemap(freearg) (int element_ids_len, int element_ids[]) {
if ($2) free($2);
}

%typemap(freearg) (int flow_values_len, double flow_values[]) {
if ($2) free($2);
}

%{
#include "pressure_resistance_flow.h"
%}

%include pressure_resistance_flow.h
