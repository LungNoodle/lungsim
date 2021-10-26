%module(package="aether") growtree
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
#include "growtree.h"
  %}

%include growtree.h
