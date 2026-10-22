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

%rename(_grow_tree_without_export) grow_tree;
%rename(_grow_tree_with_export) grow_tree_with_export;

%{
#include "growtree.h"
  %}

%include growtree.h

%pythoncode %{
def grow_tree(parent_elems, angle_max, angle_min, branch_fraction, length_limit,
              shortest_length, rotation_limit, *args, **kwargs):
    """Grow a tree, optionally exporting the seed-to-terminal element mapping.

    Supported positional endings are ``grouping``,
    ``grouping, to_export, filename``, and the historical
    ``to_export, filename, grouping``. The options can also be supplied as
    ``grouping=``, ``to_export=``, and ``filename=`` keyword arguments.
    """
    unexpected = set(kwargs) - {"grouping", "to_export", "filename"}
    if unexpected:
        name = sorted(unexpected)[0]
        raise TypeError("grow_tree() got an unexpected keyword argument %r" % name)

    grouping = kwargs.get("grouping")
    to_export = kwargs.get("to_export", False)
    filename = kwargs.get("filename", "")

    if len(args) == 1:
        if grouping is not None:
            raise TypeError("grow_tree() got multiple values for 'grouping'")
        grouping = args[0]
    elif len(args) == 3:
        if any(name in kwargs for name in ("grouping", "to_export", "filename")):
            raise TypeError("grow_tree() received duplicate export options")
        if isinstance(args[0], str):
            grouping, to_export, filename = args
        else:
            to_export, filename, grouping = args
    elif len(args) != 0:
        raise TypeError("grow_tree() expects a grouping value and optional export arguments")

    if grouping is None:
        raise TypeError("grow_tree() missing required argument: 'grouping'")
    if to_export and not filename:
        raise ValueError("filename is required when to_export is True")

    if to_export:
        return _grow_tree_with_export(parent_elems, angle_max, angle_min,
                                      branch_fraction, length_limit,
                                      shortest_length, rotation_limit, grouping,
                                      True, filename)
    return _grow_tree_without_export(parent_elems, angle_max, angle_min,
                                     branch_fraction, length_limit,
                                     shortest_length, rotation_limit, grouping)
%}
