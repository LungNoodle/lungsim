
%module(package="aether") filenames
%include symbol_export.h

// Free the returned filename
//%typemap(newfree) char * "if ($1) free($1);";

%newobject get_filename;

%{
#include "filenames.h"
%}

%include filenames.h
