
#ifndef AETHER_EXPORTS_H
#define AETHER_EXPORTS_H

#include "symbol_export.h"

SHO_PUBLIC void export_1d_elem_field(int ne_field, const char *EXELEMFILE, const char *group_name, const char *field_name );
SHO_PUBLIC void export_1d_elem_geometry(const char *EXELEMFILE, const char *name);
SHO_PUBLIC void export_elem_geometry_2d(const char *EXELEMFILE, const char *name, int offset_elem, int offset_node);
SHO_PUBLIC void export_node_field(int nj_field, const char *EXNODEFIELD, const char *name, const char *field_name);
SHO_PUBLIC void export_terminal_solution(const char *EXNODEFILE, const char *name);
SHO_PUBLIC void export_terminal_perfusion(const char *EXNODEFILE, const char *name);
SHO_PUBLIC void export_terminal_ssgexch(const char *EXNODEFILE, const char *name);
SHO_PUBLIC void export_node_geometry(const char *EXNODEFILE, const char *name);
SHO_PUBLIC void export_node_geometry_2d(const char *EXNODEFILE, const char *name, int offset);
SHO_PUBLIC void export_data_geometry(const char *EXDATAFILE, const char *name, int offset);
SHO_PUBLIC void export_elem_field(const char *EXELEMFIELD, const char *name, const char *field_name);

#endif /* AETHER_EXPORTS_H */
