
#ifndef AETHER_GROWTREE_H
#define AETHER_GROWTREE_H

#include "symbol_export.h"


SHO_PUBLIC void grow_tree(int elemlist_len, int elemlist[], int parent_ne, int supernumerary_ne, double angle_max, double angle_min, double branch_fraction, double length_limit, double shortest_length, double rotation_limit,int to_export, const char *filename, const char *grouping);
SHO_PUBLIC void smooth_1d_tree(int num_elem_start, double length_limit);

#endif /* AETHER_GROWTREE_H */
