
#ifndef AETHER_GROWTREE_H
#define AETHER_GROWTREE_H

#include "symbol_export.h"


SHO_PUBLIC void grow_tree(int elemlist_len, int elemlist[], double angle_max, double angle_min, double branch_fraction, double length_limit, double shortest_length, double rotation_limit, const char *grouping);
SHO_PUBLIC void export_terminal_seed_mapping(const char *filename);
SHO_PUBLIC void smooth_1d_tree(int num_elem_start, int n_smoothing_steps);
SHO_PUBLIC void align_terminals_with_seeds(double seed_dist);

#endif /* AETHER_GROWTREE_H */
