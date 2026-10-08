
#include "growtree.h"

#include "string.h"

#include "stdio.h"

extern void grow_tree_c(int *elemlist_len, int elemlist[], double *angle_max, double *angle_min, double *branch_fraction, double *length_limit, double *shortest_length, double *rotation_limit, const char *grouping, int *grouping_len);
extern void grow_tree_with_export_c(int *elemlist_len, int elemlist[], double *angle_max, double *angle_min, double *branch_fraction, double *length_limit, double *shortest_length, double *rotation_limit, const char *grouping, int *grouping_len, int *to_export, const char *filename, int *filename_len);

extern void smooth_1d_tree_c(int *num_elem_start, int *n_smoothing_steps);
extern void align_terminals_with_seeds_c(double *seed_dist);

void grow_tree(int elemlist_len, int elemlist[], double angle_max, double angle_min, double branch_fraction, double length_limit, double shortest_length, double rotation_limit, const char *grouping)
{
  int grouping_len = strlen(grouping);
  grow_tree_c(&elemlist_len, elemlist, &angle_max, &angle_min, &branch_fraction, &length_limit, &shortest_length, &rotation_limit, grouping, &grouping_len);
}

void grow_tree_with_export(int elemlist_len, int elemlist[], double angle_max, double angle_min, double branch_fraction, double length_limit, double shortest_length, double rotation_limit, const char *grouping, int to_export, const char *filename)
{
  int grouping_len = strlen(grouping);
  int filename_len = strlen(filename);
  grow_tree_with_export_c(&elemlist_len, elemlist, &angle_max, &angle_min, &branch_fraction, &length_limit, &shortest_length, &rotation_limit, grouping, &grouping_len, &to_export, filename, &filename_len);
}

void smooth_1d_tree(int num_elem_start, int n_smoothing_steps)
{
  smooth_1d_tree_c(&num_elem_start, &n_smoothing_steps);
}

void align_terminals_with_seeds(double seed_dist)
{
  align_terminals_with_seeds_c(&seed_dist);
}
