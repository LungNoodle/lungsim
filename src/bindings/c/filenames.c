
#include "filenames.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>


void read_geometry_main_c();
void read_geometry_evaluate_flow_c();
void get_filename_c(const char *label, int *label_len, char *filename);

void read_geometry_main()
{
  read_geometry_main_c();
}

void read_geometry_evaluate_flow()
{
  read_geometry_evaluate_flow_c();
}

char *get_filename(const char *label)
{
  int label_len = strlen(label);

  char *tmp = (char *)malloc(255 + 1);
  get_filename_c(label, &label_len, tmp);
  return tmp;
}
