
#include "filenames.h"

#include "string.h"

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

#include <stdlib.h>
#include <stdio.h>

char *get_filename(const char *label)
{
  int label_len = strlen(label);

  printf("get filename: '%s', %d\n", label, label_len);
  char *tmp = (char *)malloc(255 + 1);
  get_filename_c(label, &label_len, tmp);
  //tmp = "bob";
  printf("done: %x, %d\n", tmp, strlen(tmp));
  printf("'%s'\n",tmp);
  return tmp;//get_filename_c(label, &label_len);
}
