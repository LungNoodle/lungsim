
#include "indices.h"
#include "utils.h"
#include <string.h>

void define_problem_type_c(const char *PROBLEMTYPE,int *PROBLEMTYPE_LEN);
void ventilation_indices_c();
void perfusion_indices_c();
int get_ne_radius_c();
int get_nj_conc1_c();



void define_problem_type(const char *PROBLEMTYPE)
{
	int filename_len = strlen(PROBLEMTYPE);
	define_problem_type_c(PROBLEMTYPE, &filename_len);
}


void ventilation_indices()
{
  ventilation_indices_c();
}

void perfusion_indices()
{
  perfusion_indices_c();
}

int get_ne_radius()
{
  return get_ne_radius_c();
}

int get_nj_conc1()
{
  return get_nj_conc1_c();
}
