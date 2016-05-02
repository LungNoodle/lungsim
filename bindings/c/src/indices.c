
#include "indices.h"

void ventilation_indices_c();
void perfusion_indices_c();
int get_ne_radius_c();
int get_nj_conc1_c();

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
