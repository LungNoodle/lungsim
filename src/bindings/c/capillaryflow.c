
#include "capillaryflow.h"

#include "string.h"

void calc_cap_imped_c(double *ha, double *hv, double *omega);


void calc_cap_imped(double ha, double hv, double omega)
{
  calc_cap_imped_c(&ha, &hv, &omega);
}
