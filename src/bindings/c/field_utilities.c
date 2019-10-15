
#include "field_utilities.h"
#include "utils.h"

#include <string.h>

void scale_flow_to_inlet_c(double *INLET_FLOW,const char *VORQ,int *VQLEN);

void scale_flow_to_inlet(double INLET_FLOW,const char *VORQ)
{
	int vq_len = strlen(VORQ);
	scale_flow_to_inlet_c(&INLET_FLOW, VORQ, &vq_len);
}
