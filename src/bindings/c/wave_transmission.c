#include "wave_transmission.h"

void evaluate_wave_transmission_c(int *grav_dirn, double *grav_factor, int *n_time, double *heartrate, double *a0, int *no_freq, double *a, double *b, int *n_bcparams, double *admittance_param,int *n_model, double *model_definition, int *cap_model);

void evaluate_wave_transmission(int grav_dirn, double grav_factor,int n_time, double heartrate, double a0, int no_freq, double *a, double *b, int n_bcparams, double *admittance_param, int n_model, double *model_definition, int cap_model)
{

    evaluate_wave_transmission_c(&grav_dirn, &grav_factor,&n_time, &heartrate, &a0, &no_freq, a, b, &n_bcparams, admittance_param, &n_model, model_definition, &cap_model);
}
