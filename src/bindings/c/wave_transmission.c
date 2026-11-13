#include "wave_transmission.h"
#include <string.h>

void evaluate_wave_transmission_c(int *grav_dirn, double *grav_factor, int *n_time, double *heartrate, double *a0, int *no_freq, double *a, double *b, int *n_bcparams, double *admittance_param,int *n_model, double *model_definition, int *cap_model, double *remodeling_grade, const char *bc_type, int *bc_type_len, const char *lobe_imped, int *lobe_imped_len);
void export_time_dependent_perfusion_c(const char *flow_filename, int *flow_filename_len, const char *pressure_filename, int *pressure_filename_len);

void evaluate_wave_transmission(int grav_dirn, double grav_factor,int n_time, double heartrate, double a0, int no_freq, double *a, double *b, int n_bcparams, double *admittance_param, int n_model, double *model_definition, int cap_model, double remodeling_grade, const char *bc_type, const char *lobe_imped)
{
    int bc_type_len = strlen(bc_type);
    int lobe_imped_len = strlen(lobe_imped);
    evaluate_wave_transmission_c(&grav_dirn, &grav_factor,&n_time, &heartrate, &a0, &no_freq, a, b, &n_bcparams, admittance_param, &n_model, model_definition, &cap_model, &remodeling_grade, bc_type, &bc_type_len, lobe_imped, &lobe_imped_len);
}

void export_time_dependent_perfusion(const char *flow_filename, const char *pressure_filename)
{
    int flow_filename_len = strlen(flow_filename);
    int pressure_filename_len = strlen(pressure_filename);
    export_time_dependent_perfusion_c(flow_filename, &flow_filename_len, pressure_filename, &pressure_filename_len);
}
