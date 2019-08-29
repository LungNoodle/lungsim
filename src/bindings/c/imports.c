
#include "imports.h"
#include "utils.h"

#include <string.h>

void import_ventilation_c(const char *FLOWFILE,int *FLOWFILE_LEN);
void import_perfusion_c(const char *FLOWFILE,int *FLOWFILE_LEN);

void import_ventilation(const char *FLOWFILE)
{
	int filename_len = strlen(FLOWFILE);
	import_ventilation_c(FLOWFILE, &filename_len);
}
void import_perfusion(const char *FLOWFILE)
{
	int filename_len = strlen(FLOWFILE);
	import_perfusion_c(FLOWFILE, &filename_len);
}
