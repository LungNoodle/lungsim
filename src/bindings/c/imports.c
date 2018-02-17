
#include "imports.h"
#include "utils.h"

#include <string.h>

void import_terminal_ventilation_c(const char *FLOWFILE,int *FLOWFILE_LEN);

void import_terminal_ventilation(const char *FLOWFILE)
{
	int filename_len = strlen(FLOWFILE);
	import_terminal_ventilation_c(FLOWFILE, &filename_len);
}
