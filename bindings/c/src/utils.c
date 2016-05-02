
#include "utils.h"

#include "memory.h"

char *allocate_fortran_string(size_t size)
{
  char *string = (char *)malloc(size);
  memset(string, ' ', size);
  string[size] = '\0';
  return string;
}
