/* This is a dummy routine for unix */

#include <stdlib.h>
#include <stdio.h>

int 
DebugStr (char *string)
{
  printf ("Debug string: %s\n", string);
  exit (0);
}
