


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:

  This is roughly what the Makefile should look like. Uncomment the mv 
  statement to store binary in ~/bin

CC = gcc
CFLAGS = -c -g -I$$HOME/include -Wall
LDFLAGS= -L$$HOME/lib -lm -lkpar
BIN = $$HOME/bin

simple: simple.o
	gcc simple.o  $(LDFLAGS) -o simple
#	mv $@ $(BIN)

 

  History:
2004	ksl	Coded as better.c

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "log.h"



#define LINELENGTH 132

int
main (argc, argv)
     int argc;
     char *argv[];
{

  FILE *fptr, *fopen ();
  char parfile[LINELENGTH];
  char infile[LINELENGTH];
  char line[LINELENGTH];
  char firstword[LINELENGTH];
  double f, g;
  double ff, gg;
  int bilen ();
  int i;

  double x00[3], x01[3], x10[3], x11[3], x[3];
  double xx[3];
/* Parse the command line */
  if (argc == 1)
    strcpy (parfile, "simple.pf");
  else if (argc == 2)
    {
      strcpy (parfile, argv[1]);
      if ((strstr (parfile, ".pf")) == NULL)
	strcat (parfile, ".pf");
    }
  else
    {
      printf ("Usage: simple [file.pf]\n");
      exit (0);
    }

/* Initialize variables that will be read by rdpar */

  strcpy (infile, "file.ls");

  x00[0] = -1.0;
  x00[1] = -1.0;
  x01[0] = -1.0;
  x01[1] = 1.0;
  x10[0] = 1.0;
  x10[1] = -1.0;
  x11[0] = 1.0;
  x11[1] = 1.0;

  x[0] = 0.5;
  x[1] = 0.5;
/* Get data from a parameter file */


  opar (parfile);


  rddoub ("x00[0]", &x00[0]);
  rddoub ("x00[2]", &x00[2]);

  rddoub ("x01[0]", &x01[0]);
  rddoub ("x01[2]", &x01[2]);

  rddoub ("x10[0]", &x10[0]);
  rddoub ("x10[2]", &x10[2]);

  rddoub ("x11[0]", &x11[0]);
  rddoub ("x11[2]", &x11[2]);




  for (ff = -0.3; ff <= 1.3; ff += 0.1)
    {
      for (gg = -0.3; gg <= 1.3; gg += 0.1)
	{


	  x[0] =
	    (1. - gg) * ((1 - ff) * x00[0] + ff * x10[0]) +
	    gg * ((1 - ff) * x01[0] + ff * x11[0]);
	  x[2] =
	    (1. - gg) * ((1 - ff) * x00[2] + ff * x10[2]) +
	    gg * ((1 - ff) * x01[2] + ff * x11[2]);

/* End of input section */

	  i = bilin (x, x00, x01, x10, x11, &f, &g);

	  xx[0] =
	    (1. - g) * ((1 - f) * x00[0] + f * x10[0]) +
	    g * ((1 - f) * x01[0] + f * x11[0]);
	  xx[2] =
	    (1. - g) * ((1 - f) * x00[2] + f * x10[2]) +
	    g * ((1 - f) * x01[2] + f * x11[2]);

	  if (0.0 <= f && f <= 1 && 0.0 <= g && g <= 1 && i != 0)
	    Error ("Next answer unexpected == Supposed to be in grid \n");
	  if (i == 0 && (f < 0.0 || f > 1 || g < 0.0 || g > 1.0))
	    Error ("Next answer unexpected -- Not supposed to be in grid\n");

	  printf
	    ("%8.3f %8.3f --> %9.2e %9.2e --> %9.3f %9.3f (Returned %2d) --> %9.2e %9.2e (Diff %9.2e %9.2e) \n",
	     ff, gg, x[0], x[2], f, g, i, xx[0], xx[2], xx[0] - x[0],
	     xx[2] - x[2]);

	}
    }




/* OK now end program and go on to the next important task. */
b:
  cpar ("t_bilinear.pf");
  exit (0);

}
