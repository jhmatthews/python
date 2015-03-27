


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
#include "atomic.h"

#define NWAVES 1000

#include "log.h"

#define LINELENGTH 132
#define NTOTAL  10000000

int
main (argc, argv)
     int argc;
     char *argv[];
{

  FILE *fptr, *fopen ();
  char parfile[LINELENGTH];
//OLD  char infile[LINELENGTH];
//OLD  char line[LINELENGTH];
//OLD  char firstword[LINELENGTH];
  double lambda_min, lambda_max, t;
  double freq_min, freq_max;
  double freq;
  double planck (), timer ();
  double dlambda, lambda, nphot[NWAVES];
  int i, n, out_of_bounds, in_bounds;

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

  out_of_bounds = 0;
  in_bounds = 0;
  for (n = 0; n < NWAVES; n++)
    nphot[n] = 0.0;


/* Get data from a parameter file */


  opar (parfile);


  rddoub ("lambda_min", &lambda_min);
  rddoub ("lambda_max", &lambda_max);
  lambda_min *= 1e-8;		// Convert to Angstroms
  lambda_max *= 1e-8;		// Convert to Angstroms
  rddoub ("t", &t);

  cpar (parfile);

/* End of input section */
  freq_min = C / lambda_max;
  freq_max = C / lambda_min;
  dlambda = (lambda_max - lambda_min) / NWAVES;


  timer ();
  for (n = 0; n < NTOTAL; n++)
#define ntotal
#define ntotal
    {
      if (n % 10000 == 0)
	printf ("%10d  %6.2f   %6.2e\n", n, n * 100. / NTOTAL, timer () / n);
      t = t + 1;
      freq = planck (t, freq_min, freq_max);
      if (freq < freq_min || freq_max < freq)
	{
	  out_of_bounds++;
	  if (out_of_bounds < 10)
	    {
	      Error ("Freq returned is out of bounds %.2e    %.2e %.2e\n",
		     freq, freq_min, freq_max);
	    }
	}
      else
	{
	  in_bounds++;
	  lambda = C / freq;
	  i = (lambda - lambda_min) / dlambda;
	  if (i < 0 || i >= NWAVES)
	    {
	      Error ("i does not fit\n");
	    }
	  else
	    nphot[i] += 1;
	}
    }

  fptr = fopen ("foo.out", "w");
  for (n = 0; n < NWAVES; n++)
    {
      fprintf (fptr, "%f %f\n", (lambda_min + n * dlambda) * 1e8, nphot[n]);
    }


  printf ("Summary: in_bounds %d out_of_bounds %d \n", in_bounds,
	  out_of_bounds);
  printf ("Time to exeucte: %f\n", timer ());
/* OK now end program and go on to the next important task. */
  exit (0);

}
