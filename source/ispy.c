#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "atomic.h"

#include "python.h"

/***********************************************************
                                       Space Telescope Science Institute

 Synopsis:

 int ispy(p,n)  collects more detailed spectral information 
	about certain hardwired cells

Arguments:		

	PhotPtr	p	A photon
	int	n	The current cell number

Returns:
 
 
Description:	

This routine identifies photons as they enter a certain cell (or cells) and writes information
about these photons to a file.
		
Notes:

History:
 	98nov	ksl	Coded and debugged as part of Python effort to try to understand more
			about spectra in individual cells 
	99oct	ksl	Modified to give error message and quit if ispy file is not open.
**************************************************************/

int ispy_start = 0;
int ispy_cycle;
FILE *ispy_ptr;
int ispy_cell[] = { 91, 92, 93, 121, 122, 123 };

int ispy_ncells = 6;

int
ispy_init (filename, icycle)
     char filename[];
     int icycle;
{
  FILE *fopen ();
  char wholename[30];

  strcpy (wholename, filename);
  strcat (wholename, ".ispy");

  ispy_cycle = icycle;

  if (ispy_start == 0)
    {
      ispy_ptr = fopen (wholename, "w");
      ispy_start = 1;
    }
  return (0);
}

int
ispy_close ()
{
  if (ispy_start == 0)
    {
      Log ("ispy file is not open, why are you trying to close it?\n");
    }
  else
    {
      fclose (ispy_ptr);
      ispy_start = 0;
    }
  return (0);
}

int ispy_grid_old = -1;		// The grid cell in which a photon was laast

int ispy_phot_old = -1;		// The photon number of the last photon

int
ispy (p, n)
     PhotPtr p;			//The photon

     int n;			//The photon number

{
  int i, m;

  if (ispy_start == 0)
    {
      Log ("ispy file is not open, why are you trying to write to it?\n");
      exit (0);
    }

  if (p->grid == ispy_grid_old && n == ispy_phot_old)
    return (0);			// this photon has been seen in this cell previously

  i = ispy_grid_old = p->grid;
  ispy_phot_old = n;

  for (m = 0; m < ispy_ncells; m++)
    {
      if (i == ispy_cell[m])
	{
	  fprintf (ispy_ptr, "%2d %3d %8.2e %8.2e\n", ispy_cycle, i, p->w,
		   p->freq);
	  return (0);
	}
    }
  return (0);
}
